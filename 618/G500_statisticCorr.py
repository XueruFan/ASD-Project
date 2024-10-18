# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 09:29:39 2024
使用最小二乘法（OLS）进行线性回归来分析ASD症状与一系列脑指标之间的关系
@author: xueru
"""

# 清除所有变量
from IPython import get_ipython
get_ipython().magic('reset -sf')

import statsmodels.api as sm
import pandas as pd
import os
import numpy as np
from statsmodels.stats.multitest import multipletests

abideDir = 'E:/PhDproject/ABIDE'  # Windows
phenoDir = os.path.join(abideDir, "Preprocessed")
resDir = os.path.join(abideDir, "Analysis", "Statistic", "Gmm618")
plotDir = os.path.join(abideDir, "Plot", "Cluster", "Gmm618", "Corr")
clustDir = os.path.join(abideDir, "Analysis", "Cluster", "Gmm618")
newDate = '240610'
resDate = '240315'

pheno = pd.read_csv(os.path.join(phenoDir, f"abide_A_all_{resDate}.csv"))
pheno.columns.values[0] = "participant"

cluster = pd.read_csv(os.path.join(clustDir, f"Cluster_{newDate}.csv"))
cluster.columns.values[2:] = [col + '_centile' for col in cluster.columns[2:]]

All = pd.merge(cluster, pheno, on="participant", how="left")
All.loc[All['clusterID'] == 1, 'clusterID'] = "L"
All.loc[All['clusterID'] == 2, 'clusterID'] = "H"

All['SITE_ID'] = All['SITE_ID'].replace({
    "ABIDEII-NYU_1": "NYU",
    "ABIDEII-NYU_2": "NYU",
    "ABIDEII-KKI_1": "KKI",
    "ABIDEII-SDSU_1": "SDSU",
    "ABIDEII-UCLA_1": "UCLA",
    "UCLA_1": "UCLA",
    "UCLA_2": "UCLA",
    "ABIDEII-GU_1": "GU",
    "ABIDEII-UCD_1": "UCD",
    "ABIDEII-EMC_1": "EMC",
    "TRINITY": "TCD",
    "ABIDEII-TCD_1": "TCD",
    "ABIDEII-USM_1": "USM",
    "ABIDEII-IU_1": "IU",
    "ABIDEII-U_MIA_1": "UMIA",
    "ABIDEII-ETH_1": "ETH",
    "UM_1": "UM",
    "UM_2": "UM",
    "ABIDEII-OHSU_1": "OHSU",
    "STANFORD": "SU1",
    "ABIDEII-SU_2": "SU2",
    "LEUVEN_2": "KUL",
    "CALTECH": "CALT"
})

# 获取从第3列到最后一列的列名
names_brain = cluster.columns[2:].tolist()

# 定义认知数据列名的列表
names_cog = ["ADOS_2_SOCAFFECT", "ADOS_2_RRB", "ADOS_2_TOTAL", "ADOS_2_SEVERITY_TOTAL"]

# 将 names_brain 和 names_cog 合并到 names_col 列表
names_col = ['clusterID', 'SITE_ID'] + names_brain + names_cog

# 按照 names_col 中的列选择 DataFrame
All = All[names_col]

All[names_cog] = All[names_cog].apply(pd.to_numeric, errors='coerce')
All[names_cog] = All[names_cog].applymap(lambda x: np.nan if x < 0 else x)

L = All[All['clusterID'] == "L"]
L = L.iloc[:, 1:]
H = All[All['clusterID'] == "H"]
H = H.iloc[:, 1:]

# 定义空列表用于存储回归结果
regression_results_L = []

# 开始循环：对 names_cog 中的每个 ASD 测量变量
for brain_var in names_brain:
    for cog_var in names_cog:
        # 准备自变量和因变量
        X_site = pd.get_dummies(L[['SITE_ID']], drop_first=True)  # 只包含 SITE_ID 的哑变量
        X_site = X_site.astype(int)  # 将哑变量转换为整数类型（0 或 1）

        # 然后加入当前的 brain_var 作为自变量
        X = X_site.copy()  # 创建一个新的 DataFrame 包含哑变量
        X[brain_var] = L[brain_var]  # 加入 brain_var
        X = sm.add_constant(X)
        y = L[cog_var]

        # 确保因变量和自变量都是数值类型，并且没有 NaN
        X = X.apply(pd.to_numeric, errors='coerce')
        y = pd.to_numeric(y, errors='coerce')

        # 移除 NaN 值
        X = X.dropna()
        y = y.loc[X.index]
        
        y = y.dropna()
        X = X.loc[y.index]

        # 检查是否有足够的数据点进行回归
        if len(X) > 10 and  y.var() > 0:  # 确保有足够的数据和方差
            model = sm.OLS(y, X).fit()
            result = {
                'brain_var': brain_var,
                'cog_var': cog_var,
                'p_value': model.pvalues[brain_var],  # 当前 brain_var 的 p 值
                'coef': model.params[brain_var],      # 当前 brain_var 的回归系数
                'R-squared': model.rsquared,          # 模型的 R-squared 值
                'adj_R-squared': model.rsquared_adj   # 调整后的 R-squared 值
            }
            regression_results_L.append(result)

# 将结果转为 DataFrame 进行查看
results_L = pd.DataFrame(regression_results_L)

# FDR 校正
p_values = results_L['p_value']
reject, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

# 将校正后的 p 值添加到结果 DataFrame
results_L['FDR_corrected_p_value'] = pvals_corrected
results_L['Significant (FDR)'] = reject

