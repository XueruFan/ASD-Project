# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 17:30:08 2024

@author: xueru
"""

import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from sklearn.svm import SVR
import shap
import matplotlib.pyplot as plt
import numpy as np  # 添加 numpy 用于剪裁 SHAP 值

# 定义路径
abideDir = 'E:/PhDproject/ABIDE'  
phenoDir = os.path.join(abideDir, "Preprocessed")
resDir = os.path.join(abideDir, "Analysis", "Cluster", "Spect513")
plotDir = os.path.join(abideDir, "Plot", "Cluster", "Spect513", "SHAP")
newDate = '240610'
resDate = '240315'

pheno = pd.read_csv(os.path.join(phenoDir, f"abide_A_all_{resDate}.csv"))
pheno.columns.values[0] = "participant"

cluster_file = os.path.join(resDir, f"Cluster_{newDate}.csv")
cluster = pd.read_csv(cluster_file)

All = pd.merge(cluster, pheno, on="participant", how="left")

# 定义目标变量列表
target_variables = [
    "FIQ", "VIQ", "PIQ", "ADOS_2_SEVERITY_TOTAL", "ADOS_2_TOTAL", 
    "ADOS_2_SOCAFFECT", "ADOS_2_RRB", "SRS_AWARENESS_RAW", 
    "SRS_COGNITION_RAW", "SRS_COMMUNICATION_RAW", "SRS_MOTIVATION_RAW", 
    "SRS_MANNERISMS_RAW", "SRS_TOTAL_RAW", "ADI_R_SOCIAL_TOTAL_A", 
    "ADI_R_VERBAL_TOTAL_BV", "ADI_R_NONVERBAL_TOTAL_BV", "ADI_R_RRB_TOTAL_C"
]

# 列出要测试的模型及其参数网格
models = {
    'LinearRegression': (LinearRegression(), {}),
    'Ridge': (Ridge(), {'alpha': [0.01, 0.1, 1, 10]}),
    'RandomForest': (RandomForestRegressor(random_state=941205), {'n_estimators': [50, 100, 200], 'max_depth': [5, 10, 20]}),
    'XGBoost': (XGBRegressor(random_state=941205), {'n_estimators': [50, 100], 'learning_rate': [0.01, 0.1, 0.2]}),
    'SVR': (SVR(), {'C': [0.1, 1, 10], 'epsilon': [0.01, 0.1, 1]})
}

# 通用函数来处理数据集（L 或 H）
def process_cluster_data(data, label):
    results = []
    
    for target in target_variables:
        T = data.dropna(subset=[target])  # 删除含有 NaN 的行
        
        continuous_features = T.loc[:, 'bankssts':'insula']  # 连续的特征列
        age_feature = T[['AGE_AT_SCAN']]  # 非连续的列（AGE_AT_SCAN）
        X = pd.concat([continuous_features, age_feature], axis=1)  # 合并特征矩阵
        
        y = T[target]  # 目标变量
        
        # 划分训练集和测试集
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=941205)
        
        best_model = None
        best_score = float('-inf')  # 用于追踪最优R²值
        best_model_name = ""
        
        # 循环训练每个模型，并进行超参数调优
        for name, (model, param_grid) in models.items():
            grid_search = GridSearchCV(estimator=model, param_grid=param_grid, scoring='r2', cv=5)
            grid_search.fit(X_train, y_train)
            
            y_pred = grid_search.best_estimator_.predict(X_test)  # 使用测试集评估模型性能
            r2 = r2_score(y_test, y_pred)
            mse = mean_squared_error(y_test, y_pred)
            
            results.append({
                'Target': target,
                'Model': name,
                'Best Parameters': grid_search.best_params_,
                'R2 Score': r2,
                'Mean Squared Error': mse
            })
            
            if r2 > best_score:  # 更新最佳模型
                best_score = r2
                best_model = grid_search.best_estimator_
                best_model_name = name
        
        # 输出每个目标的最优模型信息
        print(f"Best Model for {target}: {best_model_name}")
        print(f"Best R² Score: {best_score}\n")
        
        # 计算并绘制最佳模型的SHAP图
        try:
            if best_model_name in ['RandomForest', 'XGBoost']:
                explainer = shap.TreeExplainer(best_model)
            else:
                explainer = shap.KernelExplainer(best_model.predict, X_train)
            
            shap_values = explainer.shap_values(X_test)
            
            # 剪裁掉前1%和后1%的极端 SHAP 值
            q1 = np.percentile(shap_values, 1)
            q99 = np.percentile(shap_values, 99)
            shap_values_clipped = np.clip(shap_values, q1, q99)
            
            # 绘制并保存SHAP值图，去掉极端点并只显示前10个特征
            plt.figure()
            shap.summary_plot(shap_values_clipped, X_test, show=False, max_display=10, plot_type="dot")
            plt.title(f"SHAP Summary Plot for {target} ({best_model_name})")
            plt.savefig(os.path.join(plotDir, f"{label}_{target}_{best_model_name}.png"))
            plt.close()
            
        except Exception as e:
            print(f"Error generating SHAP plot for {target} with model {best_model_name}: {str(e)}")
    
    # 将结果保存为 CSV 文件
    results_df = pd.DataFrame(results)
    results_df.to_csv(os.path.join(resDir, f"{label}_Model_Results_{newDate}.csv"), index=False)

# 处理L和H数据集
L = All[All['clusterID'] == 1]
process_cluster_data(L, 'L')

H = All[All['clusterID'] == 2]  # 假设你要处理 clusterID 为 2 的数据
process_cluster_data(H, 'H')
