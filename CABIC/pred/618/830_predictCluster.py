# -*- coding: utf-8 -*-
"""
用于用之前训练的非线性核模型对新数据进行预测
@author: xueru
"""
import os
import pandas as pd
import joblib

# 定义路径
abideDir = 'E:/PhDproject/ABIDE/Analysis/Cluster/Spect618'  # Windows
cabicDir = 'E:/PhDproject/CABIC/result'
newDate = "240928"

# 1. 加载之前保存的模型
model_file = os.path.join(abideDir, 'trained_svm_model.pkl')
loaded_model = joblib.load(model_file)

# 2. 读取新的数据集
new_data = pd.read_csv(os.path.join(cabicDir, f"cabic_centile_{newDate}.csv"))
new_data = new_data[(new_data['sex'] == 'Male') & (new_data['dx'] == 'ASD')]
new_data = new_data.drop(columns=['sex', 'dx'])

# 3. 检查新数据集中的特征是否与训练模型的特征匹配
# 获取模型期望的特征列表
model_features = loaded_model.best_estimator_.feature_names_in_

# 检查新数据集是否包含模型所需的特征
if set(model_features).issubset(new_data.columns):
    new_data_selected = new_data[model_features]
else:
    raise ValueError("新数据集缺少部分特征，请确保新数据与训练数据包含相同的特征")

# 4. 使用加载的模型进行预测
y_pred_new_data = loaded_model.predict(new_data_selected)

# 输出预测结果
new_data['predicted_cluster'] = y_pred_new_data

# 保存预测结果到 CSV 文件
output_file = os.path.join(cabicDir, f'cabic_cluster_predictions_{newDate}.csv')
new_data.to_csv(output_file, index=False)

print(f"预测结果已保存到: {output_file}")

