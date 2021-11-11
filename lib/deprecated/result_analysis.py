# 结果正确性验证完成
# python /home/lqh/Codes/Python/RNA-SSNV/lib/result_analysis.py \
# --explain_type dataset \
# --data_info /home/lqh/Codes/Python/RNA-SSNV/output/LUAD.table \
# --model_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.model \
# --one_hot_encoder_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.one_hot_encoder \
# --training_columns_path /home/lqh/Codes/Python/RNA-SSNV/model/exon_RNA_analysis_newer.training_data_col \
# --explain_plot_path /home/lqh/Codes/Python/RNA-SSNV/results/fig.3.LUAD.svg


import pandas as pd
import shap
import joblib
import matplotlib.pyplot as plt

import argparse

parser=argparse.ArgumentParser(description="A discriminate model construction pipeline for RNA-SSNV.")
# parser.add_argument('--raw_RNA_mutations', '-r' ,choices=[5,10,20],default=5,type=int,help='Number of epochs.')
# Generic parameter
parser.add_argument('--DNA_only_info', help='Type of explanation')
parser.add_argument('--RNA_DNA_overlap_info', help='Data info waiting for explanation')
parser.add_argument("--model_path", help="Path for constructed model.")
parser.add_argument("--one_hot_encoder_path", help="Path for one-hot encoder.")
parser.add_argument("--training_columns_path", help="Path for constructed model.")
parser.add_argument("--explain_plot_path", help="Path for explanation plot.")
# specific parameter
parser.add_argument('--specified_features', nargs='+', help='Specified features for plotting.')
parser.add_argument('--explain_row_index', type=int, help='Specified row index for explaination (0-based).')

args = parser.parse_args()

# 工具函数6：根据给定机器学习模型，对给定数据行进行预测，获得预测信息
def data_row_interpret(model, dataset, row_index):
    data_row = dataset.loc[row_index,]
    shap_explainer = shap.TreeExplainer(model)
    shap_value = shap_explainer.shap_values(data_row)
    return shap.force_plot(shap_explainer.expected_value[1], shap_value[1],
                    feature_names=data_row.index, matplotlib=True, show=False)

# 工具函数7：根据给定机器学习模型，对给定数据集合进行预测，获得预测信息
def dataset_interpret(model, dataset):
    shap_explainer = shap.TreeExplainer(model)
    # 计算数据集所对应的所有SHAP值
    shap_values = shap_explainer.shap_values(dataset, approximate=True)
    # 为数据集中每个样本绘制其对应特征的SHAP值，这可以更好地理解整体模式，并允许发现预测异常值
    return shap.summary_plot(shap_values[1], dataset, show=False)

# 工具函数8：根据给定机器学习模型，对给定数据集合和进行预测，获得预测信息
def features_dataset_interpret(model, dataset, features):
    shap_explainer = shap.TreeExplainer(model)
    # 获取features所对应的index情况
    dataset_columns_list = dataset.columns.tolist()
    index_list = []
    for feature in features:
        index_list.append(dataset_columns_list.index(feature))
    # 计算数据集所对应的所有SHAP值
    shap_values = shap_explainer.shap_values(dataset, approximate=True)
    # 为数据集中每个样本绘制其对应特征的SHAP值，这可以更好地理解整体模式，并允许发现预测异常值
    return shap.summary_plot(shap_values[1][:, index_list], dataset.iloc[:, index_list], show=False)

# 无关函数
def heart_disease_risk_factors(model, patient):

    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(patient)
    shap.initjs()

    return shap.force_plot(explainer.expected_value[1],shap_values[1],\
        patient,matplotlib=True,show=False)
# plt.clf()
# heart_disease_risk_factors(model, data_for_prediction)
# plt.savefig("gg.png",dpi=150, bbox_inches='tight')

if __name__ == '__main__':
    # 初始化绘图
    plt.clf()
    # 读取相应数据
    training_col = pd.read_table(args.training_columns_path)["0"]
    data_info = pd.read_table(args.data_info)[training_col]
    # 读取one-hot编码模型、机器学习模型、训练数据列名
    f1, f2 = open(args.one_hot_encoder_path, 'rb'), open(args.model_path, 'rb')
    enc, rf_gcv = joblib.load(f1), joblib.load(f2)
    # 根据情况不同进行分析
    if args.explain_type == "dataset":
        if args.specified_features is not None:
            features = args.specified_features
            features_dataset_interpret(rf_gcv, data_info, features)
            plt.savefig(args.explain_plot_path, format='svg', dpi=1000, bbox_inches='tight')
        else:
            dataset_interpret(rf_gcv, data_info)
            plt.savefig(args.explain_plot_path, format='svg', dpi=1000, bbox_inches='tight')
    elif args.explain_type == "datarow":
        data_row_interpret(rf_gcv, data_info, args.explain_row_index)
        plt.savefig(args.explain_plot_path, format='svg', dpi=1000, bbox_inches='tight')