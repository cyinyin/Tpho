from predicate_phosphatase_temperature import Tpho
from unit.path import Args
import prettytable


def result_validation():

    '''
    Obtain data from PNAS articles and use the final model for prediction
    '''

    # ================================#
    #  Obtaining  data
    predict_validation_data_tsv = Args().test_data_tsv
    predict_validation_data_fasta = Args().test_data_fasta
    all_information = []
    experimental_value = dict()
    uniprot_id = []
    predicted_value = dict()

    with open(predict_validation_data_tsv, 'r') as contrast:
        text = contrast.read()
        l = text.split('\n')
        contrast.close()
    for n in l[1:]:
        all_information.append(n.split('\t'))
    for item in all_information:
        uniprot_id.append(item[2])
        experimental_value[item[2]] = item[6]

    for item, result in enumerate(Tpho.pred_phosphatase_topt(predict_validation_data_fasta).topt_pred):
        predicted_value[uniprot_id[item]] = result

    # 创建表格对象
    table = prettytable.PrettyTable()

    # 添加表头
    table.field_names = ["uniprot id", "T opt", "predict value"]  # experimental value

    for item in uniprot_id:
        # 添加数据行
        table.add_row([item, experimental_value[item], predicted_value[item]])
    # 输出表格
    print(table)


if __name__ == "__main__":
    result_validation()