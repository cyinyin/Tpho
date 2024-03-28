import os
import sys


class Args:
    def __init__(self):
        # 获取当前文件路径
        self._curr_dir, _curr_file_name = os.path.split(os.path.abspath(__file__))
        # 根目录
        self._root = os.path.abspath(os.path.join(self._curr_dir, ".."))
        # print(self._curr_dir)
        # 源文件
        self._uniprot_sequence_temperature_tsv = 'uniprot_sequence_temperature.tsv'
        self._ec_number_tsv = 'ec_number.tsv'
        self._allEnzymeInformation_tsv = 'allEnzymeInformation.tsv'
        self._PhosphataseEnzymeInformation_tsv = 'phosphatase_information.tsv'
        self._sequence_alignment_results_intermediate_file = 'sequence_alignment_results_intermediate_file.tsv'
        self._sequence_alignment_results = 'sequence_alignment_results.tsv'
        self._predict_validation_data_fasta = 'Topt_rcaa10.fasta'
        self._predict_validation_data_tsv = 'Topt_rcaa10.tsv'
        self._predict_validation_data_complete_sequence_tsv = 'Topt10.tsv'
        self._predict_validation_data_complete_sequence_fasta = 'Topt10.fasta'
        self._PNAS_validation_data_xlsx = 'PNAS_validation_data.xlsx'
        self._test_data_fasta = 'test_data.fasta'
        self._test_data_tsv = 'test_data.tsv'
        self._complete_experimental_data_fasta = 'complete_experimental_data.fasta'
        self._complete_experimental_data_tsv = 'complete_experimental_data.tsv'
        self._allEcHtml = 'allEcHtml'
        self._phosphataseEcHtml = 'Phosphataseechtml'
        self.data_path = os.path.join(self._root, 'data')
        self.Experimentalvalidationdata_path = os.path.join(self.data_path, 'Experimentalvalidationdata')
        self.test_path = os.path.join(self._root, 'test')
        self.testdata_path = os.path.join(self.test_path, 'test_data')
        self.unit_path = os.path.join(self._root, 'unit')
        self.uniprot_sequence_temperature_tsv = os.path.join(self.data_path, self._uniprot_sequence_temperature_tsv)
        self.ec_number_tsv = os.path.join(self.data_path, self._ec_number_tsv)
        self.allEnzymeInformation_tsv = os.path.join(self.data_path, self._allEnzymeInformation_tsv)
        self.PhosphataseEnzymeInformation_tsv = os.path.join(self.data_path, self._PhosphataseEnzymeInformation_tsv)
        self.sequence_alignment_results_intermediate_file = os.path.join(self.data_path, self._sequence_alignment_results_intermediate_file)
        self.sequence_alignment_results = os.path.join(self.data_path, self._sequence_alignment_results)
        self.predict_validation_data_fasta = os.path.join(self.Experimentalvalidationdata_path, self._predict_validation_data_fasta)
        self.predict_validation_data_tsv = os.path.join(self.Experimentalvalidationdata_path, self._predict_validation_data_tsv)
        self.predict_validation_data_complete_sequence_tsv = os.path.join(self.Experimentalvalidationdata_path, self._predict_validation_data_complete_sequence_tsv)
        self.predict_validation_data_complete_sequence_fasta = os.path.join(self.Experimentalvalidationdata_path, self._predict_validation_data_complete_sequence_fasta)
        self.PNAS_validation_data_xlsx = os.path.join(self.Experimentalvalidationdata_path, self._PNAS_validation_data_xlsx)
        self.test_data_fasta = os.path.join(self.testdata_path, self._test_data_fasta)
        self.test_data_tsv = os.path.join(self.testdata_path, self._test_data_tsv)
        self.complete_experimental_data_fasta = os.path.join(self.testdata_path, self._complete_experimental_data_fasta)
        self.complete_experimental_data_tsv = os.path.join(self.testdata_path, self._complete_experimental_data_tsv)
        self.allEcHtml = os.path.join(self.data_path, self._allEcHtml)
        self.phosphataseEcHtml = os.path.join(self.data_path, self._phosphataseEcHtml)
if __name__ =='__main__':
    Args()