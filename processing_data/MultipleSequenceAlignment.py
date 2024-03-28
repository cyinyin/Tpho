"""
Multiple sequence alignment to identify conserved amino acids in the sequence
"""

from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
from unit.path import Args


def read_fasta(fasta):
    with open(fasta, 'r') as fast:
        headers = []
        for line in fast:
            if line.startswith('>'):
                substr = ''
                for i in range(line.find('|') + 1, line.rfind('|')):
                    substr += line[i]
                headers.append(substr)
    return headers


sequence_alignment_results_intermediate_file = Args().sequence_alignment_results_intermediate_file
# 读取给定序列
# given_seq = "MMKMKLARFLLIFIIFSSMFPFANASPSGVRNVIILIGDGMGFSQLQLTKLVYGHLNMEDFPYTGIELTDSLSGEVTDSAAAGTAIATGVKTYNRMISTTNVTGKLVNLTTLLEIAQMLGKATGLVTTTRITHATPAVFASHVPDRDMEEEIARQLILHNVTVLMGGGREKFSEEVLKLAEDYGYSIVYTREDLEKVKDGKVLGLFAEGHLPYVLDRSEEDVSLLEMTKKAIEILEKNPNGFFLMIEGGRIDHACHANDVASIVAETKEFDDVVGYVLDYARRRGDTLVIVLADHETGGLGIGLNYGHSVDIDSIRRIDASIEEMSKEIKSGGDIRDVIRRHTGLELTDEEVKEIEEAKNSTNKYALGNIIGEIISKKLGVGFVSHKHTGEPVPLLAYGPGAENFVGFKHHVDTAKVIAKLMIFGDRSISFTIKGVSKIKGDVTGDYRVDERDAYATLMLLLGDLVDTELENIADMDNNGIIDLLDVMAILQASS"
# seq_record = SeqRecord(Seq(given_seq), id="given_seq")

# 从FASTA文件中读取参考蛋白序列
# with open("../data_fasta/new_3.1.3.fasta", "a") as ref_file:
# SeqIO.write(seq_record, ref_file, "fasta")


clustalo_exe = 'D:/predicate_phosphatase_temperature/processing_data/clustal-omega-1.2.2-win64/clustalo.exe'
cline = ClustalOmegaCommandline(clustalo_exe)
cline.infile = "../data/compare_phosphatase.fasta"
cline.outfile = "../data/compared_phosphatase.fasta"
# cline.outfmt = "clustal"

# # 使用Clustal Omega进行多序列比对
# clustalomega_cline = ClustalOmegaCommandline(
#     infile="../data_fasta/new_3.1.3.fasta",
#     outfile="../data_fasta/aligned.fasta",
#     verbose=True,
#     auto=True
# )
_, _ = cline()

# 读取比对结果
alignment = AlignIO.read("../data/compared_phosphatase.fasta", "fasta")

# print('alignment[:, 1]')
# print(len(alignment[:, 2]))
# print(len(alignment[:774]))
# 找到相似性高的区域
threshold = 1  # 设定阈值，根据需求调整
Similar_areas = []
Similar_amino_acid = []
seq_name = []
for i in range(alignment.get_alignment_length()):
    column = alignment[:, i]  # 对每一列做比对

    given_seq_char = column[-1]
    # print('given_seq_char')
    # print(given_seq_char)
    if given_seq_char != "-":  # 跳过空白区域
        same_char_count = column.count(given_seq_char)
        similarity = same_char_count / len(column)
        if similarity >= threshold:
            Similar_areas.append(i)
            Similar_amino_acid.append(given_seq_char)
            print(f"高相似性区域：位置 {i + 1}, 氨基酸 {given_seq_char}")
file_handle = open(sequence_alignment_results_intermediate_file, mode='w', encoding='utf-8')
file_handle.write('Uniprot ID\t\t')
file_handle.write('Similar_index\t\t')
file_handle.write('Similar_areas\n')
for name in read_fasta("../data/compare_phosphatase.fasta"):
    seq_name.append(name)
for len_alignment in range(len(alignment)):
    given_seq = list(alignment[len_alignment])
    final_given_seq = []
    Similar_index = []
    for index in Similar_areas:
        given_seq[index] = 0
    for item in given_seq:
        if item != "-":
            final_given_seq.append(item)
    for index, item in enumerate(final_given_seq):
        if item == 0:
            Similar_index.append(index)
    file_handle.write(f'{seq_name[len_alignment]}\t\t')
    file_handle.write(f'{Similar_index}\t\t')
    file_handle.write(f'{Similar_amino_acid}\n')

file_handle.close()
# 清理临时文件
# os.remove("reference.fasta")
# os.remove("new_3.1.3tt.fasta")
