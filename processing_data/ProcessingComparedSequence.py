"""
Remove the conserved amino acids from the sequence based on their position in the original sequence
"""

from unit.path import Args
import os
from tqdm import tqdm

items = []
all_similar = []
primary_all_similar = []
final_similar = []
primary_final_similar = []
index_items = 1
uniprot_id = []
uniprot_id_similar_index = {}
primary_uniprot_id_similar_index = {}
sequence_alignment_results = Args().sequence_alignment_results
final_info_dict = dict()
with open(sequence_alignment_results, 'r') as contrast:
    text = contrast.read()
    l = text.split('\n')
    contrast.close()
for n in l:
    items.append(n.split('\t'))

for item in items[1:]:
    uniprot_id.append(item[0])
# print(uniprot_id)
for item in items[1:]:
    primary_all_similar.append(item[2])
for item_2 in primary_all_similar:
    length = 1
    str_item = ''
    list_item = []
    while length < len(item_2)-1:
        if item_2[length] != ',':
            str_item += item_2[length]
            length += 1
            if length == len(item_2)-1:
                list_item.append(int(str_item))
        else:
            list_item.append(int(str_item))
            length += 1
            str_item = ''
    primary_final_similar.append(list_item)
print(len(primary_final_similar))
k = 0
while k < len(primary_final_similar):
    primary_uniprot_id_similar_index[uniprot_id[k]] = primary_final_similar[k]
    k += 1

# 找到序列比对后多个百分百相似的区域坐标
# 取原序列信息
items = []
ite = []
it = []
ind = []
uniprot_ec = {}
uniprot_seq = {}
uniprot_topt = {}
this_dir, this_filename = os.path.split(__file__)
# print(this_dir)
tsv_data = os.path.join(this_dir, 'Topt249.tsv')
alignmented_information_sequnence = os.path.join(this_dir, 'Topt_rcaa218.tsv')

with open(tsv_data, 'r', encoding='utf-8') as contrast:
    text = contrast.read()
    line = text.split('\n')
    contrast.close()
for item in line:
    items.append(item.split('\t'))
for i in items:
    ite.append(i)
for j in ite[1:]:
    uniprot_ec[j[2]] = j[0]
    uniprot_seq[j[2]] = j[4]
    uniprot_topt[j[2]] = j[6]
# print(uniprot_seq)
# 在原序列中去除保守位点
select_seq_result = {}
for select_seq_key, select_seq_value in uniprot_seq.items():
    seq_str = ''
    if select_seq_key in primary_uniprot_id_similar_index.keys():
        for select_seq_index, select_seq_item in enumerate(select_seq_value):
            if select_seq_index in primary_uniprot_id_similar_index[select_seq_key]:
                continue
            else:
                seq_str += select_seq_item
        select_seq_result[select_seq_key] = seq_str
    else:
        continue


file_handle = open(alignmented_information_sequnence, mode='w', encoding='utf-8')
file_handle.write('EC number\t\t')
file_handle.write('Uniprot ID\t\t')
file_handle.write('Sequence\t\t')
file_handle.write('Topt\n')
for select_seq_result_key, select_seq_result_value in select_seq_result.items():
    file_handle.write(f'{uniprot_ec[select_seq_result_key]}\t\t')
    file_handle.write(f'{select_seq_result_key}\t\t')
    file_handle.write(f'{select_seq_result_value}\t\t')
    file_handle.write(f'{uniprot_topt[select_seq_result_key]}\n')
file_handle.close()
