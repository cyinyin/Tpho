"""
Obtain a complete sequence dataset with only sequence differences from the non conserved amino acid dataset
"""

import os
from tqdm import tqdm


items = []
itemss = []
alignmented_information_sequnence_information = []
original_complete_information_sequence_information = []
uniprot_ids = []
this_dir, this_filename = os.path.split(__file__)
alignmented_information_sequnence = os.path.join(this_dir, '../data/Topt_rcaa218.tsv')
original_complete_information_sequence = os.path.join(this_dir, '../data/Topt249.tsv')
complete_information_sequence = os.path.join(this_dir, '../data/Topt218.tsv')
with open(alignmented_information_sequnence, 'r', encoding='utf-8') as contrast:
    text = contrast.read()
    line = text.split('\n')
    contrast.close()
for item in line:
    itemss.append(item.split('\t'))
for i in itemss:
    alignmented_information_sequnence_information.append(i)


bar = tqdm(total=len(alignmented_information_sequnence_information), desc='Process completion progress1')
for jj in alignmented_information_sequnence_information[1:]:
    bar.update(1)
    uniprot_ids.append(jj[2])
with open(original_complete_information_sequence, 'r', encoding='utf-8') as contrast:
    text = contrast.read()
    line = text.split('\n')
    contrast.close()
for item in line:
    items.append(item.split('\t'))
for i in items:
    original_complete_information_sequence_information.append(i)


file_handle = open(complete_information_sequence, mode='w', encoding='utf-8')
file_handle.write('EC number\t\t')
file_handle.write('Uniprot ID\t\t')
file_handle.write('Sequence\t\t')
file_handle.write('Topt\n')
bar = tqdm(total=len(original_complete_information_sequence_information), desc='Process completion progress1')
for uniprot_id in uniprot_ids:
    for j in original_complete_information_sequence_information[1:]:
        if uniprot_id == j[2]:
            file_handle.write(f'{j[0]}\t\t')
            file_handle.write(f'{j[2]}\t\t')
            file_handle.write(f'{j[4]}\t\t')
            file_handle.write(f'{j[6]}\n')

