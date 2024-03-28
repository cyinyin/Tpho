from unit.path import Args
import os
from tqdm import tqdm

datalist = []
alldatalist = []
idlist = []
seqlist = []
this_dir, this_filename = os.path.split(__file__)
predict_validation_data_tsv = os.path.join(this_dir, '../data/Experimentalvalidationdata/Topt_rcaa10.tsv')
with open(predict_validation_data_tsv, 'r', encoding='utf-8') as contrast:
    text = contrast.read()
    line = text.split('\n')
    contrast.close()
for item in line:
     datalist.append(item.split('\t'))
for i in  datalist:
    alldatalist.append(i)

for j in alldatalist[1:]:
    idlist.append(j[2])
    seqlist.append(j[4])

predict_validation_data_fasta_path = os.path.join(this_dir,
                                                  '../data/Experimentalvalidationdata/Topt_rcaa10.fasta')
predict_validation_data_fasta = open(predict_validation_data_fasta_path, mode='w', encoding='utf-8')
for id_index, id_item in enumerate(idlist):
    predict_validation_data_fasta.write(f'>{id_item}\n')
    predict_validation_data_fasta.write(f'{seqlist[id_index]}\n')
