import urllib.request
from unit.path import Args
import os
from tqdm import tqdm

count = 0
items = []
ite = []
url = "https://www.brenda-enzymes.org/enzyme.php?ecno="
ec_number_path = Args().ec_number_tsv
phosphataseEcHtml = Args().phosphataseEcHtml

with open(ec_number_path, 'r') as contrast:
    text = contrast.read()
    line = text.split('\n')
    contrast.close()
for item in line:
    items.append(item.split('\t'))
for i in items[1:len(items)-1]:
    ite.append(i[0])


def getHtml(url):
    html = urllib.request.urlopen(url).read()
    return html


def saveHtml(file_name, file_content):
    filepath = os.path.join(phosphataseEcHtml, file_name.replace('/', '_') + ".html")
    # 注意windows文件命名的禁用符，比如 /
    with open(filepath, "wb") as f:
        print(file_name.replace('/', '_'))
        # 写文件用bytes而不是str，所以要转码
        f.write(file_content)

pbar = tqdm(total=len(ite), desc="process")
for i in ite:
    search_url = url+i
    pbar.update(1)
    html = getHtml(search_url)
    saveHtml(f"{i}", html)
pbar.close()