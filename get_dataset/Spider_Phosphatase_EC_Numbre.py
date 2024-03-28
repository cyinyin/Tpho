import os
import requests
from bs4 import BeautifulSoup

headers = {
        'authority': 'cn.bing.com',
        'sec-ch-ua': '" Not;A Brand";v="99", "Google Chrome";v="97", "Chromium";v="97"',
        'sec-ch-ua-mobile': '?0',
        'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/97.0.4692.71 Safari/537.36',
        'sec-ch-ua-arch': '"x86"',
        'sec-ch-ua-full-version': '"97.0.4692.71"',
        'sec-ch-ua-platform-version': '"10.0.0"',
        'sec-ch-ua-bitness': '"64"',
        'sec-ch-ua-model': '',
        'sec-ch-ua-platform': '"Windows"',
        'accept': 'image/avif,image/webp,image/apng,image/svg+xml,image/*,*/*;q=0.8',
        'sec-fetch-site': 'same-origin',
        'sec-fetch-mode': 'no-cors',
        'sec-fetch-dest': 'image',
        'referer': 'https://cn.bing.com/?mkt=zh-CN',
        'accept-language': 'zh-CN,zh;q=0.9',
        'cookie': 'MUID=078F787B71F46FFD09CA68D070B76E50; MUIDB=078F787B71F46FFD09CA68D070B76E50; _EDGE_V=1; SRCHD=AF=NOFORM; SRCHUID=V=2&GUID=6ABD9AEC7BE248E2962C2E4E94E03E71&dmnchg=1; TTRSL=en; _tarLang=default=zh-Hans; _TTSS_OUT=hist=WyJlbiIsInpoLUhhbnMiXQ==; _TTSS_IN=hist=WyJpdCIsInJvIiwiZW4iLCJhdXRvLWRldGVjdCJd; BCP=AD=1&AL=1&SM=1; HOOKBLOCKINDICATOR=TRUE; ABDEF=V=13&ABDV=11&MRNB=1644823879884&MRB=0; _EDGE_S=SID=320E72F8C7BF611D2D9363B4C6B760CC&mkt=zh-cn; _SS=SID=320E72F8C7BF611D2D9363B4C6B760CC; SUID=M; ZHCHATSTRONGATTRACT=TRUE; ZHCHATWEAKATTRACT=TRUE; SRCHUSR=DOB=20210907&T=1644986902000&TPC=1644976209000; ipv6=hit=1644990503642&t=4; _UR=OMD=13289460522; SNRHOP=I=&TS=; _HPVN=CS=eyJQbiI6eyJDbiI6MTcsIlN0IjoyLCJRcyI6MCwiUHJvZCI6IlAifSwiU2MiOnsiQ24iOjE3LCJTdCI6MCwiUXMiOjAsIlByb2QiOiJIIn0sIlF6Ijp7IkNuIjoxNywiU3QiOjEsIlFzIjowLCJQcm9kIjoiVCJ9LCJBcCI6dHJ1ZSwiTXV0ZSI6dHJ1ZSwiTGFkIjoiMjAyMi0wMi0xNlQwMDowMDowMFoiLCJJb3RkIjowLCJHd2IiOjAsIkRmdCI6bnVsbCwiTXZzIjowLCJGbHQiOjAsIkltcCI6MjgyfQ==; SRCHHPGUSR=SRCHLANG=zh-Hans&BRW=NOTP&BRH=M&CW=332&CH=730&SW=1536&SH=864&DPR=1.25&UTC=480&DM=0&WTS=63780583702&HV=1644988141&BZA=0',
    }
url= 'https://www.brenda-enzymes.org/all_enzymes.php'

sp_dir = os.path.abspath('..')
sp_data_dir = os.path.join(sp_dir, 'data')
ec_number_tsv = os.path.join(sp_data_dir, 'ec_number.tsv')
file_handle = open(ec_number_tsv, mode='w', encoding='utf-8')

try:
    resp = requests.get(url=url, headers=headers, timeout=(6, 6))
    resp.encoding = resp.apparent_encoding
    if resp.status_code == 200:
        print(resp)
    else:
        print(f'url: {url} error code: {resp.status_code}')

except Exception as e:
    error_info = f'{e}'
    if error_info.find('port=443') != -1:
        print(f'{url} {e}')

    else:
        print(f'url: {url} error: {e}')

# 核心爬取代码
params = {"show_ram": 1}
response = requests.get(url=url, params=params, headers=headers)  # 访问url
EcnumberData=[]  # 定义数组
soup = BeautifulSoup(response.text, 'html.parser')  # 获取网页源代码
for row in soup.find_all('div', class_='row'):
    if row.select('div')[0].text[0:6] == '3.1.3.':
        EC_Number = row.select('div')[0].text  # EC_Number
        EcnumberData.append(EC_Number)

for EcnumberData_item in EcnumberData:
    file_handle.write(f'{EcnumberData_item}\n')
file_handle.close()