from unit.path import Args
from get_dataset import parser
from pkg_resources import resource_stream
import requests
import os
from tqdm import tqdm


def get_html(url):
    """
    func: call url and return html，html.text可以获取内容
    note: 如果状态码错误返回None, 如果是服务器访问过于频繁，拒绝访问，返回443，其他错误返回None
    """
    header = {
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
    try:
        resp = requests.get(url=url, headers=header, timeout=(6, 6))
        resp.encoding = resp.apparent_encoding
        if resp.status_code == 200:
            return resp
        else:
            # print(f'url: {url} error code: {resp.status_code}')
            return None
    except Exception as e:
        error_info = f'{e}'
        if error_info.find('port=443') != -1:
            # print(f'{url} {e}')
            return 443
        else:
            # print(f'url: {url} error: {e}')
            return None

PhosphataseEnzymeInformation_tsv = Args().PhosphataseEnzymeInformation_tsv
phosphataseEcHtml = Args().phosphataseEcHtml
file_handle = open(PhosphataseEnzymeInformation_tsv, mode='w', encoding='utf-8')
file_handle.write('EC number\t\t')
file_handle.write('Uniprot ID\t\t')
file_handle.write('Sequence\t\t')
file_handle.write('Topt\n')
_phosphataseEcHtml_path = Args().phosphataseEcHtml
for root, dirs, files in os.walk(_phosphataseEcHtml_path):
    pbar = tqdm(total=len(files), desc='process')
    for file in files:
        pbar.update(1)
        # html_path = os.path.join(allEcHtml, f"{file}")
        # print(html_path)
        EXAMPLE_PAGE = resource_stream(__name__, f'../data/phosphataseEcHtml/{file}').name
        # filepath = parser.EXAMPLE_PAGE
        soup_obj = parser.open_ec(EXAMPLE_PAGE)
        if parser.TemperatureOptimum(soup_obj).get_data() is not None:
            for items in parser.TemperatureOptimum(soup_obj).get_data().values():
                for item_key in items.keys():
                    text = ""
                    url = f'https://www.uniprot.org/uniprot/{item_key}.fasta'
                    resp = get_html(url)
                    if resp is not None:
                        if resp == 443:
                # 服务器拒绝访问，返回443
                            print(f'{item_key},{url},error: return 443')
                        else:
                     # 返回正常
                            text = resp.text
                        # 返回信息是否有效
                            if len(text) == 0:
                                print(f'{item_key},{url},error: return html is Null')
                                continue
                # 去掉第一行提示信息
                            text = text[text.find('\n') + 1:]
                            text_line = text.replace('\n', '')
                            name = file[:file.find('h')-1]
                            file_handle.write(f'{item_key}\t\t')
                            file_handle.write(f"{text_line}\t\t")
                            file_handle.write(f'{items[item_key][0]}\n')
                    else:
                        print(f'{item_key},{url},error: resp.status_code or exception')
pbar.close()
file_handle.close()