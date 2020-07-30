# hugoify.py
# Author: Nathan Grant

import requests
from bs4 import BeautifulSoup
import sys
import pandas as pd
import numpy as np
import re
from urllib.parse import urlparse
from multiprocess import Pool



class Hugoify:
    def __init__(self, data, num_workers=7):
        ## File destination must be of type excel
        self.num_workers = 7
        if isinstance(data, str):
            self.df = pd.read_excel(data)
        elif isinstance(data, pd.DataFrame):
            self.df = data
        else:
            raise Exception("Must be a DataFrame object or a path to an excel file")
    def get_missing_ids(self):
        ## First argument: Make function that inputs a enterez id and spits out a hugo id
        ## Second argument: Make a list of all of the missing
        missing = []
        values = self.df.values
        for i in range(len(values)):  #
            if type(values[i][0]) is type(3.6) or values[i][0] == "nan":
                enterez = values[i][1]
                missing.append({'index':i, 'enterez':enterez})

        p = Pool(self.num_workers)
        hugo_ids = p.map(get_id, missing)
        p.terminate()
        p.join()

        ## Change Hugo Ids to new ids
        number_found = 0
        for i in range(len(hugo_ids)):
            self.df.iloc[missing[i]['index'],0] = hugo_ids[i]
            if hugo_ids[i] != "nan":
                number_found+=1
        print("Found ", number_found, "/", len(missing), "of the missing values")
        return self.df

def hugo(df, num_workers=7):
    h = Hugoify(df, num_workers=num_workers)
    return h.get_missing_ids()

def get_id(missing):
    from urllib.parse import urlparse
    import requests
    from bs4 import BeautifulSoup
    import numpy as np
    URL = 'https://www.ncbi.nlm.nih.gov/gene/?term='

    def get_text(page):
        soup = BeautifulSoup(page.content, 'html.parser')

        # kill all script and style elements
        spans = soup.findAll("dd", {"class": "noline"})
        if len(spans) ==0:
            return 'nan'
        hugo = str(spans[0].contents[0])
        if hugo == "" or hugo == None or len(hugo) == 0:
            return 'nan'

        return hugo

    def get_hugo(URL):
        headers = requests.utils.default_headers()
        headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/83.0.4103.61 Safari/537.36',
            'referer': 'https://www.google.com/search?q=genecards+729&rlz=1C1CHBF_enUS776US776&oq=genecards+729&aqs=chrome..69i57j69i64j69i60l3.9462j0j7&sourceid=chrome&ie=UTF-8'
        })
        page = requests.get(URL, headers=headers)
        if page.status_code == 200:
            return get_text(page)
        else:
            return 'nan'

    newURL = newURL = urlparse(URL+str(missing['enterez']))
    value=get_hugo(newURL.geturl())
    return value

if __name__ == "__main__":
    pass
