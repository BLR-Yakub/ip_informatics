# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 11:46:29 2022

@author: EHOUT
"""

from selenium import webdriver
import time
from selenium.webdriver.common.keys import Keys
from datetime import datetime
from selenium.webdriver.support.ui import Select

from selenium.webdriver.chrome.options import Options
import os
import zipfile
import shutil
from Bio.Seq import Seq

chrome_options = webdriver.ChromeOptions()
chrome_options.add_argument("--incognito")

driver = webdriver.Chrome("C://project_IP//chromedriver.exe",chrome_options=chrome_options)



dicts={'Ala':'a', 'Asx':'b', 'Cys':'c',
       'Asp':'d', 'Glu':'e', 'Phe':'f',
       'Gly':'g', 'His':'h', 'Ile':'i',
       'Lys':'k', 'Leu':'l', 'Met':'m',
       'Asn':'n', 'Pro':'p', 'Gln':'q',
       'Arg':'r', 'Ser':'s', 'Thr':'t',
       'Sec':'u', 'Val':'v', 'Trp':'w',
       'Xaa':'x', 'Tyr':'y', 'Glx':'z'}

lis_amin=[]
for key in dicts:
    lis_amin.append(key)
    
lis_chk=['a', 't', 'g', 'c']

def get_us_patents_1(patent_no):
    url="https://appft1.uspto.gov/netahtml/PTO/srchnum.html"
    driver.get(url)

    search_box=driver.find_element_by_xpath('//*[@id="qry"]')
    search_box.send_keys(patent_no)
    #20090023654
    search_box.send_keys(Keys.ENTER)

    lnk=driver.find_element_by_xpath('/html/body/table/tbody/tr[2]/td[3]/a')
    driver.get(lnk.get_attribute('href'))

    el = driver.find_elements_by_tag_name('p')
    ans=-1
    for i in range(len(el)):
        tsxt=el[i].text
        if(tsxt.find('Sequence CWU')!=-1):
            ans=i
            break

    import re
    data=el[ans].text    


    name=[]
    new_data=[]
    res = re.split('DNA|PRT', data)
    for i in range(1,len(res)):
        for j in range(len(res[i])-2):
            if(res[i][j].isdigit()==True and res[i][j+1].isalpha()==True):
                
                if(res[i][j+1].isupper()):
                    cnt=0
                    if res[i][j+1:j+4] in lis_amin:
                        name.append(res[i][:j])
                        new_data.append(res[i][j:])
                        break
                    
                else:
                    cnt=0
                    for k in range(j+1,j+10):
                        if(res[i][k]!=' '):
                            if res[i][k] not in lis_chk:
                                cnt+=1
                                break
                    if(cnt==0):
                        name.append(res[i][:j])
                        new_data.append(res[i][j:])
                        break
            elif(res[i][j].isdigit()==True and res[i][j+2].isalpha()==True):
                
                
                if(res[i][j+2].isupper()):
                    cnt=0
                    if res[i][j+2:j+5] in lis_amin:
                        name.append(res[i][:j])
                        new_data.append(res[i][j:])
                        break
                    
                else:
                    cnt=0
                    for k in range(j+2,j+10):
                        if(res[i][k]!=' '):
                            if res[i][k] not in lis_chk:
                                cnt+=1
                                break
                    if(cnt==0):
                        name.append(res[i][:j])
                        new_data.append(res[i][j:])
                        break
                                
                

    fin_data=[]
    for i in range(len(new_data)):
        fin_data.append(re.sub(r'[^a-zA-Z ]+', '', new_data[i]))

    seq=[]
    for i in range(len(fin_data)):
        ix=fin_data[i].find("SEQ ID NO")
        if(ix!=-1):
            fin_data[i]=fin_data[i][:ix]
        x=fin_data[i].split()
        temp=""
        if(fin_data[i].islower()):
            for j in range(len(x)):
                temp+=x[j]
        else:
            for j in range(len(x)):
                try:    
                    
                    temp+=exp_to_short(x[j]).upper()
                    
                except:
                      pass
                  
        seq.append(temp)
        
    dicts={"Name": name, "Sequence": seq}
    df=pd.DataFrame(dicts)
    return df

def get_us_patents_2(patent_no):
    url="https://patft.uspto.gov/netahtml/PTO/srchnum.htm"
    driver.get(url)

    time.sleep(3)
    search_box=driver.find_element_by_xpath('//*[@id="qry"]')
    search_box.send_keys(patent_no)
    #20090023654
    search_box.send_keys(Keys.ENTER)

    time.sleep(3)
    el = driver.find_element_by_tag_name('body')
    el2=el.text
    ans=el2.find('SEQUENCE LISTINGS')
    '''for i in range(len(el2)):
        if(el2[i]=='SEQUENCE LISTINGS'):
            ans=i+1
            break
        
    for i in range(ans,len(el2)):
        if(len(el2[i])>10):
            ans=i
            break'''


    data=el2[ans:] 


    name=[]
    new_data=[]
    res = re.split('DNA|PRT', data)
    for i in range(1,len(res)):
        for j in range(len(res[i])-2):
            if(res[i][j].isdigit()==True and res[i][j+1].isalpha()==True):
                
                if(res[i][j+1].isupper()):
                    cnt=0
                    if res[i][j+1:j+4] in lis_amin:
                        name.append(res[i][:j])
                        new_data.append(res[i][j:])
                        break
                    
                else:
                    cnt=0
                    for k in range(j+1,j+10):
                        if(res[i][k]!=' '):
                            if res[i][k] not in lis_chk:
                                cnt+=1
                                break
                    if(cnt==0):
                        name.append(res[i][:j])
                        new_data.append(res[i][j:])
                        break
            elif(res[i][j].isdigit()==True and res[i][j+2].isalpha()==True):
                
                
                if(res[i][j+2].isupper()):
                    cnt=0
                    if res[i][j+2:j+5] in lis_amin:
                        name.append(res[i][:j])
                        new_data.append(res[i][j:])
                        break
                    
                else:
                    cnt=0
                    for k in range(j+2,j+10):
                        if(res[i][k]!=' '):
                            if res[i][k] not in lis_chk:
                                cnt+=1
                                break
                    if(cnt==0):
                        name.append(res[i][:j])
                        new_data.append(res[i][j:])
                        break
                                
                

    fin_data=[]
    for i in range(len(new_data)):
        fin_data.append(re.sub(r'[^a-zA-Z ]+', '', new_data[i]))

    seq=[]
    for i in range(len(fin_data)):
        #ix=fin_data[i].find("SEQ ID NO")
        #fin_data[i]=fin_data[i][:ix]
        x=fin_data[i].split()
        temp=""
        if(fin_data[i].islower()):
            for j in range(len(x)):
                temp+=x[j]
        else:
            for j in range(len(x)):
                try:    
                    
                    temp+=exp_to_short(x[j]).upper()
                    
                except:
                      pass
                  
        seq.append(temp)
        
    dicts={"Name": name, "Sequence": seq}
    df=pd.DataFrame(dicts)
    return df
    







    
 


