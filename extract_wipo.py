# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 15:11:03 2021

@author: EHOUT

"""

from selenium import webdriver
import time
import os
from selenium.webdriver.common.keys import Keys
from datetime import datetime
from selenium.webdriver.support.ui import Select
from Bio.Seq import Seq

from selenium.webdriver.chrome.options import Options

chrome_options = webdriver.ChromeOptions()
chrome_options.add_argument("--incognito")

driver = webdriver.Chrome("C://project_IP//chromedriver.exe",chrome_options=chrome_options)

# get sequences from wipo database
def get_seq_wipo(pat_no):
    url="https://patentscope.wipo.int/search/en/search.jsf"
    driver.get(url)
    
    time.sleep(5)
    
    search_box=driver.find_element_by_id('simpleSearchForm:fpSearch:input')
    search_box.send_keys(pat_no)
    
    #WO2018081194
    search_box.send_keys(Keys.ENTER)
    
    time.sleep(5)
    
    # find the column having documents
    m=0
    ts=driver.find_elements_by_xpath('//*[@id="detailMainForm:MyTabViewId"]/ul/li')
    for i in range(len(ts)):
        if(ts[i].text=='Documents'):
            m=i+1
            
    try:
        #select the documents
        documents=driver.find_element_by_xpath('//*[@id="detailMainForm:MyTabViewId"]/ul/li['+str(m)+']/a')
        documents.send_keys(Keys.ENTER)
    except Exception:
        pass
    
    time.sleep(5)
    
    x=driver.find_elements_by_class_name('ui-datatable')
    iD=x[1].get_attribute('id') #get the second table
    
    p='//*[@id="'+iD+'"]/div[2]/table/tbody/tr' #rows of second table
    table=driver.find_elements_by_xpath(p)
    
    #find the row in which sequence listing text file is present
    fnd=0
    for i in range(len(table)):
        s=table[i].text
        if(s.find('Sequence Listing')!=-1 and s.find('TXT')!=-1):
            fnd=i+1
            break
    
    time.sleep(5)
    
    #download the text file
    sind=driver.find_element_by_xpath('//*[@id="'+iD+'"]/div[2]/table/tbody/'+'/tr['+str(fnd)+']/td[3]/div/span[1]/a')
    sind.send_keys(Keys.ENTER)
    
    tesxt=driver.find_element_by_tag_name('body')
    file=tesxt.text.split('\n')
    df=get_data(file)
    
    return df





