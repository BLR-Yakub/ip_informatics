# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 15:48:19 2022

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
import pandas as pd
import pickle

def dwnload_sequence(patent_number):
    
    url="https://patentscope.wipo.int/search/en/search.jsf"
    driver.get(url)

    time.sleep(5)

    search_box=driver.find_element_by_id('simpleSearchForm:fpSearch:input')
    search_box.send_keys(patent_number)

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
    fnd=-1
    for i in range(len(table)):
        s=table[i].text
        if(s.find('Sequence Listing')!=-1 and s.find('TXT')!=-1):
            fnd=i+1
            break

    time.sleep(5)
    
    #download the text file
    sind=driver.find_elements_by_xpath('//*[@id="'+iD+'"]/div[2]/table/tbody/'+'/tr['+str(fnd)+']/td[4]/div/span[1]/a')
    lnk=sind[0].get_attribute("href")
    if(lnk[-3:]!='app'):
        driver.get(lnk)
    else:
        raise ValueError("Can't find the file, stopping")
    
    
#p=driver.find_element_by_xpath('//*[@id="detailMainForm:MyTabViewId:j_idt9003_data"]/tr[1]/td[4]/div/span[1]/a')        
#dwnload_sequence('WO2009046288')

def exp_to_short(amino_acid):
    ''' 
      A function which can convert expanded forms
      like Ala to A
    '''

    dicts={'Ala':'a', 'Asx':'b', 'Cys':'c',
           'Asp':'d', 'Glu':'e', 'Phe':'f',
           'Gly':'g', 'His':'h', 'Ile':'i',
           'Lys':'k', 'Leu':'l', 'Met':'m',
           'Asn':'n', 'Pro':'p', 'Gln':'q',
           'Arg':'r', 'Ser':'s', 'Thr':'t',
           'Sec':'u', 'Val':'v', 'Trp':'w',
           'Xaa':'x', 'Tyr':'y', 'Glx':'z'}
    
    return dicts[amino_acid]



def translate(sequence):
    coding_dna = Seq(sequence)
    prt=str(coding_dna.translate()) #translate dna to protein
    return prt



def get_data(file, upper=True):
    
    '''
      Take the text file as input 
      and return the dataframe consisting
      of data extracted from the text file
      
      parameters
      ------------
    
      file: text file 
      
      returns
      ------------
      
      Dataframe
    '''
    
    file_content=file.readlines()

        
    seq_no=[]
    seq_len=[]
    seq_type=[]
    organism=[]
    other_info=[]
    trnsl_seq=[]

    chk="<400>" #to check the previous label
    for i in range(len(file_content)):
        x=file_content[i].split()
        
        #get sequence number
        if(len(x)>0 and x[0]=='<210>' and chk=='<400>'):
            seq_no.append(file_content[i][6:len(file_content[i])-1])
            chk=x[0]
            
        #get sequence length
        if(len(x)>0 and x[0]=='<211>' and chk=='<210>'):
            seq_len.append(file_content[i][6:len(file_content[i])-1])
            chk=x[0]
        
        #get sequence type
        if(len(x)>0 and x[0]=='<212>' and chk=='<211>'):
            seq_type.append(file_content[i][6:len(file_content[i])-1])
            chk=x[0]
            
            
        #get organism name
        if(len(x)>0 and x[0]=='<213>' and chk=='<212>'):
            organism.append(file_content[i][6:len(file_content[i])-1])
            chk=x[0]
            
            
        #get other information
        if(len(x)>0 and x[0]=='<223>' and chk=='<213>'):
            other_info.append(file_content[i][6:len(file_content[i])-1])
            chk=x[0]
            
        # if any attribute is missing add none to the list
        if(len(x)>0 and x[0]=='<400>'):
            if(len(other_info)!=len(seq_no)):
                #print(len(seq_no), len(other_info))
                
                other_info.append(" ")
            
            if(len(seq_len)!=len(seq_no)):
                seq_len.append(0)
                
            if(len(seq_type)!=len(seq_no)):
                seq_type.append("None")
                
            if(len(organism)!=len(seq_no)):
                organism.append("None")
                
            chk=x[0]
            

    seq_code=[] #getting the sequence code    
    i=0    
    while(i<len(file_content)):
        x=file_content[i].split()
        if(len(x)>0 and x[0]=='<400>'):
            i+=1
            
            temp=""
            
            while(x[0]!='<210>' and i<len(file_content)): #append all sequence between 400 and 210
                x=file_content[i].split()
                if(len(x)>0):
                    for j in range(len(x)):
                        #convert expanded forms like 'Ala' to 'a'
                        try:    
                            if(upper==True):
                                temp+=exp_to_short(x[j]).upper()
                            else:
                                temp+=exp_to_short(x[j])
                        except:
                            # don't add the digits in sequence 
                            if(x[j].isdigit()==False and x[0]!='<210>'):
                               temp+=x[j]
                else:
                    x.append('0') #to nullify the exception
                i+=1
                
            sz=len(seq_code)
            if(seq_type[sz][-3:]!="PRT"): #if sequence type is not protein
                prt=translate(temp) 
                trnsl_seq.append(prt)
            else:
                trnsl_seq.append(temp)
                
            seq_code.append(temp) 
            
        else:
            i+=1
               
             
    df_dict={'Sequence_ID_number': seq_no, 'Sequence_length':seq_len,
             'Sequence_type':seq_type, 'Organism':organism,
             'Other_info': other_info, 'Sequence':seq_code,
             'Protein': trnsl_seq}

    
    df=pd.DataFrame(df_dict)
    return df

def dwnld_seq_uspto(patent_no):
    
    
    
    url="https://seqdata.uspto.gov/"
    driver.get(url)
    
    search_box=driver.find_element_by_id('docNum')
    search_box.send_keys(patent_no)
    search_box.send_keys(Keys.ENTER)
    
    time.sleep(3)
    db=driver.find_element_by_xpath('/html/body/table[4]/tbody/tr[2]/td/center/form[2]/input[4]')
    db.send_keys(Keys.ENTER)
    
    
def get_sequences_from_zip(path):
        
    x=os.listdir(path)
    paths = [os.path.join(path, basename) for basename in x]
    latest_file1 = max(paths, key=os.path.getctime)
    
    with zipfile.ZipFile(latest_file1, 'r') as zip_ref:
        zip_ref.extractall(path)
        
    x=os.listdir(path)
    paths = [os.path.join(path, basename) for basename in x]
    latest_file2 = max(paths, key=os.path.getctime)
    
    for root, dirs, files in os.walk(latest_file2):
        for name in files:
            print(name)
            print(root)
            this_file=root+'\\'+name
            

    file1 = open(this_file,"r") 
    df=get_data(file1)
    file1.close()
    
    os.remove(latest_file1)
    shutil.rmtree(latest_file2)
    
    return df

def latest_download_file():
      path = "C:\\Users\\EHOUT\\Downloads"
      os.chdir(path)
      files = sorted(os.listdir(os.getcwd()), key=os.path.getmtime)
      newest = files[-1]

      return newest
   

data_file=pd.read_excel("patent_data.xlsx")
data_file['Publication_No'][0][0]



chrome_options = webdriver.ChromeOptions()
chrome_options.add_argument("--incognito")

driver = webdriver.Chrome("C://project_IP//chromedriver.exe",chrome_options=chrome_options)

df=pd.DataFrame()
dint_wrk=[]
wrk=[]
dn=0
#loaded_list=['WO2017125040','WO2020131413', 'WO2018205732']
#fin=0
for i in range(data_file.shape[0]):
#for i in range(len(loaded_list)):
    '''if(data_file['Publication_No'][i][0]=='U'):
        
        try:
            
            dwnld_seq_uspto(data_file['Publication_No'][i])
        
        
            dt=0
            fileends = "crdownload"
            while "crdownload" == fileends:
                time.sleep(1)
                dt+=1
                newest_file = latest_download_file()
                if "crdownload" in newest_file:
                    fileends = "crdownload"
                else:
                    fileends = "none"
            #time.sleep(10)
            path="C:\\Users\\EHOUT\\Downloads"
            df_data=get_sequences_from_zip(path)
            df=pd.concat(df,df_data)
        except:
            dint_wrk.append(data_file['Publication_No'][i])'''
     
    if(data_file['Publication_No'][i][0]=='W'):
    
        try:
            '''if(loaded_list[i][2:].lower().islower()):
                dwnload_sequence(loaded_list[i][:-2])
            else:
                dwnload_sequence(loaded_list[i])'''
            
            if(data_file['Publication_No'][i][2:].lower().islower()):
                dwnload_sequence(data_file['Publication_No'][i][:-2])
            else:
                dwnload_sequence(data_file['Publication_No'][i])
            
            
            
            dt=0
            fileends = "crdownload"
            while "crdownload" == fileends:
                time.sleep(1)
                dt+=1
                newest_file = latest_download_file()
                if "crdownload" in newest_file:
                    fileends = "crdownload"
                else:
                    fileends = "none"
    
    
            path="C:\\Users\\EHOUT\\Downloads"
            x=os.listdir(path)
            paths = [os.path.join(path, basename) for basename in x]
            latest_file = max(paths, key=os.path.getctime)
    
            file1 = open(latest_file,"r")    
            df_data=get_data(file1)
            file1.close()
            
            os.remove(latest_file)
            df=pd.concat((df,df_data))
            wrk.append(data_file['Publication_No'][i])
            dn+=1
            #fin+=1
            '''if(dn==100):
                df.to_csv("fetched_wipo"+str(i)+".csv")
                dn=0
                
                open_file = open("worked"+str(i), "wb")
                pickle.dump(wrk, open_file)
                open_file.close()
                
                
                open_file = open("dintworked"+str(i), "wb")
                pickle.dump(dint_wrk, open_file)
                open_file.close()'''
                
              
            #if(fin>4):
                #break
        except:
            dint_wrk.append(data_file['Publication_No'][i])
            
df.to_csv("fetched_wipo_final3.csv")
            
'''path="C:\\Users\\EHOUT\\Downloads"

for i in range(len(paths)):
    if(paths[i][-3:]=='txt'):
        try:
            file = open(paths[i],"r")  
            df_data=get_data(file)
            df=pd.concat((df,df_data))
            wrk.append(data_file['Publication_No'][i])
        except:
            dint_wrk.append(data_file['Publication_No'][i])
            
df.to_csv("fetched_wipo1.csv")'''

'''t2=open('tesxt.txt', "rb")
p=t2.readlines()

import PyPDF2
 
# creating a pdf file object
tes_file=open('test.pdf',"rb")
 
# creating a pdf reader object
pdfReader = PyPDF2.PdfFileReader(tes_file)
ang=[]
for i in range(pdfReader.numPages):
    pageObj = pdfReader.getPage(43)
 
# extracting text from page
    ang.append(pageObj.extractText())
    
tes_file.close()
get_data(tes_file)

sample_list=[1,2,3,4]
open_file = open("worked"+str(i), "wb")
pickle.dump(sample_list, open_file)
open_file.close()

open_file = open("dintworked2", "rb")
loaded_list = pickle.load(open_file)
open_file.close()'''

