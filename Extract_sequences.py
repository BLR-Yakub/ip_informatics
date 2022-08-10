# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 11:43:43 2021

@author: EHOUT
"""

from selenium import webdriver
import time
from selenium.webdriver.common.keys import Keys
from datetime import datetime
from selenium.webdriver.support.ui import Select
import re
from selenium.webdriver.chrome.options import Options
import os
import zipfile
import shutil
from Bio.Seq import Seq
import pandas as pd
import pickle

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

# get sequences from sequence listing text file
def get_data(file, upper=True):
    file_content=file
    
        
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
            seq_no.append(file_content[i][6:len(file_content[i])])
            chk=x[0]
            
        #get sequence length
        if(len(x)>0 and x[0]=='<211>' and chk=='<210>'):
            seq_len.append(file_content[i][6:len(file_content[i])])
            chk=x[0]
        
        #get sequence type
        if(len(x)>0 and x[0]=='<212>' and chk=='<211>'):
            seq_type.append(file_content[i][6:len(file_content[i])])
            chk=x[0]
            
            
        #get organism name
        if(len(x)>0 and x[0]=='<213>' and chk=='<212>'):
            organism.append(file_content[i][6:len(file_content[i])])
            chk=x[0]
            
            
        #get other information
        if(len(x)>0 and x[0]=='<223>' and chk=='<213>'):
            other_info.append(file_content[i][6:len(file_content[i])])
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
