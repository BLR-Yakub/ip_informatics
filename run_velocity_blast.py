# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 00:27:48 2022

@author: EHOUT
"""

import requests
import os
import json
import time
import pandas as pd
import datetime


def run_vel_blast(sequences, access_token, algo="blastn", collection="nr", evalue=10, nhits=500, word_size=11):
    url="https://velocity.ag/ncbi-blast/services/v1/blast/"
    headers = {
        'accept': 'application/json',
        'Content-Type': 'application/json',
        'authorization': 'Bearer ' + access_token
    }
    user=os.getlogin()
    job=user+"_"+str(datetime.datetime.now()).replace(" ","_")
    param_str="-evalue "+str(float(evalue))+" -num_alignments "+str(nhits)+" -word_size "+str(word_size)
    body={ "userId": user,
           "jobName": job,
           "jobDescription": job+"_"+collection,
           "blastCompletionQueue": "pd-genomics-ncbi-blast-completion-out-srgrod",
           "blastDetails": [ { "algorithm": algo,
                               "parameters": param_str,
                               "collectionName": collection, 
                               "blastFormatters": [ 0,10 ], 
                               "sequences": sequences } ]}
    
    response = requests.post('https://velocity.ag/ncbi-blast/services/v1/blast/', headers=headers, data=json.dumps(body))
    if response.status_code != 200:
        print("Couldn't submit Velocity BLAST job. Error code: "+str(response.status_code))
        return None
    
    Obj = json.loads(response.content)
    batch_id = Obj['blastExecutionBatchId']
    return batch_id


sequence=[{"name":"test",
           "text":"MQKLINSVQNYAWGSKTALTELYGMENPSSQPMAELWMGAHPKSSSRVQNAAGDIVSLRDVIESDKSTLLGEAVAKRFGELPFLFKVLCAAQPLSIQVHPNKHNSEIGFAKENAAGIPMDAAERNYKDPNHKPELVFALTPFLAMNAFREFSEIVSLLQPVAGAHPAIAHFLQQPDAERLSELFASLLNMQGEEKSRALAILKSALDSQQGEPWQTIRLISEFYPEDSGLFSPLLLNVVKLNPGEAMFLFAETPHAYLQGVALEVMANSDNVLRAGLTPKYIDIPELVANVKFEAKPANQLLTQPVKQGAELDFPIPVDDFAFSLHDLSDKETTISQQSAAILFCVEGDATLWKGSQQLQLKPGESAFIAANESPVTVKGHGRLARVYNKL*"}
          ]
bid=run_vel_blast(sequence,token,algo="blastp",collection="IC_TOXIN-All",
                  evalue=10, nhits=250, word_size=3)






def get_exec_results(batch_id, access_token, poll_time):
    url="https://velocity.ag/ncbi-blast/services/v1/blast/blast-execution-batch/"
    headers = {'accept': 'application/json', 'authorization': 'Bearer ' + access_token}
    #url_details=requests.get(url+str(batch_id),headers=headers)
    err_count = 0
    MAX_ERRS = 3
    while True:
        time.sleep(poll_time)
        response = requests.get(url+str(batch_id),headers=headers)
        # print(response.text)
        if response.ok:
            response_json = response.json()
            response_set = set(map(lambda x: x["status"], response_json['blastExecutionDetails']))
            print(response_set)
            if response_set <= {'Archive', 'Error'}:
                return response_json
        else:
            print(response.text)
            err_count += 1
            if err_count > MAX_ERRS:
                raise ValueError("Error count greater than max errs, stopping")

        
def get_result(batch_id, access_token, poll_time):
    url="https://velocity.ag/ncbi-blast/services/v1/blast/result?resultUrl="
    headers = {'accept': 'application/json', 'authorization': 'Bearer ' + access_token}
    exec_res=get_exec_results(batch_id, access_token, poll_time)
    
    results=[]
    warnings=[]
    for result_detail in exec_res['blastExecutionDetails']:
        if result_detail["status"] == 'Archive':
            result_url = result_detail['blastFormatResults'][1]['resultUrl']
            
            response = requests.get(url+result_url,headers=headers)
            return response
            #blast_out = list(filter(lambda x: len(x.rstrip()) > 0, response.text.split("\n")))
            #results.append(blast_out)
        else:
            #warnings.warn("Sequence %s has non-archive status" % result_detail["sequenceName"])
            warnings.append("Sequence %s has non-archive status:  %s"
                                        % (result_detail["sequenceName"], result_detail["errorText"]))
            
    return warnings

# Get blast results for multiple sequences
def get_blast_res(batch_id, access_token):
    try:
        
        headers = {'accept': 'application/json', 'authorization': 'Bearer ' + access_token}
        url='https://velocity.ag/ncbi-blast/services/v1/blast/blast-execution-batch/'+str(batch_id)+'/result/csv'
        res=requests.get(url, headers=headers)
        with open('blas_res6.csv', 'w') as file:
            file.write(res.text)
    
        res_file=pd.read_csv("blas_res6.csv", header=[0])  
        os.remove('blas_res6.csv')
        return res_file
    except:
        print("Couldn't find blast results")
        
df=get_blast_res(bid,token)

res=get_result(bid,token,30)

with open('blas_res2.csv', 'w') as file:
    file.write(res.text)

res_file=pd.read_csv("blas_res2.csv", header=None)  
res_file.to_csv('blas_res2.csv', header=['qacc', 'sacc', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], index=False)
    
    


















