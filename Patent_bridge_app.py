from selenium import webdriver
import time
import os
from selenium.webdriver.common.keys import Keys
from datetime import datetime
from selenium.webdriver.support.ui import Select
from Bio.Seq import Seq

from selenium.webdriver.chrome.options import Options
import requests
import os
import json
import time
import pandas as pd
import datetime
from flask import Flask, request, render_template, Response
from flask import session



def get_token(client_id, client_secret):
    '''
    Contacts the azure_url and returns an access_token
    '''
    # get credentials from a ENV
    
    access_url = 'https://login.microsoftonline.com/fcb2b37b-5da0-466b-9b83-0014b67a7c78/oauth2/v2.0/token'
    # apiURL = creds['KD_SCORE_DATA_URL']
    # Build the call to Oauth client from creds file
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    
    payload = {
                'grant_type' : 'client_credentials',
                'client_id' : client_id,
                'client_secret' : client_secret,
                'scope' : client_id + "/.default"
                }
    response = requests.post(access_url, headers=headers,data=payload)
    # check for a successful response
    if response.status_code == 200:
        accessObj = json.loads(response.content)
        accessToken = accessObj['access_token']
        return (accessToken)
    
    return (None)

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

def get_blast_res(batch_id, access_token):
    try:
        
        headers = {'accept': 'application/json', 'authorization': 'Bearer ' + access_token}
        url='https://velocity.ag/ncbi-blast/services/v1/blast/blast-execution-batch/'+str(batch_id)+'/result/csv'
        res=requests.get(url, headers=headers)
        
        return res
    except:
        print("Couldn't find blast results")
        
fl_wipo=pd.read_csv('C:\\project_IP\\IP_quality_Informatics\\fetched_wipo_all_final.csv')
fl_uspto=pd.read_csv('C:\\project_IP\\IP_quality_Informatics\\fetched_uspto_all_final.csv')

client_id = '6212df18-f0d8-49e6-a8fc-ea98aae348ad'
client_secret = '6.w_n6rrE.-zrEkPAAd2JjM0vrjktj.4E2'

token=get_token(client_id, client_secret)

app = Flask(__name__)
app.secret_key = "super secret key"

@app.route("/", methods=['POST','GET'])
def load_page():
    
    
        
    patent_no=''
    if request.method == 'POST':
        patent_no=request.form.get('pat_no')
        print(patent_no)
        
        
        session['num']=patent_no
    
    
        #df_data=pd.read_csv('C:\\Users\\EHOUT\\Downloads\\read8.csv') #comment it out
        seq_list=[]
        if(patent_no[0]=='W'):
            
            seq_list=fl_wipo[fl_wipo['patent_number']==patent_no]['Sequence'].to_list()
        else:
            seq_list=fl_uspto[fl_uspto['patent_number']==patent_no]['Sequence'].to_list()
        
        print(seq_list[0])
    
        seq_dict = [dict(zip(['seq'],[seq_list[i]])) for i in range(len(seq_list))]
        
        return render_template("index.html", data=lits2, data2=seq_dict,pat=patent_no)
    return render_template("index.html", data=lits2)
    

@app.route("/download",methods=['POST','GET'])
def downloadFile():
    #For windows you need to use drive name [ex: F:/Example.pdf]
    patent_no=session.get('num')
    df=pd.DataFrame()
    
    if(patent_no[0]=='W'):
        df=fl_wipo[fl_wipo['patent_number']==patent_no]
        
        
    else:
        df=fl_uspto[fl_uspto['patent_number']==patent_no]
    
    
    return Response(
       df.to_csv(),
       mimetype="text/csv",
       headers={"Content-disposition":
       "attachment; filename=filename.csv"})
    

@app.route('/convert', methods=['POST', 'GET'])
def get_the_data():
    
    
    #Now going to blast
 
    seq_obt=request.form.get('cod')
    print(seq_obt)
    sequence=[{"name":"test",
               "text":seq_obt}
              ]
    bid=run_vel_blast(sequence,token,algo="blastp",collection="IC_TOXIN-All",
                      evalue=10, nhits=250, word_size=3)
    
    res=get_result(bid,token,30)

    with open('C:\\Users\\EHOUT\\Downloads\\blas_res4.csv', 'w') as file:
        file.write(res.text)

    res_file=pd.read_csv("C:\\Users\\EHOUT\\Downloads\\blas_res4.csv", header=None)  
    res_file.to_csv('C:\\Users\\EHOUT\\Downloads\\blas_res4.csv', header=['qacc', 'sacc', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], index=False)
    
    
    df= pd.read_csv("C:\\Users\\EHOUT\\Downloads\\blas_res4.csv") 
    os.remove('C:\\Users\\EHOUT\\Downloads\\blas_res4.csv')
    
    
    return render_template("index.html", column_names=df.columns.values, row_data=list(df.values.tolist()),
                           zip=zip)
    
    
if __name__=="__main__":
    
    app.run(debug=True, port=8000)
