# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 19:13:39 2022

@author: EHOUT
"""

import json
import pandas as pd
import requests


    
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

client_id = '6212df18-f0d8-49e6-a8fc-ea98aae348ad'
client_secret = '6.w_n6rrE.-zrEkPAAd2JjM0vrjktj.4E2'

token=get_token(client_id, client_secret)
            





    
    










