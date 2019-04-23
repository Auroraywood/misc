# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 08:47:03 2019

@author: zdx
@usage: python3 compute_auc.py 
"""
import os
import re
import sys
import pandas as pd
from pandas import Series,DataFrame
global str
# task_folder = "/home/zdx/Downloads/topoII_0422_all_dude"
# score_file = "Merck_Kinase_0422_all_dude.csv"
# os.chdir(task_folder)

def extract_number(str):
    return(re.findall('\d+', str)[0])

#def compute_auc(task_folder, score_file):
    os.chdir(task_folder)
    df = pd.read_csv(score_file, header = None, names = ['ID', 'score'])
    df['ID'] = df['ID'].apply(extract_number) # normalize ID name
    ## process label data
    df_label = pd.read_excel('topoII.true_label.xlsx', sheetname='MERCK_KINASE') # read label file
    df_label = df_label[['Cmpd_ID', 'Active']] # remain two related columns
    df_label.columns = ['ID', 'label'] # modify column name
    df_label = df_label.reset_index(drop=True) # reset row index
    df_label['ID'] = df_label['ID'].apply(str)
    ## add label
    df = pd.merge(df, df_label, on = 'ID', how = 'left')
    
    active_num = len(df[df['label']==1])
    decoy_num = len(df[df['label']==0])
    df.sort_values(by = ['score'], inplace = True, ascending = False)
    df = df.reset_index(drop=True)
    
    act_loc = df[df['label']==1].index.values.astype(int) # location of active
    
    sum = 0
    for i in act_loc:
        df_t = df[:i]
        sum += len(df_t[df_t['label']==0].index.values.astype(int))
    auc = 1 - sum/(active_num * decoy_num)
    print(auc)

if __name__=="__main__":
    argv = sys.argv[1:]
    task_folder = argv[0]
    score_file = argv[1]
    compute_auc(task_folder, score_file)