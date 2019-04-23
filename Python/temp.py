#!/usr/bin/env python3

#######################################################
#                                                     
#Project: cnnscore                                  
#Script: cnnscore.py                                   
#
#Description: To predict binding probability from given 
#             receptor and given ligand pose(s)
#
#Reference: None
#                                                      
#1st version: Apr 19, 2019                           
#Author: Zhou dixin @ DeepDrug                                    
#                                                        
#2nd version: To_be_entered_here                      
#Author: To_be_entered_here                           
#
####################################################### 
import os
import sys
from os import listdir
from os.path import isfile, join
from subprocess import call
import glob
from pandas import Series,DataFrame
import pandas as pd

# task_folder = "/home/zdx/Documents/topoII"
# cmpd_lib = "active"

def Replace(x):
    return(x.replace(".pdbqt", ""))

def main(task_folder,
         cmpd_lib,
         model = "/home/zdx/gnina/models/refmodel3/refmodel3.model",
         weights = "/home/zdx/gnina/models/refmodel3/weights/depth3_dude_allfolds.caffemodel",
         work_dir = "/home/zdx/gnina/models/refmodel3/"
         ):
    os.chdir(work_dir)
    if not os.path.exists(task_folder):
       print("No such task {0} found, quit.".format(task_folder))
       return
    rec_pdbqt = "{0}/target/receptor.pdbqt".format(task_folder)
    if not os.path.exists(rec_pdbqt):
       print("No such receptor found, quit.")
       return
    mypath = "{0}/docking_pdbqt/{1}/tmp/*.pdbqt".format(task_folder, cmpd_lib)
    ligand_files = glob.glob(mypath)
    path = "{0}/docking_pdbqt/{1}/tmp/".format(task_folder, cmpd_lib)
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    ligand_name = list(map(Replace, onlyfiles)) # get ligand name
    cnnscore = [] # list of cnnscore
    for ligand in ligand_files:
        call(['gnina','--score_only','-r',rec_pdbqt,'-l',ligand,'--cnn_model',model,'--cnn_weights',weights,'--cnn_scoring','-o','result','--gpu'])
        with open('result', 'r') as f:
            c = f.read() # file content
            s = c.find('CNNscore') + 10# start location to extract
            e = s + 6# end location to extract
            score = c[s:e] # CNNscore
        cnnscore.append(score)
    # cnnscore = list(range(20))
    data = {'ID':ligand_name,'cnnscore':cnnscore}  
    df = DataFrame(data)
    df =  df.sort_values(by='cnnscore', ascending=False)
    f_t = "/home/zdx/{0}_score".format(cmpd_lib)
    df.to_csv(f_t, sep='\t', index=False, encoding='utf-8')
    
#    with open('active_score', 'r') as f:
        
if __name__=="__main__":
    argv = sys.argv[1:]
    task_folder = argv[0]
    cmpd_lib = argv[1]
    main(task_folder, cmpd_lib)
   
    
    
    
    
    
    
    
    
    
    
    
    
