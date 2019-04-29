# -*- coding: utf-8 -*-
'''
Description: train gnina
Script: train_gnina.py
@Author: zdx
1st version: Apr 29, 2019
'''




import os
import sys
import glob





def main(task_folder,
         cmpd_lib,
         ):
    if not os.path.exists(task_folder):
       print("No such task {0} found, quit.".format(task_folder))
       return
    rec_pdbqt = "{0}/target/receptor.pdbqt".format(task_folder)
    if not os.path.exists(rec_pdbqt):
       print("No receptor found, quit.")
       return
    mypath = "{0}/docking_pdbqt/{1}/*.pdbqt".format(task_folder, cmpd_lib)
    ligand_files = glob.glob(mypath)














if __name__=="__main__":
    argv = sys.argv[1:]
    task_folder = argv[0]
    cmpd_lib = argv[1]
    main(task_folder, cmpd_lib)