library(stringr)
library(magrittr)
library(dplyr)
library(tidyverse)

dat <- readLines("file:///C:/Users/F/Desktop/ATP subsequence.txt") %>%
  paste(collapse = " ") %>% # implode a vector of strings into a string
  str_split(">", simplify = T)#split by ">"
dat1 <- grep('PROTEIN_KINASE_ATP', dat, value = T) %>% # retain vectors have "PROTEIN_KINASE_ATP"
  trimws(which = "both") # remove leading and/or trailing whitespace from character strings
dat2 <- data.frame(v1=dat1) # separate one column into multiple columns
dat3 <- separate(dat2, col = "v1", into = c("pro_id", "ATP_bind_pocket"),  sep = " ")
bind_sq <- dat3

extract_pro_id <- function(x){
  t1 <- x[1]
  t2 <- str_extract(t1, "\\|......\\|")
  t3 <- sub("^\\|", "", t2)
  t4 <- sub("\\|$", "", t3)
  return(t4)
}

bind_sq$pro_id <- apply(dat3, 1, extract_pro_id) # extract protein id
mkl_pro_id <- read.table("file:///E:/desktop/work/MKL/pairwiseMKL-master/DTI_data/Uniprot_IDs_226kin.txt", stringsAsFactors = F, header = F)
colnames(mkl_pro_id) <- "pro_id"
res <- plyr::join(mkl_pro_id, bind_sq, by = "pro_id", type = "left")
res$ATP_bind_pocket <- toupper(res$ATP_bind_pocket)

keep_letter <- function(x){ # function of retain capital letter 
  t1 <- x[2]
  t2 <- unlist(str_extract_all(t1, "[A-Z]"))
  t3 <- paste(t2, collapse = "")
  return(t3)
}
keep_letter(res[3,])

res$ATP_bind_pocket <- apply(res, 1, keep_letter)
test <- res[1:5, 2, drop=F]
write.csv(test, file = "test1.csv")

import pandas as pd 
import numpy as np 

from collections import OrderedDict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs

import gc
import random
import os

def smiles_to_ecfp4(smiles):
    '''Use SMILES to calculate ECFP4.'''
    try:
        mol = Chem.MolFromSmiles(smiles)
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    except:
        return None

def tanimoto(fp1, fp2):
    '''Compute Tanimoto similarity between two fingerprints'''
    return DataStructs.FingerprintSimilarity(fp1,fp2)

def load_data(csv_file, sort_keyword):
    '''Load data and sort it by a keyword'''
    df = pd.read_csv(csv_file)
    df = df.sort_values(by=sort_keyword).reset_index(drop=True)
    return df

def get_neighbor(df_main, df_side, keyword, output_file, num_neighbors=3000, buffer_size=5000):
    '''Get neighbors in side dataframe for each compound in main dataframe.
       Take the ideas of merge sort and buffer, which can reduce time complexity
       and avoid repeated calculation of ECFP4.'''
    
    fp_buffer = OrderedDict()                            #A buffer, where calculated ECFP4s will be saved.
    idx_main = 0                                         #Current index in main dataframe
    idx_side = 0                                         #Current index in side dataframe
    MAX_IDX_MAIN = df_main.shape[0]                      #Length of main dataframe
    MAX_IDX_SIDE = df_side.shape[0]                      #Length of side dataframe
    fout =  open(output_file, 'w')
    while(idx_main < MAX_IDX_MAIN): 
        #Travel on main dataframe.
        
        cid_main = df_main.ix[idx_main, 'CID']           #Get current main compound ID 
        val_main = df_main.ix[idx_main, keyword]         #Get specific property value of current main compound
        smiles_main = df_main.ix[idx_main, 'SMILES']     #Get SMILES of current main compound
        fp_main = smiles_to_ecfp4(smiles_main)           #Get ECFP4 of current main compound
        
        if fp_main is None:
            # If to get ECFP4 of main compound is failed, run for next main compound.
            idx_main += 1
            continue
        
        candidates = []                                  #Candidates, from which the decoys will be randomly chosen
        
        if idx_main % 10 == 0:
            #Print log.
            print('Now finding %s decoy for  %s.  %d / %d...' % (keyword, cid_main, idx_main, MAX_IDX_MAIN))            
        
        while(idx_side < MAX_IDX_SIDE):
            #Travel on side dataframe.
            
            if idx_side % 100000 == 0:
                 #Print log.
                print('Side idx is %d' % idx_side)
                
            val_side = df_side.ix[idx_side, keyword]     #Get specific property value of current main compound
            
            if val_side < val_main:
                #If current side value is less than main value, run for next side compound, just like merge sort.
                idx_side += 1
                continue
                
            else:
                #Get index of the neighbor whose property value is the lowest and the one highest.
                #Make sure that it keep balanced and avoid boundary error. 
                lo = max(0, idx_side - num_neighbors // 2)                                 
                hi = min(MAX_IDX_SIDE - 1, idx_side + num_neighbors - num_neighbors // 2)
                if hi - lo < (num_neighbors - 1):
                    lo = max(hi - num_neighbors + 1, 0)
                    
                for i, row in df_side.ix[lo: hi].iterrows():
                    #For each one in current interval of side dataframe,
                    #that is to say, the neighbors of current main compound, whose property values are matched.
                    
                    cid_side = row['CID']                           #Get current side compound ID 
                    smiles_side = row['SMILES']                     #Get SMILES of current side compound
                    
                    try:
                        #Try to get ECFP4 of side compound in buffer
                        fp_side = fp_buffer[cid_side]
                    except KeyError:
                        #If not hit, calculate it.
                        fp_side = smiles_to_ecfp4(smiles_side)
                        fp_buffer[cid_side] = fp_side
                        if len(fp_buffer) > buffer_size:
                            #If buffer is full, pop the item that first came in, just like FIFO queue.
                            fp_buffer.popitem(last=False)
                    if fp_side is not None:
                        #Compute similarity between main fingerprints and side fingerprints
                        similarity = tanimoto(fp_main, fp_side)
                        candidates.append((cid_side, smiles_side, similarity))                        
                        
                candidates = sorted(candidates, key=lambda x: x[2])[: num_neighbors // 4 + 1]    #Sort candidates by similarity, and get the top 25% dissimilar part
                random.shuffle(candidates)      #Shuffle it
                candidates = candidates[:10]    #Choose 10 items as decoys
                
                #Output decoys and corresponding active compound. 
                for record in candidates:
                    print(record[0], record[1], sep=',', file=fout)
                print(cid_main, file=fout)
                print('sep\n', file=fout)
                
                #Run for next main compound.
                idx_main += 1 
                break
        
        #If travel on side dataframe ends, break.
        if idx_side >= MAX_IDX_SIDE:
            break

    while(idx_main < MAX_IDX_MAIN):
        #If travel on main dataframe does not end,
        #for compound in the rest of main dataframe,
        #get neighbors in the bottom of side dataframe, workflow like above block.        
        cid_main = df_main.ix[idx_main, 'CID']
        val_main = df_main.ix[idx_main, keyword]
        smiles_main = df_main.ix[idx_main, 'SMILES']
        fp_main = smiles_to_ecfp4(smiles_main)
        if fp_main is None:
            idx_main += 1
            continue
        candidates = []
        if idx_main % 10 == 0:
            print('Now finding %s decoy for  %s.  %d / %d...' % (keyword, cid_main, idx_main, MAX_IDX_MAIN))
        lo = MAX_IDX_SIDE - 1 - num_neighbors + 1
        hi = MAX_IDX_SIDE - 1
        for i, row in df_side.ix[lo: hi].iterrows():
            cid_side = row['CID']
            smiles_side = row['SMILES']
            try:
                fp_side = fp_buffer[cid_side]
            except KeyError:
                fp_side = smiles_to_ecfp4(smiles_side)
                fp_buffer[cid_side] = fp_side
                if len(fp_buffer) > buffer_size:
                    fp_buffer.popitem(last=False)
            if fp_side is not None:
                similarity = tanimoto(fp_main, fp_side)
                candidates.append((cid_side, smiles_side, similarity))
        candidates = sorted(candidates, key=lambda x: x[2])[: num_neighbors // 4 + 1]
        random.shuffle(candidates)
        candidates = candidates[:10]
        for record in candidates:
            print(record[0], record[1], sep=',', file=fout)
        print(cid_main, file=fout)
        print('sep\n', file=fout)
        idx_main += 1
    fout.close()
    print('Done.')
   
def main(csv_file_main_list, csv_file_side, output_path, keyword_list=['MWT', 'LGP', 'RTB', 'PSA', 'HBA', 'HBD'])
    '''
    Parametes:
        csv_file_main_list: a list of main data csv filenames, that is, active compound information csv filenames.
        csv_file_side: side data csv filename, from which the decoys come.
        output_path: directory where output file will be saved.
        keyword_list: a list of keywords. A keyword means a property.
        
    -----    
    Note:
        Both csv_file_main and csv_file_side must have a column of 'CID' and a column of the corresponding keyword.
        An example is as following.
        
        CID,MWT,LGP,RTB,PSA,HBA,HBD
        ID001,100,1.1,5,50.0,1,1
        ID002,200,2.2,3,120.0,2,1
        
        Default keywords are 'MWT', 'LGP', 'RTB', 'PSA', 'HBA' and 'HBD'.
        They means molecular weight, logP, rotatable bonds, polar surface area, H-bond acceptor and H-bond donor respectively.
    '''


    for keyword in keyword_list:
        print('Loading data...')
        df_side = load_data(csv_file_side, keyword) 
        gc.collect()
        output_file = os.path.join('decoys_of_%s.csv' % keyword)
        df_main = pd.read_csv(os.path.join(output_path, csv_file_main_list[0]))[:0]
        for csv_file_main in csv_file_main_list:            
            df_current = pd.read_csv(csv_file_main)
            df_main = df_main.append(df_current, ignore_index=True)
        df_main = df_main.sort_values(by=keyword).reset_index(drop=True)            
        get_neighbor(df_main, df_side, keyword, output_file)
