#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 19:34:31 2018

@author: lurban
"""

import os
import pandas as pd
import numpy as np
import h5py as h5
import scipy as sp
import scipy.stats as st
import limix
from limix.modules import qtl 
from limix.modules.qtl import test_lmm, test_lm
from statsmodels.sandbox.stats.multicomp import multipletests
import tarfile
from itertools import chain
os.chdir('/hps/nobackup2/research/stegle/users2/lurban/final_signatures/sig_paper')

def toRanks(A):
    """
    converts the columns of A to ranks
    """
    AA=sp.zeros_like(A)
    if (AA.shape[1]==1):
        AA = st.rankdata(A)
    else:
        for i in range(A.shape[1]):
            AA[:,i] = st.rankdata(A[:,i])
    AA=sp.array(sp.around(AA),dtype="int")-1
    return AA

def gaussianize(Y):
    """
    Gaussianize X: [samples x phenotypes]
    - each phentoype is converted to ranks and transformed back to normal using the inverse CDF
    """
    N,P = Y.shape
    YY=toRanks(Y)
    quantiles=(sp.arange(N)+0.5)/N
    gauss = st.norm.isf(quantiles)
    Y_gauss=sp.zeros((N,P))
    if (Y_gauss.shape[1]==1):
        Y_gauss = gauss[YY]
    else:
        for i in range(P):
            Y_gauss[:,i] = gauss[YY[:,i]]
    Y_gauss *= -1
    return Y_gauss

def rankStandardizeNormal(Y):
    """
    Gaussianize X: [samples x phenotypes]
    - each phentoype is converted to ranks and transformed back to normal using the inverse CDF
    """
    return gaussianize(Y)

# preprocessing

# load signature beta release set, map ids, and add APOBEC sinature

mutsigo = pd.read_excel('PCAWG_sub_signatures_in_samples_beta2.xlsx')
mutsigo.index = mutsigo['Sample Name']
mutsigo = mutsigo.iloc[:,2:41]

map = pd.read_csv('map.tsv',sep='\t',header=0,index_col=0)         
mapunique = pd.read_csv('mapunique.tsv',sep='\t',header=0,index_col=0)                     
map = map[['normal_wgs_aliquot_id','tumor_wgs_aliquot_id','tumor_rna_seq_aliquot_id']]

mutsigo = mutsigo.merge(map, how='left', left_index=True, right_on='tumor_wgs_aliquot_id')
mutsigo.index = mutsigo.normal_wgs_aliquot_id
del mutsigo['tumor_wgs_aliquot_id']; del mutsigo['tumor_rna_seq_aliquot_id']
mutsigo.to_csv('PCAWG_sub_signatures_normal_id.tsv', sep='\t')
del mutsigo['normal_wgs_aliquot_id']

mutsigo.loc['Signature Apo'] = mutsigo.loc['Signature 2'].astype('float') + mutsigo.loc['Signature 13'].astype('float') 
mutsigo.to_csv('PCAWG_sub_signatures_processed.tsv', sep='\t')


  

# limix
    
searchDelta = False         #specify if delta should be optimized for each SNP
test="lrt"                  #specify type of statistical test

xs = pd.read_csv('PCAWG_sub_signatures_processed.tsv', sep='\t', index_col=0, header=0)
xs = xs.transpose()

ys = pd.read_csv('peer35.pan.tsv', sep='\t', index_col=0, header=0)
ys = ys.transpose()
ys = ys.loc[xs.index]

mutsigindex = pd.read_csv('mutsigindexs.tsv', sep='\t', index_col=0, header=0) # 1159x54849

kin = h5.File('kinship.hdf5', 'r')
kin = kin['K'].value
kinh = pd.read_csv('header_kinship.tmp',header=None)
kin = pd.DataFrame(kin, columns=kinh.values, index=kinh.values)
kin.index = list(kinh.iloc[:,0])
kin.columns = list(kinh.iloc[:,0])

# choose covariates
covariates = mutsigindex.columns[39:] 
# drop unnecessary covariates 
covariates = covariates.drop(['project_code_BLCA-US','project_code_LAML-US','donor_sex_male','Bone-Cart']) 

both = set(kin.index).intersection(xs.normal_wgs_aliquot_id.values)
xs.index = xs['normal_wgs_aliquot_id']
del xs['normal_wgs_aliquot_id']
ys.index = xs.index                          
mutsigindex = mutsigindex.loc[both]
kin = kin.loc[both]
kin = kin.transpose()
kin = kin.loc[both]
kin = kin.transpose()
xs = xs.loc[both]
ys = ys.loc[both]

X0 = mutsigindex[covariates]
X0 = X0.iloc[:,np.where(X0.sum() != 0)[0]]
xs = xs.drop(xs.columns[np.where(xs.var(axis=0)==0)[0]], axis=1) 
# filter for prevalency: non-zero in at least 1% of cohort
xs = xs.convert_objects(convert_numeric=True)
# check if all TRUE
xs.apply(lambda x: sum(x==0)>=(0.01*xs.shape[0])) 

snps = xs
phenotype_std = ys
sample_relatedness = kin
covs = X0
    
columnes = phenotype_std.columns.values
indexes = snps.columns.values

# transform to matrix for limix usage

snps = np.matrix(snps)

# standardisation
snps = (snps - np.mean(snps, axis=0)) / np.std(snps, axis=0)
    
phenotype_std = np.matrix(phenotype_std)
     
# if permuted values for qq-plots
#phenotype_std = np.random.permutation(phenotype_std)
    
sample_relatedness = np.matrix(sample_relatedness)
covs = np.matrix(covs)
lmm = test_lmm(snps=snps,pheno=phenotype_std,K=sample_relatedness,covs=covs,test=test)
    
pv_lmm = lmm.getPv()
b_lmm = lmm.getBetaSNP()
pvalues_lmm = pd.DataFrame(data=np.transpose(pv_lmm), index=indexes, columns=columnes)
betas_lmm = pd.DataFrame(data=np.transpose(b_lmm), index=indexes, columns=columnes)

pvalues_lmm.to_csv('pvalues.tsv', sep='\t')
betas_lmm.to_csv('betas.tsv', sep='\t')


