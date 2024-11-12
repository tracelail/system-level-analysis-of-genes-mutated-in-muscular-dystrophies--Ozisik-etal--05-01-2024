#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 14:00:20 2023

@author: Ozan Ozisik

Converts IDs of genes

"""

import pandas as pd

def getIDMappingDict(mappingFilePath, convertFrom, convertTo):
    df=pd.read_csv(mappingFilePath, sep="\t")
    df=df[[convertFrom, convertTo]]
    df=df.dropna(axis=0)
    if('NCBI Gene ID(supplied by NCBI)' in df.columns):
        df['NCBI Gene ID(supplied by NCBI)']= df['NCBI Gene ID(supplied by NCBI)'].astype(int)
        df['NCBI Gene ID(supplied by NCBI)']= df['NCBI Gene ID(supplied by NCBI)'].astype(str)
    
    df=df.set_index(convertFrom)
    
    idMappingDict=df.to_dict()[convertTo]
    
    return idMappingDict


def mapIDsInGeneGroupsDict(geneGroupsDict, idMappingDict):
    
    idMappedGeneGroupsDict=dict()
    for group, genesList in geneGroupsDict.items():
        idMappedGeneGroupsDict[group]=[ idMappingDict[g] for g in genesList if g in idMappingDict.keys()]
        
    return idMappedGeneGroupsDict