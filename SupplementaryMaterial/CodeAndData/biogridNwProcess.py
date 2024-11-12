#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 01:14:38 2022

@author: Ozan Ozisik
"""

import networkx as nx
import pandas as pd

biogridFile='Data/Networks/BIOGRID-ORGANISM-Homo_sapiens-4.4.222.tab3.txt'
outputFile='Data/Networks/Biogrid-HSapiens-physical-4.4.222.txt'

dfBiogrid=pd.read_csv(biogridFile, sep='\t', dtype='str', na_values='-')

## Only physical interactions are kept
dfBiogrid=dfBiogrid[dfBiogrid['Experimental System Type']=='physical']

## Only human interactions will be kept
dfBiogrid=dfBiogrid[(dfBiogrid['Organism ID Interactor A']=='9606') & (dfBiogrid['Organism ID Interactor B']=='9606')]

dfBiogrid=dfBiogrid[['Official Symbol Interactor A','Official Symbol Interactor B']]

dfBiogrid=dfBiogrid.dropna()


## Using networkx to remove redundant edges (creating networkx object handles this) and self loops
G=nx.from_pandas_edgelist(dfBiogrid, source='Official Symbol Interactor A', target='Official Symbol Interactor B', create_using=nx.Graph)
G.remove_edges_from(nx.selfloop_edges(G))
nx.write_edgelist(G, outputFile, delimiter="\t", data=False)
