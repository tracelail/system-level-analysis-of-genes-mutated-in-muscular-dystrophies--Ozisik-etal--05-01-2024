#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 13:48:34 2023

@author: Ozan Ozisik

"""

import os
from enrichmentAnalysis import performEnrichmentAnalysis, applyOrsum
from NetworkAnalysis import NetworkAnalysis
from idMapping import getIDMappingDict, mapIDsInGeneGroupsDict



#### Parameters


## Choose the analyses you want to perform
runEnrichment=True
runAverageDistanceBetweenGroups=True
runRWROverlap=True


runDSD=False


rwrTopNumber=50
rwrOverlapBootstrapNumber=5 #set to a small number for people to try the code, 1000 takes time
orsumNumberOfTermsToPlot=20
averageDistanceBootstrapNumber=2000



resultFolder='Results'
enrichmentFolder=resultFolder+ os.sep+ 'EnrichmentResults'
summaryFolder=enrichmentFolder+os.sep+ 'Summary'

networkAnalysisFolder=resultFolder+ os.sep+ 'NetworkResults'

geneGroups=['DoGs', 'CoGs', 'PoGs']

pathsForGeneGroups={'DoGs':'./Data/DiseaseGenes/DoGs.txt',
                    'CoGs':'./Data/DiseaseGenes/CoGs.txt',
                    'PoGs':'./Data/DiseaseGenes/PoGs.txt'}

networkFiles=['Data/Networks/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv',
              'Data/Networks/Pathways_reactome_gene_names_190123.tsv',
              'Data/Networks/Complexes_gene_names_190123.tsv']

#networkFiles=['Data/Networks/Biogrid-HSapiens-physical-4.4.222.txt']


geneIDMappingFile='Data/HGNC 230616.txt'


## The URL for the g:Profiler request.
## We performed the analysis with e109_eg56_p17_1d3191d (database updated on 29/03/2023) 
gProfilerBaseUrl='https://biit.cs.ut.ee/gprofiler_archive3/e109_eg56_p17'

## If up-to-date data will be used,
## gProfilerBaseUrl='https://biit.cs.ut.ee/gprofiler'
## can be used. In that case the gmt files for orsum should also be updated.







#### Initialization

if not os.path.exists(resultFolder):
    os.makedirs(resultFolder)
    print('Result folder \"'+resultFolder+'\" is created.')
    


#### Reading DoGs, CoGs and PoGs
geneGroupsDict=dict()
allgenes=[]
for groupName, path in pathsForGeneGroups.items():
    with open(path) as f:
        geneGroupsDict[groupName]=[line.rstrip() for line in f]
        allgenes=allgenes+geneGroupsDict[groupName]




#### Enrichment analysis
if runEnrichment:
    
    if not os.path.exists(enrichmentFolder):
        os.makedirs(enrichmentFolder)
        print('Enrichment result folder \"'+enrichmentFolder+'\" is created.')
    
    
    ## Mapping gene symbols to ensembl IDs because g:Profiler discards genes 
    ## which have ambiguity - which are mapped to different ensembl IDs
    ## by HGNC and NCBI. Here we use HGNC to do the mapping.
    idMappingDict=getIDMappingDict(geneIDMappingFile,'Approved symbol', 'Ensembl ID(supplied by Ensembl)')
    idMappedGeneGroupsDict=mapIDsInGeneGroupsDict(geneGroupsDict, idMappingDict)
    
    
    sources=['GO:BP', 'GO:CC']
    
    performEnrichmentAnalysis(idMappedGeneGroupsDict, enrichmentFolder, summaryFolder, sources=sources, gProfiler_base_url=gProfilerBaseUrl)
    
    ## The gmt file should be changed according to the data used for enrichment analysis.
    for source in sources:
        applyOrsum(gmt='Data/Annotation/hsapiens.'+source+'.name.gmt', 
                   source=source, 
                   geneGroupsDict=idMappedGeneGroupsDict,
                   summaryFolder=summaryFolder,
                   outputFolder='orsum'+source.replace(':',''),
                   numberOfTermsToPlot=orsumNumberOfTermsToPlot
                   )


#### Network analysis

if runAverageDistanceBetweenGroups or runDSD or runRWROverlap:
    networkAnalysis=NetworkAnalysis(networkFiles, geneGroupsDict, networkAnalysisFolder)

    if runAverageDistanceBetweenGroups:
        networkAnalysis.calculateAverageDistanceBetweenGroups(averageDistanceBootstrapNumber)
    
    if runDSD:
        dfRWRScoresDict=networkAnalysis.calculateDiffusionStateDistance()
    
    if runRWROverlap:
        networkAnalysis.calculateRWRTopNodesOverlap(rwrTopNumber, orsumNumberOfTermsToPlot, rwrOverlapBootstrapNumber)

    networkAnalysis.fLog.close()

