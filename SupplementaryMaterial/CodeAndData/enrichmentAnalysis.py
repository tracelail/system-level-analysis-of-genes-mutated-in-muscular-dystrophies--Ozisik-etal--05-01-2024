#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 13:49:08 2023

@author: Ozan Ozisik
"""

from gprofiler import GProfiler
import os


def performEnrichmentAnalysis(geneGroupsDict, enrichmentFolder, summaryFolder, organism='hsapiens', sources=['GO:BP', 'GO:CC'], gProfiler_base_url = 'https://biit.cs.ut.ee/gprofiler', createOrsumInput=True):
    '''
    Parameters
    ----------
    geneGroupsDict : dict
        Dictionary mapping gene group names to gene lists.
    enrichmentFolder : string
        Enrichment folder path.
    summaryFolder : string
        Summary folder path, this will be used for orsum.
    organism : string, optional
        Organism name for g:Profiler. The default is 'hsapiens'.
    sources : list, optional
        List of source names for g:Profiler. The default is ['GO:BP', 'GO:CC'].
    gProfiler_base_url : string, optional
        The URL for the g:Profiler request. Should be changed if an older g:Profiler version is used. The default is the most recent version.
    createOrsumInput : bool, optional
        Whether enrichment ID files to use with orsum will be created. The default is True.

    Returns
    -------
    None.

    '''
    
    if createOrsumInput:
        if not os.path.exists(summaryFolder):
            os.makedirs(summaryFolder)
            print('Enrichment summary folder \"'+summaryFolder+'\" is created.')
    
    
    gp = GProfiler(return_dataframe=True, base_url=gProfiler_base_url)
    
    dfEnrichmentAll=gp.profile(organism=organism,
                            sources=sources,
                            query=geneGroupsDict,
                            no_evidences=False,#adds the intersecions and evidences columns
                            all_results=True
                            )
    
    dfEnrichmentAll.to_csv(enrichmentFolder+os.sep+'enrichmentAllIncludingInsignificant.csv', sep='\t', index=False)
    
    dfEnrichmentSignificant=dfEnrichmentAll[dfEnrichmentAll['significant']==True]
    
    for group in geneGroupsDict.keys():
        dfEnrichmentGroup=dfEnrichmentSignificant[dfEnrichmentSignificant['query']==group]
        
        dfEnrichmentGroup=dfEnrichmentGroup.sort_values(by=['source', 'p_value'], ascending=[True, True])
        dfEnrichmentGroup.to_csv(enrichmentFolder+os.sep+'enrichment'+ group +'.csv', sep='\t', index=False)
        
        for source in sources:
            dfEnrichmentGroupSource=dfEnrichmentGroup[dfEnrichmentGroup['source']==source]
            dfEnrichmentGroupSource=dfEnrichmentGroupSource.sort_values(by='p_value', ascending=True)
            dfEnrichmentGroupSource.to_csv(enrichmentFolder+os.sep+'enrichment'+ group + '-'+ source.replace(':','')+'.csv', sep='\t', index=False)
            
            if createOrsumInput:
                dfEnrichmentGroupSource['native'].to_csv(summaryFolder+os.sep+'enrichmentIDs'+ group + '-'+ source.replace(':','')+'.csv', index=False, header=None)



def applyOrsum(gmt, source, geneGroupsDict, summaryFolder, outputFolder=None, maxRepSize=None, maxTermSize=None, minTermSize=None, numberOfTermsToPlot=None):
    '''
    

    Parameters
    ----------
    gmt : string
        GMT file path.
    source : string
        Source name for g:Profiler, e.g. GO:CC.
    geneGroupsDict : dict
        Dictionary mapping gene group names to gene lists.
    summaryFolder : string
        Summary folder path, contains enrichment ID files for orsum.
    outputFolder : string, optional
        Output folder for orsum, will be created under summary folder. See orsum documentation.
    maxRepSize : integer, optional
        See orsum documentation. The default is None, orsum default will be used.
    maxTermSize : integer, optional
        See orsum documentation. The default is None, orsum default will be used.
    minTermSize : integer, optional
        See orsum documentation. The default is None, orsum default will be used.
    numberOfTermsToPlot : integer, optional
        See orsum documentation. The default is None, orsum default will be used.

    Returns
    -------
    None.

    '''

    
    noEnrichmentGroup=set()
    
    command='orsum.py '
    command=command+'--gmt \"'+gmt+'\" '
    
    command=command+'--files '
    for group in geneGroupsDict.keys():
        filePath=summaryFolder+os.sep+'enrichmentIDs'+ group + '-'+ source.replace(':','')+'.csv'
        if os.stat(filePath).st_size == 0: # in a previous version orsum used to exit when one enrichment file was empty, now it is fixed but we can keep this code.
            noEnrichmentGroup.add(group)
        else:
            command=command+'\"'+filePath+'\" '
    
    command=command+'--fileAliases '
    for group in geneGroupsDict.keys():
        if group not in noEnrichmentGroup:
            command=command+'\"'+group+'\" '
    
    if outputFolder is not None:
        command=command+'--outputFolder \"'+summaryFolder+os.sep+outputFolder+'\" '
    
    if maxRepSize is not None:
        command=command+'--maxRepSize ' + str(maxRepSize) + ' '
    
    if maxTermSize is not None:
        command=command+'--maxTermSize ' + str(maxTermSize) + ' '
    
    if minTermSize is not None:
        command=command+'--minTermSize ' + str(minTermSize) + ' '
    
    if numberOfTermsToPlot is not None:
        command=command+'--numberOfTermsToPlot '+str(numberOfTermsToPlot)
    
    
    print('orsum command:')
    print(command)
    os.system(command)
    
    
    '''
    ## This is an example command that the code above creates
    os.system('orsum.py --gmt "Data/Other/hsapiens.GO:CC.name.gmt" --files '+
              '"Results/Enrichment/Summary/enrichmentIDsDoGs-GO:CC.csv" '+
              '"Results/Enrichment/Summary/enrichmentIDsCoGs-GO:CC.csv" '+
              '"Results/Enrichment/Summary/enrichmentIDsPoGs-GO:CC.csv" '+
              '--fileAliases DoGs CoGs PoGs '+
              '--outputFolder "Results/Enrichment/Summary/orsumGOCC" '+
              ''--numberOfTermsToPlot 50)
    '''
