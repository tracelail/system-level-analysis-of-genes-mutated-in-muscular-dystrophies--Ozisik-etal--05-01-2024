#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 18:02:57 2023

@author: Ozan Ozisik
"""

import networkx as nx
import random
from sys import exit
import numpy as np
import pandas as pd
import os
import multixrank
#from scipy.stats import gmean
from enrichmentAnalysis import performEnrichmentAnalysis, applyOrsum


class NetworkAnalysis:
    def __init__(self, networkPaths, geneGroupsDict, networkResultFolder):
        '''
        Parameters
        ----------
        networkPaths : list
            List of paths for the interactions files.
        geneGroupsDict : dict
            Dictionary mapping gene group names to gene lists.
        networkResultFolder : string
            Network result folder

        Returns
        -------
        None.

        '''
        try:
            
            self.networkResultFolder=networkResultFolder
            if not os.path.exists(self.networkResultFolder):
                os.makedirs(self.networkResultFolder)
                print('Network result folder \"'+self.networkResultFolder+'\" is created.')
            
            self.fLog=open(networkResultFolder+os.path.sep+'log.txt', 'w')
            
            self.networkPaths=networkPaths
            self.singleNetwork=nx.empty_graph()#This is a single layer network, when multiple network files are given, they are merged.
            
            for networkPath in networkPaths:
                self.singleNetwork=nx.compose(self.singleNetwork, nx.read_edgelist(networkPath))
                self.printAndLogText('Number of nodes {}'.format(self.singleNetwork.number_of_nodes()) )
                self.printAndLogText('Number of edges {}'.format(self.singleNetwork.number_of_edges()) )
            
            self.singleNetwork.remove_edges_from(nx.selfloop_edges(self.singleNetwork))
            self.printAndLogText('Removing self edges.')
            self.printAndLogText('Number of nodes {}'.format(self.singleNetwork.number_of_nodes()) )
            self.printAndLogText('Number of edges {}'.format(self.singleNetwork.number_of_edges()) )
            
            self.singleNetwork=self.singleNetwork.subgraph(max(nx.connected_components(self.singleNetwork), key=len)).copy()
            self.printAndLogText('Taking the largest connected component in the network.')
            self.printAndLogText('Number of nodes {}'.format(self.singleNetwork.number_of_nodes()) )
            self.printAndLogText('Number of edges {}'.format(self.singleNetwork.number_of_edges()) )
            self.printAndLogText('')
            
            nx.write_edgelist(self.singleNetwork, self.networkResultFolder+os.sep+'aggregatedNetwork.tsv', delimiter='\t', data=False)
            
            self.geneGroupsDict=dict()
            self.allGenesInGroups=[]
            for group, genes in geneGroupsDict.items():
                genesInNetwork=[g for g in genes if g in self.singleNetwork.nodes()]
                self.printAndLogText('From ' + group + ' the following genes are not in the network and will not be used:')
                self.printAndLogText(str(set(genes).difference(set(genesInNetwork))))
                
                self.geneGroupsDict[group]=genesInNetwork
                self.allGenesInGroups=self.allGenesInGroups+genesInNetwork
            
            
            os.fsync(self.fLog.fileno())
            
            
            
        except FileNotFoundError as fnfe:
            print('Network file could not be found.')
            print(fnfe)
            exit()
        except IOError as ioe:
            print('IO Error in network analysis.')
            print(ioe)
            exit()

    def printAndLogText(self, text):
        print(text)
        self.fLog.write(text)
        self.fLog.write('\n')


    def calculateAverageDistanceBetweenTwoGeneLists(self, network, genesList1, genesList2, discardDistanceToSelf=True):
        
        avgDist=0.0
        s=0

        for i in range(len(genesList1)):
            node1=genesList1[i]

            for j in range(len(genesList2)):
                node2=genesList2[j]

                if (not discardDistanceToSelf) or (node1!=node2):
                    shortestPath=nx.bidirectional_shortest_path(network, node1, node2)
                    dist=len(shortestPath)-1
                    avgDist=avgDist+dist
                    s=s+1
        avgDist=avgDist/s
        return avgDist
    

    
    def calculateAverageDistanceBetweenGroups(self, bootstrapNumber=2000):
        '''
        Calculates the intra-group and inter-group average distances.
        This method uses the single layer network (merged network if multiple 
        networks are given to as network files)
        '''
        
        
        print('\n\n')
        print('Calculating average distance between groups.')
        
        
        groups=[k for k in self.geneGroupsDict.keys()]
        
        
        avgDistMtr=np.empty([len(groups),len(groups)],dtype=(float))
        
        # Stores distances to random gene groups, of same size with the groups and also with random sizes
        avgDistRandomMtr=np.empty([len(groups)+1,len(groups)+1],dtype=(float)) 
        avgDistRandomMtr.fill(np.nan)
        
        avgDistSignificanceMtr=np.empty([len(groups),len(groups)],dtype=(float))
        
        
        ## Intra-group and inter-group distances are calculated
        ## Significance of these distances is calculated by comparing to 
        ## distances to same size random gene groups
        ## also average distances to same size random gene groups are calculated
        for i in range(len(groups)):
            group1=groups[i]
            geneList1=self.geneGroupsDict[group1]
            
            for j in range(len(groups)):
                group2=groups[j]
                geneList2=self.geneGroupsDict[group2]
                
                distance=self.calculateAverageDistanceBetweenTwoGeneLists(self.singleNetwork, geneList1, geneList2)
                print('Average distance between '+group1+' and '+group2+' is '+str(round(distance, 2)))
                avgDistMtr[i,j]=distance
                #avgDistMtr[j,i]=distance
                
                ## In bootstrap, we check whether distance is significant 
                ## compared to random gene groups
                significanceCloser=0
                distanceToRandomGroupSum=0
                for bs in range(bootstrapNumber):
                    randomGroup=random.sample(list(nx.nodes(self.singleNetwork)),len(geneList2))
                    distanceToRandomGroup=self.calculateAverageDistanceBetweenTwoGeneLists(self.singleNetwork, geneList1, randomGroup)
                    distanceToRandomGroupSum=distanceToRandomGroupSum+distanceToRandomGroup
                    if distanceToRandomGroup<distance:
                        significanceCloser=significanceCloser+1
                significanceCloser=significanceCloser/bootstrapNumber
                avgDistSignificanceMtr[i,j]=significanceCloser
                distanceToRandomGroupAvg=(distanceToRandomGroupSum/bootstrapNumber)
                avgDistRandomMtr[i,j]=distanceToRandomGroupAvg
        
        dfAvgDist=pd.DataFrame(data=avgDistMtr, index=groups, columns=groups)
        dfAvgDist=dfAvgDist.round(2)
        dfAvgDist.to_csv(self.networkResultFolder+os.sep+'averageDistances.csv', sep='\t')
        
        dfAvgDistSignificance=pd.DataFrame(data=avgDistSignificanceMtr, index=groups, columns=groups)
        dfAvgDistSignificance.to_csv(self.networkResultFolder+os.sep+'averageDistanceSignificances.csv', sep='\t')    
        
        
        ## Finding the smallest and biggest group size to be used in random gene group generation
        groupSizes=[]
        for group in groups:
            groupSizes.append(len(self.geneGroupsDict[group]))
        smallestGroupSize=min(groupSizes)
        biggestGroupSize=max(groupSizes)
        
        
        
        ## Average distance between groups and random gene groups of random size are calculated
        for i in range(len(groups)):
            group1=groups[i]
            geneList1=self.geneGroupsDict[group1]
            
            distanceToRandomGroupSum=0
            
            for bs in range(bootstrapNumber):
                randomGroupSize=random.randint(smallestGroupSize, biggestGroupSize) #both boundaries are included in randint
                randomGroup=random.sample(list(nx.nodes(self.singleNetwork)), randomGroupSize)
                distanceToRandomGroup=self.calculateAverageDistanceBetweenTwoGeneLists(self.singleNetwork, geneList1, randomGroup)
                distanceToRandomGroupSum=distanceToRandomGroupSum+distanceToRandomGroup

            distanceToRandomGroupAvg=(distanceToRandomGroupSum/bootstrapNumber)
            avgDistRandomMtr[i,len(groups)]=distanceToRandomGroupAvg
        
        ## Average distance between two random groups are calculated
        distanceToRandomGroupSum=0
        for bs in range(bootstrapNumber):
            randomGroupSize=random.randint(smallestGroupSize, biggestGroupSize)
            randomGroupSize2=random.randint(smallestGroupSize, biggestGroupSize)
            randomGroup=random.sample(list(nx.nodes(self.singleNetwork)), randomGroupSize)
            randomGroup2=random.sample(list(nx.nodes(self.singleNetwork)), randomGroupSize2)
            distanceToRandomGroup=self.calculateAverageDistanceBetweenTwoGeneLists(self.singleNetwork, randomGroup, randomGroup2)
            distanceToRandomGroupSum=distanceToRandomGroupSum+distanceToRandomGroup

        distanceToRandomGroupAvg=(distanceToRandomGroupSum/bootstrapNumber)
        avgDistRandomMtr[len(groups),len(groups)]=distanceToRandomGroupAvg
        
        
        
        dfAvgDistRandom=pd.DataFrame(data=avgDistRandomMtr, index=groups+['Random'], columns=[group+'-Sized' for group in groups]+['Random-Sized'])
        dfAvgDistRandom=dfAvgDistRandom.round(2)
        dfAvgDistRandom.to_csv(self.networkResultFolder+os.sep+'averageDistancesRandom.csv', sep='\t')
        
        
    
    
        
    
    def calculateDiffusionStateDistance(self):
        '''
        This method calculates diffusion state distance.
        It runs RWR starting from each gene and finds the RWR score differences
        between them.
        multiXrank returns a score per multiplex per layer per node.
        We use one multiplex but we might have multiple layers.
        The scores are compared as they are, per layer and per node.
        
        
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3806810/

        Returns
        -------
        dfRWRScoresDict : TYPE
            DESCRIPTION.

        '''
        
        print('\n\n')
        print('Calculating diffusion state distances')
        
        
        rwrFolder=self.networkResultFolder+os.sep+'rwrFolder'
        seedsPath='seeds.txt'
        
        if not os.path.isdir(rwrFolder):
            os.mkdir(rwrFolder)
        
        ## Networks are located relative to rwrFolder, so we need to go to
        ## main project directory level
        directoryMoveUpText='..'+os.sep+'..'+os.sep+'..'+os.sep
        
        f=open(rwrFolder+os.sep+'multixrankConfig.yml','w')
        f.write('multiplex:\n')
        f.write('    1:\n')
        f.write('        layers:\n')
        for nwp in self.networkPaths:
            f.write('            - '+directoryMoveUpText+nwp+'\n')
        f.write('\n')
        f.write('seed:\n')
        f.write('    '+seedsPath+'\n')
        f.close()


        dfRWRScoresList=[]
        dfRWRScoresDict=dict()

        for g in self.allGenesInGroups:
            
            print('RWR running with seed',g)
            
            f=open(rwrFolder+'/'+seedsPath, 'w')
            f.write(g+'\n')
            f.close()
            
            multixrank_obj = multixrank.Multixrank(config=rwrFolder+os.sep+'multixrankConfig.yml', wdir=rwrFolder)
            dfRWRScores = multixrank_obj.random_walk_rank()
            
            
            ## We use one multiplex network with one or more layers
            ## Node and layer are set as indices together
            ## Nodes' scores are compared for different layers
            ## This is required for subtraction between two score dataframes
            dfRWRScores=dfRWRScores.drop(columns=['multiplex'])
            
            '''
            
            #Taking the geometric means of the scores across layers
            #dfRWRScoresGM=dfRWRScores.loc[dfRWRScores.score > 0].groupby(['node']).agg({'score': gmean})
            dfRWRScores=dfRWRScores.groupby(['node']).prod() ## For speed.
            '''
            
            
            dfRWRScores=dfRWRScores.set_index(['node', 'layer'])
            
            
            ## RWR score of the seed is removed. The other seeds could also be 
            ## removed but they do not cause any problems. See below.
            dfRWRScores=dfRWRScores.drop([g], level='node')
            
            ## To make probabilities sum up to 1
            dfRWRScores=(dfRWRScores/dfRWRScores.sum())
            
            #dfRWRScores=dfRWRScores.sort_values(by='score', ascending=False)
            dfRWRScores=dfRWRScores.sort_index()
            
            dfRWRScoresList.append(dfRWRScores)
            dfRWRScoresDict[g]=dfRWRScores


        #diffusion state distance (DSD)
        dsdMtr=np.empty([len(self.allGenesInGroups),len(self.allGenesInGroups)],dtype=(float))

        for i in range(len(dfRWRScoresList)-1):
            dsdMtr[i,i]=np.nan
            for j in range(i+1, len(dfRWRScoresList)):
                ## seeds get NaN in subtraction because they are missing, and sum() ignores NaNs.
                
                dsdMtr[i,j]=(dfRWRScoresList[i]-dfRWRScoresList[j]).abs().sum()
                dsdMtr[j,i]=dsdMtr[i,j]
        dsdMtr[len(dfRWRScoresList)-1,len(dfRWRScoresList)-1]=np.nan
        
        dfDSDMtr=pd.DataFrame(data=dsdMtr, index=self.allGenesInGroups, columns=self.allGenesInGroups)
        dfDSDMtr.to_csv(self.networkResultFolder+os.sep+'dsdDistances.csv', sep='\t')
        
        meanDSD=np.nanmean(dfDSDMtr.values.flatten())
        stdDSD=np.nanstd(dfDSDMtr.values.flatten())
        dfDSDMtrZ=(dfDSDMtr-meanDSD)/stdDSD
        dfDSDMtrZ.to_csv(self.networkResultFolder+os.sep+'dsdDistancesZScore.csv', sep='\t')

        return dfRWRScoresDict
    
    
    
    
    
    def runRWR(self, seeds, rwrFolder, seedsPath, topN):
        
        f=open(rwrFolder+'/'+seedsPath, 'w')
        for g in seeds:
            f.write(g+'\n')
        f.close()
        
        multixrank_obj = multixrank.Multixrank(config=rwrFolder+os.sep+'multixrankConfig.yml', wdir=rwrFolder)
        dfRWRScores = multixrank_obj.random_walk_rank()
        
        
        #We have one multiplex, we drop multiplex column
        dfRWRScores=dfRWRScores.drop(columns=['multiplex'])
        
        
        ## Taking the geometric means of the scores across layers
        #dfRWRScoresGM=dfRWRScores.loc[dfRWRScores.score > 0].groupby(['node']).agg({'score': gmean})
        
        ## As we are only comparing one to another, for speed, I use product of 
        ## scores across layers instead of geometric mean
        #dfRWRScores=dfRWRScores.astype({'score':float})
        dfRWRScoresGM=dfRWRScores.loc[dfRWRScores.score > 0].groupby(['node']).prod('score') ## For speed.
        
        
        #Dropping the seed genes
        dfRWRScoresGM=dfRWRScoresGM.drop(list(seeds), axis='index')

        #Sorting by score
        dfRWRScoresGM=dfRWRScoresGM.sort_values(by='score', ascending=False)
        
        topNgenes=list(dfRWRScoresGM.iloc[0:topN].index)
        
        
        return topNgenes
    
    
    
    
    
    def calculateRWRTopNodesOverlap(self, topN, orsumNumberOfTermsToPlot, rwrOverlapBootstrapNumber, runEnrichmentForRWRResult=False, gProfilerBaseUrl='https://biit.cs.ut.ee/gprofiler', sources=['GO:BP', 'GO:CC']):
        '''
        multiXrank returns a score per multiplex per layer per node.
        We use one multiplex but we might have multiple layers.
        For taking the top nodes we group the scores by nodes
        and we take the minimum.
        '''
        
        
        print('\n\n')
        print('Calculating overlap of top nodes from random walk with restart seeded by different groups\n')
        
        rwrFolder=self.networkResultFolder+os.sep+'rwrFolder'
        seedsPath='seeds.txt'
        
        if not os.path.isdir(rwrFolder):
            os.mkdir(rwrFolder)
        
        ## Networks are located relative to rwrFolder, so we need to go to
        ## main project directory level
        directoryMoveUpText='..'+os.sep+'..'+os.sep+'..'+os.sep
        
        f=open(rwrFolder+os.sep+'multixrankConfig.yml','w')
        f.write('multiplex:\n')
        f.write('    1:\n')
        f.write('        layers:\n')
        for nwp in self.networkPaths:
            f.write('            - '+directoryMoveUpText+nwp+'\n')
        f.write('\n')
        f.write('seed:\n')
        f.write('    '+seedsPath+'\n')
        f.close()
        
        
        rwrTopGenesDict=dict()
        
        for group, genes in self.geneGroupsDict.items():
            
            rwrTopGenesDict[group]=self.runRWR(genes, rwrFolder, seedsPath, topN)
            
        
        groups=[k for k in self.geneGroupsDict.keys()]
        
        
        
        
        ## Writing the topN genes (does not include the seed genes)
        for group in groups:
            topG=rwrTopGenesDict[group]
            with open(self.networkResultFolder+os.sep+'RWRTop'+str(topN)+'-'+group+'.txt', 'w') as f:
                for g in topG:
                    f.write(g+'\n')
        
        
        ## Creating a network file, containing edges between topN genes.
        ## Seed genes are also included in the network.
        for group in groups:
            topG=rwrTopGenesDict[group]
            topG=topG+self.geneGroupsDict[group]
            with open(self.networkResultFolder+os.sep+'RWRTop'+str(topN)+'-Network-'+group+'.txt', 'w') as f:
                for i in range(len(topG)-1):
                    g1=topG[i]
                    for j in range(i+1, len(topG)):
                        g2=topG[j]
                        if self.singleNetwork.has_edge(g1, g2):
                            f.write(g1+'\t'+g2+'\n')
        
        
        ## Writing the overlap between topN genes obtained from different groups
        f=open(self.networkResultFolder+os.sep+'RWRTop'+str(topN)+'.txt', 'w')
        for i in range(len(groups)-1):
            group1=groups[i]
            topG1=set(rwrTopGenesDict[group1])
            for j in range(i+1, len(groups)):
                group2=groups[j]
                topG2=set(rwrTopGenesDict[group2])
                
                f.write(group1 +' ' + group2+'\n\n')
                
                f.write(group1+' in the top RWR genes of '+group2+'\n')
                for g in topG2.intersection(set(self.geneGroupsDict[group1])):
                    f.write(g+', ')
                f.write('\n')
                f.write(group2+' in the top RWR nodes of '+group1+'\n')
                for g in topG1.intersection(set(self.geneGroupsDict[group2])):
                    f.write(g+', ')
                f.write('\n')
                
                f.write('Overlap of top genes between '+group1+' and ' + group2+'\n')
                for g in topG1.intersection(topG2):
                    f.write(g+', ')
                f.write('\n')
                
                #startTime=time.perf_counter()
                
                significanceSmall, significanceLarge=self.rwrOverlapBootstrap(rwrFolder, seedsPath, topN, self.allGenesInGroups, len(self.geneGroupsDict[group1]), len(self.geneGroupsDict[group2]), len(topG1.intersection(topG2)), bootstrapNumber=rwrOverlapBootstrapNumber)
                
                #print('time', time.perf_counter()-startTime)
                
                f.write('significance for smaller overlap than expected '+str(significanceSmall)+'/'+ str(rwrOverlapBootstrapNumber)+'\n')
                f.write('significance for larger overlap than expected '+str(significanceLarge)+'/'+ str(rwrOverlapBootstrapNumber)+'\n')
                
                f.write('\n\n\n')
                
                
        f.close()
        
        
        if runEnrichmentForRWRResult:
            enrichmentFolder=self.networkResultFolder+os.sep+'RWREnrichmentTop'+str(topN)
            summaryFolder=enrichmentFolder+os.sep+ 'Summary'
            
            #I used gProfiler version e108_eg55_p17_9f356ae. If g:Profiler updates the data, to reproduce the results
            #gProfiler_base_url should be set to 'https://biit.cs.ut.ee/gprofiler_archive3/e108_eg55_p17'.
            performEnrichmentAnalysis(rwrTopGenesDict, enrichmentFolder, summaryFolder, sources=sources, gProfiler_base_url=gProfilerBaseUrl)
            
            for source in sources:
                applyOrsum(gmt='Data/Annotation/hsapiens.'+source+'.name.gmt', 
                           source=source, 
                           geneGroupsDict=rwrTopGenesDict,
                           summaryFolder=summaryFolder,
                           outputFolder='orsum'+source.replace(':',''),
                           numberOfTermsToPlot=orsumNumberOfTermsToPlot
                           )
        
    
    
    def rwrOverlapBootstrap(self, rwrFolder, seedsPath, topN, genes, size1, size2, realOverlapSize, bootstrapNumber=2000):
        '''
        From the given list of genes, creates two gene groups of size1 and size2
        Runs RWR and finds their overlap
        '''
        
        genesShuffled=genes.copy()
        
        significanceOverlapLargerThanExpected=0
        significanceOverlapSmallerThanExpected=0
        
        for bs in range(bootstrapNumber):
            
            if bs%10==0:
                print('rwrOverlapBootstrap', bs, '\n')
            
            random.shuffle(genesShuffled)
            group1genes=genesShuffled[0:size1]
            group2genes=genesShuffled[size1:(size1+size2)]
            
            rwrTopGenes1=self.runRWR(group1genes, rwrFolder, seedsPath, topN)
            rwrTopGenes2=self.runRWR(group2genes, rwrFolder, seedsPath, topN)
            
            overlapSize=len(set(rwrTopGenes1).intersection(rwrTopGenes2))
            
            #print(size1, size2, overlapSize)
            
            if overlapSize<=realOverlapSize:
                significanceOverlapSmallerThanExpected=significanceOverlapSmallerThanExpected+1 #Significance drops
                
            if overlapSize>=realOverlapSize:
                significanceOverlapLargerThanExpected=significanceOverlapLargerThanExpected+1 #Significance drops
        
        
        return significanceOverlapSmallerThanExpected, significanceOverlapLargerThanExpected