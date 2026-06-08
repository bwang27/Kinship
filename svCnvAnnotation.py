#!/usr/bin/env python
# coding: utf-8

# ### Annotate SVs/CNVs

import os
import sys
import argparse, textwrap
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

chrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

def checkSvFileExists(filename, varColumns):
    """
    Check if SV input file exists
    """
    if os.path.isfile(filename):
        svMetaList = []
        nCols = len(varColumns)
        if nCols == 1:
            varColumn = varColumns[0]
            svs = list(pd.read_csv(filename, sep = '\t', usecols=[varColumn])[varColumn])
        elif nCols == 6:
            varColumn = 'curated_SV'
            svTmpDf = pd.read_csv(filename, sep = '\t', usecols=varColumns)
            if 'chr' not in list(svTmpDf[varColumns[0]].astype('str'))[0]:
                svTmpDf[varColumns[0]] = svTmpDf[varColumns[0]].astype('str')
                svTmpDf[varColumns[0]] = 'chr' + svTmpDf[varColumns[0]]
                svTmpDf[varColumns[0]] = [i.replace('.', '_') for i in svTmpDf[varColumns[0]]]
            else:
                svTmpDf[varColumns[0]] = svTmpDf[varColumns[0]].astype('str')
                svTmpDf[varColumns[0]] = [i.replace('.', '_') for i in svTmpDf[varColumns[0]]]
                
            svTmpDf[varColumns[1]] = svTmpDf[varColumns[1]].astype('str')    
            svTmpDf[varColumns[2]] = svTmpDf[varColumns[2]].astype('str')
            if 'chr' not in list(svTmpDf[varColumns[3]].astype('str'))[0]:
                svTmpDf[varColumns[3]] = svTmpDf[varColumns[3]].astype('str')
                svTmpDf[varColumns[3]] = 'chr' + svTmpDf[varColumns[3]]
                svTmpDf[varColumns[3]] = [i.replace('.', '_') for i in svTmpDf[varColumns[3]]]
            else:
                svTmpDf[varColumns[3]] = svTmpDf[varColumns[3]].astype('str')
                svTmpDf[varColumns[3]] = [i.replace('.', '_') for i in svTmpDf[varColumns[3]]]
            
            svTmpDf[varColumns[4]] = svTmpDf[varColumns[4]].astype('str')
            svTmpDf[varColumns[5]] = svTmpDf[varColumns[5]].astype('str')
            svTmpDf['curated_SV'] = svTmpDf[varColumns[0]] + "." + svTmpDf[varColumns[1]] + "." + svTmpDf[varColumns[2]] + "." + svTmpDf[varColumns[3]] + "." + svTmpDf[varColumns[4]] + "." + svTmpDf[varColumns[5]]
            svs = list(svTmpDf['curated_SV'])
        for sv in svs:
            svMetaList.append(sv.split('.')[0:6])
        svFrame = pd.DataFrame(svMetaList, columns=['ChrA', 'PosA', 'StrandA', 'ChrB', 'PosB', 'StrandB'])
        svFrame['PosA'] = svFrame['PosA'].astype('int')
        svFrame['PosB'] = svFrame['PosB'].astype('int')
        svFrame[varColumn] = svs
        subDf = svFrame[svFrame['ChrA'].isin(chrs) & svFrame['ChrA'].isin(chrs)].drop_duplicates() 
        #print(set(subDf['ChrA']), set(subDf['ChrB']))
        return subDf
    else:
        print(filename, " : File does not exist!")
        sys.exit()


def getSVinDf(filename, varColumns):
    """
    Format input SV file based on number of input columns
    """
    nCols = len(varColumns)
    if nCols == 1:
        varColumn = varColumns[0]
        svTmpDf = pd.read_csv(filename, sep = '\t')
        svTmpDf['curated_SV'] = svTmpDf[varColumns[0]].astype('str')
        if 'chr' not in list(svTmpDf['curated_SV'])[0]:
            svTmpDf['curated_SV'] = 'chr' + svTmpDf['curated_SV']
        return svTmpDf            
    elif nCols == 6:
        varColumn = 'curated_SV'
        svTmpDf = pd.read_csv(filename, sep = '\t')
        if 'chr' not in list(svTmpDf[varColumns[0]].astype('str'))[0]:
            svTmpDf[varColumns[0]] = svTmpDf[varColumns[0]].astype('str')
            svTmpDf[varColumns[0]] = 'chr' + svTmpDf[varColumns[0]]
            svTmpDf[varColumns[0]] = [i.replace('.', '_') for i in svTmpDf[varColumns[0]]]
        else:
            svTmpDf[varColumns[0]] = svTmpDf[varColumns[0]].astype('str')
            svTmpDf[varColumns[0]] = [i.replace('.', '_') for i in svTmpDf[varColumns[0]]]

        svTmpDf[varColumns[1]] = svTmpDf[varColumns[1]].astype('str')
        svTmpDf[varColumns[2]] = svTmpDf[varColumns[2]].astype('str')
        
        if 'chr' not in list(svTmpDf[varColumns[3]].astype('str'))[0]:
            svTmpDf[varColumns[3]] = svTmpDf[varColumns[3]].astype('str')
            svTmpDf[varColumns[3]] = 'chr' + svTmpDf[varColumns[3]]
            svTmpDf[varColumns[3]] = [i.replace('.', '_') for i in svTmpDf[varColumns[3]]]
        else:
            svTmpDf[varColumns[3]] = svTmpDf[varColumns[3]].astype('str')
            svTmpDf[varColumns[3]] = [i.replace('.', '_') for i in svTmpDf[varColumns[3]]]

        svTmpDf[varColumns[4]] = svTmpDf[varColumns[4]].astype('str')
        svTmpDf[varColumns[5]] = svTmpDf[varColumns[5]].astype('str')
        svTmpDf['curated_SV'] = svTmpDf[varColumns[0]] + "." + svTmpDf[varColumns[1]] + "." + svTmpDf[varColumns[2]] + "." + svTmpDf[varColumns[3]] + "." + svTmpDf[varColumns[4]] + "." + svTmpDf[varColumns[5]]
        return svTmpDf[svTmpDf[varColumns[0]].isin(chrs) & svTmpDf[varColumns[3]].isin(chrs)]
        
        
def checkCnvFileExists(filename, varColumns=['Chr', 'loc.start', 'loc.end']):
    """
    Check if CNV input file exists
    """
    if os.path.isfile(filename):
        cnvs = pd.read_csv(filename, sep = '\t', usecols=varColumns)
        return cnvs.drop_duplicates(ignore_index=True)
    else:
        print(filename, " : File does not exist!")
        sys.exit()

           
def getLargestIsoforms(geneRegionsDf):
    """
    Compute the largest Isoform for each gene
    """
    geneRegionsGrp = geneRegionsDf.reset_index().groupby(['Gene', 'RefSeq']).aggregate({'Chr':'first','Start':'min', 'End':'max', 'Strand':'first'}).reset_index()
    geneRegionsGrp['Length'] = geneRegionsGrp['End'] - geneRegionsGrp['Start']
    maxIsoforms = geneRegionsGrp.sort_values(by='Length', ascending=False).groupby('Gene').aggregate({'RefSeq':'first', 'Strand':'first','Chr':'first', 'Start':'first', 'End':'first', 'Length':'first'}).reset_index()
    largeIsoDict = dict(zip(list(maxIsoforms['Gene']), list(maxIsoforms['RefSeq'])))
    largeIsoDict['ND'] = 'ND'
    return largeIsoDict


def getDefaultIsoforms():
    """
    Annotate position with default gene isoform
    """
    isoFn = "/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/xmagrp/CohortTracker/data/defaultIsoform.human"
    iso = pd.read_csv(isoFn, sep='\t', names=['Gene', 'RefSeq'])
    defIsoDict = dict(zip(list(iso['Gene']), list(iso['RefSeq'])))
    defIsoDict['ND'] = 'ND'
    return defIsoDict


def searchUpstreamGenes(chrom, pos, fusionGenes, strand, genes, startPad=0, partner='A', location='upstream', searchStep=1e6):
    """
    Search upstream genes of the given position beyond specified padding
    """
    geneCount = 0
    searchPad = startPad
    setIndex = set(fusionGenes)
    tempDf = genes.reset_index()
    maxPos = max(tempDf['End']) + 2e6
    minPos = max(min(tempDf['Start']) - 2e6, 1)
    newPos = pos
    
    
    def searchLeft(geneCount, searchPad, searchStep, minPos, newPos, genes):
        """Search in left direction on the genome"""
        geneList = []
        while geneCount == 0 and newPos > minPos:
            searchPad += searchStep
            newPos = max(pos - searchPad, 1)
            subDf = genes.loc[(genes.index.get_level_values('Chr') == chrom) & (genes.index.get_level_values('Start') >= newPos) & (genes.index.get_level_values('End') <= pos)].reset_index()
            #groups = genes.loc[(genes.index.get_level_values('Chr') == chrom) & (genes.index.get_level_values('Start') >= newPos) & (genes.index.get_level_values('End') <= pos)].reset_index().groupby('Gene')
            subDf.sort_values('End', ascending=False, inplace=True)
            subDf['Gene_strand'] = subDf['Gene'] + '_(' + subDf['Strand'] + ')'
            geneList = list(subDf['Gene'].unique())
            geneStrandList = list(subDf['Gene_strand'].unique())
            geneDict = dict(zip(geneList, geneStrandList))
            rmFusGenes = [geneList.remove(i) for i in fusionGenes if i in geneList]
            #setAll = set(groups.groups.keys())
            #nGenes = setAll.difference(setIndex)
            if len(geneList) > 0:
                geneCount += 1

        if len(geneList) > 0 and len(geneList) <= 5:
            retGenes = ','.join([geneDict[i] for i in geneList])
        elif len(geneList) > 5:
            retGenes = ','.join([geneDict[i] for i in geneList[0:5]])
        else:
            retGenes = 'UD'
            
        gDist = str(searchPad/1e6) + 'Mb'
        return [retGenes, gDist]

    def searchRight(geneCount, searchPad, searchStep, maxPos, newPos, genes):
        """Search in Right direction on the genome"""
        geneList = []
        while geneCount == 0 and newPos <= maxPos:
            searchPad += searchStep
            newPos = pos + searchPad
            subDf = genes.loc[(genes.index.get_level_values('Chr') == chrom) & (genes.index.get_level_values('Start') >= pos) & (genes.index.get_level_values('End') <= newPos)].reset_index()
            subDf.sort_values('Start', ascending=True, inplace=True)
            subDf['Gene_strand'] = subDf['Gene'] + '_(' + subDf['Strand'] + ')'
            geneList = list(subDf['Gene'].unique())
            geneStrandList = list(subDf['Gene_strand'].unique())
            geneDict = dict(zip(geneList, geneStrandList))
            rmFusGenes = [geneList.remove(i) for i in fusionGenes if i in geneList]
            
            #groups = genes.loc[(genes.index.get_level_values('Chr') == chrom) & (genes.index.get_level_values('Start') >= pos) & (genes.index.get_level_values('End') <= newPos)].reset_index().groupby('Gene')
            #setAll = set(groups.groups.keys())
            #nGenes = setAll.difference(setIndex)
            if len(geneList) > 0:
                geneCount += 1
        
        if len(geneList) > 0 and len(geneList) <= 5:
            retGenes = ','.join([geneDict[i] for i in geneList])
        elif len(geneList) > 5:
            retGenes = ','.join([geneDict[i] for i in geneList[0:5]])
        else:
            retGenes = 'UD'
            
        gDist = str(searchPad/1e6) + 'Mb'
        return [retGenes, gDist]
    
    if strand == '+' and partner == 'A' and location == 'upstream': ## Strand + and partner A
        return searchLeft(geneCount, searchPad, searchStep, minPos, newPos, genes)
    elif strand == '+' and partner == 'A' and location == 'downstream': ## Strand + and partner A
        return searchRight(geneCount, searchPad, searchStep, maxPos, newPos, genes)
    elif strand == '-' and partner == 'A' and location == 'upstream': ## Strand - and partner A
        return searchRight(geneCount, searchPad, searchStep, maxPos, newPos, genes)
    elif strand == '-' and partner == 'A' and location == 'downstream': ## Strand - and partner A
        return searchLeft(geneCount, searchPad, searchStep, minPos, newPos, genes)
    elif strand == '+' and partner == 'B' and location == 'downstream': ## Strand + and partner B
        return searchRight(geneCount, searchPad, searchStep, maxPos, newPos, genes)
    elif strand == '+' and partner == 'B' and location == 'upstream': ## Strand + and partner B
        return searchLeft(geneCount, searchPad, searchStep, minPos, newPos, genes)
    elif strand == '-' and partner == 'B' and location == 'downstream': ## Strand - and partner B
        return searchLeft(geneCount, searchPad, searchStep, minPos, newPos, genes)
    elif strand == '-' and partner == 'B' and location == 'upstream': ## Strand - and partner B
        return searchRight(geneCount, searchPad, searchStep, maxPos, newPos, genes)


def getUpstreamGenes(chrom, pos, fusionGenes, strand, genes, padding=2e6):
    """
    Get genes upstream of the given position with specified padding
    """
    
    if strand == '+':
        newPos = max(pos - padding, 1)
        try:
            groups = genes.loc[(genes.index.get_level_values('Chr') == chrom) & (genes.index.get_level_values('Start') >= newPos) & (genes.index.get_level_values('End') <= pos)].reset_index().groupby('Gene')
            setAll = set(groups.groups.keys())
            setIndex = set(fusionGenes)
            nGenes = setAll.difference(setIndex)
            if len(nGenes) > 0:
                retGenes = ','.join(nGenes)
            else:
                retGenes = 'ND'
        except:
            retGenes = 'ND'
        return retGenes
    elif strand == '-':
        newPos = pos + padding
        try:
            groups = genes.loc[(genes.index.get_level_values('Chr') == 'chrom') & (genes.index.get_level_values('Start') >= pos) & (genes.index.get_level_values('End') <= newPos)].reset_index().groupby('Gene')
            setAll = set(groups.groups.keys())
            setIndex = {geneName}
            nGenes = setAll.difference(setIndex)
            if len(nGenes) > 0:
                retGenes = ','.join(nGenes)
            else:
                retGenes = 'ND'
        except:
            retGenes = 'ND'
        return retGenes
       
    
def getGeneName(chrom, pos, genes):
    """
    Get gene name based on position
    """
    try:
        geneName = genes.loc[(genes.index.get_level_values('Chr') == chrom) & (genes.index.get_level_values('Start') <= pos) & (genes.index.get_level_values('End') >= pos)].reset_index().groupby('Gene').aggregate({'Gene':'first'}).iloc[0]['Gene']
    except:
        geneName = 'ND'
    return geneName


def getGenesInRegion(chrom, start, end, geneRegions, maxNum='all', loc='within', padding=10e6):
    """
    Get gene name based on position
    """
    print(chrom, start, end)
    if loc == 'left':
        end = start - 1
        start = start - padding
        
    if loc == 'right':
        start = end + 1
        end = end + padding
        
    if loc == 'withFlank':
        start = start - padding
        end = end + padding
        
    try:
        geneName = list(geneRegions.loc[(geneRegions['Chr'] == chrom) & (geneRegions['Start'] >= start) & (geneRegions['Start'] <= end), 'Gene'])
        #print(start, end)
        #print(len(geneName))
        if maxNum != 'all' and loc == 'left':
            geneName = geneName[-(maxNum):]
        if maxNum != 'all' and loc == 'right':
            geneName = geneName[:(maxNum)]    
    except:
        geneName = ['ND']
    if len(geneName) == 0:
        geneName = ['ND']
    return geneName


def annotateCnvRegion(chrom, start, end, geneRegions, maxGenes=4, size='small'):
    """
    Annotate large CNV regions, 
    Large: Give 4 genes in left boundary and right boundary each. (total 8)
    Small: Give 4 genes in left boundary and right boundary each (total 8) and genes within given region
    """
    leftGenes = getGenesInRegion(chrom, start, end, geneRegions, maxNum=maxGenes, loc='left')
    rightGenes = getGenesInRegion(chrom, start, end, geneRegions, maxNum=maxGenes, loc='right')
    if size == 'small':
        withinGenes = getGenesInRegion(chrom, start, end, geneRegions)
        allG = leftGenes + withinGenes + rightGenes
        allGenes = list(dict.fromkeys(allG))
    elif size == 'large':
        allG = leftGenes + rightGenes
        withinGenes = ['ND']
        allGenes = list(dict.fromkeys(allG))
    return [', '.join(leftGenes), ', '.join(withinGenes), ', '.join(rightGenes)] 
    

def annotatePosWithIsoform(chrom, pos, isoform, genes):
    """
    Annotate position with default isoform of a gene
    """
    try:
        annotDf = genes.loc[(genes.index.get_level_values('RefSeq') == isoform) & (genes.index.get_level_values('Start') <= pos) & (genes.index.get_level_values('End') >= pos)].reset_index()
        retDf = annotDf[['Gene', 'ExonIntron', 'Strand']]
        if len(retDf) < 1:
            retDf = pd.DataFrame({'Gene':['ND'], 'ExonIntron':['ND'], 'Strand':['ND']})
    except:
        retDf = pd.DataFrame({'Gene':['ND'], 'ExonIntron':['ND'], 'Strand':['ND']})
    return retDf

        
def annotateSVListWithUpstreamGenes(svAnnoDf, geneDf, paddings):
    """
    Annotate SVs with neighboring upstream genes
    """
    #svUpDf = svAnnoDf.copy()
    
    for pad in paddings:
        suffix = str(pad/1e6) + 'M'
        for index, row in svAnnoDf.iterrows():
            chrA = row['ChrA_x']
            chrB = row['ChrB_x']
            posA = row['PosA_x']
            posB = row['PosB_x']
            strandA = row['SVIndex'].split('.')[2]
            strandB = row['SVIndex'].split('.')[5]
            geneA = row['Gene_default_A']
            geneB = row['Gene_default_B']
            geneList = [geneA, geneB]
            svAnnoDf.loc[index, 'upGene_'  + suffix + '_A'] = getUpstreamGenes(chrA, posA, geneList, strandA, geneDf, padding=pad)
            svAnnoDf.loc[index, 'upGene_'  + suffix + '_B'] = getUpstreamGenes(chrB, posB, geneList, strandB, geneDf, padding=pad)
            if pad == paddings[-1] and svAnnoDf.loc[index, 'upGene_'  + suffix + '_A'] == 'ND':
                upVal = searchUpstreamGenes(chrA, posA, geneList, strandA, geneDf, startPad=pad, partner='A')
                svAnnoDf.loc[index, 'upGene_Nearest_A'] = upVal[0]
                svAnnoDf.loc[index, 'upGenePos_Nearest_A'] = upVal[1]
            elif pad == paddings[-1] and svAnnoDf.loc[index, 'upGene_'  + suffix + '_A'] != 'ND':
                svAnnoDf.loc[index, 'upGene_Nearest_A'] = svAnnoDf.loc[index, 'upGene_'  + suffix + '_A']
                svAnnoDf.loc[index, 'upGenePos_Nearest_A'] = suffix
                
            if pad == paddings[-1] and svAnnoDf.loc[index, 'upGene_'  + suffix + '_B'] == 'ND':
                upVal = searchUpstreamGenes(chrB, posB, geneList, strandB, geneDf, startPad=pad, partner='B')
                svAnnoDf.loc[index, 'upGene_Nearest_B'] = upVal[0]
                svAnnoDf.loc[index, 'upGenePos_Nearest_B'] = upVal[1]
            elif pad == paddings[-1] and svAnnoDf.loc[index, 'upGene_'  + suffix + '_B'] != 'ND':
                svAnnoDf.loc[index, 'upGene_Nearest_B'] = svAnnoDf.loc[index, 'upGene_'  + suffix + '_B']
                svAnnoDf.loc[index, 'upGenePos_Nearest_B'] = suffix
            
    return svAnnoDf            

                   
def annotateSVListWithNeighbourGenes(svAnnoDf, geneDf):
    """
    Annotate SVs with neighboring upstream genes
    """
    #svUpDf = svAnnoDf.copy()
    #print(svAnnoDf.shape, geneDf.shape)
    for index, row in svAnnoDf.iterrows():
        chrA = row['ChrA_x']
        chrB = row['ChrB_x']
        posA = row['PosA_x']
        posB = row['PosB_x']
        print(index, chrA, posA, chrB, posB)
        strandA = row['SVIndex'].split('.')[2]
        strandB = row['SVIndex'].split('.')[5]
        geneA = row['Gene_default_A']
        geneB = row['Gene_default_B']
        geneList = [geneA, geneB]
        #print(geneList)
        upValA = searchUpstreamGenes(chrA, posA, geneList, strandA, geneDf, startPad=0, partner='A', location='upstream')
        svAnnoDf.loc[index, 'upGene_Nearest_A'] = upValA[0]
        svAnnoDf.loc[index, 'upGenePos_Nearest_A'] = upValA[1]
        dnValA = searchUpstreamGenes(chrA, posA, geneList, strandA, geneDf, startPad=0, partner='A', location='downstream')
        svAnnoDf.loc[index, 'downGene_Nearest_A'] = dnValA[0]
        svAnnoDf.loc[index, 'downGenePos_Nearest_A'] = dnValA[1]
        upValB = searchUpstreamGenes(chrB, posB, geneList, strandB, geneDf, startPad=0, partner='B', location='upstream')
        svAnnoDf.loc[index, 'upGene_Nearest_B'] = upValB[0]
        svAnnoDf.loc[index, 'upGenePos_Nearest_B'] = upValB[1]
        dnValB = searchUpstreamGenes(chrB, posB, geneList, strandB, geneDf, startPad=0, partner='B', location='downstream')
        svAnnoDf.loc[index, 'downGene_Nearest_B'] = dnValB[0]
        svAnnoDf.loc[index, 'downGenePos_Nearest_B'] = dnValB[1]                    
    return svAnnoDf            

        
def annotateSVListWithIsoform(svDf, geneDf, default=True):
    """
    Annotate SVs with default isoforms
    """
    genes = geneDf.reset_index()
    genes.set_index(['RefSeq', 'Start', 'End'], inplace=True)
    if default == True:
        defDict = getDefaultIsoforms()
        suffix = 'default'
    else:
        defDict = getLargestIsoforms(geneDf)
        suffix = 'long'
    svDefAnnot = svDf.copy()
    for index, row in svDefAnnot.iterrows():
        geneA = getGeneName(row['ChrA'], row['PosA'], geneDf)
        geneB = getGeneName(row['ChrB'], row['PosB'], geneDf)
        if geneA not in defDict.keys():
            defDict[geneA] = 'NA'
        if geneB not in defDict.keys():
            defDict[geneB] = 'NA'
        isoA = defDict[geneA]
        isoB = defDict[geneB]
        svDefAnnot.loc[index, 'Gene_' + suffix + '_A'] = geneA
        svDefAnnot.loc[index, 'RefSeq_'  + suffix + '_A'] = isoA
        annotDfA = annotatePosWithIsoform(row['ChrA'], row['PosA'], isoA, genes)
        svDefAnnot.loc[index, 'ExonIntron_'  + suffix + '_A'] = annotDfA.loc[0, 'ExonIntron']
        svDefAnnot.loc[index, 'Strand_'  + suffix + '_A'] = annotDfA.loc[0, 'Strand']
        svDefAnnot.loc[index, 'Gene_' + suffix + '_B'] = geneB
        svDefAnnot.loc[index, 'RefSeq_' + suffix + '_B'] = isoB
        annotDfB = annotatePosWithIsoform(row['ChrB'], row['PosB'], isoB, genes)
        svDefAnnot.loc[index, 'ExonIntron_' + suffix + '_B'] = annotDfB.loc[0, 'ExonIntron']
        svDefAnnot.loc[index, 'Strand_' + suffix + '_B'] = annotDfB.loc[0, 'Strand']
    svDefAnnot['SVIndex'] = svDefAnnot['ChrA'] + '.' + svDefAnnot['PosA'].astype('str') + '.' + svDefAnnot['StrandA'] + '.' + svDefAnnot['ChrB'] + '.' + svDefAnnot['PosB'].astype('str') + '.' + svDefAnnot['StrandB']
    #svDefAnnot.drop(columns=['ChrA', 'PosA', 'StrandA', 'ChrB', 'PosB', 'StrandB'], axis=1, inplace=True)
    return svDefAnnot
    
        
def annotateSVs(svFile, varColumns, genome, outFile):
    """
    Input configuration and SV annotation
    """
    #geneFile = "/research_jude/rgs01_jude/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/svs/Gene_regions.txt"
    
    if genome == 'hg19':
        geneFile = "/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/xmagrp/CohortTracker/data/Gene_regions.txt"
    
    if genome == 'hg38':
        geneFile = "/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/xmagrp/CohortTracker/data/GeneRegion_h38.txt"
    
    geneDf = pd.read_csv(geneFile, sep = '\t', names=["Chr","Gene","RefSeq","Strand","ExonIntron","ExonIntronSlim","Start","End"])
    #print(len(varColumns))
    #sys.exit()
    geneDf.set_index(['Chr', 'Start', 'End'], inplace=True)
    svDf = checkSvFileExists(svFile, varColumns)
    inDf = getSVinDf(svFile, varColumns)
    inDf['SVIndex'] = ['.'.join(sv.split('.')[0:6]) for sv in inDf['curated_SV']]
    
    if outFile is None:
        outDir = os.path.dirname(svFile)
        outFile = os.path.join(outDir, os.path.basename(svFile) + '_svAnnot.tsv')
    
    paddings = [i*1e6 for i in [0.2, 0.5, 1]]
    
    #print(geneDf.head())
    #print(paddings)    
    annotDefault = annotateSVListWithIsoform(svDf, geneDf, default=True)
    print("default isoform done")
    annotLarge = annotateSVListWithIsoform(svDf, geneDf, default=False)
    print("large isoform done")
    joinDf = pd.merge(annotDefault, annotLarge, on='SVIndex', how='left', sort=True).drop_duplicates(ignore_index=True)
    #svDefAnnot.drop(columns=['ChrA', 'PosA', 'StrandA', 'ChrB', 'PosB', 'StrandB'], axis=1, inplace=True)
    #annotNgenes = annotateSVListWithUpstreamGenes(joinDf, geneDf, paddings)
    annotNgenes = annotateSVListWithNeighbourGenes(joinDf, geneDf)
    cols = [i for i in annotNgenes.columns if i != 'SVIndex' and '_x' not in i and '_y' not in i]
    cols.insert(0, 'SVIndex')
    joinOutDf = pd.merge(inDf, annotNgenes[cols], on='SVIndex', how='left').drop_duplicates()
    #annotNgenes[cols].to_csv(outFile, sep='\t', index=False)
    joinOutDf.to_csv(outFile, sep='\t', index=False)
    print("Output saved in", outFile)


def getCancerGenes(geneString, cancerGeneSet):
    geneSet = set(geneString.split(", "))
    return ', '.join(cancerGeneSet.intersection(geneSet))


def annotateCNVs(cnvFile, varColumns, genome, outFile):
    """
    Input configuration and CNV annotation
    """
    #geneFile = "/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/xmagrp/CohortTracker/data/Gene_regions.txt"
    #geneDf = pd.read_csv(geneFile, sep = '\t', names=["Chr","Gene","RefSeq","Strand","ExonIntron","ExonIntronSlim","Start","End"])
    #geneRegions = geneDf[['Chr', 'Start', 'End', 'Gene']].groupby('Gene').aggregate({'Chr':'first','Start':'min', 'End':'max'}).reset_index().sort_values(by = ['Chr', 'Start', 'End'], ignore_index=True)
    print(genome)
    if genome == 'hg19':
        geneFile = "/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/xmagrp/CohortTracker/data/ncbiRefSeqCurated.txt"
        
    if genome == 'hg38':
        geneFile = "/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/xmagrp/CohortTracker/data/ncbiRefSeqCurated_hg38.txt"
        
    # anfn = "/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/xmagrp/CohortTracker/data/Scott_Master_Medal_Gene_List.xlsx"
    anfn = "/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/xmagrp/CohortTracker/data/Scott_Master_Medal_Gene_List.tsv"
    # andf = pd.read_excel(anfn, usecols=["Gene"])
    andf = pd.read_csv(anfn, usecols=["Gene"])
    cancerGenes = set(andf['Gene'])

    
    geneDf = pd.read_csv(geneFile, sep = '\t', usecols=['chrom', 'txStart', 'txEnd', 'name2']).drop_duplicates()
    geneDf.rename(columns={'chrom':'Chr', 'txStart':'Start', 'txEnd':'End', 'name2':'Gene'}, inplace=True)
    geneDf = geneDf.loc[~geneDf['Chr'].str.contains("_")].reset_index(drop=True)
    geneRegions = geneDf[['Chr', 'Start', 'End', 'Gene']].groupby('Gene').aggregate({'Chr':'first','Start':'min', 'End':'max'}).reset_index().sort_values(by = ['Chr', 'Start', 'End'], ignore_index=True)
    sizeCutoff = 5e6
    cnvDf = checkCnvFileExists(cnvFile, varColumns)
    inDf = pd.read_csv(cnvFile, sep='\t')
    for index, row in cnvDf.iterrows():
        chrom = row[varColumns[0]]
        start = row[varColumns[1]]
        end = row[varColumns[2]]
        size = abs(end - start)
        if size > sizeCutoff:
            cnvDf.loc[index, ['GenesA', 'GenesWithin','GenesB']] = annotateCnvRegion(chrom, start, end, geneRegions, maxGenes=4, size='large')
        else:
            cnvDf.loc[index, ['GenesA', 'GenesWithin','GenesB']] = annotateCnvRegion(chrom, start, end, geneRegions, maxGenes=4, size='small')
    
    cnvDf['CancerGenesA'] = [getCancerGenes(cnvgenes, cancerGenes) for cnvgenes in cnvDf['GenesA']]
    cnvDf['CancerGenesWithin'] = [getCancerGenes(cnvgenes, cancerGenes) for cnvgenes in cnvDf['GenesWithin']]
    cnvDf['CancerGenesB'] = [getCancerGenes(cnvgenes, cancerGenes) for cnvgenes in cnvDf['GenesB']]

    joined = pd.merge(inDf, cnvDf, on=varColumns, how='left').drop_duplicates()
    
    if outFile is None:
        outDir = os.path.dirname(cnvFile)
        outFile = os.path.join(outDir, os.path.basename(cnvFile) + '_cnvAnnot.tsv')
        
    joined.to_csv(outFile, sep='\t', index=False)
    print("Output saved in", outFile)


def annotateVars(varFile, varType, varColumns, genome, outFile):
    """
    Input configuration for SV/CNV annotations
    """
    
    if genome is None:
        genome = 'hg19'
    
    if varType == 'SV':
        if varColumns is None:
            varColumns = ['curated_SV']
        annotateSVs(varFile, varColumns, genome, outFile)
    elif varType == 'CNV':
        if varColumns is None:
            varColumns = ['Chr', 'loc.start', 'loc.end']
        annotateCNVs(varFile, varColumns[0:3], genome, outFile)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Annotation of SVs/CNVs")
    parser.add_argument("-varFile",
                    help="File with SVs/CNVs, one marker per line. ",
                    required=True)

    parser.add_argument("-varType",
                    help="Type of variant column(s): SV or CNV. ",
                    choices=['SV', 'CNV'],
                    required=True)

    parser.add_argument("-varColumns",
                    help=textwrap.dedent('''\
                        Name(s) of the variant column(s).
                        
                        #SV: Single column. Default: curated_SV. Format: chrA.PosA.strandA.chrB.PosB.strandB. OR
                        
                        #SV: Six columns. Format: chrA, PosA, strandA, chrB, PosB, strandB 
                        
                        #CNV: Three columns. Default: Chr, loc.start, loc.end. 
                        
                        #Note: chromosome values should have 'chr' prefix. 
                        
                        Only chr1 to chr22 and chrX and chrY will be processed.'''),
                    nargs='+',
                    required=False)
    
    parser.add_argument("-genome",
                    help="Genome version. hg19 or hg38. Default: hg19. ",
                    choices=['hg19', 'hg38'],
                    required=False)

    
    parser.add_argument("-outFile",
                    help="Output file name. "
                        "Default: File_varTypeAnnot.tsv",
                    required=False)
    
    args = parser.parse_args()
    varFile = args.varFile
    varType = args.varType
    varColumns = args.varColumns
    genome = args.genome
    outFile = args.outFile
    annotateVars(varFile, varType, varColumns, genome, outFile)
