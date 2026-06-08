import sys
import os
import datetime
import math
import re
from pyfaidx import Fasta
from bx.intervals.intersection import Intersecter, Interval

# this corresponds to indelannotation_v5.py under /research_jude/rgs01_jude/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/mRNA/ZL/demo/GATA1/indel_RNA
PATH = "/research_jude/rgs01_jude/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/mRNA/ZL/reference/"
refGenome_hg38 = PATH+"GRCh38_no_alt_chr.fa"
refGenome_hg19 = PATH+"GRCh37-lite-chr.fa"
geneDef_hg38 = PATH+"ncbiRefSeqCurated_hg38_default.txt"
geneDef_hg19 = PATH+"ncbiRefSeqCurated_hg19_default.txt"
SAMTOOLS = PATH+"samtools"

def revcmp_str (input_str):
    rcs = input_str[::-1]
    rcs = rcs.replace('A','X')
    rcs = rcs.replace('T','A')
    rcs = rcs.replace('X','T')
    rcs = rcs.replace('C','X')
    rcs = rcs.replace('G','C')
    rcs = rcs.replace('X','G')
    return rcs

def cmp_str (input_str):
    rcs = input_str
    rcs = rcs.replace('A','X')
    rcs = rcs.replace('T','A')
    rcs = rcs.replace('X','T')
    rcs = rcs.replace('C','X')
    rcs = rcs.replace('G','C')
    rcs = rcs.replace('X','G')
    return rcs

def get_acceptor_perc():
    acceptor_A=[111695, 100818, 91352, 81706, 75624, 69584, 64955, 58813, 53264, 48575, 50186, 56540, 60683, 65579, 49327, 53104, 137682, 33840, 562158, 0, 148617, 141450, 148802]
    acceptor_C=[138393, 142458, 146355, 146903, 146531, 152898, 150751, 148776, 150871, 141003, 151636, 160283, 175855, 180529, 183429, 156709, 146230, 356964, 0, 0, 80084, 105384, 128380]
    acceptor_G=[82215, 81296, 77751, 75478, 75078, 70635, 66407, 63238, 58487, 55987, 58834, 62651, 57034, 50278, 34624, 34457, 113304, 1461,  0, 562158, 268878, 107114, 130244]
    acceptor_T=[229855, 237586, 246700, 258071, 264925, 269041, 280045, 291331, 299536, 316593, 301502, 282684, 268586, 265772, 294778, 317888, 164942, 169893, 0, 0, 64579, 208210, 154732]
    acceptor_perc={
        "A":[],
        "C":[],
        "G":[],
        "T":[]
    }
    for i in range(len(acceptor_A)):
        total_cov=sum([acceptor_A[i],acceptor_C[i],acceptor_G[i],acceptor_T[i]])
        acceptor_perc["A"].append(float(acceptor_A[i]+1)/total_cov)
        acceptor_perc["C"].append(float(acceptor_C[i]+1)/total_cov)
        acceptor_perc["G"].append(float(acceptor_G[i]+1)/total_cov)
        acceptor_perc["T"].append(float(acceptor_T[i]+1)/total_cov)
    return acceptor_perc

def get_donor_perc():
    donor_A=[188557, 359679, 56144, 0, 0, 348773, 389990, 51154, 100836, 170028, 130168, 126921, 129363, 127164, 124677]
    donor_C=[199093, 58923, 15716, 0, 0, 14288, 40181, 31189, 81071, 103565, 134547, 142430, 129249, 132394, 129537]
    donor_G=[102467, 61614, 447479, 557358, 0, 178292, 62563, 430785, 103629, 158599, 127660, 130854, 136803, 138701, 137035]
    donor_T=[67241, 77142, 38019, 0, 557358, 16005, 64624, 44230, 271822, 125166, 164983, 157153, 161943, 159099, 166109,]
    donor_perc={
        "A":[],
        "C":[],
        "G":[],
        "T":[]
    }
    for i in range(len(donor_A)):
        total_cov=sum([donor_A[i],donor_C[i],donor_G[i],donor_T[i]])
        donor_perc["A"].append(float(donor_A[i]+1)/total_cov)
        donor_perc["C"].append(float(donor_C[i]+1)/total_cov)
        donor_perc["G"].append(float(donor_G[i]+1)/total_cov)
        donor_perc["T"].append(float(donor_T[i]+1)/total_cov)
    return donor_perc

def get_acceptor_score(seq):
    seq=seq.upper()
    score=0
    for j in range(len(seq)):
        nn=seq[j]
        if nn.upper() in ['A', 'C', 'G', 'T']: # modified 12/06/2024
            bg_frac=bg_perc[nn]
            fg_frac=acceptor_perc[nn][j]
            score+=math.log(fg_frac/bg_frac)
            score=round(score,1)
        # else: print ('Problematic acceptor score sequence: '+seq)
    return score

#def get_acceptor_score(seq):
#    seq=seq.upper()
#    score=0
#    has_n = False
#    for j in range(len(seq)):
#        nn=seq[j]
#        if nn.upper() in ['A', 'C', 'G', 'T']: # modified 12/06/2024
#            bg_frac=bg_perc[nn]
#            fg_frac=acceptor_perc[nn][j]
#            score+=math.log(fg_frac/bg_frac)
#            score=round(score,1)
#        else: has_n = True
#    return score, has_n

def get_donor_score(seq):
    seq=seq.upper()
    score=0
    for j in range(len(seq)):
        nn=seq[j]
        if nn.upper() in ['A', 'C', 'G', 'T']: # modified 12/06/2024
            bg_frac=bg_perc[nn]
            fg_frac=donor_perc[nn][j]
            score+=math.log(fg_frac/bg_frac)
            score=round(score,1)
        # else: print ('Problematic donor score sequence: '+seq)
    return score   

def getRefDefSNV4(build):
#1 bin 2 nmID 3 chrom 4 strand 5 txS 6 txE 7 cdsS 8 cdsE 9 exonN 10 exonS 11 exonE 12 score
#13 gName 14 cdsStartStat 15 cdsEndStat 16 exonFrames 17 isDefault 18 inAnother 19 gNameInAnnother
    if build == "hg19":
        refseq = geneDef_hg19
    elif build == "hg38":
        refseq = geneDef_hg38
    isodict = {}
    genedict = {}
    intersect_dict = {}
    # promoter_dict = {}
    # utr_dict = {}
    region_dict = {}
    with open(refseq,'r') as fh:
        fh.readline()
        for line in fh:
            cdsLen = 0
            new = line.strip().split('\t')
            isoform = new[1]
            old_chrom = 'chr'+new[2]
            chrom = old_chrom.replace('chrchr','chr')
            mychrom = chrom
            if mychrom not in intersect_dict:
                intersect_dict[mychrom] = Intersecter()
                region_dict[mychrom] = Intersecter()
            strand = new[3]
            txS = int(new[4])
            txE = int(new[5])
            cdsS = int(new[6])
            cdsE = int(new[7])
            nExon = int(new[8])
            exonS = new[9].split(',')[:-1] #sort to increasing
            exonE = new[10].split(',')[:-1] #sort to increasing
            genename = new[12]
            isdefault = new[16]
            anothergene = new[18]
            baselist = []
            if isdefault == 'Yes':
                canonical_splice_sites=[]
                canonical_acceptor_sites=[]
                canonical_donor_sites=[]
                exonSpos = [int(i)+1 for i in exonS] #sort to increasing
                exonEpos = [int(i) for i in exonE] #sort to increasing
                exons = zip(exonSpos, exonEpos)
                if strand == "+":
                    exons=sorted(exons,key=lambda x:x[0])
                    for i in range(len(exons)):
                        acceptor="ND"
                        donor="ND"
                        if i==0:
                            donor=mychrom+"."+str(exons[i][1]+1)
                            canonical_donor_sites+=[mychrom+"."+str(exons[i][1]+1),mychrom+"."+str(exons[i][1]+2)] ## YF updated 9/11/2024
                        elif i==len(exons)-1:
                            acceptor=mychrom+"."+str(exons[i][0]-2)
                            canonical_acceptor_sites+=[mychrom+"."+str(exons[i][0]-2),mychrom+"."+str(exons[i][0]-1)]
                        else:
                            donor=mychrom+"."+str(exons[i][1]+1)
                            acceptor=mychrom+"."+str(exons[i][0]-2)
                            canonical_acceptor_sites+=[mychrom+"."+str(exons[i][0]-2),mychrom+"."+str(exons[i][0]-1)]
                            canonical_donor_sites+=[mychrom+"."+str(exons[i][1]+1),mychrom+"."+str(exons[i][1]+2)]
                        canonical_splice_sites.append([acceptor,donor])
                else:
                    exons=sorted(exons,key=lambda x:x[0],reverse=True)
                    for i in range(len(exons)):
                        acceptor="ND"
                        donor="ND"
                        if i==0:
                            donor=mychrom+"."+str(exons[i][0]-1)
                            canonical_donor_sites+=[mychrom+"."+str(exons[i][0]-2),mychrom+"."+str(exons[i][0]-1)] ## YF updated 9/11/2024
                        elif i==len(exons)-1:
                            acceptor=mychrom+"."+str(exons[i][1]+2)
                            canonical_acceptor_sites+=[mychrom+"."+str(exons[i][1]+1),mychrom+"."+str(exons[i][1]+2)]
                        else:
                            donor=mychrom+"."+str(exons[i][0]-1)
                            acceptor=mychrom+"."+str(exons[i][1]+2)
                            canonical_donor_sites+=[mychrom+"."+str(exons[i][0]-2),mychrom+"."+str(exons[i][0]-1)]
                            canonical_acceptor_sites+=[mychrom+"."+str(exons[i][1]+1),mychrom+"."+str(exons[i][1]+2)]
                        canonical_splice_sites.append([acceptor,donor])
                baselist = [mychrom,strand,cdsS,cdsE,txS,txE,exonS,exonE,strand,isoform,genename,nExon,canonical_acceptor_sites,canonical_donor_sites,canonical_splice_sites]
                isodict[isoform] = baselist
                genedict[genename] = isoform
                intersect_dict[mychrom].add_interval(Interval(txS, txE, value=isoform))
                if txS == cdsS and txE == cdsE:
                    region_dict[mychrom].add_interval(Interval(txS, txE, value=genename+':'+isoform+':nc_Exon'))
                else:
                    if txS != cdsS:
                        region_dict[mychrom].add_interval(Interval(txS, cdsS-1, value=genename+':'+isoform+':UTR'))
                    if txE != cdsE:
                        region_dict[mychrom].add_interval(Interval(cdsE+1, txE, value=genename+':'+isoform+':UTR'))
                if strand == "+":
                    region_dict[mychrom].add_interval(Interval(txS-2000, txS+1, value=genename+':'+isoform+':promoter'))
                else:
                    region_dict[mychrom].add_interval(Interval(txE-1, txE+2000, value=genename+':'+isoform+':promoter'))
                #need to add a hash where key is chr and positions to allow binary search for genes
    return isodict, genedict, intersect_dict, region_dict

def nCodonFinder(exonS, exonE, cdsS, cdsE):
    nCodon = set([])
    cdsPos = 0
    for start, end in zip(exonS, exonE):
        start = int(start)
        end = int(end)
        for i in range(start+1, end+1):
            if i > cdsS and i <= cdsE:
                cdsPos += 1
                nCodon.add(cdsPos)
    return max(nCodon) # returning the total number of nFrame

def getBedFafiles(build, genelist):
#baselist = [mychrom,strand,cdsS,cdsE,txS,txE,exonS,exonE,strand,isoform,gene,nExon]
    if build == "hg19":
        refGenome = Fasta(refGenome_hg19)
    elif build == "hg38":
        refGenome = Fasta(refGenome_hg38)
    mychrom = genelist[0]
    cdsS = genelist[2]
    cdsE = genelist[3]
    exonS = sorted([int(i) for i in genelist[6]])
    exonE = sorted([int(i) for i in genelist[7]])
    strand = genelist[1]
    isoform = genelist[9]
    gene = genelist[10]
    nExon = genelist[11]
    
    mRNA_seq = ""
    for start, end in zip(exonS, exonE):
        seqq = str(refGenome[mychrom][start:end])
        mRNA_seq += seqq

    dnaDict = {}
    rnaDict = {}
    if strand == "+":
        codonPos = 0
    else:
        nCodon = nCodonFinder(exonS, exonE, cdsS, cdsE)
        codonPos = nCodon
    rnaPos = 0
    exon_num = 1
    spliceDict = {}
    if nExon > 1:
        for start, end in zip(exonS, exonE):
            start = int(start)
            end = int(end)
            if exon_num == 1: #need to make sure >=2 exons
                if strand == "+":
                    for ssp in range(end-5, end+1):
                        splice_pos = end - ssp + 1
                        splice = "splice_D_E_"+str(splice_pos)
                        spliceDict[ssp] = splice
                    for ssp in range(end+1, end+21):
                        splice_pos = ssp - end
                        splice = "splice_D_I_"+str(splice_pos)
                        spliceDict[ssp] = splice
                elif strand == "-":
                    for ssp in range(end-5, end+1):
                        splice_pos = end - ssp + 1
                        splice = "splice_A_E_"+str(splice_pos)
                        spliceDict[ssp] = splice
                    for ssp in range(end+1, end+21):
                        splice_pos = ssp - end
                        splice = "splice_A_I_"+str(splice_pos)
                        spliceDict[ssp] = splice
            elif exon_num > 1 and exon_num < len(exonE):
                if strand == "+":
                    for ssp in range(end-5, end+1):
                        splice_pos = end - ssp + 1
                        splice = "splice_D_E_"+str(splice_pos)
                        spliceDict[ssp] = splice
                    for ssp in range(end+1, end+21):
                        splice_pos = ssp - end
                        splice = "splice_D_I_"+str(splice_pos)
                        spliceDict[ssp] = splice
                    for ssp in range(start-19, start+1):
                        splice_pos = start - ssp + 1
                        splice = "splice_A_I_"+str(splice_pos)
                        spliceDict[ssp] = splice
                    for ssp in range(start+1, start+6):
                        splice_pos = ssp - start
                        splice = "splice_A_E_"+str(splice_pos)
                        spliceDict[ssp] = splice
                elif strand == "-":
                    for ssp in range(end-5, end+1):
                        splice_pos = end - ssp + 1
                        splice = "splice_A_E_"+str(splice_pos)
                        spliceDict[ssp] = splice
                    for ssp in range(end+1, end+21):
                        splice_pos = ssp - end
                        splice = "splice_A_I_"+str(splice_pos)
                        spliceDict[ssp] = splice
                    for ssp in range(start-19, start+1):
                        splice_pos = start - ssp + 1
                        splice = "splice_D_I_"+str(splice_pos)
                        spliceDict[ssp] = splice
                    for ssp in range(start+1, start+6):
                        splice_pos = ssp - start
                        splice = "splice_D_E_"+str(splice_pos)
                        spliceDict[ssp] = splice
            elif exon_num == len(exonE):
                if strand == "+":
                    for ssp in range(start-19, start+1):
                        splice_pos = start - ssp + 1
                        splice = "splice_A_I_"+str(splice_pos)
                        spliceDict[ssp] = splice
                    for ssp in range(start+1, start+6):
                        splice_pos = ssp - start
                        splice = "splice_A_E_"+str(splice_pos)
                        spliceDict[ssp] = splice
                elif strand == "-":
                    for ssp in range(start-19, start+1):
                        splice_pos = start - ssp + 1
                        splice = "splice_D_I_"+str(splice_pos)
                        spliceDict[ssp] = splice
                    for ssp in range(start+1, start+6):
                        splice_pos = ssp - start
                        splice = "splice_D_E_"+str(splice_pos)
                        spliceDict[ssp] = splice
            if strand == "+":
                for i in range(start+1, end+1):
                    rnaPos += 1
                    dnaPos = i
                    base = mRNA_seq[rnaPos-1]
                    if i > cdsS and i <= cdsE:
                        codonPos += 1
                        cdsPos = codonPos
                        frame = (cdsPos - 1) % 3
                    else:
                        cdsPos = -1
                        frame = -1
                    output = [cdsPos, frame, base, cmp_str(base), rnaPos, dnaPos, codonPos]
                    dnaDict[dnaPos] = output
                    rnaDict[rnaPos] = dnaPos
            else:
                for i in range(start+1, end+1):
                    rnaPos += 1
                    dnaPos = i
                    base = mRNA_seq[rnaPos-1]
                    if i > cdsS and i <= cdsE:
                        codonPos -= 1
                        cdsPos = codonPos
                        frame = (cdsPos - 1) % 3
                    else:
                        cdsPos = -1
                        frame = -1
                    output = [cdsPos, frame, base, cmp_str(base), rnaPos, dnaPos, codonPos]
                    dnaDict[dnaPos] = output
                    rnaDict[rnaPos] = dnaPos
            exon_num += 1
    else:
        if strand == "+":
            for i in range(start+1, end+1):
                rnaPos += 1
                dnaPos = i
                base = mRNA_seq[rnaPos-1]
                if i > cdsS and i <= cdsE:
                    codonPos += 1
                    cdsPos = codonPos
                    frame = (cdsPos - 1) % 3
                else:
                    cdsPos = -1
                    frame = -1
                output = [cdsPos, frame, base, cmp_str(base), rnaPos, dnaPos, codonPos]
                dnaDict[dnaPos] = output
                rnaDict[rnaPos] = dnaPos
        else:
            for i in range(start+1, end+1):
                rnaPos += 1
                dnaPos = i
                base = mRNA_seq[rnaPos-1]
                if i > cdsS and i <= cdsE:
                    codonPos -= 1
                    cdsPos = codonPos
                    frame = (cdsPos - 1) % 3
                else:
                    cdsPos = -1
                    frame = -1
                output = [cdsPos, frame, base, cmp_str(base), rnaPos, dnaPos, codonPos]
                dnaDict[dnaPos] = output
                rnaDict[rnaPos] = dnaPos
    return dnaDict, rnaDict, spliceDict

def translate(seq):      
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein =""
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if len(codon) < 3:continue
        protein+= table[codon]
    return protein

def featureScan(gmst, gmsp, txS, txE, cdsS, cdsE, dnaDict, spliceDict, strand):
    fset = set(["e_Invalid"])
    txS = int(txS)
    txE = int(txE)
    for p in range(gmst, gmsp+1):
        if p <= txS or p > txE:
            fset.add("e_Invalid")
        else:
            if p not in dnaDict:
                fset.add("d_Intron")
                if p in spliceDict:
                    fset.add("b_Splice")
            else:
                if p in spliceDict:
                    fset.add("b_Splice")
                if strand == "+":
                    if p >= int(cdsS+1) and p<= int(cdsE):
                        fset.add("a_CDS")
                    elif p < int(cdsS+1):
                        fset.add("c_UTR5") # utr region
                    elif p > int(cdsE):
                        fset.add("c_UTR3")
                elif strand == "-":
                    if p >= int(cdsS+1) and p<= int(cdsE):
                        fset.add("a_CDS")
                    elif p < int(cdsS+1):
                        fset.add("c_UTR3") # utr region
                    elif p > int(cdsE):
                        fset.add("c_UTR5")
    return fset

def spliceScan(gmst, gmsp, spliceDict):
    #output = [cdsPos, frame, base, cmp_str(base), rnaPos, dnaPos, codonPos]
    spdict = {}
    ret = "BugInSpliceScan"
    for p in range(gmst, gmsp+1):
        if p in spliceDict:
            # feature: splice_A_E_20
            feature = spliceDict[p]
            feature_region = '_'.join(feature.split('_')[1:3])
            feature_number = int(feature.split('_')[3])
            if feature_region not in spdict:
                spdict[feature_region] = set([])
            spdict[feature_region].add(feature_number)

    if len(spdict) == 0:
        ret = ""
    elif len(spdict) == 1:
        feature_region = list(spdict.keys())[0]
        min_number = min(spdict[feature_region])
        if min_number <= 2 and "I" in feature_region:
            ret = "splicesite_"+feature_region+str(min_number)
        else:
            ret = "spliceregion_"+feature_region+str(min_number)
    else:
        if ("splice_D_I" in set(spdict) and "splice_D_E" in set(spdict)): 
            ret = "splicespan_D"
        elif ("splice_A_I" in set(spdict) and "splice_A_E" in set(spdict)):
            ret = "splicespan_A"
        else:
            ret = "splicespan"
    return ret

def sequenceGenerator(gmst, gmsp, refbase, altbase, dnaDict, strand):
    realcdsindel = "" 
    mutcdsindel = ""
    refbase_check = ""
    if len(refbase) > 0:
        for p in range(gmst, gmsp+1):
            index = p-gmst
            refbase_check += refbase[index]
            if p not in dnaDict:continue
            cdsPos = dnaDict[p][0]
            if cdsPos == -1:continue # this position is not in cdsPos
            realcdsindel += refbase[index]
            if index < len(altbase):
                mutcdsindel += altbase[index]
        if len(altbase) > len(refbase):
            for i in range(len(refbase), len(altbase)):
                mutcdsindel += altbase[i]
        else:
            mutcdsindel = mutcdsindel
    else:
        realcdsindel = ""
        mutcdsindel = altbase
    pset = set(list(range(gmst, gmsp+1)))
    leftseq = ""
    midseq = ""
    rightseq = ""
    for p in sorted(dnaDict): #[cdsPos, frame, base, cmp_str(base), rnaPos, dnaPos, codonPos]
        if dnaDict[p][0] == -1:continue
        if p < gmst:
            leftseq += dnaDict[p][2]
        elif p >= gmst and p <= gmsp:
            midseq += dnaDict[p][2]
        elif p > gmsp:
            rightseq += dnaDict[p][2]
    if refbase_check != midseq:
        warning = "warning"
    else:
        warning = ""
    if realcdsindel == "":
        # for insertion
        realcdsseq = leftseq+realcdsindel+rightseq
        mutcdsseq = leftseq+mutcdsindel+rightseq
        diff = abs(len(mutcdsindel)-len(realcdsindel))
    else:
        realcdsseq = leftseq+midseq+rightseq
        mutcdsseq = leftseq+mutcdsindel+rightseq
        diff = abs(len(mutcdsindel)-len(midseq))
    if strand == "+":
        realaa = translate(realcdsseq)
        mutaa = translate(mutcdsseq)
    elif strand == "-":
        realaa = translate(revcmp_str(realcdsseq))
        mutaa = translate(revcmp_str(mutcdsseq))
    return realcdsseq, mutcdsseq, realaa, mutaa, diff, warning

def annotation(refaa, mutaa, diff, dnaDict, mutst, mutsp, strand):
    if refaa[0] != mutaa[0]:
        muttype = "startloss"
        change = "-1"
        m2 = 0
        for x in range(1,len(refaa)):
            refcodon = refaa[x]
            if refcodon == "M":
                m2 = x+1
                break
        change = "M.fs_{}".format(str(m2))
    else:
        if diff % 3 != 0:
            muttype, change = frameshiftCodonChange(refaa, mutaa)
        else:
            if diff == 0:
                codonpos = "bad"
                fullcodon = len(refaa)-1 # 0-based index
                for pos in range(mutst, mutsp+1):
                    if pos in dnaDict:
                        if dnaDict[pos][0] != -1:
                            # codonpos = int(dnaDict[pos][0])//3 # the first codon is 0
                            codonpos = int(dnaDict[pos][0])//3 + 1 if strand == "-" else int(dnaDict[pos][0])//3 # modified 09192024
                            break
                # if strand == "+":                           # modified 09192024
                #     codonpos = codonpos                     # modified 09192024
                # else:                                       # modified 09192024
                #     codonpos = fullcodon-codonpos           # modified 09192024
                muttype, change = mnvCodonChange(refaa, mutaa, codonpos)
            else:
                muttype, change = inframeCodonChange(refaa, mutaa)
    return muttype, change

def frameshiftCodonChange(refaa, mutaa):
    mutpos = 0
    for ref, mut in zip(refaa, mutaa):
        if ref == mut:
            mutpos += 1
        else:
            mutpos = mutpos+1
            break

    stopmutpos = 0
    for m in range(mutpos-1, len(mutaa)):
        if mutaa[m] != "*":
            stopmutpos += 1
        elif mutaa[m] == "*":
            stopmutpos = stopmutpos+1
            break

    muttype = "frameshift"
    beforemutpos = mutpos-1
    aachange = str(beforemutpos)+'.'+refaa[beforemutpos-1]+".fs_"+str(stopmutpos)
    return muttype, aachange

def mnvCodonChange(refaa, mutaa, codonpos):
    mutpos = 0
    aachange = "-1"
    muttype = "BugInmnvCodonChange"
    if refaa == mutaa:
        muttype = "silent"
        if codonpos == "bad":
            aachange = "badannotation"
        else:
            refAa = refaa[codonpos-1] #changed
            aachange = refAa+str(codonpos)+refAa
    else:
        mutposlst = []
        for ref, mut in zip(refaa, mutaa):
            mutpos += 1 # this is the start of mutation
            if ref != mut:
                # print (str(mutpos),ref, mut)
                mutposlst.append(mutpos)
        mutposst = min(mutposlst)
        mutpossp = max(mutposlst)
        refAa = ""
        altAa = ""
        if mutposst == mutpossp:
            refAa = refaa[mutposst-1]
            altAa = mutaa[mutposst-1]
        else:
            for mutpos in range(mutposst, mutpossp+1):
                refAa += refaa[mutpos-1]
                altAa += mutaa[mutpos-1]
        aachange = refAa+str(mutposst)+altAa
        if '*' in altAa:
            muttype = 'nonsense'
        else:
            if refAa == "*":
                muttype = "stoploss"
                aachange = "*.fs_nostopfound" # TO-DO/needs to pull mrna sequence
                for x in range(mutpos-1, len(altAa)):
                    if altAa[x] == "*":
                        aachange = "*.fs_"+str(x+1)
            else:
                muttype = "missense"
    return muttype, aachange

def inframeCodonChange(refaa, mutaa):
    leftsq = 0
    point = 0
    reflength = len(refaa)
    mutlength = len(mutaa)
    
    # squeezing from the left
    for ref, mut in zip(refaa, mutaa):
        if ref == mut:
            leftsq += 1
        else:
            leftsq += 1
            break

    refup = leftsq
    mutup = leftsq

    # squeezing from the right
    refr = 0
    mutr = 0
    refdown = 0
    mutdown = 0
    for ref, mut in zip(refaa[::-1], mutaa[::-1]):
        if ref == mut:
            point += 1
            refdown = len(refaa)-point+1
            mutdown = len(mutaa)-point+1
            if refdown <= refup or mutdown <= mutup:
                break
        else:
            break

    refaa_mid = ''.join(refaa[refup-2:refdown])
    mutaa_mid = ''.join(mutaa[mutup-2:mutdown])
    if len(refaa_mid) > len(mutaa_mid):
        muttype = "proteinDel"
        aachange = str(refup-1)+'.'+refaa_mid+'>'+mutaa_mid
        if "*" in mutaa_mid:
            muttype = "proteinDel_stopgain"
        return muttype, aachange
    else:
        muttype = "proteinIns"
        aachange = str(refup-1)+'.'+refaa_mid+'>'+mutaa_mid
        if "*" in mutaa_mid:
            muttype = "proteinIns_stopgain"
        return muttype, aachange

def reclassifyaachange(muttype, aachange):
    if "frameshift" in muttype:
        aa = aachange.split('.')
        aachange = aa[0]+'.'+aa[1]+'.unknown'
    else:
        if aachange != "-1":
            aachange = aachange
            # newaachange = ""
            # aapos = 0
            # for aa in aachange:
            #     if aa.isalpha():
            #         if aapos == 0:
            #             newaachange += aa
            #         aapos += 1
            #     elif aa.isdigit() and aapos == 1:
            #         newaachange += aa
            # aachange = newaachange
    return aachange

def preannotate(gene, isoform, genedict, isodict, build):
    baselist = isodict[isoform]
    realisoform = baselist[-3]
    realgene = baselist[-2]
    strand = baselist[1]
    cdsS = baselist[2]
    cdsE = baselist[3]
    txS = baselist[4]
    txE = baselist[5]
    exonS = baselist[6]
    exonE = baselist[7]
    dnaDict, rnaDict, spliceDict = getBedFafiles(build, baselist)
    return dnaDict, rnaDict, spliceDict

def geneassign(snv4lst, isodict, isointerval):
    isosnv4dict = {}
    for snv4 in snv4lst:
        n = snv4.split('.')
        chrom = n[0]
        gnommutPos = int(n[1])
        refbase = n[2].replace('-','')
        altbase = n[3].replace('-','')
        reflen = len(refbase)
        if reflen == 0:
            gnommutSt = gnommutPos
            gnommutSp = gnommutPos
        else:
            gnommutSt = gnommutPos
            gnommutSp = gnommutSt + reflen - 1
        a = isointerval[chrom].find(gnommutSt, gnommutSp)
        if len(a) == 0:
            gene = "ND"
            isoform = "ND"
        else:
            b = a[0]
            myisoform = b.value
            baselist = isodict[myisoform]
            realgene = baselist[-2]
            if myisoform not in isosnv4dict:
                isosnv4dict[myisoform] = []
            isosnv4dict[myisoform].append(snv4)
    return isosnv4dict

def findClosestAA(pos, isodict, isoform, dnaDict): 
    baselist = isodict[isoform]
    strand = baselist[1]
    plist = list(dnaDict)
    cdspos = 0
    fullaaseq = ""
    #output = [cdsPos, frame, base, cmp_str(base), rnaPos, dnaPos, codonPos]
    for i in range(len(plist)):
        ipos = plist[i]
        if dnaDict[ipos][0] == -1:continue
        fullaaseq += dnaDict[ipos][2]
    if strand == '+':
        fullaa = translate(fullaaseq)
    else:
        fullaa = translate(revcmp_str(fullaaseq))
    ###
    for current in range(1, len(plist)):
        pre = current - 1
        currentpos = plist[current]
        prepos = plist[pre]
        if pos >= plist[pre] and pos < plist[current]:
            cdspos = dnaDict[prepos][-1]
    if strand == "-":
        aa_pos = math.floor(cdspos/3)
    else:
        aa_pos = math.ceil(cdspos/3)
    aa = fullaa[aa_pos-1]
    aa_comment = aa+str(aa_pos)
    return aa_comment

def annotateOneSNV4withgene(snv4, gene, isoform, isodict, build, dnaDict, rnaDict, spliceDict):
    n = snv4.split('.')
    chrom = n[0]
    gnommutPos = int(n[1])
    refbase = n[2].replace('-','')
    altbase = n[3].replace('-','')
    reflen = len(refbase)
    if reflen == 0:
        gnommutSt = gnommutPos
        gnommutSp = gnommutPos
    else:
        gnommutSt = gnommutPos
        gnommutSp = gnommutSt + reflen - 1
    baselist = isodict[isoform]
    realisoform = baselist[-3]
    realgene = baselist[-2]
    strand = baselist[1]
    cdsS = baselist[2]
    cdsE = baselist[3]
    txS = baselist[4]
    txE = baselist[5]
    exonS = baselist[6]
    exonE = baselist[7]
    feature = featureScan(gnommutSt, gnommutSp, txS, txE, cdsS, cdsE, dnaDict, spliceDict, strand)
    if "a_CDS" in feature or "b_Splice" in feature:
        comment = ""
        warning = ""
        splice_feature = spliceScan(gnommutSt, gnommutSp, spliceDict)
        oricdsseq, mutcdsseq, oriaa, mutaa, seqdiff, warning = sequenceGenerator(gnommutSt, gnommutSp, refbase, altbase, dnaDict, strand)
        if "a_CDS" in feature:
            muttype, aachange = annotation(list(oriaa), list(mutaa), seqdiff, dnaDict, gnommutSt, gnommutSp, strand)
        else:
            muttype = "-1"
            aachange = "-1"
        if len(splice_feature) == 0:
            muttype = muttype
        else:
            if muttype == "-1":
                muttype = splice_feature
            else:
                muttype = muttype+':'+splice_feature
        ffeature = set([])
        for iind in list(feature):
            if "a_CDS" in feature:
                ffeature.add("CDS")
            elif "b_Splice" in feature:
                ffeature.add("Splice")
        ff = "+".join(list(ffeature))
        if "splicespan" in muttype:
            mypos = aachange.split('.')[0]
            aachange = reclassifyaachange(muttype, aachange)
            interested_splice_site = determineSpliceSite(snv4, isoform, isodict)
            # if len(aachange) == 1:
            if interested_splice_site == "D":
                if len(aachange) != 1:
                    aachange = aachange.replace(".unknown",'.donor')
                else:
                    aachange = aachange+'.donor'
            elif interested_splice_site == "A":
                if len(aachange) != 1:
                    aachange = aachange.replace(".unknown",'.acceptor')
                else:
                    aachange = aachange+'.acceptor'
                # aachange = mypos+'.'+aachange+'.acceptor'
    else:
        ff = sorted(list(feature))[0].split('_')[1]
        muttype = "-1"
        aachange = "-1"
    return ff, muttype, aachange

def determineSpliceSite(snv4, isoform, isodict, intron_flanking=20, exon_flanking=3):
    acceptor = 0
    donor = 0
    ret = "ND"
    stpos = int(snv4.split('.')[1])
    sppos = stpos+len(snv4.split('.')[2])-1
    baselist = isodict[isoform]
    strand = baselist[1]
    canonical_splice_sites = baselist[14]
    canonical_splice_acceptors=sorted([int(i[0].split('.')[1]) for i in canonical_splice_sites if i[0] != "ND"])
    canonical_splice_donors=sorted([int(i[1].split('.')[1]) for i in canonical_splice_sites if i[1] != "ND"])
    if strand == "+":
        for x in canonical_splice_acceptors:
            if stpos >= x-intron_flanking and stpos <= x+exon_flanking:
                acceptor += 1
            if sppos >= x-intron_flanking and sppos <= x+exon_flanking:
                acceptor += 1
        for x in canonical_splice_donors:
            if stpos >= x-exon_flanking and stpos <= x+intron_flanking:
                donor += 1
            if sppos >= x-exon_flanking and sppos <= x+intron_flanking:
                donor += 1
    elif strand == "-":
        for x in canonical_splice_acceptors:
            if stpos >= x-exon_flanking and stpos <= x+intron_flanking:
                acceptor += 1
            if sppos >= x-exon_flanking and sppos <= x+intron_flanking:
                acceptor += 1
        for x in canonical_splice_donors:
            if stpos >= x-intron_flanking and stpos <= x+exon_flanking:
                donor += 1
            if sppos >= x-intron_flanking and sppos <= x+exon_flanking:
                donor += 1
    if acceptor > 0 and donor == 0:
        ret = "A"
    elif acceptor == 0 and donor > 0:
        ret = "D"
    return ret

def extract_seq(refGenome, chrom, start, stop):
    seq = str(refGenome[chrom][start:stop])
    return seq

def annotateHLPA(snv4, isoform, isodict, interested_splice_site, refGenome, flank=50, intron_flank=20):
    #baselist = [mychrom,strand,cdsS,cdsE,txS,txE,exonS,exonE,strand,isoform,genename,nExon,canonical_splice_sites]
    baselist = isodict[isoform]
    mychrom = baselist[0]
    cdsS = baselist[2]
    cdsE = baselist[3]
    txS = baselist[4]
    txE = baselist[5]
    strand = baselist[1]
    isoform = baselist[9]
    gene = baselist[10]
    nExon = baselist[11]
    essencial_splice_acceptors_pos = baselist[12]
    essencial_splice_donors_pos = baselist[13]
    canonical_splice_sites = baselist[14]
    canonical_splice_acceptors=[i[0] for i in canonical_splice_sites]
    canonical_splice_donors=[i[1] for i in canonical_splice_sites]
    ###
    chro,pos,ref_allele,mut_allele=snv4.split(".")
    ref_contig=extract_seq(refGenome, mychrom, int(pos)-flank-1, int(pos)+len(ref_allele)+flank) ## modified added "+len(ref_allele)" 08122024
    ref_contig_pos=[chro+"."+str(i) for i in range((int(pos)-flank),(int(pos)+len(ref_allele)+flank+1))] ## modified added "+len(ref_allele)" 08122024
    mut_contig=ref_contig[0:flank]+mut_allele+ref_contig[(flank+len(ref_allele)):]
    mut_contig_pos=ref_contig_pos[0:flank]+["mut."+str(i) for i in range(len(mut_allele))]+ref_contig_pos[(flank+len(ref_allele)):]
    if strand=="-":
        ref_contig=revcmp_str(ref_contig)
        ref_contig_pos=ref_contig_pos[::-1]
        mut_contig=revcmp_str(mut_contig)
        mut_contig_pos=mut_contig_pos[::-1]
    ###
    all_detected_acceptors=[i.start() for i in re.finditer("AG",mut_contig)]
    all_detected_donors=[i.start() for i in re.finditer("GT",mut_contig)]
    all_detected_acceptors_update=[]
    all_detected_donors_update=[]
    for i in all_detected_acceptors:
        if i<flank-intron_flank:
            continue
        if i>flank+len(mut_allele)+intron_flank:
            continue
        all_detected_acceptors_update.append(i)

    for i in all_detected_donors:
        if i<flank-intron_flank:
            continue
        if i>flank+len(mut_allele)+intron_flank:
            continue
        all_detected_donors_update.append(i)

    # print (f"Marker: {snv4}")
    # print (f"Isoform: {isoform}")
    # print (f"all_detected_acceptors: {all_detected_acceptors}")
    # print (f"all_detected_acceptors_update: {all_detected_acceptors_update}")
    # print (f'ref_contig_pos: {ref_contig_pos}')
    # print (f'mut_contig_pos: {mut_contig_pos}')
                    
    all_detected_acceptors=all_detected_acceptors_update
    all_detected_donors=all_detected_donors_update

    if interested_splice_site=="A":
        canonical_comment="L"
        neo_comment="A"
        # if chro+"."+pos in essencial_splice_acceptors_pos:
        #     canonical_comment="H"
        for jj in range(len(ref_allele)):
            if chro+"."+str(jj+int(pos)) in essencial_splice_acceptors_pos:
                canonical_comment="H"
        neo_acceptors=[]
        for acceptor_index in all_detected_acceptors:
            acceptor=mut_contig_pos[(acceptor_index)]
            acceptor_seq=mut_contig[(acceptor_index-18):(acceptor_index+5)]
            if acceptor in canonical_splice_acceptors:
                acceptor_type="canonical"
                canonical_acceptor_index=ref_contig_pos.index(acceptor)
                canonical_acceptor_seq=ref_contig[(canonical_acceptor_index-18):(canonical_acceptor_index+5)]
                if acceptor_seq !=canonical_acceptor_seq:
                    acceptor_seq_core=acceptor_seq[17:22]
                    canonical_acceptor_seq_core=canonical_acceptor_seq[17:22]
                    c0=canonical_acceptor_seq_core[0]+">"+acceptor_seq_core[0]
                    if c0 in core_acceptor["H"]["c0"]:
                        canonical_comment="H"
                    c1=canonical_acceptor_seq_core[1]+">"+acceptor_seq_core[1]
                    if c1 in core_acceptor["H"]["c1"]:
                        canonical_comment="H"
                    c2=canonical_acceptor_seq_core[2]+">"+acceptor_seq_core[2]
                    if c2 in core_acceptor["H"]["c2"]:
                        canonical_comment="H"
            else:
                 acceptor="neo"
                 acceptor_seq_score=get_acceptor_score(acceptor_seq)
                 # acceptor_seq_score, has_n=get_acceptor_score(acceptor_seq)
                 # if has_n: print (snv4 + '\tacceptor: ' + acceptor_seq) # modified
                 neo_acceptors.append(acceptor_seq_score)
        if len(neo_acceptors)>0 and max(neo_acceptors)>2:
            neo_comment="P"
    elif interested_splice_site=="D":
        neo_donors=[]
        canonical_comment="L"
        neo_comment="A"
        # if (chro+"."+pos in essencial_splice_donors_pos):
        #     canonical_comment="H"
        for jj in range(len(ref_allele)):
            if chro+"."+str(jj+int(pos)) in essencial_splice_donors_pos:
                canonical_comment = "H" 
        for donor_index in all_detected_donors:
            donor=mut_contig_pos[donor_index]
            donor_seq=mut_contig[(donor_index-3):(donor_index+12)]
            if donor in canonical_splice_donors:
                donor_type="canonical"
                canonical_donor_index=ref_contig_pos.index(donor)
                canonical_donor_seq=ref_contig[(canonical_donor_index-3):(canonical_donor_index+12)]
                if donor_seq !=canonical_donor_seq:
                    donor_seq_core=donor_seq[0:8]
                    canonical_donor_seq_core=canonical_donor_seq[0:8]
                    
                    # print (f"Marker: {snv4}")
                    # print (f"Isoform: {isoform}")
                    # print (f"all_detected_acceptors: {all_detected_acceptors}")
                    # print (f'canonical_donor_seq_core: {canonical_donor_seq_core}')
                    # print (f"donor_seq_core: {donor_seq_core}")
                    # print (f"mut_contig: {mut_contig}")
                    # print (f"ref_contig: {ref_contig}")
                    # print (ref_contig[0:flank] + '_' + mut_allele + '_' + ref_contig[(flank+len(ref_allele)):])
                    
                    c0=canonical_donor_seq_core[0]+">"+donor_seq_core[0]
                    if c0 in core_donor["H"]["c0"]:
                        canonical_comment="H"
                        continue
                    c1=canonical_donor_seq_core[1]+">"+donor_seq_core[1]
                    if c1 in core_donor["H"]["c1"]:
                        canonical_comment="H"
                        continue
                    c2=canonical_donor_seq_core[2]+">"+donor_seq_core[2]
                    if c2 in core_donor["H"]["c2"]:
                        canonical_comment="H"
                        continue
                    c3=canonical_donor_seq_core[3]+">"+donor_seq_core[3]
                    if c0 in core_donor["H"]["c3"]:
                        canonical_comment="H"
                        continue
                    c4=canonical_donor_seq_core[4]+">"+donor_seq_core[4]
                    if c4 in core_donor["H"]["c4"]:
                        canonical_comment="H"
                        continue
                    c5=canonical_donor_seq_core[5]+">"+donor_seq_core[5]
                    if c5 in core_donor["H"]["c5"]:
                        canonical_comment="H"
                        continue
                    c6=canonical_donor_seq_core[6]+">"+donor_seq_core[6]
                    if c6 in core_donor["H"]["c6"]:
                        canonical_comment="H"
                        continue
                    c7=canonical_donor_seq_core[7]+">"+donor_seq_core[7]
                    if c7 in core_donor["H"]["c7"]:
                        canonical_comment="H"
            else:
                 donor_type="neo"
                 donor_seq_score=get_donor_score(donor_seq)
                 neo_donors.append(donor_seq_score)
        if len(neo_donors)>0 and max(neo_donors)>2:
            neo_comment="P"
    else:
        canonical_comment = "ND"
        neo_comment = "ND"
    return canonical_comment, neo_comment

def getdefaultisoform(build):
    if build == "hg38":
        genedef = geneDef_hg38
    elif build == "hg19":
        genedef = geneDef_hg19
    default_dict = {}
    with open(genedef,'r') as fh:
        fh.readline()
        for line in fh:
            new = line.strip().split('\t')
            if new[-3] != 'Yes':continue
            default_dict[new[1]] = new[12]
    return default_dict

if __name__ == "__main__":

    build = sys.argv[1]
    if build == "hg19":
        refGenome = Fasta(refGenome_hg19)
    elif build == "hg38":
        refGenome = Fasta(refGenome_hg38)

    print ('\t'.join(['Sample', 'Curated.Marker', 'GeneName', 'RefSeq', 'Region', 'Class', 'AAChange']))
    bg_perc={"A":0.3,"C":0.2,"G":0.2,"T":0.3}
    acceptor_perc=get_acceptor_perc()
    donor_perc=get_donor_perc()
    core_donor={"seq":"CAGGTAAG",
    "H":{
        "c0":["C>T","A>T","G>T"],
        "c1":["A>T","A>G","A>C"],
        #"c2":["G>A","G>T","A>C"],
        "c2":["G>A","G>T","G>C"],
        "c3":["G>A","G>T","G>C"],
        "c4":["T>A","T>G","T>C"],
        "c5":["A>T","A>C","G>T","G>C"],
        "c6":["A>T","A>G","A>C"],
        "c7":["G>A","G>C","G>T"]
        }
    }
    core_acceptor={"seq":"CAGGT",
        "H":{
            "c0":["C>A","C>G","T>A","T>G"],
            "c1":["A>T","A>G","A>C"],
            "c2":["G>T","G>C","G>A"]
        }
    }
    isodict, genedict, isointerval, regioninterval = getRefDefSNV4(build)
    default = getdefaultisoform(build) # default[gene] = isoform
    snv4file = sys.argv[2] # variant in snv4 format per line
    snv4lst = []
    allvariant = []
    snv4dict = {}
    # svlst = []
    with open(snv4file,'r') as fh:
        # print (fh.readline().strip())
        for line in fh:
            new = line.strip().split('\t')
            marker = new[1]
            sample = new[0]
            index = marker+':'+sample
            snv4lst.append(marker)
            allvariant.append(index)
            snv4dict[index] = line.strip()
    isosnv4dict = geneassign(snv4lst, isodict, isointerval) # each isoform contains multiple snv4
    snv4anndict = {}
    for isoform in sorted(isosnv4dict):
        mygene = default[isoform]
        myisoform = isoform
        for mysnv4 in isosnv4dict[myisoform]:
            dnaDict, rnaDict, spliceDict = preannotate(mygene, myisoform, genedict, isodict, build)
            feature, muttype, aachange = annotateOneSNV4withgene(mysnv4, mygene, myisoform, isodict, build, dnaDict, rnaDict, spliceDict)
            interested_splice_site = determineSpliceSite(mysnv4, myisoform, isodict)
            n = mysnv4.split('.')
            pos = int(n[1])
            if len(n[2]) == 1 and len(n[3]) == 1:
                canonical_comment, neo_comment = annotateHLPA(mysnv4, myisoform, isodict, interested_splice_site, refGenome)
            elif len(n[2]) > 1 and len(n[3]) > 1 and len(n[2]) != len(n[3]):
                canonical_comment, neo_comment = annotateHLPA(mysnv4, myisoform, isodict, interested_splice_site, refGenome)
                neo_comment = neo_comment + "_indelspan"
                # if interested_splice_site == "D":
                    # canonical_comment = ""
                #     neo_comment = "donorindel"
                # elif interested_splice_site == "A":
                    # canonical_comment = ""
                #     neo_comment = "acceptorindel"
            else:
                canonical_comment = "ND"
                neo_comment = "ND"
            if mysnv4 not in snv4anndict:
                snv4anndict[mysnv4] = []
            if "splice" in muttype or "Splice" in muttype:
                aacomment = findClosestAA(pos, isodict, isoform, dnaDict)
                if aacomment == "*0":
                    aacomment = "UTR5"
                if 'span' not in muttype:
                    ss = muttype.split('_')
                    splicepos = ss[-2]+'_'+ss[-1]
                    feature = ss[0].replace(':','_')
                    if aacomment == "*0":
                        feature = "splice"
                    if "splice" in feature or "Splice" in feature:
                        newfeature = "splice"
                    else:
                        newfeature = feature
                    if aachange == "-1":
                        aachange = aacomment
                    fullcomment = aachange+'_'+splicepos+'_'+canonical_comment+'_'+neo_comment
                    snv4anndict[mysnv4] = [mygene, myisoform, newfeature, feature, fullcomment.replace('__','_')]
                else:
                    feature = muttype.replace(':','_')
                    if "splice" in feature or "Splice" in feature:
                        newfeature = "splice"
                    else:
                        newfeature = feature
                    if aachange == "-1":
                        aachange = aacomment
                    fullcomment = aachange+'_'+canonical_comment+'_'+neo_comment
                    snv4anndict[mysnv4] = [mygene, myisoform, newfeature, feature, fullcomment.replace('__','_')]
            else:
                snv4anndict[mysnv4] = [mygene, myisoform, feature, muttype, aachange]
    for index in allvariant:
        snv4 = index.split(':')[0]
        if snv4 not in set(snv4lst):
            print (snv4dict[index]+'\t'+'\t'.join(["ND"]*2))
        else:
            if snv4 not in snv4anndict:
                n = snv4.split('.')
                chrom = n[0]
                gnommutPos = int(n[1])
                refbase = n[2].replace('-','')
                altbase = n[3].replace('-','')
                reflen = len(refbase)
                if reflen == 0:
                    gnommutSt = gnommutPos
                    gnommutSp = gnommutPos
                else:
                    gnommutSt = gnommutPos
                    gnommutSp = gnommutSt + reflen - 1
                a = regioninterval[chrom].find(gnommutSt, gnommutSp)
                if len(a) > 0:
                    if len(a) > 1:
                        print (snv4dict[index]+'\tND\tND\twarning\tND\tND')
                    else:
                        b = a[0]
                        mygenename = b.value.split(':')[0]
                        myisoform = b.value.split(':')[1]
                        myregion = b.value.split(':')[2]
                        print (snv4dict[index]+'\t'+mygenename+'\t'+myisoform+'\t'+myregion+'\tND\tND')
                else:
                    print (snv4dict[index]+'\tND\tND\tintergenic\tND\tND')
            else:
                print (snv4dict[index]+'\t'+'\t'.join(snv4anndict[snv4]))
