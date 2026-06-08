from pyfaidx import Fasta
import sys
import argparse
import os
import subprocess
import glob
import re

PATH = "/research_jude/rgs01_jude/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/mRNA/ZL/reference/"
refGenome_hg38 = PATH+"GRCh38_no_alt_chr.fa"
refGenome_hg19 = PATH+"GRCh37-lite-chr.fa"
BEDTOOLS = PATH+"bedtools"
SAMTOOLS = PATH+"samtools"
BWA = PATH+"bwa"

# Input format: a list of Indel Markers
# Output format: Indel Markers: [SampleID, marker, nM, nT, nRa, nRb]

def IndelSorting(indel_lst):
    indeldict = {}
    sortedindeldict = {}
    sorted_indel_lst = []
    for indel in indel_lst:
        chrom = indel.split('.')[0]
        pos = int(indel.split('.')[1])
        name = chrom+"."+str(pos)
        if name not in sortedindeldict:
            sortedindeldict[name] = set([])
        sortedindeldict[name].add(indel)
        if chrom not in indeldict:
            indeldict[chrom] = set([])
        indeldict[chrom].add(pos)
    for chrom in sorted(indeldict):
        sortpos = sorted(list(indeldict[chrom]))
        for mypos in sortpos:
            name = chrom+'.'+str(mypos)
            for x in sortedindeldict[name]:
                sorted_indel_lst.append(x)
    return sorted_indel_lst

def IndelBatching(indel_lst, filename, flank=5000):
    batchdict = {}
    working_indel_lst = indel_lst
    num = 1
    while len(working_indel_lst) > 0:
        batchID = filename+'_'+str(num)
        subset_indel_lst = []
        for ind in range(0, len(working_indel_lst)-1):
            nowindel = working_indel_lst[ind]
            postindel = working_indel_lst[ind+1]
            nowindelpos = int(nowindel.split('.')[1])
            nowindelchrom = nowindel.split('.')[0]
            postindelpos = int(postindel.split('.')[1])
            postindelchrom = postindel.split('.')[0]
            dist = abs(nowindelpos-postindelpos)
            if batchID not in batchdict:
                batchdict[batchID] = []
            if (nowindelchrom == postindelchrom and dist > flank) or (nowindelchrom != postindelchrom):
                batchdict[batchID].append(nowindel)
            else:
                subset_indel_lst.append(nowindel)
        prelastindel = working_indel_lst[len(working_indel_lst)-2]
        lastindel = working_indel_lst[len(working_indel_lst)-1]
        prelastindelchrom = prelastindel.split('.')[0]
        lastindelchrom = lastindel.split('.')[0]
        prelastindelpos = int(prelastindel.split('.')[1])
        lastindelpos = int(lastindel.split('.')[1])
        dist = abs(prelastindelpos-lastindelpos)
        if (prelastindelchrom == lastindelchrom and dist > flank) or (prelastindelchrom != lastindelchrom):
            batchdict[batchID].append(lastindel)
        else:
            subset_indel_lst.append(lastindel)
        if len(subset_indel_lst) == len(working_indel_lst):
            break
        else:
            working_indel_lst = subset_indel_lst
            num += 1
    ###
    if len(subset_indel_lst) > 0:
        for indel in subset_indel_lst:
            batchID = filename+'_'+str(num)
            if batchID not in batchdict:
                batchdict[batchID] = []
                batchdict[batchID].append(indel)
            else:
                batchdict[batchID].append(indel)
            num += 1

    variantdict = {}
    for batch in sorted(batchdict):
        for vv in batchdict[batch]:
            if vv not in variantdict:
                variantdict[vv] = batch
    return batchdict, variantdict

def revcmp_str (input_str):
    rcs = input_str[::-1]
    rcs = rcs.replace('A','X')
    rcs = rcs.replace('T','A')
    rcs = rcs.replace('X','T')
    rcs = rcs.replace('C','X')
    rcs = rcs.replace('G','C')
    rcs = rcs.replace('X','G')
    return rcs 

def run_bwa_genome_index(batch, label):
    faPath = "{}_{}.fa".format(batch, label)
    faIndex = "{}_{}".format(batch, label)
    bwa_index = "{} index {} -p {}".format(BWA, faPath, faIndex)
    os.system(bwa_index)

def BatchIndexing(build, batchdict, flank=1000):
    if build == "hg19":
        refGenome = Fasta(refGenome_hg19)
    elif build == "hg38":
        refGenome = Fasta(refGenome_hg38)
    batchlst = []
    for batch in sorted(batchdict):
        refref = str(batch)+'_refref.fa'
        mutref = str(batch)+'_mutref.fa'
        batchlst.append(str(batch))
        with open(refref,'w') as out1, open(mutref,'w') as out2:
            for indel in sorted(batchdict[batch]):
                chrom = 'chr'+indel.split('.')[0]
                chrom = chrom.replace('chrchr','chr')
                mutpos = int(indel.split('.')[1])
                refbase = indel.split('.')[2]
                altbase = indel.split('.')[3]
                start = mutpos-flank
                if start < 0:
                    start = 0
                stop = mutpos+flank
                seq = str(refGenome[chrom][start:stop])
                # seq = refGenome[kdict[chrom]][start:stop]
                variant = indel

                left = ""
                mutmid = ""
                refmid = ""
                right = ""
                mutst = mutpos-1
                mutsp = mutst+len(refbase)

                refposdict = {}
                mutposdict = {}
                refstop = start
                mutstop = start
                fragpos = 1

                for i in range(start, stop):
                    seqpos = i-start
                    refposdict[fragpos] = i
                    if i < mutst:
                        refstop += 1
                        mutstop += 1
                        left += seq[seqpos]
                        mutposdict[fragpos] = i
                    elif i >= mutst and i < mutsp:
                        refmid += seq[seqpos]
                        mutmid = altbase
                        refstop += 1
                        for tt in range(len(altbase)):
                            mutposdict[fragpos] = "X"
                            mutstop += 1
                    elif i >= mutsp:
                        refstop += 1
                        mutstop += 1
                        right += seq[seqpos]
                        mutposdict[fragpos] = i
                    fragpos += 1
                if len(refbase) > 0:
                    mutseq = left+mutmid+right
                else:
                    mutseq = left+altbase+right
                refseq = left+refmid+right
                out1.write(">"+variant+"\n"+refseq+'\n')
                out2.write(">"+variant+"\n"+mutseq+'\n')
        run_bwa_genome_index(str(batch), "refref")
        run_bwa_genome_index(str(batch), "mutref")

    return batchlst

def faExtract(snv, refGenome, flank=100):
    n = snv.split('.')
    chrom = n[0]
    pos = int(n[1])
    start = pos-flank
    if start < 0:
        start = 0
    stop = pos+flank
    seq = refGenome[chrom][(start-1):stop]
    return str(start), str(stop), seq

def translate_cigar (input_str):
    length = re.findall(r'[0-9]+',input_str)
    op = re.findall(r'[A-Z]+',input_str)
    rlt = {}
    rlt['op'] = op
    rlt['length'] = length
    return rlt

def call_cigar (df):
    cond1 = ["M", "I", "D"]
    cond2 = ["M", "I", "D", "N"]
    cond3 = ["M", "I", "D", "N", "S"]
    if len(df['op']) == 1:
        mut = "a1SNV"
    else:
        uniq_op = set(df['op'])
        if len(uniq_op.difference(cond1)) == 0:
            mut = "a2Indel"
        elif len(uniq_op.difference(cond2)) == 0:
            mut = "a3Intron"
        elif len(uniq_op.difference(cond3)) == 0:
            mut = "a4SV"
        else:
            mut = "nocall"
    return mut

def SeqFilterSingleTest(fq, sample_lst, out1):
    with open(fq,'r') as fqfile:
        lines = []
        readlines = []
        readdict = {}
        readset = set([])
        for line in fqfile:
            new = line.strip()
            if not new.startswith("@") and new in set(sample_lst):
                sample = new
            else:
                lines.append(new)
                if len(lines) == 4:
                    myread = lines[0].replace('@','')
                    if myread not in readset:
                        pair = 1
                        out1.write('@'+sample+':'+lines[0].replace('@','')+'_read{}'.format(str(pair))+'\n')
                        out1.write(lines[1]+'\n')
                        out1.write(lines[2]+'\n')
                        out1.write(lines[3]+'\n')
                        readset.add(myread)
                        lines = []
                    else:
                        pair += 1
                        out1.write('@'+sample+':'+lines[0].replace('@','')+'_read{}'.format(str(pair))+'\n')
                        out1.write(lines[1]+'\n')
                        out1.write(lines[2]+'\n')
                        out1.write(lines[3]+'\n')
                        lines = []

def SeqFilterSingle(fq, sample, out1):
    with open(fq,'r') as fqfile:
        lines = []
        readlines = []
        readdict = {}
        readset = set([])
        for line in fqfile:
            new = line.strip()
            lines.append(new)
            if len(lines) == 4:
                myread = lines[0].replace('@','')
                if myread not in readset:
                    pair = 1
                    out1.write('@'+sample+':'+lines[0].replace('@','')+'_read{}'.format(str(pair))+'\n')
                    out1.write(lines[1]+'\n')
                    out1.write(lines[2]+'\n')
                    out1.write(lines[3]+'\n')
                    readset.add(myread)
                    lines = []
                else:
                    pair += 1
                    out1.write('@'+sample+':'+lines[0].replace('@','')+'_read{}'.format(str(pair))+'\n')
                    out1.write(lines[1]+'\n')
                    out1.write(lines[2]+'\n')
                    out1.write(lines[3]+'\n')
                    lines = []

def run_bwa_genome_index_sv(sv, label):
    faPath = "{}_{}.fa".format(sv, label)
    faIndex = "{}_{}".format(sv, label)
    bwa_index = "{} index {} -p {}".format(BWA, faPath, faIndex)
    os.system(bwa_index)

def SamParsing(sv, sample, pcut=0.95):
    Rd = {}
    for label in ['seqA','seqB','join']:
        sam = sv+'_'+label+'.'+sample+'.sam'
        with open(sam,'r') as fh:
            for line in fh:
                if line.startswith("@"):continue
                new = line.strip().split('\t')
                read = new[0]
                qual = new[9]
                score = phred_read(qual)
                variant = new[2]
                cigar = new[5]
                seq = new[8]
                if score < pcut*len(seq):continue
                if read not in Rd:
                    Rd[read] = {}
                if label not in Rd[read]:
                    Rd[read][label] = []
                Rd[read][label].append(cigar)
    return Rd

def translatecigar(cigar):
    length = re.findall(r'[0-9]+',cigar)
    op = re.findall(r'[A-Z]+',cigar)
    match = 0
    for mylen, myop in zip(length, op):
        if myop == "M":
            match += int(mylen)
    return match

# def determineChr(bam):
#     chrstring = subprocess.check_output(["{} view -H {} | head -2 | tail -1".format(SAMTOOLS, bam)], shell=True).decode('utf-8').strip().split('\t')[1].split(':')[1]
#     return chrstring

def determineChr(bam: str) -> bool:
    """ parse bam file to get chromosome format and size; tested """
    ret = False
    cmd = f"{SAMTOOLS} view -H {bam} | grep SN | grep LN | head -1"
    rlt = subprocess.check_output(cmd, shell=True).decode('utf-8').strip().split('\t')[1]
    #rlt = subprocess.check_output(["samtools view -H {} | grep SN | grep LN | head -1".format(bam)], shell=True).decode('utf-8').strip().split('\t')[1]
    chrom = rlt.split(':')[1]
    return (chrom)

def SeqExtractor(indel, sample_lst, bam_lst, chrom_lst, build, filename, flank=200):
    fq1handle = open("{}_genome_single.fq".format(indel),'w')
    allprefq = "{}_allsamples.fq".format(indel)
    for sample, bam, mychr in zip(sample_lst, bam_lst, chrom_lst): # chrom_lst can be Chr, chr, or ''
        chrom = indel.split('.')[0]
        mutpos = int(indel.split('.')[1])
        refbase = indel.split('.')[2]
        altbase = indel.split('.')[3]
        variant = indel
        fq = '{}_{}_pre.fq'.format(indel, sample)
        mybam = "{0}_{1}.bam".format(indel, sample)
        start = mutpos-flank
        if start < 0:
            start = 0
        stop = mutpos+flank
        if mychr == '':
            region = chrom.replace('Chr','').replace('chr','')+':'+str(start)+'-'+str(stop)
        elif mychr == "chr":
            region = chrom.replace('Chr','chr')+':'+str(start)+'-'+str(stop)
        elif mychr == "Chr":
            region = chrom.replace('chr','Chr')+':'+str(start)+'-'+str(stop)
        cmd = '{0} view -h -bS {1} {2} | {0} sort -n > {3}'.format(SAMTOOLS, bam, region, mybam)
        os.system(cmd)
        cmd = '{} bamtofastq -i {} -fq {}'.format(BEDTOOLS, mybam, fq)
        os.system(cmd)
        cmd = "echo {} >> {}".format(sample, allprefq)
        os.system(cmd)
        cmd = "cat {} >> {}".format(fq, allprefq)
        os.system(cmd)
        # os.system("rm {}".format(mybam))
        # os.system("rm {}".format(fq))
    SeqFilterSingleTest(allprefq, sample_lst, fq1handle)
    os.system("rm {}".format(allprefq))
    fq1handle.close()

def Indexbuilding_OnetoOne(indel, refGenome, flank=1000):
    n = indel.split('.')
    indelindex = n[0]+'.'+n[1]+'.'+n[2]+'.'+n[3] # yfeng27 updated (prev: indelindex = n[0]+'.'+n[1]+'.'+n[2])
    refref = str(indelindex)+'_refref.fa'
    mutref = str(indelindex)+'_mutref.fa'
    with open(refref,'w') as out1, open(mutref,'w') as out2:
        chrom = 'chr'+indel.split('.')[0]
        chrom = chrom.replace('chrchr','chr')
        mutpos = int(indel.split('.')[1])
        refbase = indel.split('.')[2]
        altbase = indel.split('.')[3]
        start = mutpos-flank
        if start < 0:
            start = 0
        stop = mutpos+flank
        seq = str(refGenome[chrom][start:stop])
        variant = indel

        left = ""
        mutmid = ""
        refmid = ""
        right = ""
        mutst = mutpos-1
        mutsp = mutst+len(refbase)

        refposdict = {}
        mutposdict = {}
        refstop = start
        mutstop = start
        fragpos = 1

        for i in range(start, stop):
            seqpos = i-start
            refposdict[fragpos] = i
            if i < mutst:
                refstop += 1
                mutstop += 1
                left += seq[seqpos]
                mutposdict[fragpos] = i
            elif i >= mutst and i < mutsp:
                refmid += seq[seqpos]
                mutmid = altbase
                refstop += 1
                for tt in range(len(altbase)):
                    mutposdict[fragpos] = "X"
                    mutstop += 1
            elif i >= mutsp:
                refstop += 1
                mutstop += 1
                right += seq[seqpos]
                mutposdict[fragpos] = i
            fragpos += 1
        if len(refbase) > 0:
            mutseq = left+mutmid+right
        else:
            mutseq = left+altbase+right
        refseq = left+refmid+right
        out1.write(">"+variant+"\n"+refseq+'\n')
        out2.write(">"+variant+"\n"+mutseq+'\n')
    run_bwa_genome_index_OnetoOne(indelindex, "refref")
    run_bwa_genome_index_OnetoOne(indelindex, "mutref")

def run_bwa_genome_index_OnetoOne(indelindex, label):
    faPath = "{}_{}.fa".format(indelindex, label)
    faIndex = "{}_{}".format(indelindex, label)
    bwa_index = "{} index {} -p {}".format(BWA, faPath, faIndex)
    os.system(bwa_index)

def run_bwa_OnetoOne(variant):
    fq1 = '{}_genome_single.fq'.format(variant)
    n = variant.split('.')
    variantindex = n[0]+'.'+n[1]+'.'+n[2]+'.'+n[3] # yfeng27 updated (prev: variantindex = n[0]+'.'+n[1]+'.'+n[2])
    for label in ['refref','mutref']:
        faIndex = "{}_{}".format(variantindex, label)
        sam = label+'_'+variant+'.sam'
        sai = label+'_'+variant+'.sai'
        os.system("{} aln -t 4 {} {} > {}".format(BWA, faIndex, fq1, sai))
        os.system("{} samse {} {} {} > {}".format(BWA, faIndex, sai, fq1, sam))
    os.system("rm {}_genome_single.fq".format(variant)) # yfeng27 previously commented by zk

def run_bwa(mybatch, variant, filename):
    fq1 = '{}_genome_single.fq'.format(variant)
    for label in ['refref','mutref']:
        faIndex = "{}_{}".format(mybatch, label)
        sam = mybatch+'_'+label+'_'+variant+'.sam'
        sai = mybatch+'_'+label+'_'+variant+'.sai'
        os.system("{} aln -t 4 {} {} > {}".format(BWA, faIndex, fq1, sai))
        os.system("{} samse {} {} {} > {}".format(BWA, faIndex, sai, fq1, sam))
    os.system("rm {}_genome_single.fq".format(variant))

def run_bwa_sv(sv, sample):
    fq = sv+'.refAB.fq'
    for label in ['seqA','seqB','join']:
        faIndex = "{}_{}".format(sv, label)
        sam = sv+'_'+label+'.'+sample+'.sam'
        cmd = "{} mem {} {} > {}".format(BWA, faIndex, fq, sam)
        os.system(cmd)

def phred_read(input_str, qcut=20):
    good = 0
    for i in range(len(input_str)):
        score = ord(input_str[i])-33
        if score > qcut:
            good += 1
    return float(good)

def TableJoining(scanvariant, mybatch, outhandle, flank=1000, pericut=3):
    refsam = "{0}.{1}.tmp.sam".format(scanvariant, 'refref')
    mutsam = "{0}.{1}.tmp.sam".format(scanvariant, 'mutref')
    n = scanvariant.split('.')
    refbase = len(n[2])
    mutbase = len(n[3])
    indelindex = n[0]+'.'+n[1]+'.'+n[2]+'.'+n[3] # yfeng27 updated (prev: indelindex = n[0]+'.'+n[1]+'.'+n[2])
    # needs to define core vs peri regions for each variant
    refstart = flank
    refstop = flank+refbase
    mutstart = flank
    mutstop = flank+mutbase
    ###
    # get fasta index
    if mybatch == "":    
        refseqindex = Fasta(indelindex+'_refref.fa')
        mutseqindex = Fasta(indelindex+'_mutref.fa')
    else:
        refseqindex = Fasta(mybatch+'_refref.fa')
        mutseqindex = Fasta(mybatch+'_mutref.fa')
    ###
    readict = {}
    with open(refsam,'r') as fh, open(mutsam,'r') as sh:
        for line in fh:
            if line.startswith('@SQ'):continue
            new = line.strip().split('\t') # sample:readid_r1
            readid = '_'.join(new[0].split('_')[:2])
            pair = new[0].split('_')[2]
            variant = new[2]
            mycigar = new[5]
            myseqq = new[9]
            myst = int(new[3])
            mysp = myst+len(myseqq)-1
            myquaa = new[10]
            if readid not in readict:
                readict[readid] = {}
            if 'Refref' not in readict[readid]:
                readict[readid]['Refref'] = {}
            if variant == "*":
                if pair == 'read1':
                    readict[readid]['Refref']['R1'] = "*"
                    readict[readid]['Refref']['Qual1'] = "*"
                    readict[readid]['Refref']['Cigar1'] = "*"
                    readict[readid]['Refref']['St1'] = "*"
                    readict[readid]['Refref']['Sp1'] = "*"
                elif pair == 'read2':
                    readict[readid]['Refref']['R2'] = "*"
                    readict[readid]['Refref']['Qual2'] = "*"
                    readict[readid]['Refref']['Cigar2'] = "*"
                    readict[readid]['Refref']['St2'] = "*"
                    readict[readid]['Refref']['Sp2'] = "*"
            else:
                if pair == 'read1':
                    readict[readid]['Refref']['R1'] = myseqq
                    readict[readid]['Refref']['Qual1'] = myquaa
                    readict[readid]['Refref']['Cigar1'] = mycigar
                    readict[readid]['Refref']['St1'] = myst
                    readict[readid]['Refref']['Sp1'] = mysp
                elif pair == 'read2':
                    readict[readid]['Refref']['R2'] = myseqq
                    readict[readid]['Refref']['Qual2'] = myquaa
                    readict[readid]['Refref']['Cigar2'] = mycigar
                    readict[readid]['Refref']['St2'] = myst
                    readict[readid]['Refref']['Sp2'] = mysp
        
        for line in sh:
            if line.startswith('@SQ'):continue
            new = line.strip().split('\t')
            readid = '_'.join(new[0].split('_')[:2])
            pair = new[0].split('_')[2]
            variant = new[2]
            mycigar = new[5]
            myseqq = new[9]
            myquaa = new[10]
            myst = int(new[3])
            mysp = myst+len(myseqq)-1
            if readid not in readict:
                readict[readid] = {}
            if 'Mutref' not in readict[readid]:
                readict[readid]['Mutref'] = {}
            if variant == "*":
                if pair == 'read1':
                    readict[readid]['Mutref']['R1'] = "*"
                    readict[readid]['Mutref']['Qual1'] = "*"
                    readict[readid]['Mutref']['Cigar1'] = "*"
                    readict[readid]['Mutref']['St1'] = "*"
                    readict[readid]['Mutref']['Sp1'] = "*"
                elif pair == 'read2':
                    readict[readid]['Mutref']['R2'] = "*"
                    readict[readid]['Mutref']['Qual2'] = "*"
                    readict[readid]['Mutref']['Cigar2'] = "*"
                    readict[readid]['Mutref']['St2'] = "*"
                    readict[readid]['Mutref']['Sp2'] = "*"
            else:
                if pair == 'read1':
                    readict[readid]['Mutref']['R1'] = myseqq
                    readict[readid]['Mutref']['Qual1'] = myquaa
                    readict[readid]['Mutref']['Cigar1'] = mycigar
                    readict[readid]['Mutref']['St1'] = myst
                    readict[readid]['Mutref']['Sp1'] = mysp
                elif pair == 'read2':
                    readict[readid]['Mutref']['R2'] = myseqq
                    readict[readid]['Mutref']['Qual2'] = myquaa
                    readict[readid]['Mutref']['Cigar2'] = mycigar
                    readict[readid]['Mutref']['St2'] = myst
                    readict[readid]['Mutref']['Sp2'] = mysp

    samplecount = {}
    for read in sorted(readict):
        sample = read.split(':')[0]
        if sample not in samplecount:
            samplecount[sample] = {'nR':set([]),'nM':set([])}
        readqual = []
        if 'Refref' in readict[read]:
            if 'R1' in readict[read]['Refref']:
                if readict[read]['Refref']['R1'] == "*":
                    R1ref = "*"
                    Qual1ref = "*"
                    Cigar1ref = "*"
                    St1ref = "*"
                    Sp1ref = "*"
                else:
                    R1ref = readict[read]['Refref']['R1']
                    Qual1ref = phred_read(readict[read]['Refref']['Qual1'])/len(readict[read]['Refref']['Qual1'])
                    Cigar1ref = readict[read]['Refref']['Cigar1']
                    St1ref = readict[read]['Refref']['St1']
                    Sp1ref = readict[read]['Refref']['Sp1']
                    readqual.append(Qual1ref)
            else:
                R1ref = 'ND'
                Qual1ref = 'ND'
                Cigar1ref = 'ND'
                St1ref = 'ND'
                Sp1ref = 'ND'
            if 'R2' in readict[read]['Refref']:
                if readict[read]['Refref']['R2'] == "*":
                    R2ref = "*"
                    Qual2ref = "*"
                    Cigar2ref = "*"
                    St2ref = "*"
                    Sp2ref = "*"
                else:
                    R2ref = readict[read]['Refref']['R2']
                    Qual2ref = phred_read(readict[read]['Refref']['Qual2'])/len(readict[read]['Refref']['Qual2'])
                    Cigar2ref = readict[read]['Refref']['Cigar2']
                    St2ref = readict[read]['Refref']['St2']
                    Sp2ref = readict[read]['Refref']['Sp2']
                    readqual.append(Qual2ref)
            else:
                R2ref = 'ND'
                Qual2ref = 'ND'
                Cigar2ref = 'ND'
                St2ref = 'ND'
                Sp2ref = 'ND'
        else:
            R1ref = 'ND'
            Qual1ref = 'ND'
            Cigar1ref = 'ND'
            St1ref = 'ND'
            Sp1ref = 'ND'
            R2ref = 'ND'
            Qual2ref = 'ND'
            Cigar2ref = 'ND'
            St2ref = 'ND'
            Sp2ref = 'ND'
            ###
        if 'Mutref' in readict[read]:
            if 'R1' in readict[read]['Mutref']:
                if readict[read]['Mutref']['R1'] == "*":
                    R1mut = "*"
                    Qual1mut = "*"
                    Cigar1mut = "*"
                    St1mut = "*"
                    Sp1mut = "*"
                else:
                    R1mut = readict[read]['Mutref']['R1']
                    Qual1mut = phred_read(readict[read]['Mutref']['Qual1'])/len(readict[read]['Mutref']['Qual1'])
                    Cigar1mut = readict[read]['Mutref']['Cigar1']
                    St1mut = readict[read]['Mutref']['St1']
                    Sp1mut = readict[read]['Mutref']['Sp1']
                    readqual.append(Qual1mut)
            else:
                R1mut = 'ND'
                Qual1mut = 'ND'
                Cigar1mut = 'ND'
                St1mut = 'ND'
                Sp1mut = 'ND'
            if 'R2' in readict[read]['Mutref']:
                if readict[read]['Mutref']['R2'] == "*":
                    R2mut = "*"
                    Qual2mut = "*"
                    Cigar2mut = "*"
                    St2mut = "*"
                    Sp2mut = "*"
                else:
                    R2mut = readict[read]['Mutref']['R2']
                    Qual2mut = phred_read(readict[read]['Mutref']['Qual2'])/len(readict[read]['Mutref']['Qual2'])
                    Cigar2mut = readict[read]['Mutref']['Cigar2']
                    St2mut = readict[read]['Mutref']['St2']
                    Sp2mut = readict[read]['Mutref']['Sp2']
                    readqual.append(Qual2mut)
            else:
                R2mut = 'ND'
                Qual2mut = 'ND'
                Cigar2mut = 'ND'
                St2mut = 'ND'
                Sp2mut = 'ND'
        else:
            R1mut = 'ND'
            Qual1mut = 'ND'
            Cigar1mut = 'ND'
            St1mut = 'ND'
            Sp1mut = 'ND'
            R2mut = 'ND'
            Qual2mut = 'ND'
            Cigar2mut = 'ND'
            St2mut = 'ND'
            Sp2mut = 'ND'
        if len(readqual) == 0 or max(readqual) < 0.95:continue
        # print (read, readict[read])
        if Cigar1ref != "ND" and Cigar1mut != "ND":
            if Cigar1ref != Cigar1mut:
                if Cigar1ref != '*' and Cigar1mut != '*':
                    opref = set(re.findall(r'[A-Z]+',Cigar1ref))
                    opmut = set(re.findall(r'[A-Z]+',Cigar1mut))
                    if (len(opref) == 1 and list(opref)[0] == 'M'):
                        samplecount[sample]['nR'].add(read)
                    elif (len(opmut) == 1 and list(opmut)[0] == 'M'):
                        samplecount[sample]['nM'].add(read)
                else:
                    if Cigar1ref != '*' and Cigar1mut == '*':
                        opref = set(re.findall(r'[A-Z]+',Cigar1ref))
                        if (len(opref) == 1 and list(opref)[0] == 'M'):
                            # samplecount[sample]['nR'].add(read) # ***simpler
                            if ((St1ref < refstart and Sp1ref > refstart) or (St1ref < refstop and Sp1ref > refstop)):
                                refcoremiscount, refperimiscount = basebybaseCounting(St1ref, Sp1ref, R1ref, scanvariant, refstart, refstop, refseqindex)
                                if refcoremiscount == 0 and refperimiscount < pericut:
                                    samplecount[sample]['nR'].add(read)
                    elif Cigar1ref == '*' and Cigar1mut != '*':
                        opref = set(re.findall(r'[A-Z]+',Cigar1mut))
                        if (len(opref) == 1 and list(opref)[0] == 'M'):
                            # samplecount[sample]['nM'].add(read) # ***simpler
                            if ((St1mut < mutstart and Sp1mut > mutstart) or (St1mut < mutstop and Sp1mut > mutstop)):
                                mutcoremiscount, mutperimiscount = basebybaseCounting(St1mut, Sp1mut, R1mut, scanvariant, mutstart, mutstop, mutseqindex)
                                if mutcoremiscount == 0 and mutperimiscount < pericut:
                                    samplecount[sample]['nM'].add(read)
            else:
                if Cigar1ref != '*' and Cigar1mut != '*':
                    opref = set(re.findall(r'[A-Z]+',Cigar1ref))
                    if (len(opref) == 1 and list(opref)[0] == 'M') and (St1ref == St1mut):
                        if ((St1ref < refstart and Sp1ref > refstart) or (St1ref < refstop and Sp1ref > refstop)):
                            refcoremiscount, refperimiscount = basebybaseCounting(St1ref, Sp1ref, R1ref, scanvariant, refstart, refstop, refseqindex)
                            mutcoremiscount, mutperimiscount = basebybaseCounting(St1mut, Sp1mut, R1mut, scanvariant, mutstart, mutstop, mutseqindex)
                            if refcoremiscount == 0 and refcoremiscount < mutcoremiscount and refperimiscount < pericut and mutcoremiscount >= pericut:
                                samplecount[sample]['nR'].add(read)
                            elif mutcoremiscount == 0 and mutcoremiscount < refcoremiscount and mutperimiscount < pericut and refcoremiscount >= pericut:
                                samplecount[sample]['nM'].add(read)

        if Cigar2ref != "ND" and Cigar2mut != "ND":
            if Cigar2ref != Cigar2mut:
                if Cigar2ref != '*' and Cigar2mut != '*':
                    opref = set(re.findall(r'[A-Z]+',Cigar2ref))
                    opmut = set(re.findall(r'[A-Z]+',Cigar2mut))
                    if (len(opref) == 1 and list(opref)[0] == 'M'):
                        samplecount[sample]['nR'].add(read)
                    elif (len(opmut) == 1 and list(opmut)[0] == 'M'):
                        samplecount[sample]['nM'].add(read)
                else:
                    if Cigar2ref != '*' and Cigar2mut == '*':
                        opref = set(re.findall(r'[A-Z]+',Cigar2ref))
                        if (len(opref) == 1 and list(opref)[0] == 'M'):
                            # samplecount[sample]['nR'].add(read) # ***simpler
                            if ((St2ref < refstart and Sp2ref > refstart) or (St2ref < refstop and Sp2ref > refstop)):
                                refcoremiscount, refperimiscount = basebybaseCounting(St2ref, Sp2ref, R2ref, scanvariant, refstart, refstop, refseqindex)
                                if refcoremiscount == 0 and refperimiscount < pericut:
                                    samplecount[sample]['nR'].add(read)
                    elif Cigar2ref == '*' and Cigar2mut != '*':
                        opref = set(re.findall(r'[A-Z]+',Cigar2mut))
                        if (len(opref) == 1 and list(opref)[0] == 'M'):
                            # samplecount[sample]['nM'].add(read) # ***simpler
                            if ((St2mut < mutstart and Sp2mut > mutstart) or (St2mut < mutstop and Sp2mut > mutstop)):
                                mutcoremiscount, mutperimiscount = basebybaseCounting(St2mut, Sp2mut, R2mut, scanvariant, mutstart, mutstop, mutseqindex)
                                if mutcoremiscount == 0 and mutperimiscount < pericut:
                                    samplecount[sample]['nM'].add(read)
            else:
                if Cigar2ref != '*' and Cigar2mut != '*':
                    opref = set(re.findall(r'[A-Z]+',Cigar2ref))
                    if (len(opref) == 1 and list(opref)[0] == 'M') and (St2ref == St2mut):
                        if ((St2ref < refstart and Sp2ref > refstart) or (St2ref < refstop and Sp2ref > refstop)):
                            refcoremiscount, refperimiscount = basebybaseCounting(St2ref, Sp2ref, R2ref, scanvariant, refstart, refstop, refseqindex)
                            mutcoremiscount, mutperimiscount = basebybaseCounting(St2mut, Sp2mut, R2mut, scanvariant, mutstart, mutstop, mutseqindex)
                            if refcoremiscount == 0 and refcoremiscount < mutcoremiscount and refperimiscount < pericut and mutcoremiscount >= pericut:
                                samplecount[sample]['nR'].add(read)
                            elif mutcoremiscount == 0 and mutcoremiscount < refcoremiscount and mutperimiscount < pericut and refcoremiscount >= pericut:
                                samplecount[sample]['nM'].add(read)

    for sample in sorted(samplecount):
        nT = len(samplecount[sample]['nR'])+len(samplecount[sample]['nM'])
        output = [sample, scanvariant, str(len(samplecount[sample]['nM'])), str(nT)]
        outhandle.write('\t'.join(output)+'\n')

def regionCounting(startlst, stoplst, queryseq, querySt, querySp, seqindex, variant):
    miscount = 0
    for start, stop in zip(startlst, stoplst):
        seq = seqindex[variant][start-1:stop-1]
        mystart = start-querySt
        mystop = stop-querySp
        qseq = queryseq[mystart:mystop-1]
        for i,j in zip(qseq, seq):
            if i != j:
                miscount += 1
    return miscount

def delIndex():
    os.system("rm *_genome_single.fq *.amb *.ann *.bwt *.fa *.fai *.pac *.sa")

def basebybaseCounting(readSt, readSp, readseq, variant, variantSt, variantSp, seqindex):
    coreSt = []
    coreSp = []
    periSt = []
    periSp = []
    if variantSt >= readSt and variantSp <= readSp:
        coreSt = [variantSt]
        coreSp = [variantSp]
        periSt = [readSt, variantSp+1]
        periSp = [variantSt-1, readSp]
    elif variantSt < readSt and variantSp >= readSt:
        coreSt = [readSt]
        coreSp = [variantSp]
        periSt = [variantSp+1]
        periSp = [readSp]
    elif variantSt <= readSp and variantSp > readSp:
        coreSt = [variantSt]
        coreSp = [readSp]
        periSt = [readSt]
        periSp = [variantSp+1]
    coremiscount = regionCounting(coreSt, coreSp, readseq, readSt, readSp, seqindex, variant)
    perimiscount = regionCounting(periSt, periSp, readseq, readSt, readSp, seqindex, variant)
    return coremiscount, perimiscount

def IndelParsing_OnetoOne(mybatch, variant, outhandle, pcut=0.95):
    for label in ['refref','mutref']:
        mfiles = glob.glob('{0}_{1}.sam'.format(label, variant))
        with open("{0}.{1}.tmp.sam".format(variant, label),'w') as out:
            for sfile in mfiles:
                with open(sfile,'r') as fh:
                    for line in fh:
                        new = line.strip().split('\t')
                        if new[2] == variant or new[2] == '*':
                            out.write(line)
    mybatch = ""
    TableJoining(variant, mybatch, outhandle)
    for label in ['refref','mutref']:
        cmd = "rm {0}.{1}.tmp.sam".format(variant, label)
        os.system(cmd)
        cmd = "rm {0}*{1}_{2}*sam {0}*{1}_{2}*sai".format(mybatch, label, variant)
        os.system(cmd)

def IndelParsing(mybatch, variant, outhandle, pcut=0.95):
    for label in ['refref','mutref']:
        mfiles = glob.glob('{0}_{1}_{2}.sam'.format(mybatch, label, variant))
        with open("{0}.{1}.tmp.sam".format(variant, label),'w') as out:
            for sfile in mfiles:
                with open(sfile,'r') as fh:
                    for line in fh:
                        new = line.strip().split('\t')
                        if new[2] == variant or new[2] == '*':
                            out.write(line)
    TableJoining(variant, mybatch, outhandle)
    for label in ['refref','mutref']:
        cmd = "rm {0}.{1}.tmp.sam".format(variant, label)
        ## os.system(cmd)
        cmd = "rm {0}*{1}_{2}*sam {0}*{1}_{2}*sai".format(mybatch, label, variant)
        ## os.system(cmd)

def bamExtract(snv, bam, mychr, flank=200):
    n = snv.split('.')
    chrom = n[0]
    pos = int(n[1])
    start = pos-flank
    if start < 0:
        start = 0
    stop = pos+flank
    if mychr == 'chr':
        newchr = chrom.replace('Chr','chr')
        region = newchr+':'+str(start)+'-'+str(stop)
    elif mychr == 'Chr':
        newchr = chrom.replace('chr','Chr')
        region = newchr+':'+str(start)+'-'+str(stop)
    elif mychr == '':
        newchr = chrom.replace('chr','').replace('Chr','')
        region = newchr+':'+str(start)+'-'+str(stop)
    sam = snv+'.sam'
    cmd = "{} view {} {} > {}".format(SAMTOOLS, bam, region, sam)
    os.system(cmd)

def countfromSam(snv, seq, qcut=20, pcut=0.95):
    sam = snv+'.sam'
    mypos = int(snv.split('.')[1])
    good_nucl = set(["A","G","T","C"])
    valid_nucl = ["A","C","G","T","N"]

    holder = {}

    used = set([])
    readholder = {}
    with open(sam,'r') as fh:
        for line in fh:
            new = line.strip().split('\t')
            read = new[0]
            start = int(new[3])
            cigar = new[5]
            seq = new[9]
            seqlen = len(seq)
            qual = new[10]
            good = phred_read(str(seq))
            seqqual = good/len(seq)
            if seqqual < pcut:continue
            mut = call_cigar(translate_cigar(cigar))
            if 'SNV' not in mut:continue
            if read not in readholder:
                readholder[read] = {}
            base1 = list(seq)
            qual1 = list(qual)
            pys1 = start-1
            for ii in range(pys1, pys1+seqlen):
                myqual = ord(qual1[ii-pys1])-qcut
                mybase = base1[ii-pys1]
                if mybase not in good_nucl:continue
                if ii not in readholder[read]:
                    readholder[read][ii] = []
                readholder[read][ii].append([mybase, myqual])
    
    for read in sorted(readholder):
        for pos in sorted(readholder[read]):
            if pos not in holder:
                holder[pos] = {}
            if len(readholder[read][pos]) == 1:
                mystring = readholder[read][pos][0]
                mybase = mystring[0]
                myqual = mystring[1]
                if myqual >= 33:
                    if mybase not in holder[pos]:
                        holder[pos][mybase] = 0
                    holder[pos][mybase] += 1
            elif len(readholder[read][pos]) == 2:
                mystring1 = readholder[read][pos][0]
                mybase1 = mystring1[0]
                myqual1 = mystring1[1]
                ###
                mystring2 = readholder[read][pos][1]
                mybase2 = mystring2[0]
                myqual2 = mystring2[1]
                if myqual1 >= 33 or myqual2 >= 33: # checking for single nucleotide quality
                    if myqual1 >= myqual2:
                        if mybase1 not in holder[pos]:
                            holder[pos][mybase1] = 0
                        holder[pos][mybase1] += 1
                    else:
                        if mybase2 not in holder[pos]:
                            holder[pos][mybase2] = 0
                        holder[pos][mybase2] += 1
    ## os.system("rm {}*".format(snv))
    return holder[mypos-1]

def countfromSamSV(Rd):
    nRa = set([])
    nRb = set([])
    nM = set([])
    nT = set([])
    for read in sorted(Rd):
        if ('seqA' not in Rd[read] and 'seqB' not in Rd[read]):continue
        if ('seqA' in Rd[read]) and ('seqB' in Rd[read]):
            for refA, refB in zip(Rd[read]['seqA'],Rd[read]['seqB']):
                if refA == refB:
                    label = 'seqA'
                else:
                    refAlen = translatecigar(refA)
                    refBlen = translatecigar(refB)
                    if refAlen >= refBlen:
                        label = "seqA"
                    else:
                        label = "seqB"
        elif ('seqA' in Rd[read]) and ('seqB' not in Rd[read]):
            label = 'seqA'
        elif ('seqA' not in Rd[read]) and ('seqB' in Rd[read]):
            label = 'seqB'
        for ref, mut in zip(Rd[read][label],Rd[read]['join']):
            if ref == mut:continue
            opref = set(re.findall(r'[A-Z]+',ref))
            opmut = set(re.findall(r'[A-Z]+',mut))
            if (len(opref) == 1 and list(opref)[0] == 'M') and (len(opmut) > 1):
                r = read.split(':')
                realread1 = ':'.join(r[1:])
                realread = '_'.join(realread1.split('_')[:-1])
                if len(realread) == 0:continue
                if label == 'seqA':
                    nRa.add(realread)
                elif label == 'seqB':
                    nRb.add(realread)
                # print ("Ref", realread, read)
                nT.add(realread)
            elif (len(opmut) == 1 and list(opmut)[0] == 'M') and (len(opref) > 1):
                r = read.split(':')
                realread1 = ':'.join(r[1:])
                realread = '_'.join(realread1.split('_')[:-1])
                if len(realread) == 0:continue
                nM.add(realread)
                nT.add(realread)
    return len(nM),len(nRa),len(nRb),len(nT)

def RegionExt(chrA, chrB, posA, posB, strandA, strandB, onebam, flank):
    startA = posA-flank
    endA = posA+flank
    startB = posB-flank
    endB = posB+flank
    string = determineChr(onebam)
    if 'chr' in string:
        newchrom = 'chr'+chrA
        chromA = newchrom.replace('chrchr','chr')
        regionA = chromA+':'+str(startA)+'-'+str(endA)
        newchrom = 'chr'+chrB
        chromB = newchrom.replace('chrchr','chr')
        regionB = chromB+':'+str(startB)+'-'+str(endB)
    else:
        newchrom = chrA.replace('chr','')
        regionA = newchrom+':'+str(startA)+'-'+str(endA)
        newchrom = chrB.replace('chr','')
        regionB = newchrom+':'+str(startB)+'-'+str(endB)
    return regionA, regionB

def IndexBuildingSV(onesv, chrA, chrB, posA, posB, strandA, strandB, insertion, refGenome, flank):
    
    startA = posA-flank
    endA = posA+flank
    startB = posB-flank
    endB = posB+flank
    if strandA == '+' and strandB == '+':
        refseqA = str(refGenome[chrA][startA:endA])
        refseqB = str(refGenome[chrB][startB:endB])
        refseqjoin = refseqA[:flank]+insertion+refseqB[flank-1:]
    elif strandA == '+' and strandB == '-':
        refseqA = str(refGenome[chrA][startA:endA])
        refseqB = revcmp_str(str(refGenome[chrB][startB:endB]))
        refseqjoin = refseqA[:flank]+insertion+refseqB[flank:]
    elif strandA == '-' and strandB == '+':
        refseqA = revcmp_str(str(refGenome[chrA][startA:endA]))
        refseqB = str(refGenome[chrB][startB:endB])
        refseqjoin = refseqA[:flank+1]+insertion+refseqB[flank-1:]
    elif strandA == '-' and strandB == '-':
        refseqA = revcmp_str(str(refGenome[chrA][startA:endA]))
        refseqB = revcmp_str(str(refGenome[chrB][startB:endB]))
        refseqjoin = refseqA[:flank+1]+insertion+refseqB[flank:]
    refseqAfile = onesv+'_seqA.fa'
    refseqBfile = onesv+'_seqB.fa'
    refseqjoinfile = onesv+'_join.fa'
    with open(refseqAfile,'w') as fileA, open(refseqBfile,'w') as fileB, open(refseqjoinfile,'w') as filejoin:
        fileA.write(">"+onesv+'\n'+refseqA)
        fileB.write(">"+onesv+'\n'+refseqB)
        filejoin.write(">"+onesv+'\n'+refseqjoin)
    run_bwa_genome_index_sv(onesv, "seqA")
    run_bwa_genome_index_sv(onesv, "seqB")
    run_bwa_genome_index_sv(onesv, "join")
    # return regionA, regionB

def MNVParsing(variant, sample_lst, bam_lst, chrom_lst, refGenome, outhandle, qcut=20, pcut=0.95):
    refbase = variant.split('.')[2]
    altbase = variant.split('.')[3]
    start, stop, seq = faExtract(variant, refGenome)

    for sample, bam, mychr in zip(sample_lst, bam_lst, chrom_lst):
        bamExtract(variant, bam, mychr)
        sam = variant+'.sam'
        mypos = int(variant.split('.')[1])
        good_nucl = set(["A","G","T","C"])
        valid_nucl = ["A","C","G","T","N"]

        holder = {}

        used = set([])
        readholder = {}
        with open(sam,'r') as fh:
            for line in fh:
                new = line.strip().split('\t')
                read = new[0]
                start = int(new[3])
                cigar = new[5]
                seq = new[9]
                seqlen = len(seq)
                qual = new[10]
                good = phred_read(str(seq))
                seqqual = good/len(seq)
                if seqqual < pcut:continue
                mut = call_cigar(translate_cigar(cigar))
                if 'SNV' not in mut:continue
                if read not in readholder:
                    readholder[read] = {}
                base1 = list(seq)
                qual1 = list(qual)
                pys1 = start-1
                for ii in range(pys1, pys1+seqlen):
                    myregion = range((ii-pys1),(ii-pys1+len(refbase)))
                    mybase = ''.join(base1[(ii-pys1):(ii-pys1+len(refbase))])
                    myqual = qual1[(ii-pys1):(ii-pys1+len(refbase))]
                    tqual = 0
                    for q in myqual:
                        value = ord(q)-qcut
                        tqual += value
                    if tqual/len(refbase) < 33:continue
                    if ii not in readholder[read]:
                        readholder[read][ii] = []
                    readholder[read][ii].append([mybase, tqual])
        
        for read in sorted(readholder):
            for pos in sorted(readholder[read]):
                if pos not in holder:
                    holder[pos] = {}
                if len(readholder[read][pos]) == 1:
                    mystring = readholder[read][pos][0]
                    mybase = mystring[0]
                    myqual = mystring[1]
                    if mybase not in holder[pos]:
                        holder[pos][mybase] = 0
                    holder[pos][mybase] += 1
                elif len(readholder[read][pos]) == 2:
                    mystring1 = readholder[read][pos][0]
                    mybase1 = mystring1[0]
                    myqual1 = mystring1[1]
                    ###
                    mystring2 = readholder[read][pos][1]
                    mybase2 = mystring2[0]
                    myqual2 = mystring2[1]
                    if myqual1 >= myqual2:
                        if mybase1 not in holder[pos]:
                            holder[pos][mybase1] = 0
                        holder[pos][mybase1] += 1
                    else:
                        if mybase2 not in holder[pos]:
                            holder[pos][mybase2] = 0
                        holder[pos][mybase2] += 1
        if mypos-1 not in holder:
            nR = 0
            nM = 0
        else:
            if refbase in holder[mypos-1]:
                nR = holder[mypos-1][refbase]
            else:
                nR = 0
            if altbase in holder[mypos-1]:
                nM = holder[mypos-1][altbase]
            else:
                nM = 0
        nT = nM + nR
        outhandle.write(sample+'\t'+variant+'\t'+str(nM)+'\t'+str(nT)+'\n')
    ## os.system("rm {}.sam".format(variant))

def SNVParsing(variant, sample_lst, bam_lst, chrom_lst, refGenome, outhandle):
    n = variant.split('.')
    refbase = n[2]
    altbase = n[3]
    start, stop, seq = faExtract(variant, refGenome)

    for bam, sample, mychr in zip(bam_lst, sample_lst, chrom_lst):
        bamExtract(variant, bam, mychr)
        countdict = countfromSam(variant, seq)
        nMbase = []
        nMcount = []
        nT = 0
        for x in countdict:
            if x == refbase:
                nR = countdict[x]
                nT += nR
            else:
                nMbase.append(x)
                nMcount.append(countdict[x])
                nT += countdict[x]
        if len(nMbase) > 0:
            for base, count in zip(nMbase, nMcount):
                outhandle.write(sample+'\t'+variant+'\t'+str(count)+'\t'+str(nT)+'\t'+base+'\n')
        else:
            outhandle.write(sample+'\t'+variant+'\t'+str(0)+'\t'+str(nT)+'\t'+'ND'+'\n')

def SVParsing(variant_lst, sample_lst, bam_lst, refGenome, outhandle):
    if len(set(variant_lst)) == 1:
        onesv = variant_lst[0]
        n = onesv.split('.')
        chrA = n[0]
        posA = int(n[1])
        strandA = n[2]
        chrB = n[3]
        posB = int(n[4])
        strandB = n[5]
        insertion = n[6]
        ret = "NA"
        if (chrA != chrB) or (chrA == chrB and abs(posA-posB) > 500):
            IndexBuildingSV(onesv, chrA, chrB, posA, posB, strandA, strandB, insertion, refGenome, 500)
            ret = "large"
        else:
            IndexBuildingSV(onesv, chrA, chrB, posA, posB, strandA, strandB, insertion, refGenome, 200)        
            ret = "small"
        for onesample, onebam, onesv in zip(sample_lst, bam_lst, variant_lst):
            if ret == "large":
                regionA, regionB = RegionExt(chrA, chrB, posA, posB, strandA, strandB, onebam, 500)
                mybam = onesv+'.refA.bam'
                fq = onesv+'.refA.pre.fq'
                fqsingle = onesv+'.refA.fq'
                fqhandle = open(fqsingle,'w')
                cmd = '{0} view -h -bS {1} {2} | {0} sort -n > {3}'.format(SAMTOOLS, onebam, regionA, mybam)
                os.system(cmd)
                cmd = '{} bamtofastq -i {} -fq {}'.format(BEDTOOLS, mybam, fq)
                os.system(cmd)
                ## os.system("rm {}".format(mybam))
                SeqFilterSingle(fq, onesample, fqhandle)
                fqhandle.close()
                ###
                mybam = onesv+'.refB.bam'
                fq = onesv+'.refB.pre.fq'
                fqsingle = onesv+'.refB.fq'
                fqhandle = open(fqsingle,'w')
                cmd = '{0} view -h -bS {1} {2} | {0} sort -n > {3}'.format(SAMTOOLS, onebam, regionB, mybam)
                os.system(cmd)
                cmd = '{} bamtofastq -i {} -fq {}'.format(BEDTOOLS, mybam, fq)
                os.system(cmd)
                ## os.system("rm {}".format(mybam))
                SeqFilterSingle(fq, onesample, fqhandle)
                fqhandle.close()
                ###
                cmd = "cat {0}.refA.fq {0}.refB.fq > {0}.refAB.fq".format(onesv)
                os.system(cmd)
                ###
                ## os.system("rm {0}.refA.fq {0}.refB.fq".format(onesv))
            elif ret == "small":
                regionA, regionB = RegionExt(chrA, chrB, posA, posB, strandA, strandB, onebam, 200)
                mybam = onesv+'.refA.bam'
                fq = onesv+'.refA.pre.fq'
                fqsingle = onesv+'.refA.fq'
                fqhandle = open(fqsingle,'w')
                cmd = '{0} view -h -bS {1} {2} | {0} sort -n > {3}'.format(SAMTOOLS, onebam, regionA, mybam)
                os.system(cmd)
                cmd = '{} bamtofastq -i {} -fq {}'.format(BEDTOOLS, mybam, fq)
                os.system(cmd)
                ## os.system("rm {}".format(mybam))
                SeqFilterSingle(fq, onesample, fqhandle)
                fqhandle.close()
                ###
                mybam = onesv+'.refB.bam'
                fq = onesv+'.refB.pre.fq'
                fqsingle = onesv+'.refB.fq'
                fqhandle = open(fqsingle,'w')
                cmd = '{0} view -h -bS {1} {2} | {0} sort -n > {3}'.format(SAMTOOLS, onebam, regionB, mybam)
                os.system(cmd)
                cmd = '{} bamtofastq -i {} -fq {}'.format(BEDTOOLS, mybam, fq)
                os.system(cmd)
                ## os.system("rm {}".format(mybam))
                SeqFilterSingle(fq, onesample, fqhandle)
                fqhandle.close()

                mybam = onesv+'.refAB.bam'
                fq = onesv+'.refAB.pre.fq'
                fqsingle = onesv+'.refAB.fq'
                newchrom = regionA.split(':')[0]
                startA = int(regionA.split(':')[1].split('-')[0])
                endA = int(regionA.split(':')[1].split('-')[1])
                startB = int(regionB.split(':')[1].split('-')[0])
                endB = int(regionB.split(':')[1].split('-')[1])
                regionAB = newchrom+':'+str(min([startA,endA,startB,endB]))+'-'+str(max([startA,endA,startB,endB]))
                fqhandle = open(fqsingle,'w')
                cmd = '{0} view -h -bS {1} {2} | {0} sort -n > {3}'.format(SAMTOOLS, onebam, regionAB, mybam)
                os.system(cmd)
                cmd = '{} bamtofastq -i {} -fq {}'.format(BEDTOOLS, mybam, fq)
                os.system(cmd)
                ## os.system("rm {}".format(mybam))
                SeqFilterSingle(fq, onesample, fqhandle)
                fqhandle.close()
            run_bwa_sv(onesv, onesample)
            Rd = SamParsing(onesv, onesample)
            nM, nRa, nRb, nT= countfromSamSV(Rd)
            ## os.system("rm {}*{}*".format(onesv, onesample))
            onesample1 = onebam.split('/')[-1].replace('.bam','')
            outhandle.write(onesv+'\t'+onesample1+'\t'+str(nRa)+'\t'+str(nRb)+'\t'+str(nM)+'\n')
        ## os.system("rm {}*".format(onesv))
    else:
        for onesample, onebam, onesv in zip(sample_lst, bam_lst, variant_lst):
            n = onesv.split('.')
            chrA = n[0]
            posA = int(n[1])
            strandA = n[2]
            chrB = n[3]
            posB = int(n[4])
            strandB = n[5]
            insertion = n[6]
            if (chrA != chrB) or (chrA == chrB and abs(posA-posB) > 500):
                regionA, regionB = RegionExt(chrA, chrB, posA, posB, strandA, strandB, onebam, 500)
                IndexBuildingSV(onesv, chrA, chrB, posA, posB, strandA, strandB, insertion, refGenome, 500)
                mybam = onesv+'.refA.bam'
                fq = onesv+'.refA.pre.fq'
                fqsingle = onesv+'.refA.fq'
                fqhandle = open(fqsingle,'w')
                cmd = '{0} view -h -bS {1} {2} | {0} sort -n > {3}'.format(SAMTOOLS, onebam, regionA, mybam)
                os.system(cmd)
                cmd = '{} bamtofastq -i {} -fq {}'.format(BEDTOOLS, mybam, fq)
                os.system(cmd)
                ## os.system("rm {}".format(mybam))
                SeqFilterSingle(fq, onesample, fqhandle)
                fqhandle.close()
                ###
                mybam = onesv+'.refB.bam'
                fq = onesv+'.refB.pre.fq'
                fqsingle = onesv+'.refB.fq'
                fqhandle = open(fqsingle,'w')
                cmd = '{0} view -h -bS {1} {2} | {0} sort -n > {3}'.format(SAMTOOLS, onebam, regionB, mybam)
                os.system(cmd)
                cmd = '{} bamtofastq -i {} -fq {}'.format(BEDTOOLS, mybam, fq)
                os.system(cmd)
                ## os.system("rm {}".format(mybam))
                SeqFilterSingle(fq, onesample, fqhandle)
                fqhandle.close()
                ###
                cmd = "cat {0}.refA.fq {0}.refB.fq > {0}.refAB.fq".format(onesv)
                os.system(cmd)
                ###
                ## os.system("rm {0}.refA.fq {0}.refB.fq".format(onesv))
            else:
                regionA, regionB = RegionExt(chrA, chrB, posA, posB, strandA, strandB, onebam, 200)
                IndexBuildingSV(onesv, chrA, chrB, posA, posB, strandA, strandB, insertion, refGenome, 200)
                mybam = onesv+'.refA.bam'
                fq = onesv+'.refA.pre.fq'
                fqsingle = onesv+'.refA.fq'
                fqhandle = open(fqsingle,'w')
                cmd = '{0} view -h -bS {1} {2} | {0} sort -n > {3}'.format(SAMTOOLS, onebam, regionA, mybam)
                os.system(cmd)
                cmd = '{} bamtofastq -i {} -fq {}'.format(BEDTOOLS, mybam, fq)
                os.system(cmd)
                ## os.system("rm {}".format(mybam))
                SeqFilterSingle(fq, onesample, fqhandle)
                fqhandle.close()
                ###
                mybam = onesv+'.refB.bam'
                fq = onesv+'.refB.pre.fq'
                fqsingle = onesv+'.refB.fq'
                fqhandle = open(fqsingle,'w')
                cmd = '{0} view -h -bS {1} {2} | {0} sort -n > {3}'.format(SAMTOOLS, onebam, regionB, mybam)
                os.system(cmd)
                cmd = '{} bamtofastq -i {} -fq {}'.format(BEDTOOLS, mybam, fq)
                os.system(cmd)
                ## os.system("rm {}".format(mybam))
                SeqFilterSingle(fq, onesample, fqhandle)
                fqhandle.close()

                mybam = onesv+'.refAB.bam'
                fq = onesv+'.refAB.pre.fq'
                fqsingle = onesv+'.refAB.fq'
                newchrom = regionA.split(':')[0]
                startA = int(regionA.split(':')[1].split('-')[0])
                endA = int(regionA.split(':')[1].split('-')[1])
                startB = int(regionB.split(':')[1].split('-')[0])
                endB = int(regionB.split(':')[1].split('-')[1])
                regionAB = newchrom+':'+str(min([startA,endA,startB,endB]))+'-'+str(max([startA,endA,startB,endB]))
                fqhandle = open(fqsingle,'w')
                cmd = '{0} view -h -bS {1} {2} | {0} sort -n > {3}'.format(SAMTOOLS, onebam, regionAB, mybam)
                os.system(cmd)
                cmd = '{} bamtofastq -i {} -fq {}'.format(BEDTOOLS, mybam, fq)
                os.system(cmd)
                ## os.system("rm {}".format(mybam))
                SeqFilterSingle(fq, onesample, fqhandle)
                fqhandle.close()
            run_bwa_sv(onesv, onesample)
            Rd = SamParsing(onesv, onesample)
            nM, nRa, nRb, nT= countfromSamSV(Rd)
            onesample1 = onebam.split('/')[-1].replace('.bam','')
            outhandle.write(onesv+'\t'+onesample1+'\t'+str(nRa)+'\t'+str(nRb)+'\t'+str(nM)+'\n')
            ## os.system("rm {}*".format(onesv))

def GenotypingAlltoAll(build, bamfile, filename):
    if build == "hg19":
        refGenome = Fasta(refGenome_hg19)
    elif build == "hg38":
        refGenome = Fasta(refGenome_hg38)
    sample_lst = []
    bam_lst = []
    with open(bamfile,'r') as fh:
        for line in fh:
            bam = line.strip()
            new = line.strip().split('/')
            sample = new[-1].split('/')[-1].replace('.bam','')
            sample_lst.append(sample)
            bam_lst.append(bam)
    variant_lst = []
    working_indel_lst = []
    working_mnv_lst = []
    working_snv_lst = []
    working_sv_lst = []
    with open(filename,'r') as sh:
        for line in sh:
            variant = line.strip()
            n = variant.split('.')
            variant_lst.append(variant)
            if len(n) == 4:
                if len(n[2]) > 1 or len(n[3]) > 1: # this includes both indel and MNV
                    if len(n[2]) != len(n[3]):
                        working_indel_lst.append(variant)
                    else:
                        working_mnv_lst.append(variant)
                else:
                    if len(n[2]) == 1 and len(n[3]) == 1:
                        working_snv_lst.append(variant)
            else:
                working_sv_lst.append(variant)
    ###
    sample_lst = list(sample_lst)
    variant_lst = list(variant_lst)
    bam_lst = list(bam_lst)
    print (len(sample_lst),len(variant_lst),len(bam_lst))
    ###
    if len(variant_lst) == 0:
        sys.exit("No enough valid variant for analysis!")
    else:
        chrom_lst = []
        for bam in bam_lst:
            chromstr = determineChr(bam) #
            if 'chr' in chromstr:
                chrom_lst.append('chr')
            elif 'Chr' in chromstr:
                chrom_lst.append('Chr')
            else:
                chrom_lst.append('')
        if len(working_indel_lst) > 0:
            sortindellst = IndelSorting(list(working_indel_lst))
            batchindel, batchdict = IndelBatching(sortindellst, filename)
            batchlst = BatchIndexing(build, batchindel)
        outhandle = open("counts_"+filename+'_'+bamfile+'.txt','w')
        for variant in variant_lst:
            n = variant.split('.')
            if variant in set(working_indel_lst):
                mybatch = str(batchdict[variant])
                SeqExtractor(variant, sample_lst, bam_lst, chrom_lst, build, filename)
                run_bwa(mybatch, variant, filename)
                IndelParsing(mybatch, variant, outhandle)
            elif variant in set(working_mnv_lst):
                MNVParsing(variant, sample_lst, bam_lst, chrom_lst, refGenome, outhandle)
            elif variant in set(working_snv_lst): # for SNP
                SNVParsing(variant, sample_lst, bam_lst, chrom_lst, refGenome, outhandle)
            elif variant in set(working_sv_lst): # for SNP
                full_variant_lst = [variant]*len(sample_lst)
                SVParsing(full_variant_lst, sample_lst, bam_lst, refGenome, outhandle)
        outhandle.close()
    ## if len(working_indel_lst) > 0:
        ## os.system("rm {}_*".format(filename))

def GenotypingOnetoOne(build, filename):
    if build == "hg19":
        refGenome = Fasta(refGenome_hg19)
    elif build == "hg38":
        refGenome = Fasta(refGenome_hg38)
    sample_lst = []
    bam_lst = []
    variant_lst = []
    working_indel_lst = []
    working_snv_lst = []
    working_mnv_lst = []
    working_sv_lst = []
    with open(filename,'r') as sh:
        for line in sh:
            new = line.strip().split('\t') # bam path \t variant \t sample
            bam = new[0]
            variant = new[1]
            sample = new[2]
            n = variant.split('.')
            sample_lst.append(sample)
            variant_lst.append(variant)
            bam_lst.append(bam)
            if len(n) == 4:
                if len(n[2]) > 1 or len(n[3]) > 1: # this includes both indel and MNV
                    if len(n[2]) != len(n[3]):
                        working_indel_lst.append(variant)
                    else:
                        working_mnv_lst.append(variant)
                else:
                    if len(n[2]) == 1 and len(n[3]) == 1:
                        working_snv_lst.append(variant)
            else:
                working_sv_lst.append(variant)
    ###
    if len(variant_lst) == 0:
        sys.exit("No enough valid variant for analysis!")
    else:
        chrom_lst = []
        for bam in bam_lst:
            chromstr = determineChr(bam)
            if 'chr' in chromstr:
                chrom_lst.append('chr')
            elif 'Chr' in chromstr:
                chrom_lst.append('Chr')
            else:
                chrom_lst.append('')
        outhandle = open("counts_"+filename+'.onetoone.txt','w')
        for variant, bam, sample in zip(variant_lst, bam_lst, sample_lst):
            print (variant, bam, sample)
            if variant in set(working_indel_lst):
                Indexbuilding_OnetoOne(variant, refGenome)
                SeqExtractor(variant, [sample], [bam], chrom_lst, build, filename)
                run_bwa_OnetoOne(variant)
                mybatch = ""
                IndelParsing_OnetoOne(mybatch, variant, outhandle)
            elif variant in set(working_mnv_lst): # this is a MNV
                MNVParsing(variant, [sample], [bam], chrom_lst, refGenome, outhandle)
            elif variant in set(working_snv_lst): # for SNP
                SNVParsing(variant, [sample], [bam], chrom_lst, refGenome, outhandle)
            elif variant in set(working_sv_lst): # for SNP
                SVParsing([variant], [sample], [bam], refGenome, outhandle)
        outhandle.close()
        if len(working_indel_lst) > 0:
            delIndex()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Variant Genotyping Pipeline",
                                    add_help=False)

    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                    help='-bamfile and -variantlst are used together for "all-to-all"; \n -variantfile is used for "one-to-one"')

    parser.add_argument("-action", 
        help="determining the action for genotyping, all-to-all or one-to-one",
        required=True)

    parser.add_argument("-build",
        help="the reference genome version that is used for data alignment, hg19 or hg38",
        required=True)

    parser.add_argument("-bamfile",
        help="a bamfile including the bam directory of each genotyping sample",
        required=False)

    parser.add_argument("-variantlst", 
        help="a variant file including variant per line",
        required=False)

    parser.add_argument("-variantfile", 
        help="a variant file including variant and targeted bam location per file",
        required=False)

    args = parser.parse_args()
    action = args.action
    build = args.build
    bamfile = args.bamfile
    variantlst = args.variantlst
    variantfile = args.variantfile

    if action == "all-to-all":
        GenotypingAlltoAll(build, bamfile, variantlst)
    elif action == "one-to-one":
        GenotypingOnetoOne(build, variantfile)
    else:
        sys.exit("Invalid action! Choosing all-to-all or one-to-one")
