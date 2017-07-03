#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:30:37 2017

@author: igorhut
"""

#FASTA and FASTQ formats#
#########################
#Problem 1.
#Given ucsc.hg19.fasta file determine the following:
# 1. Number of contigs (chromosomes)
# 2. Number of bases
# 3. Number of bases for each contig
#Extract chromosome 20 (chr20) in separate FASTA file and validate 
#created FASTA file by running an alignment (given the folowing reads: 
#    mate1, mate2).

## Try to code this not using python modules other then built-in ones

filename = '/Users/igorhut/Documents/Onboarding_failed_tasks/Projects/253b1ef0-de03-4a2c-95f4-4d8fc2963aaa/ucsc.hg19.fasta'

#Count contigs
def countContigs(filename):
    num_contigs = 0
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == '>':
                num_contigs += 1
    return num_contigs

num_cont = countContigs(filename)
print('There are ' + str(num_cont) + ' contigs.') # so there are 93 contigs

# Read in the whole genome
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
genome = readGenome(filename)
genome[:100]

# And now count the bases
import collections
all_bases = collections.Counter(genome)
all_bases

#Counter({'A': 425987911,
#         'C': 288802645,
#         'G': 289047188,
#         'N': 239850717,
#         'T': 426549515,
#         'a': 428975238,
#         'c': 304164079,
#         'g': 304278040,
#         'n': 85,
#         't': 429505846})
    
# Number of all bases:
sum(all_bases.values()) #3137161264
    
# Count bases by contigs (using dictionary)
def countBasesByContig(filename):
    cnt = {}
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == '>':
                cnt_line = line                
                cnt[cnt_line] = 0
                continue
            cnt[cnt_line] += len(line)                
    return cnt

dict_num_by_contig = countBasesByContig(filename)
dict_num_by_contig

# Count bases by contigs (using Counter)
def countBasesByContig(filename):
    cnt = collections.Counter()
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == '>':
                cnt_line = line                
                cnt[cnt_line] = 0
                continue
            cnt[cnt_line] += len(line)-1                
    return cnt

bases_by_contig = countBasesByContig(filename)
bases_by_contig

#Counter({'>chr1\n': 249250621,
#         '>chr10\n': 135534747,
#         '>chr11\n': 135006516,
#         '>chr11_gl000202_random\n': 40103,
#         '>chr12\n': 133851895,
#         '>chr13\n': 115169878,
#         '>chr14\n': 107349540,
#         '>chr15\n': 102531392,
#         '>chr16\n': 90354753,
#         '>chr17\n': 81195210,
#         '>chr17_ctg5_hap1\n': 1680828,
#         '>chr17_gl000203_random\n': 37498,
#         '>chr17_gl000204_random\n': 81310,
#         '>chr17_gl000205_random\n': 174588,
#         '>chr17_gl000206_random\n': 41001,
#         '>chr18\n': 78077248,
#         '>chr18_gl000207_random\n': 4262,
#         '>chr19\n': 59128983,
#         '>chr19_gl000208_random\n': 92689,
#         '>chr19_gl000209_random\n': 159169,
#         '>chr1_gl000191_random\n': 106433,
#         '>chr1_gl000192_random\n': 547496,
#         '>chr2\n': 243199373,
#         '>chr20\n': 63025520,
#         '>chr21\n': 48129895,
#         '>chr21_gl000210_random\n': 27682,
#         '>chr22\n': 51304566,
#         '>chr3\n': 198022430,
#         '>chr4\n': 191154276,
#         '>chr4_ctg9_hap1\n': 590426,
#         '>chr4_gl000193_random\n': 189789,
#         '>chr4_gl000194_random\n': 191469,
#         '>chr5\n': 180915260,
#         '>chr6\n': 171115067,
#         '>chr6_apd_hap1\n': 4622290,
#         '>chr6_cox_hap2\n': 4795371,
#         '>chr6_dbb_hap3\n': 4610396,
#         '>chr6_mann_hap4\n': 4683263,
#         '>chr6_mcf_hap5\n': 4833398,
#         '>chr6_qbl_hap6\n': 4611984,
#         '>chr6_ssto_hap7\n': 4928567,
#         '>chr7\n': 159138663,
#         '>chr7_gl000195_random\n': 182896,
#         '>chr8\n': 146364022,
#         '>chr8_gl000196_random\n': 38914,
#         '>chr8_gl000197_random\n': 37175,
#         '>chr9\n': 141213431,
#         '>chr9_gl000198_random\n': 90085,
#         '>chr9_gl000199_random\n': 169874,
#         '>chr9_gl000200_random\n': 187035,
#         '>chr9_gl000201_random\n': 36148,
#         '>chrM\n': 16571,
#         '>chrUn_gl000211\n': 166566,
#         '>chrUn_gl000212\n': 186858,
#         '>chrUn_gl000213\n': 164239,
#         '>chrUn_gl000214\n': 137718,
#         '>chrUn_gl000215\n': 172545,
#         '>chrUn_gl000216\n': 172294,
#         '>chrUn_gl000217\n': 172149,
#         '>chrUn_gl000218\n': 161147,
#         '>chrUn_gl000219\n': 179198,
#         '>chrUn_gl000220\n': 161802,
#         '>chrUn_gl000221\n': 155397,
#         '>chrUn_gl000222\n': 186861,
#         '>chrUn_gl000223\n': 180455,
#         '>chrUn_gl000224\n': 179693,
#         '>chrUn_gl000225\n': 211173,
#         '>chrUn_gl000226\n': 15008,
#         '>chrUn_gl000227\n': 128374,
#         '>chrUn_gl000228\n': 129120,
#         '>chrUn_gl000229\n': 19913,
#         '>chrUn_gl000230\n': 43691,
#         '>chrUn_gl000231\n': 27386,
#         '>chrUn_gl000232\n': 40652,
#         '>chrUn_gl000233\n': 45941,
#         '>chrUn_gl000234\n': 40531,
#         '>chrUn_gl000235\n': 34474,
#         '>chrUn_gl000236\n': 41934,
#         '>chrUn_gl000237\n': 45867,
#         '>chrUn_gl000238\n': 39939,
#         '>chrUn_gl000239\n': 33824,
#         '>chrUn_gl000240\n': 41933,
#         '>chrUn_gl000241\n': 42152,
#         '>chrUn_gl000242\n': 43523,
#         '>chrUn_gl000243\n': 43341,
#         '>chrUn_gl000244\n': 39929,
#         '>chrUn_gl000245\n': 36651,
#         '>chrUn_gl000246\n': 38154,
#         '>chrUn_gl000247\n': 36422,
#         '>chrUn_gl000248\n': 39786,
#         '>chrUn_gl000249\n': 38502,
#         '>chrX\n': 155270560,
#         '>chrY\n': 59373566})
    
sum(bases_by_contig.values()) #3137161264 - bombona!

#Extract chromosome 20 (chr20) in separate FASTA file 

parse = False
with open(filename, 'r') as f1:
    with open('chr20.fasta', 'a') as f2:
        for line in f1:
            if line.startswith('>chr20'):
                parse = True
            elif line.startswith('>chr21'):
                parse = False
            if parse:
                f2.write(line)
                
!ls               
# Count # of bases in chr20.fasta to check whether all is fine
chr_20 = readGenome('chr20.fasta')
import collections
chr_20_bases = collections.Counter(chr_20)
chr_20_bases
sum(chr_20_bases.values()) #63025520 - Valja!



# And now repeat everything with 'pysam'#
#########################################
import pysam as ps

# Fetch names of reference sequences
fasta = ps.Fastafile(filename)
fasta.references

# Tuple with the lengths of reference sequences
fasta.lengths

#Make a dictionary with keys = references and values = lengths
ref_len_dict = dict(zip(fasta.references, fasta.lengths))
ref_len_dict # This gives num of bases per contig

#{'chr1': 249250621,
# 'chr10': 135534747,
# 'chr11': 135006516,
# 'chr11_gl000202_random': 40103,
# 'chr12': 133851895,
# 'chr13': 115169878,
# 'chr14': 107349540,
# 'chr15': 102531392,
# 'chr16': 90354753,
# 'chr17': 81195210,
# 'chr17_ctg5_hap1': 1680828,
# 'chr17_gl000203_random': 37498,
# 'chr17_gl000204_random': 81310,
# 'chr17_gl000205_random': 174588,
# 'chr17_gl000206_random': 41001,
# 'chr18': 78077248,
# 'chr18_gl000207_random': 4262,
# 'chr19': 59128983,
# 'chr19_gl000208_random': 92689,
# 'chr19_gl000209_random': 159169,
# 'chr1_gl000191_random': 106433,
# 'chr1_gl000192_random': 547496,
# 'chr2': 243199373,
# 'chr20': 63025520,
# 'chr21': 48129895,
# 'chr21_gl000210_random': 27682,
# 'chr22': 51304566,
# 'chr3': 198022430,
# 'chr4': 191154276,
# 'chr4_ctg9_hap1': 590426,
# 'chr4_gl000193_random': 189789,
# 'chr4_gl000194_random': 191469,
# 'chr5': 180915260,
# 'chr6': 171115067,
# 'chr6_apd_hap1': 4622290,
# 'chr6_cox_hap2': 4795371,
# 'chr6_dbb_hap3': 4610396,
# 'chr6_mann_hap4': 4683263,
# 'chr6_mcf_hap5': 4833398,
# 'chr6_qbl_hap6': 4611984,
# 'chr6_ssto_hap7': 4928567,
# 'chr7': 159138663,
# 'chr7_gl000195_random': 182896,
# 'chr8': 146364022,
# 'chr8_gl000196_random': 38914,
# 'chr8_gl000197_random': 37175,
# 'chr9': 141213431,
# 'chr9_gl000198_random': 90085,
# 'chr9_gl000199_random': 169874,
# 'chr9_gl000200_random': 187035,
# 'chr9_gl000201_random': 36148,
# 'chrM': 16571,
# 'chrUn_gl000211': 166566,
# 'chrUn_gl000212': 186858,
# 'chrUn_gl000213': 164239,
# 'chrUn_gl000214': 137718,
# 'chrUn_gl000215': 172545,
# 'chrUn_gl000216': 172294,
# 'chrUn_gl000217': 172149,
# 'chrUn_gl000218': 161147,
# 'chrUn_gl000219': 179198,
# 'chrUn_gl000220': 161802,
# 'chrUn_gl000221': 155397,
# 'chrUn_gl000222': 186861,
# 'chrUn_gl000223': 180455,
# 'chrUn_gl000224': 179693,
# 'chrUn_gl000225': 211173,
# 'chrUn_gl000226': 15008,
# 'chrUn_gl000227': 128374,
# 'chrUn_gl000228': 129120,
# 'chrUn_gl000229': 19913,
# 'chrUn_gl000230': 43691,
# 'chrUn_gl000231': 27386,
# 'chrUn_gl000232': 40652,
# 'chrUn_gl000233': 45941,
# 'chrUn_gl000234': 40531,
# 'chrUn_gl000235': 34474,
# 'chrUn_gl000236': 41934,
# 'chrUn_gl000237': 45867,
# 'chrUn_gl000238': 39939,
# 'chrUn_gl000239': 33824,
# 'chrUn_gl000240': 41933,
# 'chrUn_gl000241': 42152,
# 'chrUn_gl000242': 43523,
# 'chrUn_gl000243': 43341,
# 'chrUn_gl000244': 39929,
# 'chrUn_gl000245': 36651,
# 'chrUn_gl000246': 38154,
# 'chrUn_gl000247': 36422,
# 'chrUn_gl000248': 39786,
# 'chrUn_gl000249': 38502,
# 'chrX': 155270560,
# 'chrY': 59373566}

# Num of contigs
len(ref_len_dict) #93

#Num of bases
sum(ref_len_dict.values()) # 3137161264 - odlican

#Extract chromosome 20 (chr20) into separate FASTA file 
f = open('chr_20_pysam.fasta', 'wt')
f.write(fasta.fetch(reference = 'chr20'))
fasta.close()



############################################
#Heng Li's function for parsing FASTA/FASTQ#
############################################
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

#Problem 4.
#Given CCLE sample (mate1, mate2), discard all reads that contain at least one base 
#with quality of 10 (Sanger quality scale). This assumes creating new FASTQ files 
#which contain only valid reads. Check if you can run alignment successfully with 
#new files.

from Bio import SeqIO

#First let's count the # of reads
count = 0
for rec in SeqIO.parse("/Users/igorhut/Documents/Onboarding_failed_tasks/Projects/253b1ef0-de03-4a2c-95f4-4d8fc2963aaa/G20479.HCC1143.2.converted.pe_1_1Mreads.fastq", "fastq"):
    count += 1
print( "%i reads" % count) #1000000 reads

# And now let's trim reads which have that contain at least one base with quality of 
# 10 (Sanger quality scale)
good_reads_1 = (rec for rec in \
              SeqIO.parse("/Users/igorhut/Documents/Onboarding_failed_tasks/Projects/253b1ef0-de03-4a2c-95f4-4d8fc2963aaa/G20479.HCC1143.2.converted.pe_1_1Mreads.fastq", "fastq") \
              if min(rec.letter_annotations["phred_quality"]) >= 10)
count_1 = SeqIO.write(good_reads_1, "good_quality_pe_1.fastq", "fastq")

print("Saved %i reads" % count_1) # Saved 738538 reads

good_reads_2 = (rec for rec in \
              SeqIO.parse("/Users/igorhut/Documents/Onboarding_failed_tasks/Projects/253b1ef0-de03-4a2c-95f4-4d8fc2963aaa/G20479.HCC1143.2.converted.pe_2_1Mreads.fastq", "fastq") \
              if min(rec.letter_annotations["phred_quality"]) >= 10)
count_2 = SeqIO.write(good_reads_2, "good_quality_pe_2.fastq", "fastq")

print("Saved %i reads" % count_2) # Saved 606352 reads 


# A bit of introductory playing with 'pysam'

import pysam as ps

fname_1 = '/Users/igorhut/Documents/GitHub/bds-files/chapter-11-alignment/NA12891_CEU_sample.bam'
bamfile_1 = ps.AlignmentFile(filename = fname_1, mode = 'rb')

type(bamfile_1)

#Fetch aligned reads from a particular region of an indexed BAM file
for read in bamfile_1.fetch('1', start=215906528, end=215906567):
    print(read.qname, "aligned at position", read.pos)
    
for read in bamfile_1:
    status = "unaligned" if read.is_unmapped else "aligned"
    print(read.qname, "is", status)
 
bamfile_1.close()

# Let's repeat this for homework assignment file
fname = '/Users/igorhut/Documents/Onboarding_failed_tasks/Projects/253b1ef0-de03-4a2c-95f4-4d8fc2963aaa/NA12878_NA12878-Garvan-Vial1_R.realigned.base_recalibrated.merged.bam'
bamfile = ps.AlignmentFile(filename = fname, mode = 'rb')

bam_dict = {}
for read in bamfile.fetch('20', start = 20327462, end = 20327573):
    print(read.qname, 'aligned at position', read.pos)
    bam_dict[read.qname] = read.pos
    
#ST-E00185:49:H5LVWCCXX:5:1201:9779:63438 aligned at position 20327324
#ST-E00185:49:H5LVWCCXX:5:2205:12439:8692 aligned at position 20327325
#ST-E00185:49:H5LVWCCXX:5:1116:23878:51501 aligned at position 20327329
#ST-E00185:49:H5LVWCCXX:5:2116:21797:46068 aligned at position 20327329
#ST-E00185:49:H5LVWCCXX:5:1209:9993:48336 aligned at position 20327330
#ST-E00185:49:H5LVWCCXX:5:2213:21553:38649 aligned at position 20327333
#ST-E00185:49:H5LVWCCXX:5:1220:14885:64651 aligned at position 20327339
#ST-E00185:49:H5LVWCCXX:5:1114:4522:48529 aligned at position 20327345
#ST-E00185:49:H5LVWCCXX:5:2217:8237:27468 aligned at position 20327349
#ST-E00185:49:H5LVWCCXX:5:1115:29876:19909 aligned at position 20327351
#ST-E00185:49:H5LVWCCXX:5:2108:15646:73072 aligned at position 20327351
#ST-E00185:49:H5LVWCCXX:5:1109:21726:22915 aligned at position 20327355
#ST-E00185:49:H5LVWCCXX:5:1109:21807:22880 aligned at position 20327355
#ST-E00185:49:H5LVWCCXX:5:1108:14763:24673 aligned at position 20327356
#ST-E00185:49:H5LVWCCXX:5:1108:14773:24691 aligned at position 20327356
#ST-E00185:49:H5LVWCCXX:5:1108:14804:24603 aligned at position 20327356
#ST-E00185:49:H5LVWCCXX:5:2110:23908:49022 aligned at position 20327357
#ST-E00185:49:H5LVWCCXX:5:1120:26233:60924 aligned at position 20327362
#ST-E00185:49:H5LVWCCXX:5:1213:15788:2469 aligned at position 20327363
#ST-E00185:49:H5LVWCCXX:5:2201:28973:13175 aligned at position 20327373
#ST-E00185:49:H5LVWCCXX:5:2212:6359:29174 aligned at position 20327378
#ST-E00185:49:H5LVWCCXX:5:2208:7222:25570 aligned at position 20327382
#ST-E00185:49:H5LVWCCXX:5:1203:31937:14019 aligned at position 20327386
#ST-E00185:49:H5LVWCCXX:5:2208:6247:50182 aligned at position 20327388
#ST-E00185:49:H5LVWCCXX:5:1102:7749:46912 aligned at position 20327389
#ST-E00185:49:H5LVWCCXX:5:1119:4075:37348 aligned at position 20327389
#ST-E00185:49:H5LVWCCXX:5:2210:19584:11224 aligned at position 20327391
#ST-E00185:49:H5LVWCCXX:5:2110:16895:5581 aligned at position 20327393
#ST-E00185:49:H5LVWCCXX:5:1203:31937:14019 aligned at position 20327393
#ST-E00185:49:H5LVWCCXX:5:1116:23878:51501 aligned at position 20327396
#ST-E00185:49:H5LVWCCXX:5:2221:12215:38157 aligned at position 20327396
#ST-E00185:49:H5LVWCCXX:5:1112:29257:10011 aligned at position 20327398
#ST-E00185:49:H5LVWCCXX:5:2107:3669:47580 aligned at position 20327400
#ST-E00185:49:H5LVWCCXX:5:1103:9840:47369 aligned at position 20327402
#ST-E00185:49:H5LVWCCXX:5:2117:31592:44328 aligned at position 20327406
#ST-E00185:49:H5LVWCCXX:5:1219:18001:58304 aligned at position 20327408
#ST-E00185:49:H5LVWCCXX:5:2123:18965:6636 aligned at position 20327408
#ST-E00185:49:H5LVWCCXX:5:2208:6613:35766 aligned at position 20327409
#ST-E00185:49:H5LVWCCXX:5:1111:7617:63877 aligned at position 20327410
#ST-E00185:49:H5LVWCCXX:5:1111:7891:63965 aligned at position 20327410
#ST-E00185:49:H5LVWCCXX:5:1111:7902:63983 aligned at position 20327410
#ST-E00185:49:H5LVWCCXX:5:1221:23472:32022 aligned at position 20327414
#ST-E00185:49:H5LVWCCXX:5:1117:14733:49303 aligned at position 20327424
#ST-E00185:49:H5LVWCCXX:5:1103:29714:61820 aligned at position 20327429
#ST-E00185:49:H5LVWCCXX:5:1115:29876:19909 aligned at position 20327433
#ST-E00185:49:H5LVWCCXX:5:2204:12916:36909 aligned at position 20327446
#ST-E00185:49:H5LVWCCXX:5:2217:29562:63596 aligned at position 20327447
#ST-E00185:49:H5LVWCCXX:5:2106:17737:34571 aligned at position 20327451
#ST-E00185:49:H5LVWCCXX:5:1105:18316:48758 aligned at position 20327457
#ST-E00185:49:H5LVWCCXX:5:2213:21553:38649 aligned at position 20327458
#ST-E00185:49:H5LVWCCXX:5:1209:6785:21966 aligned at position 20327460
#ST-E00185:49:H5LVWCCXX:5:2219:11068:32936 aligned at position 20327463
#ST-E00185:49:H5LVWCCXX:5:2201:28973:13175 aligned at position 20327466
#ST-E00185:49:H5LVWCCXX:5:2220:10490:63438 aligned at position 20327467
#ST-E00185:49:H5LVWCCXX:5:2218:28384:15742 aligned at position 20327472
#ST-E00185:49:H5LVWCCXX:5:2101:26365:13193 aligned at position 20327475
#ST-E00185:49:H5LVWCCXX:5:1203:10429:11505 aligned at position 20327476
#ST-E00185:49:H5LVWCCXX:5:2107:4633:59025 aligned at position 20327476
#ST-E00185:49:H5LVWCCXX:5:1102:22913:25288 aligned at position 20327486
#ST-E00185:49:H5LVWCCXX:5:2219:14936:53030 aligned at position 20327487
#ST-E00185:49:H5LVWCCXX:5:2219:15057:51764 aligned at position 20327487
#ST-E00185:49:H5LVWCCXX:5:1213:15788:2469 aligned at position 20327490
#ST-E00185:49:H5LVWCCXX:5:1108:28080:23829 aligned at position 20327492
#ST-E00185:49:H5LVWCCXX:5:1113:9810:13457 aligned at position 20327501
#ST-E00185:49:H5LVWCCXX:5:1205:18214:19891 aligned at position 20327504
#ST-E00185:49:H5LVWCCXX:5:1111:10246:10029 aligned at position 20327505
#ST-E00185:49:H5LVWCCXX:5:1209:25492:63227 aligned at position 20327506
#ST-E00185:49:H5LVWCCXX:5:1103:24477:23530 aligned at position 20327507
#ST-E00185:49:H5LVWCCXX:5:2122:15930:31547 aligned at position 20327509
#ST-E00185:49:H5LVWCCXX:5:1207:3243:17729 aligned at position 20327512
#ST-E00185:49:H5LVWCCXX:5:1109:6572:47264 aligned at position 20327517
#ST-E00185:49:H5LVWCCXX:5:1203:20833:26009 aligned at position 20327518
#ST-E00185:49:H5LVWCCXX:5:2210:19584:11224 aligned at position 20327519
#ST-E00185:49:H5LVWCCXX:5:2205:12439:8692 aligned at position 20327527
#ST-E00185:49:H5LVWCCXX:5:2205:12520:8798 aligned at position 20327527
#ST-E00185:49:H5LVWCCXX:5:1108:14763:24673 aligned at position 20327529
#ST-E00185:49:H5LVWCCXX:5:1108:14773:24691 aligned at position 20327529
#ST-E00185:49:H5LVWCCXX:5:1108:14804:24603 aligned at position 20327529
#ST-E00185:49:H5LVWCCXX:5:1110:19402:14072 aligned at position 20327532
#ST-E00185:49:H5LVWCCXX:5:1117:14733:49303 aligned at position 20327532
#ST-E00185:49:H5LVWCCXX:5:2219:11068:32936 aligned at position 20327532
#ST-E00185:49:H5LVWCCXX:5:1209:26212:31424 aligned at position 20327540
#ST-E00185:49:H5LVWCCXX:5:1115:19645:65231 aligned at position 20327541
#ST-E00185:49:H5LVWCCXX:5:1221:25847:39194 aligned at position 20327542
#ST-E00185:49:H5LVWCCXX:5:2109:29552:23038 aligned at position 20327543
#ST-E00185:49:H5LVWCCXX:5:2218:2816:2082 aligned at position 20327552
#ST-E00185:49:H5LVWCCXX:5:1102:7749:46912 aligned at position 20327552
#ST-E00185:49:H5LVWCCXX:5:1216:13007:47404 aligned at position 20327556
#ST-E00185:49:H5LVWCCXX:5:1223:2410:60836 aligned at position 20327556
#ST-E00185:49:H5LVWCCXX:5:1214:17443:62928 aligned at position 20327558
#ST-E00185:49:H5LVWCCXX:5:2110:14895:12947 aligned at position 20327565
#ST-E00185:49:H5LVWCCXX:5:1109:28080:50270 aligned at position 20327566
#ST-E00185:49:H5LVWCCXX:5:1109:28161:50164 aligned at position 20327566
#ST-E00185:49:H5LVWCCXX:5:1219:18001:58304 aligned at position 20327566
#ST-E00185:49:H5LVWCCXX:5:1217:22020:70224 aligned at position 20327569

print(bam_dict['ST-E00185:49:H5LVWCCXX:5:1111:10246:10029']) #20327505

bamfile.reset()

single_read = next(bamfile)
print(single_read.query_name)
print(single_read.reference_start)
print(single_read.query_sequence)
print(single_read.query_length)
print(single_read.is_paired)
print(single_read.is_proper_pair)

#ST-E00185:49:H5LVWCCXX:5:1223:30110:1397
#9997
#CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAGAGCCTAAGCCTAACCCTCACCCCCACTCCAACCCCACGCATGCCGCTAACCGCATCCCAAACGATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC
#150
#True
#False



fname = '/Users/igorhut/Documents/Onboarding_failed_tasks/Projects/253b1ef0-de03-4a2c-95f4-4d8fc2963aaa/G20479.HCC1143.2.converted.pe.Aligned.out.sorted.bam'
bamfile = ps.AlignmentFile(filename = fname, mode = 'rb')

for read in bamfile.fetch():
    if read.is_secondary:
        print(read)
        break

# C0VHUACXX120727:5:2212:10271:38999      355     0       10516   0       101M    0       10539   101     CTCTGCTCCGCCTTCGCAGTACCACCGAAATCTGTGCAGAGGAGAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAACTCC   array('B', [34, 34, 34, 37, 37, 37, 37, 37, 39, 39, 39, 39, 39, 41, 41, 41, 41, 41, 41, 38, 40, 41, 41, 41, 41, 41, 41, 41, 40, 41, 41, 41, 41, 41, 38, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 41, 40, 41, 41, 41, 39, 39, 39, 39, 37, 37, 35, 35, 35, 35, 35, 35, 35, 31, 35, 35, 34, 34, 35, 35, 35, 35, 35, 33, 27, 31, 34, 35, 34, 35, 35, 35, 35, 34, 35, 33, 34, 33, 35, 35, 35, 35, 35, 35, 33, 35, 35, 34, 34])        [('NH', 5), ('HI', 3), ('AS', 182), ('nM', 8), ('RG', '1')]

bamfile.reset()
counter = 0
for read in bamfile.fetch():
    if read.query_name == "C0VHUACXX120727:5:2212:10271:38999":
        print(read)
        counter += 1
print("Number of secondary alignments for read C0VHUACXX120727:5:2212:10271:38999 is : ", str(counter/2 - 1))


#Problem 10.
#How many soft clipped reads are there in the range from Problem 1?

fname = '/Users/igorhut/Documents/Onboarding_failed_tasks/Projects/253b1ef0-de03-4a2c-95f4-4d8fc2963aaa/NA12878_NA12878-Garvan-Vial1_R.realigned.base_recalibrated.merged.bam'
bamfile = ps.AlignmentFile(filename = fname, mode = 'rb')
counter = 0
for read in bamfile.fetch('20', start = 20327462, end = 20327573):
    if 'S' in read.cigarstring:
        counter += 1
        print(read.cigarstring)

print('There are ', str(counter), 'soft clipped reads in this region.')


# Iterate through all reads of RNA-seq BAM file in order to find a read which spans exons.
# What is the length of the intervening intron?

import pysam as ps
fname = '/Users/igorhut/Documents/Onboarding_failed_tasks/Projects/253b1ef0-de03-4a2c-95f4-4d8fc2963aaa/G20479.HCC1143.2.converted.pe.Aligned.out.sorted.bam'
bamfile = ps.AlignmentFile(filename = fname, mode = 'rb')

for read in bamfile.fetch():
    if 'N' in read.cigarstring:
        print(read)
        print(read.cigarstring)
        break
bamfile.close()
    
# C0VHUACXX120727:7:1306:16072:8082       355     0       11585   0       86M338N15M      0       12128   101     GTACTGTTCTGTATCCCACCAGCAATGTCTAGGAATACCTGTTTCTCCACAAAGTGTTTACTTTTGGATTTTTGCCAGTCTAACAGGTGTCTGACTTCCAG   array('B', [34, 31, 31, 37, 37, 37, 37, 37, 39, 39, 39, 39, 39, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 40, 41, 40, 41, 41, 40, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 40, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 37, 39, 38, 40, 40, 41, 41, 41, 41, 41, 41, 41, 40, 40, 41, 41, 41, 41, 41, 41, 41, 39, 39, 39, 39, 39, 39, 37, 37, 37, 37, 37, 35, 34, 36, 34, 36, 36, 35, 34, 35, 35, 35, 35, 35, 35, 34])        [('NH', 5), ('HI', 3), ('AS', 196), ('nM', 2), ('RG', '1')]
# 86M338N15M
# So the length of intervening intron is 338 bases

# Iterate through all reads of RNA-seq BAM file in order to find a read which spans more than 2 exons.
# What are the lengths of the introns? Can you tell what’s the length of each exon?
bamfile.reset()

import pysam as ps
fname = '/Users/igorhut/Documents/Onboarding_failed_tasks/Projects/253b1ef0-de03-4a2c-95f4-4d8fc2963aaa/G20479.HCC1143.2.converted.pe.Aligned.out.sorted.bam'
bamfile = ps.AlignmentFile(filename = fname, mode = 'rb')

for read in bamfile.fetch():
    if read.cigarstring.count('N') >= 2:
        print(read.query_name, ': ', read.cigarstring)
    break

bamfile.close()  


#Problem 8.
#How many PCR duplicates are there in the whole RNA-seq sample?
#What percentage of reads are PCR duplicates? What does this number tell you?
#How long does the execution of this command last? You can use %%time magic function to do this.

import pysam as ps
import time
fname = '/Users/igorhut/Documents/Onboarding_failed_tasks/Projects/253b1ef0-de03-4a2c-95f4-4d8fc2963aaa/G20479.HCC1143.2.converted.pe.Aligned.out.sorted.bam'
bamfile = ps.AlignmentFile(filename = fname, mode = 'rb')

counter_reads = 0
counter_dup = 0
start = time.time()
for read in bamfile.fetch():
    counter_reads += 1
    if read.is_duplicate:
        counter_dup += 1
end = time.time()

print('Counting lasted for total of ', str(end-start), 'seconds.')
print('There are ', str(counter_dup), ' PCR duplicates in the whole RNA-seq sample.')
print(str(counter_dup/counter_reads),'% of reads are PCR duplcate reads.')
bamfile.close()

#Counting lasted for total of  3760.9297931194305 seconds.
#There are  0  PCR duplicates in the whole RNA-seq sample.
#0.0 % of reads are PCR duplcate reads.

import pysam as ps
import time
fname = '/Users/igorhut/Documents/Onboarding_failed_tasks/Projects/253b1ef0-de03-4a2c-95f4-4d8fc2963aaa/G20479.HCC1143.2.converted.pe.Aligned.out.sorted.deduped.bam'
bamfile = ps.AlignmentFile(filename = fname, mode = 'rb')

counter_reads = 0
counter_dup = 0
start = time.time()
for read in bamfile:
    counter_reads += 1
    if read.is_duplicate:
        counter_dup += 1
end = time.time()

print('Counting lasted for total of ', str(end-start), 'seconds.')
print('There are ', str(counter_dup), ' PCR duplicates in the whole RNA-seq sample.')
print(str(counter_dup/counter_reads),'% of reads are PCR duplcate reads.')
bamfile.close()

#Counting lasted for total of  2711.204339027405 seconds.
#There are  34515194  PCR duplicates in the whole RNA-seq sample.
#0.22335784019599 % of reads are PCR duplcate reads.