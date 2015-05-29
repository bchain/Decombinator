####################################################################################################
# Build CDR3 files for each .freq file
# 
####################################################################################################

from __future__ import division
import os
import ast
import collections as coll
import string
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

genepath = '/Volumes/chain/TcRSeq_Processed/Unileverdata/KB_analysis_scripts/'
indir='/Volumes/chain/TcRSeq_Processed/Unileverdata/freq/'
outdir='/Volumes/chain/TcRSeq_Processed/Unileverdata/cdr3s/'

####################################################################################################
# Functions to find CDR3s
def import_genes(genepath, chain):
    """
    Imports V and J genes from genepath, as well as importing the position of the conserved
    cysteine in the translated V region

    Returns lists of Vregions, Jregions, VconservedC
    """
    vfilename = os.path.join(genepath, 'exthuman_TR'+string.upper(chain)[0]+'V_region.fasta')
    vregions = [str(string.upper(item.seq))
                for item in list(SeqIO.parse(vfilename, 'fasta'))]

    jfilename = os.path.join(genepath, 'exthuman_TR'+string.upper(chain)[0]+'J_region.fasta')
    jregions = [str(string.upper(item.seq))
                for item in list(SeqIO.parse(jfilename, 'fasta'))]

    cysteinefilename = os.path.join(genepath,
                                    'TR'+string.upper(chain)[0]+'V_ConservedC.txt')
    conserved_cysteine = [int(x) for x in
                          open(cysteinefilename).read().splitlines()]

    phefilename = os.path.join(genepath,
                               'TR'+string.upper(chain)[0]+'J_ConservedF.txt')
    conserved_phe = [int(x) for x in
                     open(phefilename).read().splitlines()]

    return vregions, jregions, conserved_cysteine, conserved_phe

def find_nt(dcr, vregions, jregions):

    elems = dcr.split(',')
    v = int(elems[0])
    j = int(elems[1])
    vdel = int(elems[2])
    jdel = int(elems[3])
    ins = elems[4]

    if vdel != 0:
        vused = vregions[v][:-vdel]
    else:
        vused = vregions[v]
    jused = jregions[j][jdel:]

    nt = ''.join([vused, ins, jused])

    return nt

def find_cdr3(nt, conserved_cysteine, conserved_phe):

    aa = str(Seq(nt, generic_dna).translate())

    if aa[conserved_cysteine-1] != 'C':
        return 0
    if aa[-conserved_phe] != 'F' and aa[-conserved_phe] != 'W':
        return 0

    cdr3 = aa[(conserved_cysteine-1):(-conserved_phe+4)]
    return cdr3

def check_functionality(nt, conserved_cysteine, conserved_phe):
    # check for in-frame
    if (len(nt)-1) % 3 != 0:
        return 0

    # check for CDR3
    elif find_cdr3(nt, conserved_cysteine, conserved_phe) == 0:
        return 0
    
    # check for stop codon in CDR3
    elif '*' in find_cdr3(nt, conserved_cysteine, conserved_phe):
        return 0

    else:
        return 1

alpha_vregions, alpha_jregions, alpha_conservedC, alpha_conservedF = import_genes(genepath, 'alpha')
beta_vregions, beta_jregions, beta_conservedC, beta_conservedF = import_genes(genepath, 'beta')
####################################################################################################

####################################################################################################
print 'read dcr data in'
infiles = [x for x in os.listdir(indir) if x.endswith('.freq')]

dcr_counters = list()
sample_ids = list()
for i, f in enumerate(infiles):
    print i, f
    sample_ids.append('_'.join(f.rstrip('.freq').split('_')[2:]))
    dcr_counters.append(coll.Counter())
    inhandle = open(indir + f, 'r')
    for line in inhandle:
        bits = line.rstrip('\n').split(',')
        dcr = ','.join(bits[:5])
        freq = int(bits[5])
        dcr_counters[i][dcr] += freq
    inhandle.close()
####################################################################################################

####################################################################################################
print 'translate dcrs to cdr3s'

cdr3_counters = list()

for i, c in enumerate(dcr_counters):
    print i, sample_ids[i]
    cdr3_dict = coll.Counter()
    if sample_ids[i].endswith('a'):
        vr = alpha_vregions
        jr = alpha_jregions
        cc = alpha_conservedC
        cf = alpha_conservedF
    elif sample_ids[i].endswith('b'):
        vr = beta_vregions
        jr = beta_jregions
        cc = beta_conservedC
        cf = beta_conservedF
    for dcr, freq in c.iteritems():
        v = int(dcr.split(',')[0])
        j = int(dcr.split(',')[1])
        nt = find_nt(dcr, vr, jr)
        if not check_functionality(nt, cc[v], cf[j]):
            continue
        cdr3 = find_cdr3(nt, cc[v], cf[j])
        cdr3_dict[cdr3] += freq
    cdr3_counters.append(cdr3_dict)
####################################################################################################

####################################################################################################
print 'save cdr3 information to file'

for i, c in enumerate(cdr3_counters):
    print i, sample_ids[i]
    outfile = open(outdir + sample_ids[i] + '_cdr3s.txt', 'w')
    for cdr3, freq in c.most_common():
        print >> outfile, ','.join([cdr3, str(freq)])
    outfile.close()




