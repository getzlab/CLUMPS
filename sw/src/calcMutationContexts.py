import sys

from gzip import GzipFile
from twobitreader import TwoBitFile


MAFNM = sys.argv[1]
if MAFNM.endswith('.gz'):
    maf = GzipFile(MAFNM)
else:
    maf = file(MAFNM)
hdr = '#'
while hdr.startswith('#'):
    hdr = maf.readline()

hdr = hdr.strip().split('\t')
iPat = hdr.index('Tumor_Sample_Barcode')
iClass = hdr.index('Variant_Type')
iChr = hdr.index('Chromosome')
iPos = hdr.index('Start_Position')
iRefA = hdr.index('Reference_Allele')
iNewbase = hdr.index('Tumor_Seq_Allele2')
iTtype = hdr.index('ttype')

idx = {}
i = 0
for b1 in 'acgt':
    for b3 in 'acgt':
        for b2 in 'ac':
            for n in 'acgt':
                if n == b2:
                    continue
                idx[(b1+b2+b3, n)] = i
                i += 1

compl = {'a':'t', 'c':'g', 'g':'c', 't':'a'}
def reverse_complement(abc):
    return compl[abc[2]] + compl[abc[1]] + compl[abc[0]]

def encode(abc, n):
    if abc[1] == 'g' or abc[1] == 't':
        abc = reverse_complement(abc)
        n = compl[n]
    return idx[(abc,n)]

hg = TwoBitFile('hg19.2bit')    

pat2vec = {}

while 1:
    l = maf.readline()
    if not l:
        break
    l = l.strip('\n').split('\t')
    if l[iClass] != 'SNP':
        continue
    ref = l[iRefA].lower()
    newb = l[iNewbase]
    if len(ref) != 1 or ref == '-':
        continue
    if len(newb) != 1 or newb == '-':
        continue
    pos = int(l[iPos])
    chr = l[iChr]
    if chr == '23':
        chr = 'X'
    elif chr == '24':
        chr = 'Y'
    elif chr == 'MT':
        chr = 'M'
    abc = hg['chr'+chr][pos-2:pos+1].lower()
    if abc[1] != ref and ref != '--':
        print abc, ref, l
        print 'non-matching reference.'
        continue
    pat = (l[iPat][:12], l[iTtype]) # l[iTtype].split('-')[0]
    if pat not in pat2vec:
        pat2vec[pat] = [0]*96
    try:
        pat2vec[pat][encode(abc,l[iNewbase].lower())] += 1
    except:  ## because of Ns
        continue

fo = file('sampleMutSpectra.txt', 'w')
hdr = idx.items()
hdr.sort(key=lambda x:x[1])
fo.write('\t'.join(['patient', 'ttype'] + [i[0][0]+'-'+i[0][1] for i in hdr]) + '\n')
for p in pat2vec:
    l = [p[0], p[1]] + map(str,pat2vec[p])
    fo.write('\t'.join(l) + '\n')

fo.close()
