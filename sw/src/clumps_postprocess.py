import os
import sys
import urllib2

from gzip import GzipFile
from statsmodels.sandbox.stats.multicomp import multipletests


print "TAKING MIN INSTEAD OF MEDIAN P-VALUE PER CLUSTER"

PARAM = sys.argv[1]
RESDIR = sys.argv[2]
CANGENES = sys.argv[3]
HUNIPROT2PDBMAP = sys.argv[4]
TTYPE = sys.argv[5]

fi = file(PARAM)
for l in fi.readlines():
    l = l.strip('\n')
    exec(l)

if TTYPE == 'PanCan':
    TTYPE = None


"""
def getFragmentAnnotation(pdb, ch):
    #Get the PDB 'fragment' annotation for pdb chain ch
    poly = parsePDBHeader('../dat/pdb/%s.pdb.gz' % pdb, 'polymers')
    for p in poly:
        if p.chid != ch:
            continue
        return p.fragment
    return ''
"""

## parse gene IDs and symbols
mapdi = {}
fi = urllib2.urlopen('http://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:9606&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length,database(geneid)&format=tab')
hdr = fi.readline().strip().split('\t')
iEntry = hdr.index('Entry')
iLength = hdr.index('Length')
iGn = hdr.index('Gene names')
iEntrez = hdr.index('Cross-reference (geneid)')
while 1:
    l = fi.readline()
    if not l:
        break
    l = l.strip('\n').split('\t')
    mapdi[l[iEntry]] = [map(lambda x:x.strip(';').strip('.'), l[iEntrez].split()), l[iGn], l[iLength]]
    

## get list of all known or candidate cancer genes
cangenes = {}
for l in file(CANGENES).readlines()[1:]:
    l = l.strip('\n').split('\t')
    cangenes[l[0]] = l[2]  # cangenes[l[1]] = l[2]  ## for gene symbols

canproteins = {}
for u1 in mapdi:
    for e in mapdi[u1][0]:
        if e in cangenes:
            canproteins[u1] = cangenes[e]

## get structure coverage start and end coordinates
struct2startend = {}
struct2covsamples = {}
prot2muts = {}
fi = GzipFile(HUNIPROT2PDBMAP)
while 1:
    l = fi.readline()
    if not l:
        break
    l = l.strip('\n').split('\t',4)
    if l[0] not in prot2muts:
        di = {}
        try:
            fmut = file('/sw/dat/splitByProtein/%s' % l[0])
        except:
            print 'WW mutation file not found - probably a bug'
            continue
        for ml in fmut.readlines():
            ml = ml.strip().split('\t')
            if ml[6] != 'M':
                continue
            if TTYPE and ml[0] != TTYPE:
                continue
            sample = (ml[1], ml[0])
            rs = int(ml[5])
            if rs not in di:
                di[rs] = set([])
            di[rs].add(sample)
        prot2muts[l[0]] = di
    muts = prot2muts[l[0]]
    covmuts = set([])  ## covered mutation sites
    covsamples = set([])  ## samples contributing covered mutations
    for r in [int(i.split(':')[0]) for i in l[4].split()]:
        if r in muts:
            covmuts.add(r)
            covsamples.update(muts[r])
    s = l[4].split(' ',1)[0].split(':')[0]
    e = l[4].rsplit(' ',1)[-1].split(':')[0]
    iden = 100.0
    if l[3] != '-':
        for i in l[3].split(' '):
            if i.startswith('pdb_identity'):
                iden = float(i.split(':')[1])
    struct2startend[(l[0],l[1],l[2],s)] = [int(s), int(e), iden, len(covmuts)]
    struct2covsamples[(l[0],l[1],l[2],s)] = covsamples
    
u1structs = {}

for fn in os.listdir(RESDIR):
    fileline,u1,u2,pdbch,rs = fn.split('_')
#    if u1 == 'P42226':
#        print '###'
    if u1 not in u1structs:
        u1structs[u1] = []
    try:
        dat = file(RESDIR + '/' + fn).readlines()
        if not dat or dat == ['#\n']:  ## empty file: not tested dur to some filter
            continue
        if dat[-1] != '#0\n':
            print 'WW: potentially insufficient number of simulations for %s' % fn
        dim = len(dat[0].split('\t'))  ## p-value array length
        P = [[0,0] for x in xrange(dim)]
        for i in xrange(0,len(dat),2):
            l = dat[i].strip().split('\t')
            for j in xrange(dim):
                e,d = map(int, l[j].split('/'))  ## enumerator, denominator
                P[j][0] += e
                P[j][1] += d
        P = [x[0]/float(x[1]) for x in P]  ## now it's an array of p-values
    except:
        print 'EE', fn
        continue
    if (u1,u2,pdbch,rs) not in struct2startend:
#        if u1 == 'P42226':
#            print 'NOT'
        continue
    s,e,iden,covmuts = struct2startend[(u1,u2,pdbch,rs)]
    if covmuts < 3:  ## filter based on number of covered mutated residues
#        if u1 == 'P42226':
#            print '##CM', covmuts
        continue
    covsamples = len(struct2covsamples[(u1,u2,pdbch,rs)])
    u1structs[u1].append([u2,pdbch,s,e,iden,P,covmuts,covsamples])

## select structures according to heuristic
u1structs_filt = {}
for u1 in u1structs:
    ss = u1structs[u1]
    ss.sort(key=lambda x:x[3]-x[2], reverse=True)  ## key is coverage length
    u1structs_filt[u1] = []
    clusters = []  ## clusters of homologous structures
    for i in ss:
        newcl = True  ## does this found a new cluster
        for cl in clusters:
            num_overlapping = 0  ## how many members of the cluster cl does this struct overlap with (at least 90% jaccard overlap)
            for j in cl:
                ol = max(0, min(i[3], j[3]) - max(i[2], j[2]))
                if float(ol)/((i[3]-i[2]) + (j[3]-j[2]) - ol) > 0.9:
                    num_overlapping += 1
            if num_overlapping/float(len(cl)) >= 0.5:  ## overlaps with at least half of the members of the cluster by at least 90%
                cl.append(i)
                newcl = False
                break
        if newcl:
            clusters.append([i])
    for cl in clusters:
        #cl.sort(key=lambda x:x[4]*(x[3]-x[2]), reverse=True)  ## sort by identity * mapped length
        cl.sort(key=lambda x:x[5][2])  ## sort by p-value (6A)
        #j = cl[int(round(len(cl)/2.0))-1] ## median p
        j = cl[0] ## min p
        ok = True
        for i in u1structs_filt[u1]:
            ol = max(0, min(i[3], j[3]) - max(i[2], j[2]))
            if float(ol)/min(i[3]-i[2], j[3]-j[2]) > 0.1:
                ok = False
        if ok:
            u1structs_filt[u1].append(j)
    #print '######################', u1
    #print u1structs[u1]
    #print '######################'
    #print u1structs_filt[u1]
    #print '######################'

## select the structs for cancer genes
#u1structs_filt_cangenes = {u1:u1structs_filt[u1] for u1 in u1structs_filt if u1 in canproteins}

## pre-format for output
outdata = []
for u1 in u1structs_filt:
    for l in u1structs_filt[u1]:
        cancerannot = ''
        if u1 in canproteins:
            cancerannot = canproteins[u1]
        pdb,ch = l[1].split('-')
        if u1 in mapdi:
            if mapdi[u1][1]:
                gename = mapdi[u1][1].split()[0]
            else:
                gename = ''
            protlen = mapdi[u1][2]
            present = 1
        else:
            gename = ''
            protlen = ''
            present = 0
        if present and not gename: # or genename.startswith('HLA-'):
            print "filtering out %s due to missing gene name (Ig?)" % u1
            continue
        if l[7] < CLUMPSMUTSAMPLEFILT:  ## filter on the number of samples
            continue
        l = {'UNIPROT_ID':u1, 
             'GENE_NAMES':gename, 
             'IN_CANCER_GENE_LISTS':cancerannot,
             'MAPPED_UNIPROT_ID':l[0],
             'PDBID-CHAIN':l[1], 
             'PDB_FRAGMENT':'NA', #getFragmentAnnotation(pdb, ch), 
             'UNIPROT_SEQ_LENGTH':str(protlen),
             'MAP_START':str(l[2]), 
             'MAP_END':str(l[3]), 
             'PERCENT_IDENTITY':'%.1f' % l[4], 
             'NSITES':str(l[6]),
             'NSAMPLES':str(l[7]),
             'CLUMPS_P':l[5][2], 
             'CLUMPS_Q_FULL':-1, 
             'CLUMPS_Q_RESTRICTED':-1}
        outdata.append(l)

## correct for multiple testing
pvals = [i['CLUMPS_P'] for i in outdata]
#qvals = list(R.r['p.adjust'](R.FloatVector(pvals), "fdr"))
qvals = multipletests(pvals, method='fdr_bh')[1]
for i in xrange(len(outdata)):
    outdata[i]['CLUMPS_Q_FULL'] = qvals[i]

## just for cancer proteins
pvals = [i['CLUMPS_P'] for i in outdata if i['IN_CANCER_GENE_LISTS']]
idx = [i for i in xrange(len(outdata)) if outdata[i]['IN_CANCER_GENE_LISTS']]
#qvals = list(R.r['p.adjust'](R.FloatVector(pvals), "fdr"))
qvals = multipletests(pvals, method='fdr_bh')[1]
for i in xrange(len(idx)):
    outdata[idx[i]]['CLUMPS_Q_RESTRICTED'] = qvals[i]

outdata.sort(key=lambda x:x['CLUMPS_P'])  ## sorted by p-value

fo = file('clumps.summary.txt', 'w')
hdr = ['UNIPROT_ID', 'GENE_NAMES', 'IN_CANCER_GENE_LISTS', 'MAPPED_UNIPROT_ID', 'PDBID-CHAIN', 'PDB_FRAGMENT', 'UNIPROT_SEQ_LENGTH', 'MAP_START', 'MAP_END', 'PERCENT_IDENTITY', 'NSITES', 'NSAMPLES', 'CLUMPS_P', 'CLUMPS_Q_FULL', 'CLUMPS_Q_RESTRICTED']
fo.write('\t'.join(hdr) + '\n')
for i in xrange(len(outdata)):
    l = outdata[i]
    ll = [l[i] for i in hdr[:-3]] + ['%g' % l[i] for i in hdr[-3:]]
    fo.write('\t'.join(ll) + '\n')

fo.close()