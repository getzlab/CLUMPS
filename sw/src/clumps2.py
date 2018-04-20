import os
import random
import scipy as sp
import sys
import time
import urllib2

from gzip import GzipFile
from lib import *
from multiprocessing import Process, Queue
from samplers.UniformSampler import *
from samplers.CoverageSampler import *
from samplers.MutspecCoverageSampler import *


PARAM = sys.argv[1]
TIMED = int(sys.argv[2])    ## whether the p-value calculations are timed: for scans through all protens/structures, TIMED should be set to a positive value (minutes per protein/structure pair); for assessing individual candidates, set TIMED=0
NTHREADS = int(sys.argv[3])
MAPSFILE = sys.argv[4]
LINEIDX = int(sys.argv[5])
TTYPE = sys.argv[6]
SAMPLER = sys.argv[7]
SAMPLEMUTFREQWEIGHT = sys.argv[8]  ## per-patient mutational load score
PATMUTSPECTRA = sys.argv[9]   ## per-patient mutational spectra (i.e. frequency of base mutations considering trinucleotide context)
COVERAGETRACK = sys.argv[10]  ## read coverage track used by some of the samplers


MODE = 1 #int(sys.argv[3])     ## if need to calculate significance of distance between two groups of mutations, expose the parameter
#if MODE not in [1,2]:
#    raise Exception('unknown mode; should be 1 or 2.')
#if abs(MODE) == 2:        ## a file containing a second group of residues whose average distance to the mutations is to be measured
#    GROUP2 = sys.argv[5]

param = file(PARAM).readlines()
for l in param:
    l = l.strip('\n')
    exec(l)

if TTYPE == 'PanCan':
    TTYPE = None

if TTYPE and PANCANFACTOR != 1.0:
    print 'WARNING: pancanfactor is not 1 althought TTYPE is set. Correcting to pancanfactor=1'
    PANCANFACTOR = 1.0

xpol = len(XPO)

## read the mapped residues
fi = GzipFile(MAPSFILE, 'r')
i = 1
while 1:
    ln = fi.readline()
    if not ln:
        raise Exception('line index out of range')
    if i == LINEIDX:
        break
    i += 1

fi.close()

SHARD = MAPSFILE.rsplit('.',1)[0].rsplit('_',1)[1]  ## chunk identifier of the mapsfile

u1,u2,pdbch,alidt,resmap = ln.strip('\n').split('\t', 4)
pdbch = pdbch.split('-')

print u1,u2,pdbch

pth = '../res/clumps/%s-%d_%s_%s_%s-%s_%s' % (SHARD, LINEIDX, u1, u2, pdbch[0], pdbch[1], resmap.split(':',1)[0])
if 0 and os.path.exists(pth):
    print 'FILE %s EXISTS ALREADY; CHANGE SCRIPT SETTINGS TO OVERRIDE. EXITING.' % pth
    exit(0)
else:
    if not os.path.exists('../res/clumps/'):
        os.makedirs('../res/clumps')
    fo = file('../res/clumps/%s-%d_%s_%s_%s-%s_%s' % (SHARD, LINEIDX, u1, u2, pdbch[0], pdbch[1], resmap.split(':',1)[0]), 'a')


ur = []  ## uniprot (u1) residues
pr = []  ## pdb residues
prd = {} ## pdb residues as dict (for containment checks)

for i in resmap.split():
    x,y = i.split(':')
    ur.append(int(x))
    pr.append(int(y))
    prd[int(y)] = True

rn = len(ur)  ## number of mapped residues (between uniprot and pdb)

if rn < 5:
    fo.write('#\n')
    fo.close()
    exit(0)

## read the residue coordinates
D,x = getResidueDistanceMatrix(pdbch, pdb_resids=pr)

## transform distance matrix
DDt = []  ## array of transformed distance matrices
for x in xrange(xpol):
    den = 2.0 * XPO[x]**2
    m = []
    for i in xrange(rn):
        mrow = sp.zeros(i, dtype=sp.float32)
        for j in xrange(i):
            mrow[j] = sp.exp(-(D[i][j]**2)/den)
        m.append(mrow)
    DDt.append(m)
    
#            s[0] += f * sp.exp(-(d**2)/18.0)   ## 18.0 = (2.0 * 3**2)
#            s[1] += f * sp.exp(-(d**2)/40.5)   ## 40.5 = (2.0 * 4.5**2)
#            s[2] += f * sp.exp(-(d**2)/72.0)   ## 72.0 = (2.0 * 6**2)
#            s[3] += f * sp.exp(-(d**2)/128.0)  ## 128.0 = (2.0 * 8**2)
#            s[4] += f * sp.exp(-(d**2)/200.0)  ## 200.0 = (2.0 * 10**2)


## get mutation data
# first get mutation frequencies of samples (for allele weighting based on that)
if SAMPLEMUTFREQWEIGHT:
    mfreq = {}
    fi = file(SAMPLEMUTFREQWEIGHT)
    l = fi.readline() ## hdr
    while 1:
        l = fi.readline()
        if not l:
            break
        l = l.strip().split('\t')
        mfreq[(l[0],l[1])] = float(l[4])

# count and weight the mutation positions
md = {}
try:
    if MUTURL.startswith('http://'):
        gm = urllib2.urlopen(MUTURL + u1).read()
    else:
        gm = file(MUTURL + u1).read()
except:
    print 'File not found: %s' % (MUTURL + u1)
    exit(0)

for l in map(lambda x:x.split('\t'), gm.split('\n')):
    if l == [''] or not l[5]:
        continue
    if l[6][0] in MUTTYPES:
        ## ttype selection
        if TTYPE and not l[0].startswith(TTYPE):
            continue
        p = int(l[5])
        if p not in md:
            md[p] = [0, set([]), set([])]
        elif USEPROVIDEDVALUES:
            print 'USEPROVIDEDVALUES is True and there are several lines for residue %d. Taking the average value.' % p
        if USEPROVIDEDVALUES:
            md[p][0] += float(l[7])
        elif SAMPLEMUTFREQWEIGHT:  ## sample-mutation-frequency based weighting of mutations
            md[p][0] += mfreq[(l[0],l[1])]
        else:  ## weight mutations equally (CLUMPS-1 approach)
            md[p][0] += 1
        md[p][1].add(l[0])
        md[p][2].add((l[1],l[0]))

if USEPROVIDEDVALUES:  ## average value
    for p in md:
        if len(md[p][2]) > 1:
            md[p][0] /= len(md[p][2])

mi = []  ## index of mutated residue
mv = []  ## normalized mutation count at each residue
mt = []  ## cancer types contributing mutations
for i in xrange(rn):
    if ur[i] in md:
        mi.append(i)
        mt.append(md[ur[i]][1])
        if USEPROVIDEDVALUES:
            mv.append(md[ur[i]][0])
        else:
            mv.append(hill(2.0, HILLEXP, md[ur[i]][0]))
        #mv.append(md[ur[i]])  ## raw
        #mv.append(sp.log(md[ur[i]]) or 0.1)

mvcorr = range(len(mv))  ## correspondence between mi and mv
cm = len(mi)   ## mutated residue count

if MODE == 1:
    Mmv = []   ## matrix that holds mv[i]*mv[j] values (sqrt or not)
    for i in xrange(cm):
        mrow = sp.zeros(cm, sp.float64)  
        for j in xrange(cm):
            #mrow[j] = sp.sqrt(mv[i]*mv[j])  ## geometric mean; actually does not perform better in most cases
            if PANCANFACTOR == 1.0:
                mrow[j] = mv[i]*mv[j]          ## product
            else:
                mrow[j] = (PANCANFACTOR + (1.0-PANCANFACTOR)*(len(mt[i] & mt[j])>0)) * mv[i]*mv[j]          ## product
        Mmv.append(mrow)

    def wap(mi,mvcorr):
        s = sp.zeros(len(DDt), sp.float64)
        for mat in xrange(xpol):
            d = DDt[mat]
            for i in xrange(cm):
                dcol = d[mi[i]]
                for j in xrange(i):
                    s[mat] += Mmv[mvcorr[i]][mvcorr[j]] * dcol[mi[j]]
        return s

elif MODE == 2:
    nd = [x.strip().split('\t') for x in file(GROUP2).readlines()]
    nd = {int(x[0]):float(x[1]) for x in nd}
    
    mi2 = []  ## index of mutated residue
    mv2 = []  ## normalized mutation count at each residue
    #mt2 = []  ## cancer types contributing mutations
    for i in xrange(rn):
        if ur[i] in nd:
            mi2.append(i)
            mv2.append(nd[ur[i]])
            #mt2.append(md[ur[i]][1])

    mvcorr = range(len(mv))  ## correspondence between mi and mv
    cm = len(mi)    ## mutated residue count
    cm2 = len(mi2)   ## group2 residue count

    Mmv = []   ## matrix that holds mv[i]*mv[j] values (sqrt or not)
    for i in xrange(cm):
        mrow = sp.zeros(cm2, sp.float64)  
        for j in xrange(cm2):
            #mrow[j] = sp.sqrt(mv[i]*mv[j])  ## geometric mean; actually does not perform better in most cases
            if PANCANFACTOR == 1.0:
                mrow[j] = mv[i]*mv2[j]          ## product
            else:
                raise 'Not implemented'
                #mrow[j] = (PANCANFACTOR + (1.0-PANCANFACTOR)*(len(mt[i] & mt[j])>0)) * mv[i]*mv[j]          ## product
        Mmv.append(mrow)

    def wap(mi,mvcorr):
        s = sp.zeros(len(DDt), sp.float64)
        for mat in xrange(xpol):
            d = DDt[mat]
            for i in xrange(cm):
                for j in xrange(cm2):
                    a,b = sorted([mi[i],mi2[j]])
                    if a == b:
                        s[mat] += Mmv[mvcorr[i]][j] * 1
                    else:
                        s[mat] += Mmv[mvcorr[i]][j] * d[b][a]
        return s

## get the observed score
wap_obs = wap(mi,mvcorr)

## create the null model
rnd = 0
P = [0]*xpol
WAP_RND = [0]*xpol
mireal = [i for i in mi]
if SAMPLER == 'UniformSampler':
    sam = UniformSampler(ur)
elif SAMPLER == 'CoverageSampler':
    sam = CoverageSampler(ur, u1, COVERAGETRACK)
elif SAMPLER == 'MutspecCoverageSampler':
    sam = MutspecCoverageSampler(ur, u1, COVERAGETRACK, PATMUTSPECTRA)
    sam.calcMutSpecProbs(md)

if TIMED:
    STARTTIME=time.time()

def rndThread(qu):
    def booster():
        """ Implements the booster algorithm by Getz et al. that saves CPU time.
            Returns False if the randomization should be stopped.
        """
        ret = False
        for i in xrange(xpol):
            s = (rnd - p[i] + 1.0) / ((rnd + 3)*(p[i] + 1))
            if s >= 0.05271:  ## =0.9/(2*1.96)]**2:
                ret = True
                break
        return ret
    sp.random.seed()
    p = [0]*xpol
    wap_rnd = [0]*xpol
    rnd = 0
    exitstatus=0  ## 0 means terminated OK, 1 means it had to abort due to timeout
    while rnd < MAXRAND/NTHREADS and (rnd%1000 or booster()):  ## booster is applied once per 1000 randomizations
        if not rnd%1000 and TIMED and (time.time()-STARTTIME)/60.0 > TIMED:
            exitstatus=1
            break
        x = None
        while x is None:
            x = sam.sample(mireal)  ## some samplers will fail to yield a sample in some (small number of) of runs due to combinatorics
        mi,mvcorr = x
        r = wap(mi,mvcorr)
        for rr in xrange(xpol):
            wap_rnd[rr] += r[rr]
            if r[rr] >= wap_obs[rr]:
                p[rr] += 1
        rnd += 1
    qu.put((rnd,p,wap_rnd,exitstatus))

queue = Queue()
pcs = []
for r in xrange(NTHREADS):
    x = Process(target=rndThread, args=(queue,))
    x.start()
    pcs.append(x)

for x in pcs:
    x.join()

time.sleep(2)

## get results from subprocesses
if queue.qsize() != NTHREADS:
    raise Exception('something went wrong with the randomizations.')

totalrnd = 0
totalexitstatus = 0
for r in xrange(NTHREADS):
    rnd,p,wap_rnd,exitstatus = queue.get()
    totalrnd += rnd
    totalexitstatus += exitstatus
    for i in xrange(xpol):
        P[i] += p[i]
        WAP_RND[i] += wap_rnd[i]

fo.write('\t'.join(['%d/%d' % (P[i], totalrnd) for i in xrange(len(P))]) + '\n')
fo.write('#%d\n' % totalexitstatus)
fo.close()
