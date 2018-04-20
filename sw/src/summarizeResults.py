import os
import sys

from statsmodels.sandbox.stats.multicomp import multipletests


RESDIR = sys.argv[1]


out = []
for fn in os.listdir(RESDIR):
    dat = file(RESDIR + '/' + fn).readlines()
    if not dat or dat == ['#\n']:  ## empty file: not tested e.g. due to 
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
    out.append((fn, P, []))

dim = len(out[0][1])  ## p-value array length
for i in xrange(dim):
    pp = [x[1][i] for x in out]
    qq = multipletests(pp, method='fdr_bh')[1]
    for j in xrange(len(out)):
        out[j][2].append(qq[j])

if dim >= 3:
    out.sort(key=lambda x:x[1][2])

fo = file('resultsSummary.txt', 'w')
hdr = ['protein_pdbchain']
for i in xrange(dim):
    hdr.append('P[%d]' % i)
    hdr.append('Q[%d]' % i)

fo.write('\t'.join(hdr) + '\n')

for r in out:
    line = [r[0]]
    for i in xrange(dim):
        line.append('%g' % r[1][i])
        line.append('%g' % r[2][i])
    fo.write('\t'.join(line) + '\n')

fo.close()