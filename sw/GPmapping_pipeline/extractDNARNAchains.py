import os
import urllib

from gzip import GzipFile
from lxml import etree

print """
## need to download from PDB the list of PDB id's that contain DNA or RNA molecules. To do this: go to PDB -> Advanced search -> Choose a query type: set to Macromolecule type -> select Contains DNA  OR  Contains RNA to Yes-> Submit query -> Download results -> copy and paste in a new file.
"""

HBPLUS = '../sw/hbplus/hbplus'

availpdb = {i.split('.',1)[0][3:].lower():True for j in os.listdir('../dat/pdb/divided/') for i in os.listdir('../dat/pdb/divided/%s' % j) if i.endswith('.ent.gz')}  ## the PDB files we need must have been downloaded by the mapping pipeline

pids_with_dna = map(lambda x:x.strip().lower(), file('../dat/pids_with_dna.txt').readlines() + file('../dat/pids_with_nahybrid.txt').readlines())
pids_with_dna = dict(zip(pids_with_dna, [True]*len(pids_with_dna)))

pids_with_rna = map(lambda x:x.strip().lower(), file('../dat/pids_with_rna.txt').readlines())
pids_with_rna = dict(zip(pids_with_rna, [True]*len(pids_with_rna)))

nachains = {}
for p in pids_with_dna:
    if p not in availpdb:
        continue
    nachains[p] = {}
for p in pids_with_rna:
    if p not in availpdb:
        continue
    nachains[p] = {}

## download unavailable files
present = {i:True for i in os.listdir('../dat/pdbml/')}
for pid in nachains:
    fn = '%s.xml.gz' % pid.upper()
    if fn not in present:
        urllib.urlretrieve('http://www.rcsb.org/pdb/files/%s' % fn, '../dat/pdbml/%s' % fn)

for pid in nachains:
    fi = GzipFile('../dat/pdbml/%s.xml.gz' % pid.upper())
    try:
        tree = etree.fromstring(fi.read())
    except:
        print('FIXME', pid)
        raise Exception('fixme')
    for i in xrange(len(tree)):
        if tree[i].tag.split('}')[-1] == 'entity_polyCategory':
            for j in xrange(len(tree[i])):
                if tree[i][j].tag.split('}')[-1] != 'entity_poly':
                    continue
                chain = None
                etp = None
                for k in xrange(len(tree[i][j])):
                    if tree[i][j][k].tag.split('}')[-1] == 'type':
                        etp = tree[i][j][k].text
                    elif tree[i][j][k].tag.split('}')[-1] == 'pdbx_strand_id':
                        chain = tree[i][j][k].text
                if etp in ['polydeoxyribonucleotide', 'polydeoxyribonucleotide/polyribonucleotide hybrid']:
                    for x in chain.split(','):
                        x = x.strip()
                        nachains[pid][x] = ('D',{})
                elif etp == 'polyribonucleotide':
                    for x in chain.split(','):
                        x = x.strip()
                        nachains[pid][x] = ('R',{})
    fi.close()


## find out which DNA chains are paired
for pid in nachains:
    di = nachains[pid]
    r = os.system('zcat ../dat/pdb/divided/%s/pdb%s.ent.gz > ../tmp/%s.pdb && %s ../tmp/%s.pdb' % (pid[1:3],pid, pid, HBPLUS, pid))
    if r:
        print 'EE cannot get PDB file for %s or cannot run HBPLUS on it. Moving on.' % pid
        continue
    fi = file('%s.hb2' % pid)
    for i in xrange(8):
        x = fi.readline()
    while 1:
        l = fi.readline()
        if not l:
            break
        if len(l) != 76:
            raise Exception('FIX PARSER[0]')
        c1,c2 = l[0],l[14]
        if c1 == c2:
            continue
        if c1 not in di or c2 not in di:
            continue
        r1,r2 = l[6:9].strip(), l[20:23].strip()
        if c2 not in di[c1][1]:
            di[c1][1][c2] = 0
        if c1 not in di[c2][1]:
            di[c2][1][c1] = 0
        di[c1][1][c2] += 1
        di[c2][1][c1] += 1
    os.system('rm -f %s.hb2' % pid)

## create NA molecule clusters out of contact pairs
for pid in nachains:
    for r in xrange(3):
        for c1 in nachains[pid]:
            for c2 in nachains[pid][c1][1].keys():
                for c3 in nachains[pid][c2][1]:
                    if c3 != c1 and c3 not in nachains[pid][c1][1]:
                        nachains[pid][c1][1][c3] = -1

#for pid in nachains:
#    for c in nachains[pid]:
#        nachains[pid][c][1].clear()

fo = file('../dat/pids_with_dnarna.txt', 'w')
for pid in nachains:
    for x in nachains[pid]:
        t,di = nachains[pid][x]
        di = di.items()
        di.sort(key=lambda x:x[1], reverse=True)
        contact = ','.join(['%s:%d' % c for c in di])
        fo.write('%s\t%s\t%sNA\t%s\n' % (pid, x, t, contact))

fo.close()
