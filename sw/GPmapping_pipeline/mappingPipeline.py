"""
Maps uniprot to PDB positions
"""

import ftplib
import os
import sys
import time

from gzip import GzipFile
from lib import *
from prody import parsePDBStream

HOMOLOGY = True  ## whether homology modeling should be applied
#BLASTED_SIFTS = '../dat/uniprot.human.sp/'
BLASTED_SIFTS = sys.argv[1]

if BLASTED_SIFTS[-1] != '/':
    BLASTED_SIFTS += '/'


allhumupids = [i.split('.',1)[0] for i in os.listdir(BLASTED_SIFTS) if i.endswith('.seq.blasted.gz')]
print 'seeking', len(allhumupids), 'maps'

###################
## GET SIFTS ID MAPS
###################
siftsIdMaps = {}
fi = file('../dat/sifts/pdb_chain_uniprot.tsv')
x = fi.readline()
x = fi.readline()
while 1:
    l = fi.readline()
    if not l:
        break
    l = l.strip('\n').split('\t')
    if l[2] not in siftsIdMaps:
        siftsIdMaps[l[2]] = {}
    siftsIdMaps[l[2]][(l[0], l[1])] = True

siftsResidMaps = {}  ## uniprot 2 pdb residue maps

###################
## GET AVAILABLE ALIGNMENTS
###################
alignedUpids = {}
if HOMOLOGY:
    for fn in os.listdir(BLASTED_SIFTS):
        if fn.endswith('.seq.blasted.gz'):
            alignedUpids[fn.split('.',1)[0]] = True

print len(alignedUpids), 'total human uniprot IDs have been aligned'

def getAlignments(uniprotId, evalue_threshold=0.001, relIdentity_threshold=0.15):
    """  """
    try:
        fi = GzipFile(BLASTED_SIFTS + '%s.seq.blasted.gz' % uniprotId)
    except:
        return {}
    ## parse the alignment xml file
    results = []
    cur = etree.fromstring(fi.read())
    qlen = 0  ## query length
    for i in cur:
        if i.tag == 'BlastOutput_iterations':
            cur = i
        elif i.tag == 'BlastOutput_query-len':
            qlen = int(i.text)
    cur = cur[-1]  ## last iteration
    for i in cur:
        if i.tag == 'Iteration_hits':
            cur = i
            break
    for hit in cur:
        hit_def = None
        hsps = None
        for i in hit:
            if i.tag == 'Hit_def':
                hit_def = i.text
            elif i.tag == 'Hit_hsps':
                hsps = i
        for hsp in hsps:
            di = {}
            di['Hit_def'] = hit_def
            for i in hsp:
                di[i.tag] = i.text
            if float(di['Hsp_evalue']) > evalue_threshold:  ### !!! FILTER !!!
                continue
            relIdentity = float(di['Hsp_identity'])/float(di['Hsp_align-len'])
            if relIdentity < relIdentity_threshold:         ### !!! FILTER !!!
                continue
            di['relIdentity'] = relIdentity
            results.append(di)
    return results

alignments = {}

###################
## GET DIRECT AND INDIRECT MAPS
###################
u1u2pdbch = {}  ## human uniprot 2 homologous uniprots 2 pdb-chain

## direct maps
for u1 in allhumupids:
    if u1 not in siftsIdMaps:
        continue
    ## otherwise, there are direct maps
    u1u2pdbch[u1] = {u1:{}}
    for k in siftsIdMaps[u1]:
        u1u2pdbch[u1][u1][k] = []

## indirect maps
if HOMOLOGY:
    for u1 in allhumupids:
        if u1 not in alignedUpids:
            continue
        alis = getAlignments(u1)
        if not alis:
            continue
        if u1 not in alignments:
            alignments[u1] = {}
        au1 = alignments[u1]
        if u1 not in u1u2pdbch:
            u1u2pdbch[u1] = {}
        for ali in alis:
            u2 = ali['Hit_def'].split('|',2)[1]
            if u2 not in au1:
                au1[u2] = []
            au1[u2].append(ali)
            if u2 not in u1u2pdbch[u1]:
                u1u2pdbch[u1][u2] = {}
            if u2 not in siftsIdMaps:
                print u2, 'not in sifts'
                continue
            for k in siftsIdMaps[u2]:
                u1u2pdbch[u1][u2][k] = []


## figure out which residues are mapped

allpdb = {}
for u1 in u1u2pdbch:
    for u2 in u1u2pdbch[u1]:
        for k in u1u2pdbch[u1][u2]:
            allpdb[k[0]] = 1

print 'total pdb files that have maps from mutated proteins:', len(allpdb)

availpdb = {i[3:7]:True for j in os.listdir('../dat/pdb/divided') for i in os.listdir('../dat/pdb/divided/%s' % j) if i.endswith('.ent.gz')}

## download pdb files if necessary
if 0:
    print 'available are %d pdb files' % len(availpdb)
    cnt = 0  ## how many new pdb files have been downloaded
    na = 0   ## how many structures (from sifts) are not present in pdb
    for i in range(200):
        try:
            ftp = ftplib.FTP("ftp.wwpdb.org")
            break
        except:
            time.sleep(5)
    ftp.login()
    ftp.cwd('/pub/pdb/data/structures/all/pdb')
    for p in allpdb:
        if p not in availpdb:
            fo = file('../dat/pdb/divided/%s/pdb%s.ent.gz' % (p[1:3], p), 'w')
            try:
                ftp.retrbinary('RETR pdb%s.ent.gz' % p, fo.write)
                fo.close()
                cnt += 1
                availpdb[p] = True
                print('downloaded', p)
            except:
                print "STRUCTURE FILE NOT AVAILABLE", p
                fo.close()
                os.system('touch ../dat/pdb/divided/%s/pdb%s.ent.gz' % (p[1:3], p))  ## touch the empty file
                os.remove('../dat/pdb/divided/%s/pdb%s.ent.gz' % (p[1:3], p))  ## delete the empty file
                na += 1
    print 'got %d new pdb files; not available: %d' % (cnt, na)
    ftp.close()

"""
## check which sifts residue map files are available
print 'CHECKME CHECKME CHECKME'
availsifts = {}
for dir in os.listdir('../dat/sifts/residmaps/'):
    dl = {i.split('.',1)[0]:True for i in os.listdir('../dat/sifts/residmaps/%s/' % dir) if i.endswith('.xml.gz')}
    availsifts.update(dl)

print 'available are %d sifts xml files' % len(availsifts)
"""

print 'len(u1u2pdbch)', len(u1u2pdbch)

## write output
present = {}
if 1:  ## overwrite!!!
    #fo = GzipFile('../res/huniprot2pdb.run18.txt.gz_%s' % time.time(), 'w')
    fo = GzipFile('../res/huniprot2pdb.run18.%s.txt.gz' % (BLASTED_SIFTS.strip('/').split('/')[-1]), 'w')
else:  ## append
    fi = GzipFile('../res/huniprot2pdb.run18.txt.gz', 'r')
    while 1:
        l = fi.readline()
        if not l:
            break
        u1,u2 = l.split('\t',2)[:2]
        if u1 not in present:
            present[u1] = {}
        present[u1][u2] = True
    fi.close()
    print 'appending to list already containing %d proteins' % len(present)
    fo = GzipFile('../res/huniprot2pdb.run18.txt.gz', 'a')

for u1 in u1u2pdbch:
    for u2 in u1u2pdbch[u1]:
        if u1 in present and u2 in present[u1]:
            continue
        for pdbch in u1u2pdbch[u1][u2]:
            if pdbch[0] not in availpdb:
                continue
            if u2 in siftsResidMaps and pdbch in siftsResidMaps[u2]:
                srm = siftsResidMaps[u2][pdbch]
            else:
                if u2 not in siftsResidMaps:
                    siftsResidMaps[u2] = {}
                srm = getSiftsUniprot2PdbResidueMap(pdbch[0], pdbch[1], u2)
                if not srm:
                    continue
                ## get the residue numbers from the pdb structure 
                ## because not all residues mapped in sifts are present
                try:
                    aa = parsePDBStream(GzipFile('../dat/pdb/divided/%s/pdb%s.ent.gz' % (pdbch[0][1:3], pdbch[0])), chain=pdbch[1])
                    d = set([])
                    resnames = aa.getResnames()
                    resnums = aa.getResnums()
                    for i in xrange(len(resnames)):
                        if resnames[i] in aamap:
                            d.add(resnums[i])
                    aa = d
                except:  ## not all structures are available in .pdb format
                    print 'could not parse %s.' % pdbch[0]
                    continue
                srm = {i:srm[i] for i in srm if srm[i] in aa}
                siftsResidMaps[u2][pdbch] = srm
            if u1 == u2:  ## native structure
                ss = srm.items()
                if not ss:
                    continue
                ss.sort()
                l = [u1,
                     u2,
                     '%s-%s' % pdbch,
                     '-',
                     ' '.join(map(lambda x:'%d:%d' % x, ss))]
                fo.write('\t'.join(l) + '\n')
            else:  ## homology structure
                for ali in alignments[u1][u2]:
                    srm_u1 = {}
                    query_ind = int(ali['Hsp_query-from'])
                    hit_ind = int(ali['Hsp_hit-from'])
                    qa = ali['Hsp_qseq']
                    ha = ali['Hsp_hseq']
                    pdbhitovl = [-1,-1]  ##
                    for i in xrange(len(ha)):
                        if qa[i] == '-':
                            hit_ind += 1
                            continue
                        if ha[i] == '-':
                            query_ind += 1
                            continue
                        else:  ## not a gap
                            if hit_ind in srm:
                                srm_u1[query_ind] = srm[hit_ind]
                                if pdbhitovl[0] < 0:
                                    pdbhitovl[0] = i
                                pdbhitovl[1] = i
                            query_ind += 1
                            hit_ind += 1
                    ss = srm_u1.items()
                    if len(ss) < 3:  ## less than 3 mapped residues
                        continue
                    ## calculate relative identity on the alignment part overlapping with pdb seq
                    pdb_iden = sum([qa[i]==ha[i] for i in xrange(pdbhitovl[0], pdbhitovl[1]+1)]) / float(pdbhitovl[1]-pdbhitovl[0]+1)
                    if pdb_iden < 0.20:
                        continue
                    ss.sort()
                    l = [u1,
                         u2,
                         '%s-%s' % pdbch,
                         'uhit_identity:%.4g evalue:%.1g pdb_identity:%.4g' % (100.0*ali['relIdentity'], float(ali['Hsp_evalue']), 100.0*pdb_iden),
                         ' '.join(map(lambda x:'%d:%d' % x, ss))]
                    fo.write('\t'.join(l) + '\n')

fo.close()
