import ftplib
import os
import random
import re
import scipy as sp
import urllib
import urllib2

from gzip import GzipFile
from lxml import etree
from prody import parsePDBStream, parsePDBHeader
from scipy.spatial.distance import euclidean


#############################################
#SIFTS_IDMAP = '../dat/sifts/pdb_chain_uniprot.tsv'
#SIFTS_RESIDMAPS = '../dat/sifts/residmaps/'
# BLASTED_SIFTS = '../dat/uniprot.human.sp/'  ## xml.gz files with alignments of all human proteins to custom SIFTS-based database
PDB_STRUCTURES = '../dat/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/'
#############################################


aamap = {
'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','PYX':'C','SEP':'S','TPO':'T','TYS':'Y','MK8':'L','M3L':'K','DPR':'P','DSN':'S','ALN':'A','DLY':'K','MVA':'V','MLE':'L','DLE':'L','CSO':'C','PTR':'Y','BMT':'T','DAL':'A','FAK':'K','MSE':'M','SMF':'A','HYP':'P'}


def hill(t,n,x):
    return x**n/(t**n + x**n)


def booster(x, r):
    """ implements the booster algorithm by Getz et al.
        returns False if the randomization should be stopped.
        x is the number of non-rejected null hypotheses; r is the number of randomizations.
    """
    s = (r - x + 1.0) / ((r + 3)*(x + 1))
    if s >= 0.05271: # =0.9/(2*1.96)]**2:
        return True
    else:
        return False


def getResidueDistanceMatrix(pdbch, point='centroid', pdb_resids=None, return_centroid_coordinates=False):
    """ Returns a distance matrix summarizing the Euclidean 
        distances between residues in structure pdbch """
    if return_centroid_coordinates and point != 'centroid':
        raise Exception('return_centroid_coordinates is True but point argument is not set to centroid')

    if pdbch[1]:
        aa = parsePDBStream(GzipFile('%s/%s/pdb%s.ent.gz' % (PDB_STRUCTURES, pdbch[0][1:3], pdbch[0])), chain=pdbch[1])  ## gets ALL atoms
    else:
        aa = parsePDBStream(GzipFile('%s/%s/pdb%s.ent.gz' % (PDB_STRUCTURES, pdbch[0][1:3], pdbch[0])))  ## gets ALL atoms
    xx = aa.getResnums()
    yy = aa.getCoords()
    zz = aa.getResnames()

    if pdb_resids is None:
        pdb_resids = {}
        for i in xrange(len(xx)):
            if zz[i] in aamap:
                pdb_resids[xx[i]] = True
        pdb_resids = sorted(pdb_resids.keys())

    coords = {}
    for i in xrange(len(xx)):
        if xx[i] not in pdb_resids:
            continue
        if xx[i] not in coords:
            coords[xx[i]] = []
        coords[xx[i]].append(yy[i])  ## add coordinates of an atom belonging to this residue

    rn = len(pdb_resids)  ## number of residues
    
    ## Euclidean distance matrix
    #D = sp.zeros((rn,rn))
    D = []
    for i in xrange(rn):
        D.append(sp.zeros(i, dtype=sp.float32))
    
    if point == 'centroid':
        ## distance between centroids
        ## calculate residue centroid positions
        centroids = {}
        for k in coords:
            cc = coords[k]
            centroids[k] = [sp.mean(map(lambda x:x[0],cc)), sp.mean(map(lambda x:x[1],cc)), sp.mean(map(lambda x:x[2],cc))]
        co = [centroids[i] for i in pdb_resids]  ## pdb residue coordinates
        for i in xrange(rn):
            for j in xrange(i):
                D[i][j] = euclidean(co[i], co[j])
    elif point == 'min':
        ## min-distance (atom pairs)
        co = [coords[i] for i in pdb_resids]     ## pdb atom coordinates
        for i in xrange(rn):
            for j in xrange(i):
                m = 10000000
                for x in co[i]:
                    for y in co[j]:
                        e = euclidean(x, y)
                        if e < m:
                            m = e
                D[i][j] = m
    else:
        raise Exception('Unknown setting for point: %s' % point)
    if return_centroid_coordinates:
        return (D, pdb_resids, co)
    return (D, pdb_resids)


def getFragmentAnnotation(pdb, ch):
    """ Get the PDB 'fragment' annotation for chain ch """
    poly = parsePDBHeader('%s/%s/pdb%s.ent.gz' % (PDB_STRUCTURES, pdb[1:3], pdb), 'polymers')
    for p in poly:
        if p.chid != ch:
            continue
        return p.fragment
    return ''


def getMutations(url, u1, ttype, muttypes=['M'], samplemutfreq=None, returnSampleIDs=False):
    md = {}
    gm = urllib2.urlopen(url + u1).read()
    for l in map(lambda x:x.split('\t'), gm.split('\n')):
        if l == [''] or not l[5]:
            continue
        if l[6][0] not in muttypes:
            continue
        ## ttype selection
        if ttype and not l[0].startswith(ttype):
            continue
        p = int(l[5])
        if returnSampleIDs:
            val = set([])
        else:
            val = 0
        if p not in md:
            md[p] = val
        if samplemutfreq:  ## sample-mutation-frequency based weighting of mutations
            md[p] += samplemutfreq[l[1]]
        elif returnSampleIDs:
            md[p].add((l[1],l[0]))
        else:  ## weight mutations equally (CLUMPS-1 approach)
            md[p] += 1
    return md


def affinityPropagation(S, itrs=300, lam=0.9):
    dim = len(S)
    R1 = sp.zeros((dim,dim),float)
    A1 = sp.zeros((dim,dim),float)
    R2 = sp.zeros((dim,dim),float)
    A2 = sp.zeros((dim,dim),float)
    for itr in xrange(itrs):
        # calc R2
        for i in xrange(dim):
            for k in xrange(dim):
                m = -sp.inf  ## the maximum
                for kp in xrange(dim):
                    if kp == k:
                        continue
                    t = A1[i,kp] + S[i,kp]
                    if t > m:
                        m = t
                R2[i,k] = (1-lam) * (S[i,k] - m)  +  lam * (R1[i,k])
        # calc A2
        for i in xrange(dim):
            for k in xrange(dim):
                if i == k:
                    nv = 0
                    for ip in xrange(dim):
                        if ip == k:
                            continue
                        nv += max(0,R2[ip,k])
                else:
                    sm = R2[k,k]
                    for ip in xrange(dim):
                        if ip == i or ip == k:
                            continue
                        if R2[ip,k] > 0:
                            sm += R2[ip,k]
                    nv = min(sm,0)
                A2[i,k] = (1-lam) * nv  +  lam * A1[i,k]
        (R1,R2) = (R2,R1)
        (A1,A2) = (A2,A1)
    res = {}
    for i in xrange(dim):
        # find the exemplar for this point
        exm = -1
        m = 0
        for k in xrange(dim):
            t = A1[i,k] + R1[i,k]
            if t > m:
                m = t
                exm = k
        res[i] = exm
    return (res,R1,A1)


def getSiftsUniprot2PdbResidueMap(pdbid, matchChain, matchUniprot_id):
    fname = pdbid.lower() + '.xml.gz'
    intre = re.compile('-?\d+')
    ## get the file
    try:
        fi = GzipFile(SIFTS_RESIDMAPS + pdbid[1:3] + '/' + fname, 'r')
    except:
        print 'getting sifts file from ftp'
        ftp = ftplib.FTP("ftp.ebi.ac.uk")
        ftp.login()
        ftp.cwd('/pub/databases/msd/sifts/xml')
        fo = file(SIFTS_RESIDMAPS + pdbid[1:3] + '/' + fname, 'w')
        ftp.retrbinary('RETR %s' % fname, fo.write)
        fo.close()
        fi = GzipFile(SIFTS_RESIDMAPS + pdbid[1:3] + '/' + fname, 'r')
    ## parse xml
    ret = {}
    tree = etree.fromstring(fi.read())
    for i in xrange(len(tree)):
        if tree[i].tag.split('}')[-1] != 'entity':
            continue
        for j in xrange(len(tree[i])):
            if tree[i][j].tag.split('}')[-1] != 'segment':
                continue
            for k in xrange(len(tree[i][j])):
                if tree[i][j][k].tag.split('}')[-1] != 'listResidue':
                    continue
                for resid in tree[i][j][k].iterchildren():
                    pdb_id = None
                    pdb_rn = None  ## residue number
                    pdb_aa = None
                    uniprot_id = None
                    uniprot_rn = None  ## residue number
                    uniprot_aa = None
                    for ref in resid.iterchildren():
                        if ref.tag.split('}')[-1] != "crossRefDb":
                            continue
                        if ref.attrib['dbSource'] == 'PDB':
                            if matchChain != ref.attrib['dbChainId']:
                                continue
                            pdb_id = ref.attrib['dbAccessionId']
                            pdb_chain = ref.attrib['dbChainId']
                            pdb_aa = ref.attrib['dbResName']
                            pdb_rn = int(intre.findall(ref.attrib['dbResNum'])[0])
                        elif ref.attrib['dbSource'] == 'UniProt':
                            if matchUniprot_id != ref.attrib['dbAccessionId']:
                                continue
                            uniprot_id = ref.attrib['dbAccessionId']
                            uniprot_aa = ref.attrib['dbResName']
                            uniprot_rn = int(ref.attrib['dbResNum'])
                    if pdb_rn is not None and uniprot_rn is not None:
                        if ret.has_key(uniprot_rn) and ret[uniprot_rn] != pdb_rn:  ## should not happen; otherwise: redundant maps
                            print pdbid, matchChain, matchUniprot_id, ret[uniprot_rn], (pdb_rn, pdb_aa)
                            #raise Exception("Fix this! 1")
                        #ret[uniprot_rn] = (pdb_rn, pdb_aa)
                        ret[uniprot_rn] = pdb_rn
    return ret


def getPdbsumPPIs(pdb, protchains={}):
    ppi = {}
    try:
        fi = GzipFile('../dat/pdbsum/%s_pp.gz' % (pdb))
    except:
        urllib.urlretrieve('http://www.ebi.ac.uk/thornton-srv/databases/PDBsum/%s/%s/igrow.out.gz' % (pdb[1:3], pdb), '../dat/pdbsum/%s_pp.gz' % pdb)
        fi = GzipFile('../dat/pdbsum/%s_pp.gz' % (pdb))
    while 1:
        try:
            l = fi.readline()  ## some of the pdbsum files are corrupt...
        except:
            continue
        if not l:
            break
        l = l.split()
        if not l[4].isdigit() or not l[12].isdigit():
            continue
        c1c2 = [(l[3], int(l[4])), (l[13],int(l[12]))]
        c1c2.sort()
        c1,c2 = map(lambda x:x[0], c1c2)
        r1,r2 = map(lambda x:int(x[1]), c1c2)
        if protchains and (c1 not in protchains or c2 not in protchains):
            continue
        if not ppi.has_key(c1):
            ppi[c1] = {}
        ## the commented part is helpful if we need residue-pair resolution. for now we don't.
        #if not ppi[c1].has_key(c2):
        #    ppi[c1][c2] = {}
        #ppi[c1][c2][(r1, r2)] = True
        if c2 not in ppi[c1]:
            ppi[c1][c2] = ({},{})
        ppi[c1][c2][0][r1] = 1
        ppi[c1][c2][1][r2] = 1
    fi.close()
    return ppi


def getPdbsumPOIs(pdb, protchains={}, dnarna={}):
    poi = {}
    try:
        fi = GzipFile('../dat/pdbsum/%s_po.gz' % (pdb))
    except:
        urllib.urlretrieve('http://www.ebi.ac.uk/thornton-srv/databases/PDBsum/%s/%s/grow.out' % (pdb[1:3], pdb), '../dat/pdbsum/%s_po' % pdb)
        if file('../dat/pdbsum/%s_po' % pdb).readline().startswith('<!DOCTYPE HTML'):  ## has not been found. create empty file
            fo = file('../dat/pdbsum/%s_po' % pdb, 'w')
            fo.close()
        os.system('gzip ../dat/pdbsum/%s_po' % pdb)
        fi = GzipFile('../dat/pdbsum/%s_po.gz' % (pdb))
    while 1:
        try:
            l = fi.readline().strip()  ## some of the pdbsum files are corrupted...
        except:
            continue
        if not l:
            break
        l = l.split()
        if len(l) < 5:
            continue
        if not l[4].isdigit():
            continue
        c1c2s = []
        if l[3] != l[13]  and  dnarna  and  l[13] in dnarna:  ## different chain that is a DNA or RNA chain
            c1c2 = [(l[3], int(l[4])), ((l[13], 0, dnarna[l[13]][0]), int(l[12]))]  ## the second element is something like (('B', 0, 'D'), 12) where 'D' is DNA
            c1c2s.append(c1c2)
            paired = [l[13]] + dnarna[l[13]][1].keys()
            paired = ','.join(sorted(paired))
            c1c2 = [(l[3], int(l[4])), ((paired, 0, 'DNA/RNA'), 0)]  ## the second element is something like (('B', 0, 'D'), 12) where 'D' is DNA
            c1c2s.append(c1c2)
        elif l[3] == l[13]:  ## same chain, should be a small molecule
            c1c2 = [(l[3], int(l[4])), ((l[13], int(l[12]), l[11]), 0)]
            c1c2s.append(c1c2)
        else:  ## probably small molecule on different chain
            #print 'DROPPING', ' '.join(l)
            continue
        for c1c2 in c1c2s:
            c1,c2 = map(lambda x:x[0], c1c2)
            r1,r2 = map(lambda x:x[1], c1c2)
            if c1 not in protchains:
                continue
            if not poi.has_key(c1):
                poi[c1] = {}
            ##  only if residue-level resolution is needed, and for compatibility reasons.
            #if not poi[c1].has_key(c2):
            #    poi[c1][c2] = {}
            #poi[c1][c2][(r1, r2)] = True
            if c2 not in poi[c1]:
                poi[c1][c2] = ({},{})
            poi[c1][c2][0][r1] = 1
    fi.close()
    return poi


def getFoldxPIs(pdb, protchains={}, dnarna={}):
    pi = {}
    try:
        fi = file('../fdx/if_resid/%s' % pdb)
    except:
        os.chdir('../fdx/tmp/')
        os.system('zcat ../../dat/pdb/%s.pdb.gz > %s.pdb' % (pdb,pdb))
        os.system('foldx4 --command=AnalyseComplex --pdb=%s.pdb' % pdb)
        os.system('mv Interface_Residues_%s_AC.fxout ../if_resid/%s' % (pdb,pdb))        
        fi = file('../if_resid/%s' % pdb)
    while 1:
        l = fi.readline()
        if not l:
            break
        if l.startswith('interface residues between '):
            chch = l[27:].strip().split(' and ')
            r1r2 = [{},{}]
            l = fi.readline().strip()
            if not l:
                continue
            allifres = l.strip().split()
            for r in allifres:
                if r[1] == chch[0]:
                    i = 0
                else:
                    i = 1
                r1r2[i][int(r[2:])] = 1
            if chch[0] > chch[1]:
                chch.reverse()
                r1r2.reverse()
            if chch[0] not in pi:
                pi[chch[0]] = {}
            if chch[1] not in pi[chch[0]]:
                pi[chch[0]][chch[1]] = tuple(r1r2)
    ppi = {}
    poi = {}
    for x1 in pi:
        for x2 in pi[x1]:
            c1,c2 = x1,x2
            if c1 not in protchains and c2 not in protchains:
                continue
            elif (c1 in protchains and c2 in dnarna) or (c2 in protchains and c1 in dnarna):  ## one of the chains is DNA/RNA
                if c2 in protchains:
                    r2,r1 = pi[c1][c2]
                    c1,c2 = c2,c1
                else:
                    r1,r2 = pi[c1][c2]
                c2 = (c2, 0, dnarna[c2][0])
                if c1 not in poi:
                    poi[c1] = {}
                poi[c1][c2] = (r1, {})
                paired = [c2[0]] + dnarna[c2[0]][1].keys()
                paired = (','.join(sorted(paired)), 0, 'DNA/RNA')
                if paired not in poi[c1]:
                    poi[c1][paired] = ({k:r1[k] for k in r1},{})
                else:
                    poi[c1][paired][0].update(r1)
            elif c1 in protchains and c2 in protchains:  ## its' a PPI
                if c1 not in ppi:
                    ppi[c1] = {}
                ppi[c1][c2] = pi[c1][c2]
    return (ppi, poi)


def calculateClusters(D, mi, mv):
    """ Calculate clusters of mutated residues.
        D is a distance matrix;
        mi and mv are the indices (relative to D) and normalized recurrences, respectively, of the mutated residues.
    """
    cm = len(mi)
    S = sp.zeros((cm,cm))
    for i in xrange(cm):
        for j in xrange(cm):
            if i == j:
                continue
            d = D[mi[i],mi[j]]
            S[i,j] = S[j,i] = -d*d
    #med = sp.median(S)
    #print 'med', med
    for i in xrange(cm):
        S[i,i] = 225 * (mv[i]/2-1) + random.random()  ## the random noise is just to prevent equal affinities/responsibilities
    cl = affinityPropagation(S)
    cl = cl[0].items()
    cc = []
    for x in cl:
        ok = False
        for y in cc:
            if y[1] == x[1]:
                y[0].append(x[0])
                ok = True
        if not ok:
            cc.append([[x[0]],x[1]])
    cc = [i[0] for i in cc]
    cc.sort(key=lambda x:len(x), reverse=True)
    change = True
    it = 0
    while change and it < 200:
        change = False
        ## remove residues that are >15A away from half of the cluster's members
        removed = []  ## residues removed from clusters by c15
        ncc = []
        for c in cc:
            l = len(c)
            nc = []
            for i in xrange(l):
                c15 = mv[c[i]]
                denom = mv[c[i]]
                for j in xrange(i+1, l):
                    denom += mv[j]
                    if D[mi[c[i]],mi[c[j]]] < 15:
                        c15 += mv[c[j]]
                if c15 < denom/2:
                    removed.append(c[i])
                else:
                    nc.append(c[i])
            ncc.append(nc)
        if removed:
            cc = ncc + [[i] for i in removed]    ## create individual clusters with removed residues
            change = True
            it += 1
            continue
        ## merge clusters that are close
        mp2s = []  ## list of ((c1,c2), s) where c1 and c2 are clusters and s is the score (normalized c15)
        for i in xrange(len(cc)):
            c1 = cc[i]
            for j in xrange(i+1, len(cc)):
                c2 = cc[j]
                c15 = 0
                denom = 0
                for x in c1:
                    for y in c2:
                        d = D[mi[x],mi[y]]
                        wd = mv[x] * mv[y] * D[mi[x],mi[y]]
                        denom += wd
                        if d < 15:
                            c15 += wd
                mp2s.append(((i,j), c15/denom))
        best = max(mp2s, key=lambda x:x[1])  ## best candidate of cluster pairs to merge
        if best[1] > 0.5:
            ncc = [cc[i] for i in xrange(len(cc)) if i not in best[0]]
            ncc.append(cc[best[0][0]] + cc[best[0][1]])
            #print 'merger', cc, ncc
            cc = ncc
            change = True
        it += 1
    return cc