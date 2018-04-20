import random
import scipy as sp

from Bio.Data.CodonTable import standard_dna_table
from CoverageSampler import CoverageSampler
from GPmapper import GPmapper


compl = {'a':'t', 'c':'g', 'g':'c', 't':'a'}
def reverse_complement(abc):
    return compl[abc[2]] + compl[abc[1]] + compl[abc[0]]

class MutspecCoverageSampler(CoverageSampler):
    def __init__(self, availUPresid, upid, covtrack, mutSpectraFn):
        CoverageSampler.__init__(self, availUPresid, upid, covtrack)
        # self.gpm
        # self.covprobs
        # self.availUPresid
        # self.availUPresidIdx
        # self.upid
        
        self.mutspecprobs = {}
        self.presamples = []
        self.presampleindex = 1000
                
        ## extract mutation spectra from file
        self.patcounts = {}
        self.conefindex = {}  ## index of context-effects
        fi = file(mutSpectraFn)
        conef = fi.readline().strip().split('\t')[2:]
        for i in xrange(len(conef)):
            self.conefindex[tuple(conef[i].split('-'))] = i
        while 1:
            l = fi.readline()
            if not l:
                break
            l = l.strip().split('\t')
            self.patcounts[tuple(l[:2])] = map(int, l[2:])

        ## codon table
        self.codonTable = {i.lower():standard_dna_table.forward_table[i] for i in standard_dna_table.forward_table}

        ## calculate possible effects
        aa2genpos = {}
        self.availUPresidSet = set(self.availUPresid)
        x = self.gpm.prot2gen[upid][0]
        gp = self.gpm.gen2prot[x[0]][x[1]]
        chr = x[0]
        trdir = gp[3]
        gex = gp[2]
        pex = gp[5]
        blacklist = gp[6]
        ig = gex[0][0]  ## genomic position
        ie = 0          ## exon index
        phase = 0
        if trdir == '+':
            forward = True
        else:
            forward = False
        if forward:
            ip = pex[0][0][0]+1  ## protein position
        else:
            ip = pex[0][1][0]    ## protein position
        while 1:
            if ip in self.availUPresidSet and ig not in blacklist:
                if ip not in aa2genpos:
                    aa2genpos[ip] = []
                aa2genpos[ip].append(ig)
            ig += 1
            phase += 1
            if phase == 3:
                phase = 0
                if forward:
                    ip += 1
                else:
                    ip -= 1
            if ig >= gex[ie][1]:
                ie += 1
                if ie == len(gex): ## we've reached the end of the last exon
                    break
                ig = gex[ie][0]    ## otherwise, correct ig
                if forward:
                    ip = pex[ie][0][0] + 1
                else:
                    ip = pex[ie][1][0] + (pex[ie][1][1] > 0)
        self.aa2conteff = {}
        for ip in aa2genpos:
            self.aa2conteff[ip] = []
            gposs = aa2genpos[ip]   ## genomic positions of the codon bases
            origcodon = ''.join([self.gpm.hg[chr][i:i+1] for i in gposs]).lower()
            if not forward:
                origcodon = reverse_complement(origcodon)
                gposs.reverse()
            origaa = self.codonTable[origcodon]
            if origaa != self.gpm.sp[upid][ip-1]:
                print origaa, self.gpm.sp[upid][ip-1], ip-1, gposs
                raise Exception('Translation does not match the reference!')
            for i in xrange(3):  ## codon position
                ## test which change will create a missense mutation
                for j in ['a','c','g','t']:
                    if origcodon[i] == j:
                        continue
                    newcodon = origcodon[:i]+j+origcodon[i+1:]
                    if newcodon in self.codonTable and self.codonTable[newcodon] != origaa:  ## this change creates a missense mutation
                        trin = ''.join(self.gpm.hg[chr][gposs[i]-1:gposs[i]+2]).lower()
                        if 'n' in trin:
                            continue
                        if not forward:
                            j = compl[j]
                        if trin[1] in ['a','c']:
                            self.aa2conteff[ip].append(self.conefindex[(trin,j)])
                        else:
                            trin = reverse_complement(trin)
                            self.aa2conteff[ip].append(self.conefindex[(trin,compl[j])])


    def calcMutSpecProbs(self, md):
        """ calculate probabilities per uniprot position """
        patprobs = {}
        for pos in md:
            for pat in md[pos][2]:
                if pat not in patprobs:
                    patprobs[pat] = []
        for pat in patprobs:
            #totalmut = float(sum(self.patcounts[pat] ))  ## total mutations in patient
            for i in self.availUPresid:
                if i not in self.aa2conteff:  ## because it was in the blacklist
                    #print 'WARNING: check if residue %s maps to blacklisted positions (it should)' % i
                    p = 0
                else:
                    p = sum([self.patcounts[pat][x] for x in self.aa2conteff[i]]) #/totalmut  ## normalization not required
                ## subtract one if the patient has a mutation at that residue (to avoid leaving the residue as the only option or the most likely option, especially in patients with low mutation counts)
                if i in md and pat in md[i][2]:
                    p -= 1
                patprobs[pat].append(p)
            denom = float(sum(patprobs[pat]))
            if denom == 0.0:  ## because we subtract 1 above, this may happen (patient has only 1 mutation)
                patprobs[pat] = [1.0/len(patprobs[pat]) for x in patprobs[pat]]
            else:
                patprobs[pat] = [x/denom for x in patprobs[pat]]
        for pos in md:
            if pos not in self.aa2conteff:  ## because it was in the blacklist
                #print 'WARNING: check if residue %s maps to blacklisted positions (it should)' % pos
                continue
            if len(md[pos][2]) == 1:
                for j in md[pos][2]:  ## there is actually one element but I don't want to use pop()
                    p = patprobs[j]
            else:
                ## prob product
                #p = [1]*len(self.availUPresid)
                #for i in xrange(len(self.availUPresid)):
                #    for j in md[pos][2]:
                #        p[i] *= patprobs[j][i]
                
                p = [[] for i in xrange(len(self.availUPresid))]
                for i in xrange(len(self.availUPresid)):
                    for j in md[pos][2]:
                        p[i].append(patprobs[j][i])
                p = map(lambda x:len(x) and sp.median(x) or 0, p)
            ## multiply by coverage vector
            p = [p[i]*self.covprobs[i] for i in xrange(len(self.availUPresid))]
            p = [i/sum(p) for i in p]
            self.mutspecprobs[pos] = p


    def presample(self, mireal):  ## mi are indices
        self.presamples = []
        self.presampleindex = 0
        for i in mireal:
            m = self.availUPresid[i]
            r = sp.random.choice(self.availUPresidIdx, 1000, p=self.mutspecprobs[m])
            self.presamples.append(r)


    def sample(self, mireal):
        if self.presampleindex == 1000:
            self.presample(mireal)
        idx = sp.random.permutation(len(mireal))
        ret = {}
        for m in idx:
            n = self.presamples[m][self.presampleindex]
            while n in ret:
                for n in sp.random.choice(self.availUPresidIdx, 1000, p=self.mutspecprobs[self.availUPresid[mireal[m]]]):
                    if n not in ret:
                        break
                if n in ret:
                    return
            ret[n] = m
        self.presampleindex += 1
        ret = ret.items()
        ret.sort()
        return ([i[0] for i in ret], [i[1] for i in ret])
