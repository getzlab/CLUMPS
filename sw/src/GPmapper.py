import os
import sys

from Bio.Seq import Seq
from gzip import GzipFile
from twobitreader import TwoBitFile


class GPmapper():
    def __init__(self, hgfile='hg19.2bit', spfile='UP000005640_9606.fasta.gz', mapfile='/sw/dat/genomeProteomeMaps.txt'):
        self.hg = TwoBitFile(hgfile)
        ##
        self.sp = {}
        fi = GzipFile(spfile)
        annot = fi.readline()
        seq = ''
        while 1:
            l = fi.readline()
            if not l or l.startswith('>'):
                an = annot.strip().lstrip('>').split('|')[1]
                self.sp[an] = seq
                if not l:
                    break
                annot = l
                seq = ''
            else:
                seq += l.strip()
        fi.close()
        ##
        self.gen2prot = {}
        self.prot2gen = {}
        fi = file(mapfile)
        while 1:
            l = fi.readline()
            if not l:
                break
            l = l.strip('\n').split('\t')
            if len(l) < 8:
                break
            if l[0] not in self.gen2prot:
                self.gen2prot[l[0]] = []
            bl = []
            if l[7]:
                bl = map(int,l[7].split(','))
            bl = set(bl)
            self.gen2prot[l[0]].append([int(l[1]), int(l[2]), [(int(i[0]),int(i[1])) for i in map(lambda x:x.split('-'), l[3].split(','))], l[4], l[5], [((int(i[0][:-2]), int(i[0][-1])), (int(i[1][:-2]), int(i[1][-1]))) for i in map(lambda x:x.split('-'), l[6].split(','))], bl])
            if l[5] not in self.prot2gen:
                self.prot2gen[l[5]] = []
            self.prot2gen[l[5]].append((l[0], len(self.gen2prot[l[0]])-1))
        self.baseComplement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
        #for chr in self.gen2prot:  MESSES UP prot2gen!!!
        #    self.gen2prot[chr].sort()
    
    
    def findExons(self, chr, pos):
        """ Finds exons that contain pos """
        ret = []
        for gi in xrange(len(self.gen2prot[chr])):
            g = self.gen2prot[chr][gi]
            if not g[0] <= pos < g[1]:
                continue
            if pos in g[6]:
                continue
            for ei in xrange(len(g[2])):
                if not g[2][ei][0] <= pos < g[2][ei][1]:
                    continue
                ret.append((gi,ei))
                break
        return ret
    
    
    def mapGPchange(self, chr, pos, gi, ei, newbase):
        """ Maps a genomic change to change(s) in protein(s).
            Either pos or both gi,ei should be None. """
        pos -= 1
        if chr == '23':
            chr = 'X'
        elif chr == '24':
            chr = 'Y'
        elif chr == 'MT':
            chr = 'M'
        chr = 'chr'+chr
        if gi is None:
            exons = self.findExons(chr, pos)
        else:
            exons = [(gi,ei)]
        ret = []
        for gi,ei in exons:
            g = self.gen2prot[chr][gi]
            if g[4] not in self.sp:
                continue  ## deprecated UniProt ID
            if g[3] == '+':
                offset = pos - g[2][ei][0]
            else:
                offset = g[2][ei][1] - pos - 1
            poffset = (offset/3, offset%3)              ## offset in protein coordinates
            qblockstart = g[5][ei][0]
            pp = (qblockstart[0] + poffset[0] + (qblockstart[1] + poffset[1])/3, (qblockstart[1] + poffset[1])%3)
            if g[3] == '+':
                r = pp[1]      ## residual
            else:
                r = 2-pp[1]    ## residual
            codon = []
            fromPrevExon = min(0, pos-r - g[2][ei][0])
            if fromPrevExon:
                codon += range(g[2][ei-1][1] + fromPrevExon, g[2][ei-1][1])
            fromCurrExon = [max(pos-r, g[2][ei][0]), min(pos+(3-r), g[2][ei][1])]
            codon += range(fromCurrExon[0], fromCurrExon[1])
            fromNextExon = max(0, pos+(3-r) - g[2][ei][1])
            if fromNextExon:
                codon += range(g[2][ei+1][0], g[2][ei+1][0]+fromNextExon)
            codonseq = ''.join([self.hg[chr][i:i+1] for i in codon])  ## (reference) codon sequence
            ncodonseq = list(codonseq)                                ## new (mutant) codon sequence
            if newbase:
                ncodonseq[r] = newbase
            ncodonseq = ''.join(ncodonseq)
            if g[3] == '+':
                try:
                    raa = self.sp[g[4]][pp[0]]            ## actual amino acid according to uniprot
                except:
                    continue  ## index error (seq has changed)
                taa = str(Seq(codonseq).translate())  ## amino acid resulting from translation
                naa = str(Seq(ncodonseq).translate()) ## new (mutant) amino acid 
            else:
                try:
                    raa = self.sp[g[4]][pp[0]]            ## actual amino acid according to uniprot
                except:
                    continue  ## index error (seq has changed)
                taa = str(Seq(codonseq).reverse_complement().translate())  ## amino acid resulting from translation
                naa = str(Seq(ncodonseq).reverse_complement().translate()) ## new (mutant) amino acid 
            if taa != raa:
                print 'WARNING: non-matching reference and translated AA', chr, pos, fromPrevExon, fromCurrExon, fromNextExon, len(codon), ei, len(g[2])
                continue
            ret.append((g[4], 1+pp[0], taa, naa))
        return ret


def makeGPfile():
    gpm = GPmapper()
    inMaf = sys.argv[1]
    if inMaf.endswith('gz'):
        fi = GzipFile(inMaf)
    else:
        fi = file(inMaf)
    fo = GzipFile('muts.gp.gz', 'w')
    h = fi.readline()
    h = h.strip().split('\t')
    iChr = h.index('Chromosome')
    iPos = h.index('Start_Position')
    iRefAllele = h.index('Reference_Allele')
    iTumAllele2 = h.index('Tumor_Seq_Allele2')
    iPatient = h.index('Tumor_Sample_Barcode')
    iTtype = h.index('ttype')
    fo.write('\t'.join(['patient', 'ttype', 'chr', 'pos', 'refbase', 'newbase', 'uniprot_change']) + '\n')
    while 1:
        l = fi.readline()
        if not l:
            break
        l = l.strip('\n').split('\t')
        if len(l[iTumAllele2]) > 1 or l[iTumAllele2] == '-':
            continue
        if len(l[iRefAllele]) > 1 or l[iRefAllele] == '-':
            continue
        newbase = l[iTumAllele2]
        r = gpm.mapGPchange(l[iChr], int(l[iPos]), None, None, newbase)
        if not r:
            continue
        r = '; '.join(map(lambda x:'%s:%s%d%s' % (x[0], x[2], x[1], x[3]), r))
        nl = [l[iPatient][:12], l[iTtype], l[iChr], l[iPos], l[iRefAllele], newbase, r]  #l[iPat] l[iRefAllele] 
        fo.write('\t'.join(nl) + '\n')
    fo.close()


def splitMutFile(gpfile):
    """ Takes the mutations file created by GPmapper and splits it to individual mutation files """
    mutTypes = set(['M']) # mutation types
    pdata = {}  ## data acessible by protein
    if gpfile.endswith('.gz'):
        fi = GzipFile(gpfile)
    else:
        fi = file(gpfile)
    hdr = fi.readline().strip('\n').split('\t')
    iUniprotChange = hdr.index('uniprot_change')
    iPatient = hdr.index('patient')
    iTtype = hdr.index('ttype')
    while 1:
        l = fi.readline()
        if not l:
            break
        l = l.strip('\n').split('\t')
        if not l[iUniprotChange]:
            continue
        uchanges = l[iUniprotChange].split('; ')
        for um in uchanges:
            u,m = um.split(':')
            if m[-1] == '*':
                mt = 'N'
            elif m[-1] == m[0]:
                mt = 'S'
            else:
                mt = 'M'
            if mt not in mutTypes:
                continue
            mi = m[1:-1]
            while not mi[-1].isdigit():  ## because some mutations are annotated as SNPs aohtough they are not
                mi = mi[:-1]
            dline = '\t'.join([l[iTtype], l[iPatient], 'na', u, 'p.' + m, mi, mt])
            if not pdata.has_key(u):
                pdata[u] = []
            pdata[u].append(dline)
    os.makedirs('splitByProtein/')
    for u in pdata:
        fo = file('splitByProtein/%s' % u, 'w')
        fo.write('\n'.join(pdata[u]) + '\n')
        fo.close()


if __name__ == '__main__':
    makeGPfile()
    splitMutFile('muts.gp.gz')
