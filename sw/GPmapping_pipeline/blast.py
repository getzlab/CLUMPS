import os
import sys

from Bio.Seq import Seq
from gzip import GzipFile
from lxml import etree
from math import ceil
from twobitreader import TwoBitFile


#################
###  SECTION 1: SIFTS TO HUMAN-UNIPROT
#################


def buildBlastDb():
    """ Build a custom BLAST+ database from the uniprot IDs that are in the sifts map file """
    sifts = sys.argv[1]
    sprot = sys.argv[2]
    trembl = sys.argv[3]
    ##
    upids = {}
    for l in file(sifts).readlines()[1:]:
        upids[l.strip().split('\t')[2]] = 0
    print len(upids), 'uniprot ids in the sifts file.'
    ##
    fo = file('../dat/uniprot.sifts.custom.db.seq', 'w')
    for fn in [sprot, trembl]:
        fi = GzipFile(fn)
        annot = fi.readline()
        seq = []
        while 1:
            l = fi.readline()
            if not l or l.startswith('>'):
                an = annot.split('|')
                if an[1] in upids:
                    fo.write(annot + (''.join(seq)))
                    upids[an[1]] = 1
                if not l:
                    break
                annot = l
                seq = []
            else:
                seq.append(l)
    fo.close()
    cnt = upids.values().count(0)
    c = 0
    for k in upids:
        if c > 5:
            break
        if not upids[k]:
            print k
            c += 1
    print cnt, 'upids from sifts were not found in the uniprot fasta file.'

#buildBlastDb()
## python blast.py ../dat/sifts/pdb_chain_uniprot.tsv ../../dat/uniprot/uniprot_sprot.fasta.gz ../../dat/uniprot/uniprot_trembl.fasta.gz
## makeblastdb -in uniprot.sifts.custom.db.seq -dbtype 'prot' -out siftsdb


def prepareBlastInput():
    """ Formats the input for the blast command """
    fi = GzipFile('../dat/uniprot/UP000005640_9606.fasta.gz')  ## uniprot.proteomes.HUMAN.fasta.15-04-28.gz
    foall = file('../dat/uniprot.human.sp.fa', 'w')
    annot = fi.readline()
    seq = []
    while 1:
        l = fi.readline()
        if not l or l.startswith('>'):
            ## run blast if the sequence qualifies (swissprot)
            if annot.startswith('>sp|'):
                annot = annot.split('|')
                if '-' not in annot[1]:
                    foall.write('>%s\n' % annot[1])
                    foall.write(''.join(seq) + '\n')
                    fo = file('../dat/uniprot.human.sp/%s.seq' % annot[1], 'w')
                    fo.write(''.join(seq) + '\n')
                    fo.close()
            if not l:
                break
            annot = l
            seq = []
        else:
            seq.append(l.strip())

#prepareBlastInput()

###############################
# for i in `\ls -1 *.seq`; do blastp -query $i -db siftsdb/siftsdb -out $i.blasted -outfmt 5 && gzip $i.blasted && rm -f $i; done  ## the backslash in flont of ls makes sure that color is disabled (otherwise leading to problems)
### DEPRECATED: I HATE PSIBLAST!!!!!!! CANNOT BLAST PROBERLY HRAS TO KRAS!!!! for i in `ls *.seq`; do psiblast -query $i -num_iterations 5 -db siftsdb/siftsdb -out $i.blasted -outfmt 5 && gzip $i.blasted && rm -f $i; done
### DEPRECATED line: get pdbaa from ftp://ftp.ncbi.nih.gov/blast/db/
### DEPRECATED line:  for i in `ls *.seq`; do bsub -q hour -P clumps2 -o blast.log psiblast -query $i -db pdbaa -out $i.blasted -outfmt 5 && gzip $i.blasted; done
###############################


##################
###  SECTION 2: HUMAN-UNIPROT TO GENOMIC POSITIONS
##################

def buildBlatQueryList(chunksize=100):
    """ Creates a fasta files with SwissProt protein sequences (per chromosome and in chunks) """
    uniprotfa = '../dat/uniprot.human.sp.fa'
    uniprot2entrezg = '../dat/uniprot2entrezg.tab'
    entrezg2chromosome = '../dat/entrezg2chromosome.tab'
    e2c = {}
    for l in file(entrezg2chromosome).readlines():
        e,c = l.strip().split()
        if c in ('-','Un') or '|' in c:
            continue
        if c == 'MT':
            c = 'M'
        e2c[e] = c
    u2e = {}
    for l in file(uniprot2entrezg).readlines():
        u,ee = l.strip('\n').split('\t')
        if not ee:
            continue
        f = [i.strip() for i in ee.strip(';').split(';') if i in e2c]
        if f:
            u2e[u] = f
    ##
    c2u = {c:[] for c in [str(i) for i in range(1,23)] + ['X','Y','M','NA']}
    ##
    upids = {}
    fi = file(uniprotfa)
    annot = fi.readline()
    seq = []
    while 1:
        l = fi.readline()
        if not l or l.startswith('>'):
            an = annot.strip()[1:]
            upids[an] = ''.join(seq)
            if not l:
                break
            annot = l
            seq = []
        else:
            seq.append(l)
    for u in upids:
        if u in u2e:
            for e in u2e[u]:
                c = e2c[e]
                c2u[c].append(u)
        else:
            c2u['NA'].append(u)
    for c in c2u:
        cups = c2u[c]
        nchunks = len(cups)/chunksize
        if len(cups)%chunksize:
            nchunks += 1
        for i in range(nchunks):
            fo = file('../dat/blat.in/upids.chr%s.%d.fa' % (c,i), 'w')
            for u in range(i*chunksize,min((i+1)*chunksize, len(cups))):
                fo.write('>%s\n' % cups[u])
                fo.write('%s' % upids[cups[u]])
            fo.close()
        
# buildBlatQueryList()

## for chr in {1..22} X Y M NA; do for i in `ls *.chr$chr.*fa`; do srun --partition=short -t 12:00:00 -o ../../logs/$i.blat.log blat ../../dat/hg19/chr$chr.2bit $i $i.psl -t=dnax -q=prot -out=psl & done; done


def parseBlatOutput():
    alifiles = '../dat/blat.in/'
    upid2bestali = {}
    for alifile in os.listdir(alifiles):
        if not alifile.endswith('psl'):
            continue
        fi = file(alifiles+alifile)
        x = fi.readline()
        while not x.startswith('-----'):
            x = fi.readline()
            if not x:
                raise Exception('parsing failed')
        while 1:
            l = fi.readline()
            if not l:
                break
            l = l.strip('\n').split('\t')
            for i in [0,1,2,3,4,5,6,7,10,11,12]:
                l[i] = int(l[i])
            if l[9] not in upid2bestali:
                upid2bestali[l[9]] = l
            else:
                best = upid2bestali[l[9]]
                if l[0] > l[10]:  ## number of matches is greater than protein length
                    continue
                if l[0]-l[1] > best[0]-best[1]:
                    upid2bestali[l[9]] = l
    fo = file('../res/blat_sp_hg.psl', 'w')            
    for u in upid2bestali:
        fo.write('\t'.join(map(str,upid2bestali[u])) + '\n')
    fo.close()

parseBlatOutput()

print('blast output parsed')

#def mapGenomic2ProteinPos():
""" Pretty messy code (even though it works); consider rewriting """
if 1:
    #psl = '../test.psl'  ## the psl file resulting from BLAT
    psl = '../res/blat_sp_hg.psl'
    hgfile = '../dat/hg19/25chr.2bit'  ## reference human genome in 2bit format
    qfile = '../dat/uniprot.human.sp.fa'  ## file containing the BLAT queries (usually human SwissProt proteome)
    upid2seq = {}
    fi = file(qfile)
    u = None
    seq = []
    while 1:
        l = fi.readline()
        if not l:
            break
        if l.startswith('>'):
            upid2seq[u] = ''.join(seq)
            u = l[1:].strip()
            seq = []
        else:
            seq.append(l.strip())
    upid2seq[u] = ''.join(seq)
    
    maxoffset = 10      ## the maximum offset to apply for splice junctions
    hg = TwoBitFile(hgfile)
    fo = file('../res/genomeProteomeMaps.run2.txt', 'w')
    fi = file(psl)
    cnt100 = 0
    cnt99 = 0
    cnt95 = 0
    cnt90 = 0
    rest = []
    while 1:
        l = fi.readline()
        if not l:
            break
        l = l.strip('\n').split('\t')  ## 0)match  1)mismatch  2)rep.match  3)N's  4)Q_gap_count  5)Q_gap_bases  6)T_gap_count  7)T_gap_bases  8)strand  9)Q_name  10)Q_size  11)Q_start  12)Q_end  13)T_name  14)T_size  15)T_start  16)T_end  17)block_count  18)blockSizes  19)qStarts  20)tStarts
        #print l[9], l[8]
        chr = hg[l[13]]  ## chromosome name
        tsize = int(l[14])
        if l[8] == '++':
            forward = True
        else:
            forward = False
        bs = map(int, l[18].strip(',').split(','))
        qs = map(int, l[19].strip(',').split(','))
        ts = map(int, l[20].strip(',').split(','))
        if 1:  ## remove small exons (<3 residues) and treat them as gaps (the splice donor/acceptor search will consider them and extend other exons)
            toRemove = set([])
            for i in xrange(len(bs)):
                if bs[i] < 3:
                    toRemove.add(i)
            bs = [bs[i] for i in range(len(bs)) if i not in toRemove]
            qs = [qs[i] for i in range(len(qs)) if i not in toRemove]
            ts = [ts[i] for i in range(len(ts)) if i not in toRemove]
        #if l[9] not in ['Q8TC05','Q8TBF8']:  ## short exons examples: - and +
        #    continue
        if 1: ## try to perfectly map short exons that blat has missed
            toAdd = []
            df = [qs[0]] + [(qs[i+1]-qs[i]-bs[i]) for i in xrange(len(qs)-1)] + [int(l[10])-(bs[-1]+qs[-1])]
            #print l[9], df
            for i in xrange(len(df)):
                if df[i] < 4 or df[i] > 50:  ## either no unmapped exon, or exon size is too small or too big
                    continue
                if i == 0:             ## look before first exon
                    region = [max(ts[0]-750000, 0), ts[0]]  ## chromosome region to search for exon
                    q = upid2seq[l[9]][0:qs[0]]
                elif i == len(df)-1:   ## look after last exon
                    region = [ts[-1]+3*bs[-1], min(ts[-1]+3*bs[-1]+750000, len(chr))]
                    q = upid2seq[l[9]][qs[-1]+bs[-1]:int(l[10])]
                else:                  ## look between exons
                    region = [ts[i-1]+3*bs[i-1], ts[i]]
                    q = upid2seq[l[9]][qs[i-1]+bs[i-1]:qs[i]]
                #print 'looking for', l[9], 'query', q, 'in region', region, l[8]
                k = None  ## genomic position of the start of the newly mapped exon
                for off in range(3):
                    if forward:
                        f = str(Seq(chr[region[0]+off:region[1]+off]).translate())
                    else:
                        f = str(Seq(chr[int(l[14])-region[1]-off:int(l[14])-region[0]-off]).reverse_complement().translate())
                    try:
                        if i == 0:
                            k = f.rindex(q)
                        else:
                            k = f.index(q)
                        #if forward:
                        k = region[0] + off + 3*k
                        #else:
                        #    k = region[0] - off + 3*k
                        #print 'k', k
                        break
                    except:
                        pass ## not found in this frame
                if k:
                    if i < len(qs):
                        b = qs[i]-df[i]
                    else:
                        b = int(l[10])-df[i]
                    toAdd.append((df[i], b, k))
            if toAdd:
                for a,b,c in toAdd:
                    ok = 0
                    for j in xrange(len(qs)):
                        if qs[j] > b:
                            ok = 1
                            break
                    if not ok:
                        j = len(qs)
                    bs = bs[:j] + [a] + bs[j:]
                    qs = qs[:j] + [b] + qs[j:]
                    ts = ts[:j] + [c] + ts[j:]
                ## correct the query/target start & stop positions
                #print 'was', l[11], l[12], l[15], l[16]
                l[11] = min(int(l[11]), qs[0])
                l[12] = max(int(l[12]), bs[-1]+qs[-1])
                if forward:
                    l[15] = min(int(l[15]), ts[0])
                    l[16] = max(int(l[16]), 3*bs[-1] + ts[-1])
                else:
                    l[15] = min(int(l[15]), int(l[14])-ts[-1]-3*bs[-1])
                    l[16] = max(int(l[16]), int(l[14])-ts[0])
                #print 'is', l[11], l[12], l[15], l[16]
            if 0: ## just for debugging
                print '#'
                x = raw_input().strip()
                while x:
                    try:
                        exec(x)
                    except:
                        print '!!!'
                    x = raw_input().strip()

        qgaps = [qs[i+1] - (qs[i] + bs[i]) for i in xrange(len(bs)-1)]
        if forward:
            sdon = lambda:chr[iend + os*offset : iend + os*offset + 2].lower() in ['gt','gc']
            sacc = lambda:chr[jstart - 3*qgap + os*offset - 2 : jstart - 3*qgap + os*offset].lower() in ['ag']
        else:
            sdon = lambda:chr[tsize-(iend + os*offset) - 2 : tsize-(iend + os*offset)].lower() in ['ac','gc']
            sacc = lambda:chr[tsize-(jstart - 3*qgap + os*offset) : tsize-(jstart - 3*qgap + os*offset) + 2].lower() in ['ct']
        exons = []
        qblocks = []
        nextStart = ts[0]   ## start of next exon in the target DB (initiated with exon1_start)
        nextStartQ = qs[0]
        resid = 0
        codon = 0
        ## in each iteration, we simultaneously find the end of one exon i and the start of the next exon j:
        for i in xrange(len(bs)-1):  ## i is an exon index
            qgap = qgaps[i]
            iend = ts[i] + 3*bs[i]
            jstart = ts[i+1]
            calcOffset = True
            if qgap and iend+3*qgap > jstart:  ## in rare cases, this happens (there is a 'deletion' in the chromosome sequence); better not try to save some residues
                qgap = qgaps[i] = 0
                calcOffset = False
            iendQ = qs[i] + bs[i]
            jstartQ = qs[i+1]
            offset = 0
            while calcOffset and offset < maxoffset+3*qgap and qgap <= 3:
                ok = 0
                for os in [1, -1]:  ## offset sign
                    if sdon() and sacc():
                        offset *= os
                        ok = 1
                        break
                if ok:
                    break
                offset += 1
            #print i, offset, qgap
            if offset == maxoffset+3*qgap:
                offset = 0
                while offset < 4:
                    ok = 0
                    for os in [1, -1]:  ## offset sign
                        if sdon() or sacc():
                            offset *= os
                            ok = 1
                            break
                    if ok:
                        break
                    offset += 1
                if offset == 4:
                    offset = 0
                #print i, offset, qgap, '###'
            iendCorr = iend+offset
            jstartCorr = jstart+offset
            iendCorrQ = iendQ+offset/3.0
            jstartCorrQ = jstartQ+offset/3.0
            if 0 < qgap <= 3:
                jstartCorr -= 3*qgap
                jstartCorrQ -= qgap
            if forward:
                exons.append((nextStart, iendCorr))
            else:
                exons.append((tsize-iendCorr, tsize-nextStart))
            qblocks.append((nextStartQ, iendCorrQ))
            nextStart = jstartCorr
            nextStartQ = jstartCorrQ
        if forward:
            exons.append((nextStart, int(l[16])))
        else:
            exons.append((int(l[15]), tsize-nextStart))
        qblocks.append((nextStartQ, int(l[12])))        
        #print 'exons', el, sum(el), sum(qbl), int(sum(el)), int(sum(el)) == sum(qbl), qbl
        if round(sum(map(lambda x:(x[1]-x[0])/3.0, exons))) != round(sum(map(lambda x:x[1]-x[0], qblocks))):
            raise 0
        if forward:
            cdna = Seq(''.join([chr[i[0]:i[1]] for i in exons]))
        else:
            exons.reverse()
            cdna = Seq(''.join([chr[i[0]:i[1]] for i in exons]))
            cdna = cdna.reverse_complement()
        transl = cdna.translate()
        orig = ''.join(upid2seq[l[9]][int(i[0]):int(i[1])] for i in qblocks)
        ###########
        blacklist = []
        upseq = upid2seq[l[9]]
        if forward:
            eqind = 0 ## exon/qblock index
            epos = exons[eqind][0]
            qpos = int(qblocks[eqind][0])
            while qpos < qblocks[-1][1]:
                origaa = upseq[qpos]
                codon = range(epos,min(epos+3,exons[eqind][1]))
                cl = len(codon)
                if cl < 3:
                    eqind += 1
                    while exons[eqind][0] == exons[eqind][1]:  ## skip zero-length exons
                        eqind += 1
                    codon += range(exons[eqind][0],exons[eqind][0]+(3-cl))
                    epos = exons[eqind][0]+(3-cl)
                    qpos = int(ceil(qblocks[eqind][0]))
                elif epos+3 == exons[eqind][1] and eqind < len(exons)-1:
                    eqind += 1
                    while exons[eqind][0] == exons[eqind][1]:  ## skip zero-length exons
                        eqind += 1
                    epos = exons[eqind][0]
                    qpos = int(qblocks[eqind][0])
                else:
                    epos += 3
                    qpos += 1
                if str(Seq(''.join([chr[i:i+1] for i in codon])).translate()) != origaa:
                    blacklist += codon
                    #print '########', str(Seq(''.join([chr[i:i+1] for i in codon])).translate()), codon, qpos, origaa
                #else:
                #    print '@@@@@@@@', str(Seq(''.join([chr[i:i+1] for i in codon])).translate()), codon, qpos, origaa
                #print codon, epos, qpos, Seq(''.join([chr[i:i+1] for i in codon])).translate(), origaa
        else:
            eqind = 0 ## exon/qblock index
            epos = exons[-eqind-1][1]
            qpos = int(qblocks[eqind][0])
            while qpos < qblocks[-1][1]:
                codon = range(max(epos-3, exons[-eqind-1][0]), epos)
                cl = len(codon)
                origaa = upseq[qpos]
                if cl < 3:
                    eqind += 1
                    while exons[-eqind-1][0] == exons[-eqind-1][1]:  ## skip zero-length exons
                        eqind += 1
                    epos = exons[-eqind-1][1]
                    qpos = int(qblocks[eqind][0]) + 1
                    codon = range(epos-(3-cl),epos) + codon
                    epos -= (3-cl)
                elif epos-3 == exons[-eqind-1][0] and eqind < len(exons)-1:
                    eqind += 1
                    while exons[-eqind-1][0] == exons[-eqind-1][1]:  ## skip zero-length exons
                        eqind += 1
                    epos = exons[-eqind-1][1]
                    qpos = int(qblocks[eqind][0])
                else:
                    epos -= 3
                    qpos += 1
                if str(Seq(''.join([chr[i:i+1] for i in codon])).reverse_complement().translate()) != origaa:
                    blacklist += codon
                #print codon, epos, qpos, Seq(''.join([chr[i:i+1] for i in codon])).translate(), origaa
        if 0 and blacklist:
            print l[9]
            print blacklist
            x = raw_input().strip()
            while x:
                try:
                    exec(x)
                except:
                    print '!!!'
                x = raw_input().strip()

        match = 100.0*sum([(orig[i] == transl[i]) for i in xrange(max(len(transl),len(orig)))]) / max(len(transl),len(orig))
        if match == 100.0 and blacklist:
            raise 8
        if match > 99.999:
            cnt100 += 1
        elif match > 99:
            cnt99 += 1
        elif match > 95:
            cnt95 += 1
        elif match > 90:
            cnt90 += 1
        else:
            rest.append(l[9])
        strexons = ','.join(['%d-%d' % i for i in exons])
        if not forward:
            qblocks.reverse()
        strqblocks = ','.join(['%.1f-%.1f' % i for i in qblocks]).replace('.3','.1').replace('.7','.2')
        lfo = [l[13], str(exons[0][0]), str(exons[-1][1]), strexons, l[8][1], l[9], strqblocks, ','.join(map(str,blacklist))]
        fo.write('\t'.join(lfo) + '\n')
        #if cnt100 > 10:
        #    break
    fo.close()
    print (cnt100, cnt99, cnt95, cnt90, len(rest))
        
