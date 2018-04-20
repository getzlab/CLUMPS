import scipy as sp
import sys

from gzip import GzipFile


print """ SHOULD BE A MAF WITH ALL CODING MUTATIONS (NOT JUST MISSENSES) """
MAF_FILE = sys.argv[1]


if MAF_FILE.endswith('.gz'):
    fi = GzipFile(MAF_FILE)
else:
    fi = file(MAF_FILE)
h = '#'
while h.startswith('#'):
    h = fi.readline()
hdr = h.strip('\n').split('\t')
iTumorSampleBarcode = hdr.index('Tumor_Sample_Barcode')
iVarClass = hdr.index('Variant_Classification')
iTtype = hdr.index('ttype')


codingMuts = set(["Translation_Start_Site", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_frame_Del", "In_Frame_Ins", "Missense_Mutation", "Missense", "Nonsense_Mutation", "Nonsense", "Nonstop_Mutation", "Silent", "Synonymous", "Splice_Site", "Start_Codon_DNP", "Start_Codon_Del", "Start_Codon_Ins", "Start_Codon_ONP", "Stop_Codon_DNP", "Stop_Codon_Del", "Stop_Codon_Ins", "RNA", "lincRNA", "De_novo_Start_InFrame", "Start_Codon_SNP", "Read-through", "Splice_Region", "Splice_Site_DNP", "Splice_Site_Del", "Splice_Site_Ins", "Splice_Site_ONP", "Splice_Site_SNP", "Splice_site", "Splice_site_SNP"])
nonCodingMuts = set(["3'UTR", "5'Flank", "3'Flank", "5'UTR", "De_novo_Start_OutOfFrame", "IGR", "Intron", "Non-coding_Transcript"])

di = {}
while 1:
    l = fi.readline()
    if not l:
        break
    l = l.strip('\n').split('\t')
    vc = l[iVarClass]
    if vc in codingMuts:
        pass
    elif vc in nonCodingMuts:
        continue
    else:
        raise Exception('mutation type %s unknown' % vc)
    ttype = l[iTtype] # l[iTtype].split('-')[0]
    sample = l[iTumorSampleBarcode][:12]
    if ttype not in di:
        di[ttype] = {}
    if sample not in di[ttype]:
        di[ttype][sample] = 0
    di[ttype][sample] += 1

fo = file('sampleMutFreq.txt', 'w')
fo.write('TTYPE\tSAMPLE\tMUT_COUNT\tTTYPE_RANK_SCORE\tZLOG_SCORE\n')
print 'ttypes:', len(di)
for t in di:
    ss = di[t].items()
    print '', t, len(ss)
    ss.sort(key=lambda x:x[1])
    lt = float(len(ss))
    logmf = [sp.log10(x[1]) for x in ss]  ## log mutation frequencies
    mean = sp.mean(logmf) # logmf[int(lt/10):int(9*lt/10)]
    std = sp.std(logmf)   # logmf[int(lt/10):int(9*lt/10)]
    for i in xrange(len(ss)):
        s = ss[i]
        zlog = 1.0/max(1, (sp.log10(s[1]) - mean) / std)
        fo.write('\t'.join([t, s[0],  ## ttype and sample
                            '%d' % s[1], ## raw count
                            '%g' % (i/lt), ## rank score
                            '%g' % zlog]) + '\n')

fo.close()
