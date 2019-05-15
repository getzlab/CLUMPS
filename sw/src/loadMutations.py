"""
Loads a MAF file into memory and, upon query, returns the mutations for specific genes.
"""

import sys

from gzip import GzipFile
from wsgiref import simple_server


#GP_FILE = '../../cl3d/dat/pancan4700.v1f.oncotator.maf.gupmap'  ## original 5k
GP_FILE = sys.argv[1]
PORT = int(sys.argv[2])
#GP_FILE='../dat/PCAWG.SNPs.maf.gp'
#GP_FILE='../dat/dlbcl.final.maf.gp'
#GP_FILE='../../bjoern/tbl1xr1_ourslit.gp'
mutTypes = set(['M']) # , 'S'  ## missense and silent

prot2gene = {'P0CI26':'TRIM49C'}
if 0:
    for l in file('../dat/uniprot.genenames.txt').readlines()[1:]:
        l = l.strip('\n').split('\t')
        if not l[1]:
            prot2gene[p] = ''
            continue
        p = l[0]
        g = l[1].split()[0]
        prot2gene[p] = g

gdata = {}  ## data acessible by gene
pdata = {}  ## data acessible by protein
if GP_FILE.endswith('.gz'):
    fi = GzipFile(GP_FILE)
else:
    fi = file(GP_FILE)
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
        if u in prot2gene:
            gene = prot2gene[u]
        else:
            gene = ''
        mi = m[1:-1]
        while not mi[-1].isdigit():  ## because some mutations are annotated as SNPs aohtough they are not
            mi = mi[:-1]
        dline = '\t'.join([l[iTtype], l[iPatient], gene, u, 'p.' + m, mi, mt])
        if not pdata.has_key(u):
            pdata[u] = []
        pdata[u].append(dline)
        if gene:
            if not gdata.has_key(gene):
                gdata[gene] = []
            gdata[gene].append(dline)


def app(env, resp):
    resp('200 OK', [('Content-type', 'text/txt')])
    pi = env['PATH_INFO']
    if pi.startswith('/g='):
        gene = pi.split('=')[-1]
        if not gdata.has_key(gene):
            return []
        return ['\n'.join(gdata[gene])]
    elif pi.startswith('/p='):
        protein = pi.split('=')[-1]
        if not pdata.has_key(protein):
            return []
        return ['\n'.join(pdata[protein])]
    else:
        return []

server=simple_server.make_server('', PORT, app)
server.serve_forever()
