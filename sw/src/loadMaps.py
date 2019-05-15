"""
Loads residue maps into memory and, upon query, returns the maps for specific proteins.
"""

from gzip import GzipFile
from wsgiref import simple_server

MAPS     = '../res/huniprot2pdb.run18.filt.txt.gz'
FILTMAPS = '../res/huniprot2pdb.run18.filt.txt.gz'

#MAPS     = '../res/jake_structs.txt.gz'
#FILTMAPS = '../res/jake_structs.txt.gz'

## identify the maps that passed the filter
filt = {}
fi = GzipFile(FILTMAPS)
while 1:
    l = fi.readline()
    if not l:
        break
    u1,u2,pdbch,alidt,resmap = l.strip('\n').split('\t',4)
    filt[(u1,u2,pdbch,resmap.split(':',1)[0])] = True
    
## extract the data from the unfiltered maps file
data = []
fi = GzipFile(MAPS)
while 1:
    l = fi.readline()
    if not l:
        break
    u1,u2,pdbch,alidt,resmap = l.strip('\n').split('\t',4)
    if (u1,u2,pdbch,resmap.split(':',1)[0]) in filt:
        ff = 1
    else:
        ff = 0
    d = {}
    for i in resmap.split():
        i = map(int,i.split(':'))
        d[i[0]] = i[1]
    resmap = d
    data.append([u1,u2,pdbch,alidt,resmap,ff])

u1_lines = {}
pdb_lines = {}
cnt = 0
for i in data:
    if i[0] not in u1_lines:
        u1_lines[i[0]] = []
    u1_lines[i[0]].append(cnt)
    p = i[2][:4]  ## pdb id
    if p not in pdb_lines:
        pdb_lines[p] = []
    pdb_lines[p].append(cnt)
    cnt += 1

def app(env, resp):
    resp('200 OK', [('Content-type', 'text/txt')])
    pi = env['PATH_INFO']
    lines = []
    if pi.startswith('/p='):
        u1 = pi.split('=')[-1]
        if u1 not in u1_lines:
            return []
        lines = u1_lines[u1]
    elif pi.startswith('/s='):
        pdb = pi.split('=')[-1]
        if pdb not in pdb_lines:
            return []
        lines = pdb_lines[pdb]
    ret = []
    for l in lines:
        dl = data[l]
        nl = dl[:4] + [' '.join(['%d:%d' % (x,dl[4][x]) for x in sorted(dl[4])])] + [str(dl[5])]
        ret.append('\t'.join(nl))
    return ['\n'.join(ret)]

server = simple_server.make_server('', 9000, app)
server.serve_forever()
