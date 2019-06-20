"""
Selects a set of minimally overlapping, representative structures
"""

from gzip import GzipFile
#from prody import parsePDBStream

## read structures and residue maps from the mappingPipeline
struct2resol = {}
u1structs = {}
fi = GzipFile('../res/huniprot2pdb.run18.txt.gz', 'r')
#fi = GzipFile('../res/q.gz', 'r')
lnc = -1  ## line counter
while 1:
    ln = fi.readline()
    if not ln:
        break
    lnc += 1
    u1,u2,pdbch,alidt,resmap = ln.strip('\n').split('\t', 4)
#    if u1 not in ['Q00653']:
#        continue

    #pdb = pdbch.split('-')[0]
    #hdr = parsePDBStream(GzipFile('../dat/pdb/%s.pdb.gz' % pdb), header=True, model=0) ## parse only pdb header
    #if 'resolution' in hdr:
    #    resol = float(hdr['resolution'])
    #else:
    #    resol = 100
    iden = 0
    if u1 == u2:
        direct = True
        iden = 100
    else:
        direct = False
        for i in alidt.split():
            if i.startswith('pdb_identity:'):
                iden = float(i.split(':')[1])
    if iden < 10:  ## WE DON'T REALLY NEED THIS IF THE INPUT IS GOOD 
        continue   ##
    #mr = map(lambda x:int(x.split(':')[0]), resmap.split())
    #mn = min(mr)
    #mx = max(mr)
    #lenmr = len(mr)
    lenmr = resmap.count(':')
    if lenmr < 3:
        continue
    mn = int(resmap.split(' ',1)[0].split(':')[0])
    mx = int(resmap.rsplit(' ',1)[1].split(':')[0])
    if lenmr < 0.5*(mx-mn):  ## half of the sequence (between start and end in pdb) is missing, e.g. residues missing from structure (e.g. 4ql6-A)
        #print '##', u1, u2, pdbch
        continue
    if u1 not in u1structs:
        u1structs[u1] = []
    u1structs[u1].append( ((u2,pdbch), direct, (mn,mx,lenmr), iden, lnc) ) # , resol

fi.close()

## select all native structures and fill the gaps with non-native ones
u1structs_filt = {}
cnt = 0
u1structslen = len(u1structs)
print u1structslen, "u1's"
for u1 in u1structs:
    cnt += 1
    if not cnt % 10:
        print cnt
    u1structs_filt[u1] = []
    ss = []  ## all homologous structures
    for i in xrange(len(u1structs[u1])):
        if u1structs[u1][i][1]:  ## native structure; directly passes the filter
            u1structs_filt[u1].append(u1structs[u1][i])
        else:
            ss.append(u1structs[u1][i])
    ss.sort(key=lambda x:x[2][2], reverse=True)  ## sorted by mapped residues, longest first
    hclusters = []  ## clusters of homologous structures
    for i in ss:
        newcl = True  ## does this found a new cluster
        for cl in hclusters:
            num_overlapping = 0  ## how many members of the cluster cl does this struct overlap with (at least 90% jaccard overlap)
            for j in cl:
                ol = max(0, min(i[2][1], j[2][1]) - max(i[2][0], j[2][0]))
                if float(ol)/((i[2][1]-i[2][0]) + (j[2][1]-j[2][0]) - ol) > 0.9:
                    num_overlapping += 1
            if num_overlapping/float(len(cl)) >= 0.5:  ## overlaps with at least half of the members of the cluster by at least 90%
                cl.append(i)
                newcl = False
                break
        if newcl:
            hclusters.append([i])
    for cl in hclusters:
        cl.sort(key=lambda x:x[3]*x[2][2], reverse=True)  ## sort by identity * mapped length
        for j in cl:
            ok = True
            for i in u1structs_filt[u1]:
                ol = max(0, min(i[2][1], j[2][1]) - max(i[2][0], j[2][0]))
                if float(ol)/min(i[2][2], j[2][2]) > 0.1:
                    ok = False
            if ok:
                u1structs_filt[u1].append(j)  ## this member of cl is selected; move on to next cluster
                break
    #u1structs_filt[u1] = {i[0]:True for i in u1structs_filt[u1]}  ## id is (u1,pdbch)
    u1structs_filt[u1] = {i[4]:True for i in u1structs_filt[u1]}   ## id is line number

fo = GzipFile('../res/huniprot2pdb.run18.filt.txt.gz', 'w')
fi = GzipFile('../res/huniprot2pdb.run18.txt.gz', 'r')
lnc = -1
while 1:
    ln = fi.readline()
    if not ln:
        break
    lnc += 1
    u1,u2,pdbch,alidt,resmap = ln.strip('\n').split('\t', 4)
    if u1 in u1structs_filt and lnc in u1structs_filt[u1]:
        fo.write(ln)


fi.close()
fo.close()
