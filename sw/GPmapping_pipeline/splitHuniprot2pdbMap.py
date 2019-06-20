import os
import string

from gzip import GzipFile


alpha = string.ascii_lowercase+string.digits

fi = GzipFile('../res/huniprot2pdb.run18.txt.gz')
ofiles = {ab:'../res/huniprot2pdb.run18.split/%s.gz' % ab for ab in [x+y for x in alpha for y in alpha]}
while 1:
    l = fi.readline()
    if not l:
        break
    m = l.split('\t',3)
    ab = m[2][1:3]
    fo = GzipFile(ofiles[ab], 'a')
    fo.write(l)
    fo.close()
