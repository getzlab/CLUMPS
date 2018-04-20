import os
import sys


RESDIR = sys.argv[1]
CHUNKPATHS = sys.argv[2].split(',')

fo1 = file('candidates_chunkpaths.txt', 'w')
fo2 = file('candidates_lineids.txt', 'w')
fo3 = file('candidates_idx.txt', 'w')

n = 0
for fn in os.listdir(RESDIR):
    dat = file(RESDIR + '/' + fn).readlines()
    if len(dat) > 0 and dat[-1] == '#1\n':
        chunk_line = fn.split('_', 1)[0]
        fo1.write(CHUNKPATHS[int(chunk_line.split('-')[0])] + '\n')
        fo2.write(chunk_line.split('-')[1] + '\n')
        fo3.write('%d\n' % n)
        n += 1

fo1.close()
fo2.close()
fo3.close()