import __main__
__main__.pymol_argv = [ 'pymol', '-qc']   ## Quiet and no GUI

import cgi
import pickle
import pymol
import sys
import urllib2

from gzip import GzipFile
from math import log
#from lib import gene2prots, getSiftsUniprotIdsForPdb, getSiftsPdbIdsForUniprot_homol

pymol.finish_launching()

#MUTURL = 'http://cga2:8200/p='
#MUTURL = 'http://cga2:8301/p='
#MAPURL = 'http://cga3:9000/s='  ## now uses the split gz file rather than MAPURL

MUTFILE = '/mounted2/dat/TRIM49C.server'

COL = {'bgcol':'white', 'focusmappedcol':'gray70', 'focusnonmappedcol':'paleyellow', 'nonfocuscol':'aquamarine'}
#COL = {'bgcol':'black', 'focusmappedcol':'white', 'focusnonmappedcol':'paleyellow', 'nonfocuscol':'gray40'}

if 1:  ## for DBG only
    pid = sys.argv[1].lower()
    chain2proteins = sys.argv[2]
    focuschain = sys.argv[3]
    ttype = sys.argv[4]
    filter = int(sys.argv[5])
    delOtherChains = int(sys.argv[6])
else:
    args = cgi.FieldStorage()
    pid = args['pid'].value
    chain2proteins = args['chain2proteins'].value
    focuschain = args['focuschain'].value
    ttype = args['ttype'].value
    filter = int(args['filter'].value)
    delOtherChains = int(args['delOtherChains'].value)

ch2p_name = chain2proteins.replace(':','-')

if ttype == 'all':
    ttype = None
    
if chain2proteins == '.':
    chain2proteins = None

if focuschain == '.':
    focuschain = None

## parse argument chain2proteins
di = {}
if chain2proteins:
    for cu in chain2proteins.split('_'):
        cu = cu.split(':')
        if cu[0] not in di:
            di[cu[0]] = []
        di[cu[0]].append(cu[1])
chain2proteins = di

try:
    pdbstring = GzipFile('../dat/pdb/%s.pdb.gz' % pid.lower()).read()
except:
    pdbstring = urllib2.urlopen("http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s" % pid).read()

pymol.cmd.delete("all")
pymol.cmd.read_pdbstr(pdbstring, pid)
pymol.cmd.disable("all")
pymol.cmd.enable(pid)
pymol.cmd.remove('solvent')
pymol.cmd.hide('lines')
pymol.cmd.show('cartoon')
pymol.cmd.color('skyblue')
pymol.cmd.bg_color(COL['bgcol'])
pymol.cmd.set('cartoon_transparency', '0.6')
#pymol.cmd.set("surface_mode", "1")
#pymol.cmd.set("surface_color", "wheat")
#pymol.cmd.set("transparency", "0.5")
#pymol.cmd.show("surface")

## get chain2up maps
chain2up = {}
#for l in urllib2.urlopen(MAPURL + pid).readlines():
fi = GzipFile('/mounted2/res/huniprot2pdb.run6.split/%s.gz' % pid[1:3])
while 1:
    l = fi.readline()
    if not l:
        break
    l = l.strip('\n').split('\t', 5)
    if not l[2].startswith(pid):
        continue
    pdbch = l[2].split('-')
    if pdbch[1] in chain2proteins and l[0] not in chain2proteins[pdbch[1]]:
        continue
    if l[3] == '-':
        identity = 100
    else:
        attr = l[3].split()
        identity = 0
        for i in attr:
            i = i.split(':')
            if i[0] == 'pdb_identity':
                identity = float(i[1])
    rmap = dict([tuple(map(int, i.split(':'))) for i in l[4].split()])
    s = len(rmap)*identity
    if pdbch[1] not in chain2up or s > chain2up[pdbch[1]][1]:
        chain2up[pdbch[1]] = (l[0], s, rmap)

print chain2up

## get mutations
chain2muts = {}
for c in chain2up:
    md = {}
    u1,iden,rm = chain2up[c]
    #for l in map(lambda x:x.split('\t'), urllib2.urlopen(MUTURL + u1).read().split('\n')):
    for l in map(lambda x:x.split('\t'), open(MUTFILE).read().split('\n')):
        if l == [''] or not l[5]:
            continue
        if l[6][0] == 'M':
            ## ttype selection
            if ttype and not l[0].startswith(ttype):
                continue
            pos = int(l[5])
            if pos not in rm:
                continue
            ppos = rm[pos]
            if ppos not in md:
                md[ppos] = {}
            md[ppos][(l[0],l[1])] = True
    chain2muts[c] = {k:len(md[k]) for k in md if len(md[k]) >= filter}

## color the mutations
maxv = 0
for vv in chain2muts.values():
    if not len(vv):
        continue
    n = log(max(vv.values()))
    if n > maxv:
        maxv = n

maxv = min(maxv,4)  ## !!!!!!!!!!!

defcols = {}  ## defined colors
for c in pymol.cmd.get_chains(pid):
    if delOtherChains and c not in chain2proteins:
        pymol.cmd.remove("chain %s" % (c))
        continue
    if not focuschain or c == focuschain:
        if c not in chain2up:
            pymol.cmd.color(COL['nonfocuscol'], "chain %s" % (c))
            continue
        pymol.cmd.color(COL['focusnonmappedcol'], "chain %s and polymer" % c)  ## the non-mapped part will remain pale yellow
        ## get residue map coverage
        p,iden,rm = chain2up[c]
        crm = rm.values()
        rc = 0
        while rc < len(crm):
            pymol.cmd.color(COL['focusmappedcol'], "chain %s and resid %s and polymer" % (c, '+'.join(map(str,crm[rc:rc+10]))))  ## only the blasted part is colored white
            rc += 10
    elif focuschain:
        pymol.cmd.color(COL['nonfocuscol'], "chain %s" % (c))
    #pymol.cmd.bg_color(color='white')
    pymol.cmd.show('sticks', "chain %s" % c)
    #pymol.cmd.hide('sticks', "chain %s and resn ALA+ARG+ASN+ASP+CYS+GLN+GLU+GLY+HIS+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL" % (c))
    pymol.cmd.hide('sticks', "chain %s and polymer" % (c))
    if not chain2muts.has_key(c):
        continue
    msp = chain2muts[c]
    for r1 in msp:
        m = msp[r1]
        if maxv == 0:
            v = 0
        else:
            v = min(log(m), 4) / maxv  ### !!!!!!!!!!!!!!!
        v = round(v*100)
        color = 'col%g' % (v)
        v /= 100.0
        if not defcols.has_key(color):
            pymol.cmd.set_color(color, (1.0, 0.5*(1-v), 0.5*(1-v)))  ## define a new color
            defcols[color] = True
        pymol.cmd.color(color, "chain %s and resid %s" % (c, r1))
        #pymol.cmd.color('marine', "chain %s and resid %s" % (c, r1))  ## blue mutations, not scaled
        if msp[r1] > 2:
            pymol.cmd.show("sticks", "chain %s and resid %s" % (c, r1))
        else:
            pymol.cmd.show("lines", "chain %s and resid %s" % (c, r1))
            if msp[r1] == 2:
                pymol.cmd.set_bond("line_width", "3", "chain %s and resid %s" % (c, r1))



pymol.cmd.reset()  ## re-centers the molecule, in case we've deleted something

FNAME = "%s_%s_%s_%s_%s%s.pse" % (pid, (ttype or 'all'), ch2p_name, focuschain, filter, delOtherChains)

pymol.cmd.save("/mounted/%s" % FNAME)

pymol.cmd.quit()

print "Content-Type: text/html"
print ''
html = """
<!DOCTYPE html>
<html>
<title>VIZ</title>
<head>
<meta charset="utf-8">
<script type="text/javascript" src="/js/jmol/JSmol.min.js"></script>
</head>
<body>
<script type="text/javascript"> 
$(document).ready(function() {
Info = {
	width: $(window).width()-30,
	height: $(window).height()-40,
	debug: false,
	color: "0xC0C0C0",
	disableJ2SLoadMonitor: true,
	disableInitialConsole: true,
	//addSelectionOptions: true,
	serverURL: "http://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php",
	use: "HTML5",
	readyFunction: null,
	script: "load 3d/%s"  // + window.location.href.split('?s=')[1]
}
$("#mydiv").html(Jmol.getAppletHtml("jmolApplet0",Info))
});
</script>
<span id=mydiv></span>
<script 
type="text/javascript">Jmol._alertNoBinary=false;Jmol._binaryTypes=[".pse"];</script>
<a href="javascript:Jmol.script(jmolApplet0, 'spin on')">spin on</a>
<a href="javascript:Jmol.script(jmolApplet0, 'spin off')">spin off</a>
<br>
TODO: CHAIN ANNOTATION
</body>
</html>

""" % FNAME

print html

#subprocess.call('sleep 20; rm ')
