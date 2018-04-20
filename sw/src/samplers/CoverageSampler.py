import jpype
import scipy as sp

from GPmapper import GPmapper


#jpype.startJVM(jpype.getDefaultJVMPath(), '-ea', '-Djava.class.path=/xchip/cga/reference/mutsig_params/FixedWidthBinary.jar')
jpype.startJVM(jpype.getDefaultJVMPath(), '-ea', '-Djava.class.path=/sw/etc/FixedWidthBinary.jar')
FixedWidthBinary = jpype.JPackage('org').broadinstitute.cga.tools.seq.FixedWidthBinary

class CoverageSampler():
    def __init__(self, availUPresid, upid, covtrack):
        self.availUPresid = availUPresid
        self.availUPresidIdx = range(len(self.availUPresid))
        self.upid = upid
        self.gpm = GPmapper()
        self.fwb = FixedWidthBinary(covtrack)
        self.covprobs = self.calcProbabilities()  ## array of probabilities (derived from coverage) for each position
    
    def __del__(self):
        jpype.shutdownJVM()
    
    def calcProbabilities(self):
        resid2coverage = {}
        genelocs = self.gpm.prot2gen[self.upid]   ## more than one gene may exist! thus, genelocs may have len>1
        aur = set(self.availUPresid)  # searchable availUPresid
        for chr, gi in genelocs:
            if chr[3:] == 'X':
                chrn = 23
            elif chr[3:] == 'Y':
                chrn = 24
            elif chr[3:] == 'M':
                chrn = 25
            else:
                chrn = int(chr[3:])
            for ei in xrange(len(self.gpm.gen2prot[chr][gi][2])):  ## for each exon
                for pos in xrange(self.gpm.gen2prot[chr][gi][2][ei][0], self.gpm.gen2prot[chr][gi][2][ei][1]):
                    res = self.gpm.mapGPchange(chr[3:], pos+1, gi, ei, None)
                    for r in res:
                        if r[0] != self.upid:
                            continue
                        if r[1] not in aur:
                            continue
                        if r[1] not in resid2coverage:
                            resid2coverage[r[1]] = []
                        cov = self.fwb.get(chrn,pos+1)
                        resid2coverage[r[1]].append(cov)
        coverage = []
        for i in self.availUPresid:
            if i in resid2coverage:
                cov = sum(resid2coverage[i]) / float(len(resid2coverage[i]))
                if cov < 0:
                    cov = 0
            else:
                cov = 0
            coverage.append(cov)
        if sum(coverage) == 0:
            coverage = [1.0]*len(coverage)
        return [i/sum(coverage) for i in coverage]
    
    def sample(self, mi):
        return (sorted(sp.random.choice(self.availUPresidIdx, len(mi), p=self.covprobs, replace=False)), sp.random.permutation(len(mi)))
   
