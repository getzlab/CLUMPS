import random
import scipy as sp


class UniformSampler:
    def __init__(self, availUPresid):
        self.availUPresid = availUPresid
        self.availUPresidIdx = range(len(self.availUPresid))
        
    def sample(self, mi):
        return (sorted(random.sample(self.availUPresidIdx, len(mi))),  sp.random.permutation(len(mi)))
