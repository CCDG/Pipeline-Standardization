#!/usr/bin/env python

import re
import sys

class Uninformative(RuntimeError):
    pass

class MIE(object):
    def __init__(self):
        pass

    def __call__(self, fields):
        #assume last 3 columns are the genotypes in KID DAD MOM order
        kid, dad, mom = [ re.split('[/|]', x) for x in fields[-3:] ]
        dad_alleles, mom_alleles, kid_alleles = map(set, (dad, mom, kid))
        if len(dad_alleles) == 2 and len(mom_alleles) == 2 and len(dad_alleles.union(mom_alleles).union(kid_alleles)) == 2:
            # both are het and not a multi-allelic site, this site is UNINFORMATIVE assuming input includes no missing genotypes
            raise Uninformative
        a1, a2 = kid
        classification = 'MIE'
        if (a1 in dad_alleles and a2 in mom_alleles) or (a2 in dad_alleles and a1 in mom_alleles):
            # NOT AN MIE
            classification = 'NONE'
        return fields[:-3] + [ classification ]

if __name__ == '__main__':
    mie_classify = MIE()
    for line in sys.stdin:
        fields = line.rstrip().split('\t')
        try:
            print '\t'.join(mie_classify(fields))
        except Uninformative:
            pass

