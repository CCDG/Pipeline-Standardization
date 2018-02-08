#!/usr/bin/env python
from __future__ import division

import sys
import argparse

#TODO what do we do with missing data? Right now is counted as discordance or match depending on if the other sample has data or not

class Result(object):
    def __init__(self):
        self.unique = 0
        self.match = 0
        self.discordant = 0
        self.match_discordant_type = 0
        self.discordant_discordant_type = 0
        self.unmatched = dict()

class Bedpe(object):
    def __init__(self, bed_list):
        self.c1 = bed_list[0]
        self.s1 = int(bed_list[1])
        self.e1 = int(bed_list[2])
        self.c2 = bed_list[3]
        self.s2 = int(bed_list[4])
        self.e2 = int(bed_list[5])
        self.name = bed_list[6]
        self.score = bed_list[7]
        self.o1 = bed_list[8]
        self.o2 = bed_list[9]
        self.svtype = bed_list[10]
        self.filter = bed_list[11]
        self.orig_name1 = bed_list[12]
        self.orig_ref1 = bed_list[13]
        self.orig_alt1 = bed_list[14]
        self.orig_name2 = bed_list[15]
        self.orig_ref2 = bed_list[16]
        self.orig_alt2 = bed_list[17]
        self.info1 = bed_list[18]
        self.info2 = bed_list[19]
        # remainder will be entries then second bedpe entry will follow
        self.format = bed_list[20]
        self.misc = bed_list[21:]

class PairToPair(object):
    def __init__(self, label):
        self.label = label
        self.results = Result()
        self.deferred = dict()

    def process(self, fields_list):
        bedpe1 = Bedpe(fields_list[0:22])
        if len(fields_list) == 22:
            # this is unique to the first file
            # process the line and record variants as unique
            for index, x in enumerate(bedpe1.misc):
                fields = x.rstrip().split(':', 1)
                if fields[0] in ('0/1', '1/1'):
                    #should be a variant
                    self.results.unique += 1
        else:
            # this MAY be a match
            # check strands and then determine class
            # for each sample record results
            bedpe2 = Bedpe(fields_list[22:])
            if strands_match(bedpe1.o1, bedpe1.o2, bedpe2.o1, bedpe2.o2):
                # matched by strand
                if bedpe1.name in self.deferred:
                    del self.deferred[bedpe1.name]
                for index in range(1):
                    discordant_type = not types_match(bedpe1.svtype, bedpe2.svtype)
                    gt1 = gt_for_index(bedpe1, index)
                    gt2 = gt_for_index(bedpe2, index)
                    if gt1 in ('0/1', '1/1') or gt2 in ('0/1', '1/1'):
                        if gt1 == gt2:
                            self.results.match += 1
                            if discordant_type:
                                self.results.match_discordant_type += 1
                        else:
                            self.results.discordant += 1
                            if discordant_type:
                                self.results.discordant_discordant_type += 1
            else:
                self.deferred[bedpe1.name] = bedpe1

    def process_shadow_unique(self):
        # count up things that overlapped but ended up being unique
        for bedpe in self.deferred.values():
            for index, x in enumerate(bedpe.misc):
                gt = gt_for_index(bedpe, index)
                if gt in ('0/1', '1/1'):
                    self.results.unique += 1

    def summarize(self):
        unique_label = '-'.join((self.label, 'only'))
        print '\t'.join((unique_label, str(self.results.unique)))
        print '\t'.join(('match', str(self.results.match)))
        print '\t'.join(('discordant', str(self.results.discordant)))
        print '\t'.join(('match_discordant_type', str(self.results.match_discordant_type)))
        print '\t'.join(('discordant_discordant_type', str(self.results.discordant_discordant_type)))

def gt_for_index(bedpe, index):
    format_fields = bedpe.misc[index].rstrip().split(':', 1)
    return format_fields[0]


def strands_match(eval_strand_a, eval_strand_b, compare_strand_a, compare_strand_b):
    if '.' in [eval_strand_a, eval_strand_b, compare_strand_a, compare_strand_b]:
        return False
    if eval_strand_a != eval_strand_b:
        return (eval_strand_a == compare_strand_a 
                and eval_strand_b == compare_strand_b)
    else:
        return compare_strand_a == compare_strand_b

def types_match(type1, type2):
    return type1 == type2

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter output of pairtopair based on strand.')
    parser.add_argument('-l', '--label', metavar='<STR>', type=str, help='Label for a file')
    args = parser.parse_args()
    
    p = PairToPair(args.label)
    try:
        for line in sys.stdin:
            fields = line.rstrip().split('\t')
            p.process(fields)
        p.process_shadow_unique()
        p.summarize()
    except ValueError as e:
        sys.exit(e)

