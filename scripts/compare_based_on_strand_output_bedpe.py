#!/usr/bin/env python
from __future__ import division

import sys
import argparse

#TODO what do we do with missing data? Right now is counted as discordance or match depending on if the other sample has data or not

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

    def toString(self, status, sample, region):
        return '\t'.join(str(v) for v in (self.c1, self.s1, self.e1, self.c2, self.s2, self.e2, self.name, self.score, self.o1, self.o2, self.svtype, self.filter, self.orig_name1, self.orig_ref1, self.orig_alt1, self.orig_name2, self.orig_ref2, self.orig_alt2, self.info1+";Regions="+region, self.info2, self.format+":RE", self.misc[sample]+":"+status))+"\n"

class PairToPair(object):
    def __init__(self, sample_header, out_dir, region):
        self._load_samples(sample_header)
        self.outFiles = [ open(out_dir+"/"+x+".bedpe", 'a') for x in self.samples ]
        self.region = region
        self.deferred = dict()

    def _load_samples(self, sample_header):
        with open(sample_header) as fh:
            for line in fh:
                fields = line.rstrip().split('\t')
                self.samples = fields[21:]
                self.num_samples = len(self.samples)
                self.index = 21 + self.num_samples

    def process(self, fields_list):
        bedpe1 = Bedpe(fields_list[0:self.index])
        if len(fields_list) == self.index:
            # this is unique to the first file
            # process the line and record variants as unique
            for index, x in enumerate(bedpe1.misc):
                fields = x.rstrip().split(':', 1)
                #self.results[index].unique += 1
                if fields[0] in ('0/1', '1/1'):
                    #should be a variant
                    self.outFiles[index].write(bedpe1.toString("0-only", index, self.region))
        else:
            # this MAY be a match
            # check strands and then determine class
            # for each sample record results
            bedpe2 = Bedpe(fields_list[self.index:])
            if strands_match(bedpe1.o1, bedpe1.o2, bedpe2.o1, bedpe2.o2):
                # matched by strand
                if bedpe1.name in self.deferred:
                    del self.deferred[bedpe1.name]
                for index in range(self.num_samples):
                    discordant_type = not types_match(bedpe1.svtype, bedpe2.svtype)
                    gt1 = gt_for_index(bedpe1, index)
                    gt2 = gt_for_index(bedpe2, index)
                    if gt1 in ('0/1', '1/1') or gt2 in ('0/1', '1/1'):
                        if gt1 == gt2:
                            self.outFiles[index].write(bedpe1.toString("match", index, self.region))
                        else:
                            self.outFiles[index].write(bedpe1.toString("discordant", index, self.region))
            else:
                self.deferred[bedpe1.name] = bedpe1

    def process_shadow_unique(self):
        # count up things that overlapped but ended up being unique
        for bedpe in self.deferred.values():
            for index, x in enumerate(bedpe.misc):
                gt = gt_for_index(bedpe, index)
                if gt in ('0/1', '1/1'):
                    self.outFiles[index].write(bedpe.toString("0-only", index, self.region))

    def finalize(self):
        for index, f in enumerate(self.outFiles):
            f.close()

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
    parser.add_argument('-s', '--sample-header', metavar='<FILE>', type=str, help='Filename of file containing sample header of the bedpe')
    parser.add_argument('-o', '--outDir', metavar='<STR>', type=str, help='Directory where output files should be written')
    parser.add_argument('-r', '--region', metavar='<STR>', type=str, help='Region to record in INFO field')
    args = parser.parse_args()
    
    p = PairToPair(args.sample_header, args.outDir, args.region)
    try:
        for line in sys.stdin:
            fields = line.rstrip().split('\t')
            p.process(fields)
        p.process_shadow_unique()
        p.finalize()
    except ValueError as e:
        sys.exit(e)

