#!/bin/bash

set -eo pipefail
BEDTOOLS=$1 #path to bedtools executable
pad1bp () {
    cat "$1" | perl -ape '$F[1] -= 1; $F[2]+=1; $F[4] -= 1; $F[5] += 1; $_ = join("\t", @F)."\n"' > $1.padded
}

OUTPUT_DIR=$2 #output directory
OUTPUT=$OUTPUT/results.txt
DIFFICULT=$3 #BED file of difficult regions
REST=$4 #BED file of remaining regions
SAMPLE_HEADER=$5 #BEDPE header line listing all samples in correct order
echo "Pipeline1	Pipeline2	Sample	Class	Count" > $OUTPUT

#file containing list of pipelines and paths to final BEDPE files from each pipeline,
#separated by tabs
BEDPE_LIST=$6
while read pipeline1 pipeline1_file
do
    while read pipeline2 pipeline2_file
    do
        if [ $pipeline1 == $pipeline2 ]; then
            continue
        fi
        out_dir=$OUTDIR/${pipeline1}_${pipeline2}
        mkdir -p $out_dir/logs

        for sample in $( cut -f 22- $SAMPLE_HEADER ); do
            grep "^##" $pipeline1_file > $out_dir/$sample.bedpe
            echo "##FORMAT=<ID=RE,Number=1,Type=String,Description=\"Matching status in pairwise comparison\">" >> $out_dir/$sample.bedpe
            echo "##INFO=<ID=Regions,Number=1,Type=String,Description=\"Genomic regions\">" >> $out_dir/$sample.bedpe
            grep "^#CHROM_A" $pipeline1_file | cut -f 1-21 | sed "s/$/\t$sample/" >> $out_dir/$sample.bedpe
        done
        $BEDTOOLS pairtopair -is -a $pipeline1_file.padded -b $pipeline2_file.padded -type both -slop 50 | $BEDTOOLS pairtobed -a - -b $DIFFICULT -type either | rev | cut -f 4- | rev | sort -u | python compare_based_on_strand_output_bedpe.py -s $SAMPLE_HEADER -o $out_dir -r hard
        $BEDTOOLS pairtopair -is -a $pipeline1_file.padded -b $pipeline2_file.padded -type both -slop 50 | $BEDTOOLS pairtobed -a - -b $DIFFICULT -type neither | $BEDTOOLS pairtobed -a - -b $REST -type either | rev | cut -f 4- | rev | sort -u | python compare_based_on_strand_output_bedpe.py -s $SAMPLE_HEADER -o $out_dir -r medium
        $BEDTOOLS pairtopair -is -a $pipeline1_file.padded -b $pipeline2_file.padded -type both -slop 50 | $BEDTOOLS pairtobed -a - -b $DIFFICULT -type neither | $BEDTOOLS pairtobed -a - -b $REST -type neither | sort -u | python compare_based_on_strand_output_bedpe.py -s $SAMPLE_HEADER -o $out_dir -r easy
        $BEDTOOLS pairtopair -is -a $pipeline1_file.padded -b $pipeline2_file.padded -type notboth -slop 50 | $BEDTOOLS pairtobed -a - -b $DIFFICULT -type either | rev | cut -f 4- | rev | sort -u | python compare_based_on_strand_output_bedpe.py -s $SAMPLE_HEADER -o $out_dir -r hard
        $BEDTOOLS pairtopair -is -a $pipeline1_file.padded -b $pipeline2_file.padded -type notboth -slop 50 | $BEDTOOLS pairtobed -a - -b $DIFFICULT -type neither | $BEDTOOLS pairtobed -a - -b $REST -type either | rev | cut -f 4- | rev | sort -u | python compare_based_on_strand_output_bedpe.py -s $SAMPLE_HEADER -o $out_dir -r medium
        $BEDTOOLS pairtopair -is -a $pipeline1_file.padded -b $pipeline2_file.padded -type notboth -slop 50 | $BEDTOOLS pairtobed -a - -b $DIFFICULT -type neither | $BEDTOOLS pairtobed -a - -b $REST -type neither | sort -u | python compare_based_on_strand_output_bedpe.py -s $SAMPLE_HEADER -o $out_dir -r easy

    done < $BEDPE_LIST
done < $BEDPE_LIST
