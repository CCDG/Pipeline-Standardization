#Alignment pipeline standards

##Reference genome version
Summary: We agreed that each center should use exactly the same reference genome.

Standard:
* GRCh38DH, 1000 Genomes Project version
* Includes the standard set of chromosomes and alternate sequences from GRCh38
* Includes the decoy sequences
* Includes additional alternate versions of the HLA locus

##Alignment
Summary: We agreed that each center should use exactly the same alignment strategy

Standard:
* Aligner: BWA-MEM
* Version: We will use 0.7.15 (https://github.com/lh3/bwa/releases/tag/v0.7.15)
* Standardized parameters:
    * Do not use `-M` since it causes split-read alignments to be marked as "secondary" rather than "supplementary" alignments, violating the BAM specification
    * Use `-K 100000000` to achieve deterministic alignment results
    * Use `-Y` to force soft-clipping rather than default hard-clipping of supplementary alignments
    * Include a `.alt` file for consumption by BWA-MEM; do not perform post-processing of alternate alignments
* Optional parameters (may be useful for convenience and not expected to alter results):
    * `-p` (for interleaved fastq)
    * `-C` (append FASTA/FASTQ comment to SAM output)
    * `-v` (logging verbosity)
    * `-t` (threading)
    * `-R` (read group header line)
* Post-alignment modification:
    * In order to reduce false positive calls due to bacterial contamination randomly aligning to the human genome, the Broad Institute has started marking (setting 0x4 bit in the SAM flag) reads (and their mates) if the following conditions apply:
        1. The primary alignment has less than 32 aligned bases
        2. The primary alignment is (soft) clipped on both sides
    * This filtering is optional
    * The original mapping information will be encoded in a tag on the marked reads using the same format as the SA tag in the BAM specification.
    * Modification of other flags after alignment will not be performed.

Notes:
* Our current understanding is that `-K 100000000` produces completely deterministic alignment results, regardless of the number of threads used, assuming that identical fastq and reference genome files are used. This should be verified.
* We verified that reference genome index files are not a source of stochastic alignment results (via md5 checksum of index files produced at different centers)
* `-Y` was suggested by the Michigan group to help ensure compatibility with SRA.
* Broad can provide some examples of problematic regions (not whole samples) where bacterial contamination is problematic

##Duplicate marking
Summary: This processing step is a source of considerable variability among centers, with three different tools being used at the beginning of this exercise: Picard (Broad, NYGC, Baylor, WashU), BamUtil (Michigan). These tools differ in their behavior at supplementary alignments, at “orphan” alignments where one of the two reads is unmapped, and based on whether they select the “best” read-pair in a set of duplicates (Picard & BamUtil), or the first read-pair (Samblaster). We agreed that it was acceptable for different centers to use different tools, so long as the same number of reads were marked duplicate and results were functionally equivalent.

Standard:
* Match Picard’s current definition of duplicates for primary alignments where both reads of a pair align to the reference genome. Both Samblaster and BamUtil already attempt to match Picard for this class of alignments.
* If a primary alignment is marked as duplicate, then all supplementary alignments for that same read should also be marked as duplicates. Both Picard and BamUtil have modified to exhibit this behavior.
* Orphans will be marked as duplicates if there’s another read with the same alignment (mated, or orphaned)
* The unmapped mate of duplicate orphan reads is required to also be marked as a duplicate.
* It is not a requirement for duplicate marking software to choose the best pair based on base quality sum, but results must be functionally equivalent.

Notes:
* The Broad has modified Picard MarkDuplicates to meet these specifications when run on a queryname sorted input file. This functionality was available beginning with version 2.4.1 of Picard.
* WashU and Baylor have moved to using Picard MarkDuplicates >2.4.1 to meet this specification.
* Michigan is modifying bamUtil to meet these specifications
* We agreed that if a primary alignment is marked as duplicate, then all secondary alignments for that read should also be marked as duplicates. However, given that no secondary alignments will exist using our proposed alignment strategy, we decided that it should be optional for different groups to incorporate this behavior into their software.
* There was a discussion about whether duplicate marking should be deterministic. We did not reach a decision on this.
* We have discussed the preferred behavior for marking duplicates in datasets with multiple sequencing libraries and have decided that this is a minor concern given that very few samples should have multiple libraries. Currently MarkDuplicates supports multiple libraries with the caveat that the term “Library” isn’t exactly defined (consider a technical replicate that starts somewhere in the middle of the LC process, how early must it be to be called a different library?)

##Indel realignment
Summary: We agreed that this computationally expensive data processing step is dispensable given the state of current Indel detection algorithms, and should be removed.

##Base quality score recalibration
Summary: Currently, this step is performed by all centers except WashU, using two different tools: GATK (Broad, NYGC & Baylor) and BamUtil (Michigan). There was discussion about dropping BQSR given evidence from WashU and Michigan that the impact on variant calling performance is minimal. However, given that this project will involve combined analysis of data from multiple centers and numerous sequencers, generated over multiple years, and that we cannot ensure the consistency of Illumina base-calling software over time, we decided that it is preferable to perform BQSR.

Standard:
* We will use the following files from the [GATK hg38 bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/) for the site list:
    * Homo_sapiens_assembly38.dbsnp138.vcf
    * Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    * Homo_sapiens_assembly38.known_indels.vcf.gz
* The recalibration table may optionally be generated using only the autosomes (chr1-chr22)
* Downsampling of the reads is allowed (but optional) to generate the recalibration table

Command line:
For users of GATK, the following command line options should be utilized for the BaseRecalibrator tool:
```
-R ${ref_fasta} \
-I ${input_bam} \
-O ${recalibration_report_filename} \
-knownSites "Homo_sapiens_assembly38.dbsnp138.vcf" \
-knownSites “Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
-knownSites “Homo_sapiens_assembly38.known_indels.vcf.gz”
```

For users of GATK, the following command line options are optional for efficiency and can be utilized for the BaseRecalibrator tool:
```
--downsample_to_fraction .1 \
    -L chr1 \
    -L chr2 \
    -L chr3 \
    -L chr4 \
    -L chr5 \
    -L chr6 \
    -L chr7 \
    -L chr8 \
    -L chr9 \
    -L chr10 \
    -L chr11 \
    -L chr12 \
    -L chr13 \
    -L chr14 \
    -L chr15 \
    -L chr16 \
    -L chr17 \
    -L chr18 \
    -L chr19 \
    -L chr20 \
    -L chr21 \
    -L chr22
```

For users of GATK, the following command line options are optional:
* `--interval_padding 200`
* `-rf BadCigar`
* `--preserve_qscores_less_than 6`
* `--disable_auto_index_creation_and_locking_when_reading_rods`
* `--disable_bam_indexing`
* `-nct`
* `--useOriginalQualities`

Notes:
* NYGC and WashU have expressed concern about using two different BQSR tools given our decision to compress base quality scores and not retain the original values (see below). Under this scenario, BQSR becomes an important and irreversible data processing step.
* We either need to use a single tool, or we need to be absolutely sure that the different tools are performing BQSR in identical or nearly identical ways.
* Since BQSR includes BAQ, which is very compute intensive, broad is evaluating whether this can be removed while retaining the quality of the calls.

##Base quality score binning scheme
Summary: We agreed that additional base quality score compression was required to reduce file size, and that it should be possible to achieve this with minimal adverse impacts on variant calling.

Standard:
* 4-bin quality score compression. The 4-bin scheme is 2-6, 10, 20, 30. The 2-6 scores correspond to Illumina error codes and will be left as-is by recalibration.
* Bin base quality scores by rounding off to the nearest bin value, in probability space. This feature is already implemented in the current version of GATK.

Command line:
For users of GATK, the following command line options should be utilized for the PrintReads (GATK3) or ApplyBQSR (GATK4) tool:
```
-R ${ref_fasta} \
-I ${input_bam} \
-O ${output_bam_basename}.bam \
-bqsr ${recalibration_report} \
-SQQ 10 -SQQ 20 -SQQ 30
```

For users of GATK, the following command line options are optional:
* `--globalQScorePrior -1.0`
* `--preserve_qscores_less_than 6`
* `--disable_indel_quals`
* `--useOriginalQualities`
* `-nct`
* `-rf BadCigar`
* `--emit_original_quals`
* `--createOutputBamMD5`
* `--addOutputSAMProgramRecord`

Notes:
* The average file size for a 30X CRAM file compressed using Samtools is 13.7 Gb using 3 bins and 16.8 Gb using 4 bins. These sizes are based on 50 WGS datasets from WashU, downsampled to exactly 30X coverage. All but two datasets are from the HiSeqX. Original quality scores were not retained.
* After some discussion of alternate binning schemes, the GATK scheme using 2-6, 10, 20, 30 was decided on 6/1.
* After additional discussion the 40 bin was added on 10/10 and subsequently removed on 11/14
* What to do with bases recalibrated to have qualities between 2-6 is still under discussion.
* Broad and Michigan have expressed a preference for 4 bins due to the greater information content and minor difference in file size.
* Variant calling performance was compared using 3 or 4 bins and 30X and 20X data.

##File format
Summary: We agreed that each center should use the same file format, while retaining flexibility to include additional information for specific centers or projects.

Standard:
* Lossless CRAM. Upon conversion to BAM, the BAM file should be valid according to Picard’s ValidateSamFile.
* Read group (@RG) tags should be present for all reads. The header for the RG should contain minimally the ID tag, PL tag, PU tag, SM tag, and LB tag. The CN tag is recommended. Other tags are optional. The ID tag must be unique within the CRAM. ID tags may be freely renamed to maintain uniqueness when merging CRAMs. No assumptions should be made about the permanence of RG IDs. The PL tag should indicate the instrument vendor name according to the SAM spec (CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT, ONT, and PACBIO); PL values are case insensitive. The PU tag is used for grouping reads for BQSR and should uniquely identify reads as belonging to a sample-library-flowcell-lane (or other appropriate recalibration unit) within the CRAM file. PU is not required to contain values for fields that are uniform across the CRAM (e.g., single sample CRAM or single library CRAM). The PU tag is not guaranteed to be sufficiently informative after merging with other CRAMs, and anyone performing a merge should consider modifying PU values appropriately. SM should contain the individual identifier for the sample (e.g., NA12878) without any other process or aliquot-specific information. The LB tag should uniquely identify the library for the sample; it must be present even if there is only a single library per sample or CRAM file. If the PM tag is used, values should conform to one of the following (for Illumina instruments): “HiSeq-X”, “HiSeq-4000”, “HiSeq-2500”, “HiSeq-2000”, “NextSeq-500”, or “MiSeq”.
* Retain original query names.
* Retain @PG records for bwa, duplicate marking, quality recalibration, and any other tools that was run on the data.
* Retain the minimal set of tags (RG, MQ, MC and SA).  NOTE: an additional tool may be needed to add the MQ and MC tags if none of the tools add these tags otherwise.  One option is to pipe the alignment through samblaster (https://github.com/GregoryFaust/samblaster) with the options `-a --addMateTags` as it comes out of BWA
* Groups can add custom tags as needed.
* Do not retain the original base quality scores (OQ tag).

Notes:
* We have discussed the minimal set of tags to include. RG, MQ, MC and SA are all required. We will get MD and NM for free by using CRAM. OQ will be explicitly excluded, but all additional tags may be included by any group.
* The use of htsjdk/picard for converting bam to cram is not currently condoned. We are working on making picard output equivalent to samtools, but are not there yet. Broad currently converts using samtools.
* Should retention of original query names be mandatory, or optional? Query names will not be retained in current cSRA format.
* dbGaP and cSRA currently accept CRAM files
*  it is recommended that users use samtools version 1.3.1 to convert from bam/sam to Cram (not picard). Users that would like to convert back from cram to bam (and want to avoid ending up with a working, but invalid bam) need to either convert to sam and then to bam (piping works) or compile samtools with HTSLib version 1.3.2. To enable this you need to: configure the build of samtools with the parameter `--with-htslib=/path/to/htslib-1.3.2`.

#Functional equivalence evaluation
Summary:  We agreed that all pipelines used for this effort need to be validated as functionally equivalent.

Resources available:
* NA12878 (WashU - 2 replicates)
* NA12878 - 2 replicates from same DNA, NA12891, NA12892 (Broad)
* NA19238 x 2 (2 different centers, but particulars unable to be shared)
* Luyha sample + 1 more non-European (NYGC)
* NA12878 replicates downsampled to ~24X (the two Broad replicates and 1 of the WashU replicates).
* NA19238 downsampled to ~30X to match the coverage of the other NA19238 sample.
* A software tool to measure functional equivalence based on two CRAM files.
* A software tool to measure functional equivalence based on two VCF files.

Proposed functional equivalence metrics:

1. CRAM: identical counts for different alignment classes (e.g., run samtools flagstat).

2. CRAM: correlation between individual base quality scores at read level. This aims to test consistency of BQSR.

3. CRAM: correlation between summed base quality scores at each reference genome coordinate. This aims to test differences arising from duplicate marking tools that retain a different set of non-duplicate reads, with potentially different base qualities and clipping information.

4. CRAM: similarity between the overall distributions of base qualities

5. VCF: variant genotype concordance (separate out common/1000 genomes vs ‘rare’)

6. VCF: variant quality score correlation

7. VCF: genotype quality correlation

8. VCF: for NA12878, sensitivity to gold standard variants (GiaB version 3.2, lifted over to b38)

9. VCF: for NA12878, FDR in GiaB high confidence regions

10. VCF: for NA12878, genotype concordance with gold standard variants

11. VCF: for NA12878-NA12891-NA12892 trio, # of Mendelian errors (may want to use jointly genotyped variants for this)

12. SV: lumpy equivalence metrics TBD (possibly Mendelian errors, genotype concordance, comparison to NA12878 validation set that needs to be lifted over)

Proposed procedure:

1. Each center generates CRAM files for each of the 148 data sets (10 original plus 4 downsampled) and shares with the group. The CRAM tab in the manifest will be updated upon sharing.

2. Generate 2 pairwise CRAM metrics (#3,4) for NA12878 run on the SAME pipeline but from different input sets that have been downsampled to ~24X.  This can be done 4 times: twice with 2 replicates from the same center, and twice with replicates from two different centers.

3. Generate 2 pairwise CRAM metrics (#3,4) for NA19238 run on the SAME pipeline but from different input sets that have been downsampled to ~30X

4. Generate pairwise CRAM metrics for the SAME input set run on different pipelines (10 input sets by 4 metrics). Downsampled data is not included in this comparison, only the original data for all samples and replicates.

5. Use the results from steps 2 and 3 to determine acceptable ranges for metrics from step 4.  The details of this are still TBD.

6. Each interested center calls variants using software of choice (GATK, vt, lumpy)

7. Generate 6 pairwise VCF metrics (#5,6,7,8,9,10) for NA12878 run on the SAME pipeline and SAME variant caller but from different input sets. This is done only with replicates at the same coverage (~24X).

8. Generate 3 pairwise VCF metrics (#5,6,7) for NA19238 run on the SAME pipeline and SAME variant caller but from different input sets. This is done only with replicates at the same coverage (~24X).

9. Generate pairwise VCF metrics for the SAME input set run on different pipelines with the SAME variant caller (10 input sets by #5,6,7 metrics by 2? variant callers). This is done for non-downsampled, original data only.

10. Use results from steps 7 and 8 to determine acceptable ranges for metrics from step 9.

Notes:
* In cases where two pipelines are not functionally equivalent, it will be necessary to determine which one is more “correct” according to our standards, and which one needs to be modified.  Comparison to GiaB should help with this, although in some cases it may be a case of different tradeoffs being made (e.g. sensitivity vs specificity)
* There are several proposed ways to compare qualities, etc.  Just looking at the correlation coefficient may not be sufficient to detect biased differences.  We may also need to look at the distribution of differences (t-test, sign test, KS test)

#Pathway for updates to this standard
Summary: We agreed that pipelines will need to be updated during the project, but that this should be a tightly controlled process given the need to reprocess vast amounts of data each time substantial pipeline modifications occur.

Draft plan:
* Initial pipeline versions should serve for project years 1-2.
* Efficiency updates passing functional equivalence tests are always allowed.
* Propose to start a review process in late 2017: invite proposals for pipeline updates that incorporate new aligners, reference genomes and data processing steps.
* If substantial improvements are achievable, implement new pipelines for project years 3-4.
* There will need to be a decision about how large the potential variant calling improvements should be to warrant pipeline modification and data reprocessing.

##Data resources to generate for mid-project pipeline improvement
Summary: We agreed that additional data resources would be needed for development and testing of future pipelines. We have not spent much time on this topic yet.

Draft plan:
* 32 genomes
* 2 hydatidiform moles (CHM1 and CHM13).
* A “pseudodiploid” mixture of CHM1 and CHM13. Will this mixing happen at the level of data or DNA?
* 10 trios with diverse ancestry from the 1000 Genomes Project, selected to maximize overlap with ongoing genome assembly projects.
* CEPH: NA12878, NA12892, NA12891. Best characterized trio. PacBio assembly for NA12878 available from the WashU Reference Genomes Improvement project.
* Yoruban: NA19240, NA19238, NA19239. PacBio assembly for NA19240 (WashU).
* Puerto Rican trio: HG00733, HG00732, HG00731. PacBio assembly for HG00733 (WashU).
* Han Chinese: HG00514, HG00513, HG00512. PacBio assembly for HG00514 (WashU).
* Columbian: HG01352, HG01351, HG01350. PacBio assembly for HG01352 (WashU)
* Gambian: HG02818, HG02816, HG02817. PacBio assembly for HG02818 (WashU)
* Vietnamese: HG02059, HG02060, HG02061. PacBio assembly for HG02059 (WashU)
* Ashkenazi trio. PacBio assembly for child (Genome in a Bottle Consortium).
* Two more?
* Use mendelian segregation as an unbiased measure in trios.
* Use false heterozygous calls as an unbiased measure in haploid mole genomes.
* Develop variant “truthsets” for fully assembled genomes. Leverage ongoing efforts such as Genome in a Bottle.

#Miscellaneous issues
Notes:
* There is the question about whether to trim adaptor sequences. Broad clips adaptor sequences. Michigan and WashU do not. Does not appear that NYGC and Baylor clip adaptor sequences either. WashU thinks that BWA-MEM is capable of doing adequate clipping during alignment.
* Will our pipelines be designed to accommodate multiple sequencing libraries? Implications for duplicate marking.
