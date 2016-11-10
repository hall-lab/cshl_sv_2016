cd ~/HALL_2014/results

#----------------------------------------------------------------------------------------------------------
# 1) HOW MANY LUMPY DELETION CALLS ARE SUPPORTED BY CNV CALLS MADE BY READ-DEPTH ANALYSIS WITH CNVNATOR?
#----------------------------------------------------------------------------------------------------------
# The CNVNator file is in the supplemental directory: ~/HALL_2014/supplemental/NA12878.cnvnator.cnvs.bed
# For this we will use intersectBed because we want to define LUMPY/CNVNator matches based on "reciprocal overlap"
# You should understand each part of this command
cat ../supplemental/NA12878.cnvnator.cnvs.bed | grep DEL \
    | bedtools intersect -a deletions.naive.bed -b stdin -r -f 0.5 -wo \
    | cut -f 4 | sort -u | wc -l
#   35
# so, 35 of 467 (7.5%) of the naive LUMPY calls have read-depth support 

cat ../supplemental/NA12878.cnvnator.cnvs.bed | grep DEL \
    | bedtools intersect -a deletions.strict.bed -b stdin -r -f 0.5 -wo \
    | cut -f 4 | sort -u | wc -l
# 24 of the 82 (29.3%) of the strict LUMPY calls have read-depth support
# QUESTION: why don't more of the strict LUMPY calls overlap with read-depth CNV calls?

#----------------------------------------------------------------------------------------------------------
# 2) HOW MANY LUMPY DELETION CALLS ARE INSTEAD LIKELY TO BE TRANSPOSON INSERTIONS IN THE REFERENCE GENOME?
#----------------------------------------------------------------------------------------------------------
# One reason is that they are not actually deletions in the "test" genome, but insertions in the reference genome
# Which can be difficult to detect by read-depth analysis because transposons are too repetitive
# Only SINEs and LINEs are active in human, so we'll look only at these
# We're using bedtools intersect for the reciprocal overlap feature
bedtools intersect -a deletions.strict.bed -b ../annotations/SINE.b37.sorted.bed -r -f 0.5 -wo | cut -f 4 | sort -u | wc -l
#   18
bedtools intersect -a deletions.strict.bed -b ../annotations/LINE.b37.sorted.bed -r -f 0.5 -wo | cut -f 4 | sort -u | wc -l
#   3

#----------------------------------------------------------------------------------------------------------
# 3) HOW MANY DELETION CAUSE EXON LOSS?
#----------------------------------------------------------------------------------------------------------
# First, go get the gene annotations file:
# Look at the remote file using curl
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz" | gunzip -c | head
# Reformat and write file, making a bed entry for each exon in the gene file
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz" | gunzip -c \
    | awk '{gsub("^chr","",$3); n=int($9); split($10,start,",");split($11,end,","); for(i=1;i<=n;++i) {print $3,start[i],end[i],$2"."i,$13,$2}}' OFS="\t" \
    | bedtools sort > ~/HALL_2014/annotations/refGene.exons.b37.bed

# Now intersect the LUMPY calls in bedpe format against the exons in bed format (bedtools pairtobed)
# (we can just use the strict calls this time)
bedtools pairtobed -a deletions.strict.bedpe -b ../annotations/refGene.exons.b37.bed -type ospan \
    | cut -f 7 | sort -u | wc -l
#       5

# How do we summarize the results in a more meaningful way?
# Enter "bedtools groupby" (this tool can work wonders........)
bedtools pairtobed -a deletions.strict.bedpe -b ../annotations/refGene.exons.b37.bed -type ospan \
    | sort -k7,7 | bedtools groupby -g 7 -c 21,20,20 -o distinct,count,collapse
# This prints the breakpoint ID, the gene ID, the number of affected exons, and the IDs of the affected exons
# If you want to see the entire intersected line with the full bedpe and bed entries (not just variant id), add "-full" to the command

#----------------------------------------------------------------------------------------------------------
# 4) OTHER CHALLENGING QUESTIONS YOU COULD WORK ON IF THERE IS TIME:
#----------------------------------------------------------------------------------------------------------
# OR, FEEL FREE TO COME UP WITH YOUR OWN

# a) Can you find any putative retrogene insertions? Hint, look for breakpoints joining exons of the same gene
# You should use the entire genome callset for this question to ensure that some are found
# ~/HALL_2014/supplemental/NA12878.lumpy.bedpe
# Don't feel bad if you can't figure this one out; Colby had trouble with it too
# You can always check the completed script to see how we did this and the following commands
# reminder, completed scripts are here: ~/HALL_2014/supplemental/completedScriptsWithAnswers/

pairToBed -type either -a ../supplemental/NA12878.lumpy.bedpe -b ../annotations/refGene.exons.b37.bed \
    | sort -k7,7 \
    | groupBy -g 7 -c 21,22,20,22,20 -o distinct,distinct,distinct,count_distinct,count_distinct \
    | awk '$5==1 && $6>1'
## ID   GENE    TRANSCRIPT  EXONS   NUM_TRANSCRIPTS NUM_EXONS
# 2829  LINC00202-1 NR_026795   NR_026795.3,NR_026795.5 1   2
# 3086  MTCH2   NM_014342   NM_014342.1,NM_014342.2 1   2
# 3087  MTCH2   NM_014342   NM_014342.1,NM_014342.3 1   2
# 3089  MTCH2   NM_014342   NM_014342.11,NM_014342.12,NM_014342.13  1   3
# 3401  TDG NM_003211   NM_003211.8,NM_003211.9 1   2
# 3402  TDG NM_003211   NM_003211.10,NM_003211.8    1   2
# 3403  TDG NM_003211   NM_003211.10,NM_003211.9    1   2
# 3455  FBRSL1  NM_001142641    NM_001142641.1,NM_001142641.2   1   2


# b) Can you find "complex variants" where deletion, duplication and/or inversion calls are clustered or overlapping each other?
# hint, use pairtopair of the variant callset file against itself, using "-rdn" to not report self-hits
bedtools pairtopair -type either -is -rdn -a breakpoints.strict.bedpe -b breakpoints.strict.bedpe 

# c) Did we capture any believable inversion calls where both breakpoints were detected?
# hint 1: this information is contained in column 13 of the breakpoints.strict.bedpe file
# which describes the different read-pair strand combinations that were clustered
# hint 2: searching for the literal '-' requires a preceding backslash as an escape character ('\-')
cat breakpoints.strict.bedpe | grep '++' | grep '\-\-'
# the backslash preceding the '-' is an escape character, used to search for the literal '-'

# d) Are there any believable loss-of-function variants, where a LUMPY call suggests an exon deletion
# that is also supported by CNVNator read-depth analysis? 
# (use the entire genome callset for this question to help chances of finding one)
cat ../supplemental/NA12878.lumpy.bedpe \
    | awk '{ if ($11=="DEL") print $1,$2,$6,$7 }' OFS="\t" \
    | bedtools intersect -wa -r -f 0.5 -a stdin -b ../supplemental/NA12878.cnvnator.cnvs.bed \
    | bedtools intersect -wo -a stdin -b ../annotations/refGene.exons.b37.bed \
    | groupBy -g 1,2,3,4 -c 9 -o distinct

# e) are any of the loss-of-function variants from (d) homozygous?
# hint: begin by retrieving the homozygous variants from ../supplemental/NA12878.lumpy.vcf, then use
# ~/HALL_2014/bin/zjoin to join them against the BEDPE file

# retrieve the IDs of homozygous variants from the VCF
cat ../supplemental/NA12878.lumpy.vcf | grep -v '^#' | awk '{ if ($10~"1/1") print $3 }' > homAlt.ids.txt

# join the homozygous variant against the command from (d).
cat ../supplemental/NA12878.lumpy.bedpe \
    | awk '{ if ($11=="DEL") print $1,$2,$6,$7 }' OFS="\t" \
    | bedtools intersect -wa -r -f 0.5 -a stdin -b ../supplemental/NA12878.cnvnator.cnvs.bed \
    | bedtools intersect -wo -a stdin -b ../annotations/refGene.exons.b37.bed \
    | groupBy -g 1,2,3,4 -c 9 -o distinct \
    | zjoin -1 4 -2 1 -a stdin -b homAlt.ids.txt
