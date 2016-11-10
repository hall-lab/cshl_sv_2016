cd ~/HALL_2014/results

#----------------------------------------------------------------------------------------------------------
# 1) EXTRACT DELETION CALLS IN BOTH BEDPE AND BED FORMAT
#----------------------------------------------------------------------------------------------------------

# For this exercise, we'll only use deletion calls; we'll want them in both BEDPE and BED format

cat breakpoints.naive.bedpe | awk '$11=="DEL"' > deletions.naive.bedpe
cat breakpoints.naive.bedpe | awk '$11=="DEL"' | cut -f 1,2,6,7 > deletions.naive.bed

cat breakpoints.strict.bedpe | awk '$11=="DEL"' > deletions.strict.bedpe
cat breakpoints.strict.bedpe | awk '$11=="DEL"' | cut -f 1,2,6,7 > deletions.strict.bed

#----------------------------------------------------------------------------------------------------------
# 2) ASSESS VARIANT DETECTION SENSITIVITY
#----------------------------------------------------------------------------------------------------------
# How many of the published 1000 Genomes (1kg) project deletion calls did we detect?
# (for this excercise we will assume that all 1kg calls are correct)
# Use bedtools pairtopair for this (see slide on bedtools commands)
# FIRST, count how many 1kg calls are on chr20:
cat ../annotations/1kg.na12878.sv.merged.bedpe | awk '$1==20' | wc -l
#      76

# SECOND, count how many 1kg calls were found in our dataset:
# Why do we need to count unique occurences of the ID field?

# the naive callset
cat ../annotations/1kg.na12878.sv.merged.bedpe | awk '$1==20' \
    | bedtools pairtopair -a stdin -b deletions.naive.bedpe -type both -is \
    | cut -f 7 | sort -u | wc -l
#      48
# Sensitivity: 48 / 76 = 0.6316

# the strict callset
cat ../annotations/1kg.na12878.sv.merged.bedpe | awk '$1==20' \
    | bedtools pairtopair -a stdin -b deletions.strict.bedpe -type both -is \
    | cut -f 7 | sort -u | wc -l
#      42
# Sensitivity: 42 / 76 = 0.5526

# If there is time, look at a few 1kg calls that we missed in IGV to help us troubleshoot LUMPY.
# And, to discern whether these might actually be false positives in 1kg

#----------------------------------------------------------------------------------------------------------
# 3) ASSESS VARIANT DETECTION FALSE DISCOVERY RATE (FDR)
#----------------------------------------------------------------------------------------------------------
# We will briefly discuss different methods to validate SVs and to estimate FDR.
# Namely, independent datasets (1000 Genomes), PCR validation, de novo assembly, long-read data, segregation patterns.
# For this excercise, assume any deletion not previously found by 1kg is a false positive.

# FIRST, count how many total calls are there in each deletion set?
wc -l deletions.*.bedpe
#     467 deletions.naive.bedpe
#      82 deletions.strict.bedpe

# SECOND, count how many were not found by 1kg.
bedtools pairtopair -a deletions.naive.bedpe -b ../annotations/1kg.na12878.sv.merged.bedpe -type notboth -is \
    | cut -f 7 | sort -u | wc -l
# 416
# False discovery rate (FDR): 416 / 467 = 0.8908

bedtools pairtopair -a deletions.strict.bedpe -b ../annotations/1kg.na12878.sv.merged.bedpe -type notboth -is \
    | cut -f 7 | sort -u | wc -l
# 40
# False discovery rate (FDR): 40 / 82 = 0.4878

# If there is time, look a randomly chosen few of the false positives from the strict dataset?
# Are these really false positives, or are they variants that were missed by 1kg?

#----------------------------------------------------------------------------------------------------------
# 4) DISCUSS HOW WE CAN FILTER A DATASET TO IMPROVE PRECISION (REDUCE FDR)
#----------------------------------------------------------------------------------------------------------
#       a) more supporting reads
#       b) variant quality score
#       c) exclude regions of genome prone to false positive (high depth, simple repeats; satellite repeats; segmental duplications)
#       d) uniqueness: mapping quality or repeat content
#       f) multiple signals (i.e, require paired-end and split-read)
#       g) independent variant detection tools (intersection, voting schemes)
#       h) for variant types with two detectable breakpoints (i.e., inversions & reciprocal translocations), 
#          we can require both to be detected.

# 5) HOW DO WE TUNE PARAMETERS/FILTERS TO ACHIEVE THE DESIRED TRADEOFF BETWEEN SENSITIVITY AND ACCURACY? 
# see the ROC curve in the practical session slides
