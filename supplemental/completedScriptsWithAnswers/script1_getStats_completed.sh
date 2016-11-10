
# RULES:
# All comments, explanations and notes begin with a "#". 
# (You should take notes within this document, following this rule.)
# All commands that you should run from the command line do not begin with a "#";
# It is important to run commands from the correct directory.
# When you need to change directories I try to include a "cd" command.
# Type "pwd" to check which directory you are in.
# Type "ls" or "ls -lh" to see the contents of the current directory.

# IMPORTANT: if you fall behind and need to skip forward, or if you want to check your work,
# we have included a completed version of these scripts – including answers – in the following directory:
# ~/HALL_2014/supplemental/completedScriptsWithAnswers/

#----------------------------------------------------------------------------------------------------------
# 1) LOOK AT THE BAM FILE TO ENSURE IT IS COMPLETE AND HAS BEEN ALIGNED CORRECTLY
#----------------------------------------------------------------------------------------------------------

# navigate to the top-most directory:
cd ~/HALL_2014

# take a look at the directory structure
ls 
ls -lh

# navigate to the directory with our primary data
cd results

# Look at the first 10 lines of the BAM file:
samtools view NA12878.20.bam | head

# Using "less" allows you to scroll down; "-S" allows you to scroll to right rather than wrapping text.
samtools view NA12878.20.bam | less -S

# What is each field in the bamfile? 

# Col   Field   Description
# 1	QNAME	Query (pair) NAME
# 2	FLAG	bitwise FLAG
# 3	RNAME	Reference sequence NAME
# 4	POS	1-based leftmost POSition/coordinate of clipped sequence
# 5	MAPQ	MAPping Quality (Phred-scaled)
# 6	CIAGR	extended CIGAR string
# 7	MRNM	Mate Reference sequence NaMe (‘=’ if same as RNAME)
# 8	MPOS	1-based Mate POSistion
# 9	ISIZE	Inferred insert SIZE
# 10	SEQ	query SEQuence on the same strand as the reference
# 11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
# 12	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE

# What is the read-length? Simply the length of the sequence in column 10
samtools view NA12878.20.bam | awk '{print length($10)}' | head

# Run "samtools flagstat" to check whether alignment was done properly.
# What information do these flags encode?
# See the first slide of the practical presentation: ~/HALL_2014/slides/HALL_CSHL_SEQ_2014_PRACTICAL.pdf
# See http://samtools.github.io/hts-specs/SAMv1.pdf for formal format specification.
# This website allows you to look up a flag as well as combinations of flags: http://broadinstitute.github.io/picard/explain-flags.html
# With the new version of samtools you can also just type "samtools flags" in the command line 
# The 0x just means that the following string is in hexidecimal
samtools flags

# run samtools flagstat and write the results to a file
samtools flagstat NA12878.20.bam > NA12878.20.bam.flagstat

cat NA12878.20.bam.flagstat
#31776951 + 0 in total (QC-passed reads + QC-failed reads)
#192207 + 0 secondary
#0 + 0 supplimentary
#178569 + 0 duplicates
#31689975 + 0 mapped (99.73%:nan%)
#31584744 + 0 paired in sequencing
#15792372 + 0 read1
#15792372 + 0 read2
#30709224 + 0 properly paired (97.23%:nan%)
#31410792 + 0 with itself and mate mapped
#86976 + 0 singletons (0.28%:nan%)
#502518 + 0 with mate mapped to a different chr
#216020 + 0 with mate mapped to a different chr (mapQ>=5)

# QUESTION: what does each line mean? 
# QUESTION: why are there are two numbers on each line (the latter of which is zero)?
# QUESTION: do these results look OK?

#----------------------------------------------------------------------------------------------------------
# 2) EXAMINE THE INSERT SIZE DISTRIBUTION OF THE SEQUENCING LIBRARY
#----------------------------------------------------------------------------------------------------------
# QUESTION: what does each line of the following command do?
# ANSWER: readpair is "concordant" (or "properly paired"); query is mapped; mate is mapped; fragment length > 0;
# cut out column 9; go 1 million lines into the file; extract 100 thousand lines.
# QUESTION: what does the "&" at the end do?
# ANSWER: the "&" above puts the command in the background, letting you use the terminal while it is running.
# type "fg" and hit return to bring it back to the foreground; type "jobs" to get a list of running jobs in the background
samtools view -u -f 0x0002 NA12878.20.bam \
    | samtools view -u -F 0x0400 - \
    |  samtools view -F 0x0100 - \
    | awk '$9>0' | cut -f 9 \
    | head -n 1000000 | tail -n 100000 > concordant.size &

# check the file to make sure it makes sense; it should be composed of positive numbers between 0 and 1000
cat concordant.size | head

# open a new terminal window (command + N) and open the R program by typing "r" in the terminal

# plot the insert size distribution using a histogram
setwd("~/HALL_2014/results")
concordant = scan("concordant.size")
hist(concordant, breaks = 1:1000)

# zoom in using xlim and ylim
# Use "dev.new()" to generate a new plot window
dev.new(); hist(concordant, breaks = 1:1000,xlim = c(0,700), ylim = c(0,100))

# characterize the insert size distribution by running these commands in the R window
mean(concordant)
sd(concordant)
median(concordant)
mad(concordant)
min(concordant)
max(concordant)

# > mean(concordant)
# [1] 319.9139
# > sd(concordant)
# [1] 73.49944
# > median(concordant)
# [1] 321
# > mad(concordant)
# [1] 65.2344
# > min(concordant)
# [1] 19
# > max(concordant)
# [1] 866

# QUESTION: how do we define "concordant" (or "properly paired") vs. "discordant" alignments
# ANSWER: usually the mean insert size +/- 4 or 5 standard deviations. We'll use 5 standard deviations.
# Now we'll calculate the minimum and maximum fragment length for concordant mappings.
# 319.9139 + (5 * 73.49944) = 687.4111; thus, 687 is the maximum insert size of a concordant readpair
# 319.9139 - (5 * 73.49944) = -47.5833; 0 is the minimum size of a concordant readpair.
# Or, we can define discordants as the median +/- 5-10 median absolute deviations. This is less sensitive to outliers.

# QUESTION: why is there one "concordant" fragment at size 866 when this is clearly outside the normal distribution?

# If you don't trust the aligner, this is how you would generate an unbiased insert size distribution:
# QUESTION: what does each line of this command do?
# ANSWER: query not unmapped; mate not unmapped; not a secondary alignment; query aligns to plus strand; 
# mate aligns to minus strand; fragment length between 0 and 1000; cut the length field out; 
# go 1 million lines into file; extract 100 thousand lines, write to file. 
samtools view -u -F 0x0004 NA12878.20.bam \
    | samtools view -u -F 0x0008 - \
    | samtools view -u -F 0x0100 - \
    | samtools view -u -F 0x0010 - \
    | samtools view -f 0x0020 - \
    | awk '$9>0 && $9<1000' | cut -f 9 \
    | head -n 1000000 | tail -n 100000 > unbiased.size &

# Navigate to the terminal window running R, or open a new window (command + N) and open the R program
# by typing "r" in the terminal

setwd("~/HALL_2014/results") 
unbiased = scan("unbiased.size")
dev.new()
hist(unbiased, breaks = 1:1000)

# QUESTION: is there a significant different between the two histogram plots?

#----------------------------------------------------------------------------------------------------------
# 3) ESTIMATE THE COVERAGE OF OUR DATASET
#----------------------------------------------------------------------------------------------------------

# "Sequence coverage" is the number of times each base in the genome is sequenced, on average
# "Physical coverage" is the mean number of times each base is contained in a fragment from the sequenced library
# The haploid human reference genome is 3.1 billion basepairs; chromosome 20 is 63025520 bp (63 Mb)

# First calculate the sequence coverage. To do this properly is actually a bit complicated. To do it approximately is easy.
# Roughly, it is this is number of mapped reads that are primary alignments, and are not duplicates, 
# multiplied by readlength, and divided by genome size (chr20 size)
(31689975 - 178569 - 192207) * 101 / 63025520
# 50.1898X coverage

# Now calculate physical coverage:
# How many read-pairs are there that are not duplicates? 
# We can't tell this from the samtools flagstat results so we'll count ourselves
# QUESTION: what does each line of this command do?
# ANSWER: not secondary alignment, query not unmapped, mate not unmapped, not a duplicate, first read in pair
samtools view -u -F 0x100 NA12878.20.bam \
    | samtools view -u -F 0x4 - \
    | samtools view -u -F 0x8 - \
    | samtools view -u -F 0x400 - \
    | samtools view -F 0x40 - \
    | wc -l
# 15622923 read-pairs of mean insert size 319.9139 bp

# Now calculate physical coverage for the haploid genome (chr20)
# This is the number of readpairs multiplied by the mean insert size, divided by genome size (chr20 size)
(15622923*319.9139)/63025520
# = 79.3011X coverage
# Physical coverage for the diploid genome is half that

# However, this calculation uses the "outer span" and does not account for the aligned reads, which are necessary to detect an SV breakpoint
# How many times is each base, on average, spanned by an aligned read-pair?  
# We refer this as "spanning coverage" to distinguish it from physical coverage
(15622923*(319.9139-202))/63025520
# 29.2288  for haploid genome, half that for diploid

# QUESTION: how many read-pairs will detect the average heterozygous breakpoint? 
# Since aligners can align small (~25-50) portions of reads to the reference genome, 
# It is somehwere between physical and spanning coverage. 
