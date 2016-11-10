cd ~/HALL_2014/results

#----------------------------------------------------------------------------------------------------------
# 1) GENERATE BREAKPOINT CALLS IN BED12 FORMAT FOR EASY IGV VISUALIZATION
#----------------------------------------------------------------------------------------------------------

~/HALL_2014/bin/bedpeToBed12 -i breakpoints.naive.bedpe -n breakpoints.naive > ../tracks/ourTracks/breakpoints.naive.bed 
~/HALL_2014/bin/bedpeToBed12 -i breakpoints.strict.bedpe -n breakpoints.strict > ../tracks/ourTracks/breakpoints.strict.bed 

# COLOR KEY: 
# Deletion = DEL = red; 
# Tandem Duplication = DUP = green; 
# Inversion = INV = blue
# Inter-chromosomal = INT = gray (this also includes intra-chromosomal calls > 1 Mb, to enable easy viewing)

# In BED12 format, the predicted breakpoint-containing intervals are shown as solid blocks.
# The thinner line just connects the blocks. 
# Note that when mappings are very far apart (e.g., >1 Mb or different chromosomes),
# we draw the two read mapping clusters as distinct features (to avoid too many lines taking over IGV screen space)
# BED12 was designed to visualize genes, with exons as blocks and introns as the thin line.
# But, since IGV and UCSC do not accommodate paired variant call formats, we've adapted BED12
# Later on we'll visualize these in VCF format as well.

#----------------------------------------------------------------------------------------------------------
# 2) LOAD TRACKS
#----------------------------------------------------------------------------------------------------------

# Start the Integrative Genomics Viewer (IGV)

# FIRST, make sure you choose the "Human (b37)" genome before doing anything else.

# SECOND, go to the settings menu at the top: view > preferences > alignments:
# Change "Visibility range threshold" to 20
# Change "Max read count" to 100
# Change "per window size" to 50

# THIRD, change "insert size options" at the bottom; set min to 0 and max to 687 (this is 5 standard deviations)

# FOURTH, load the bamfile: ~/HALL_2014/results/NA12878.20.bam

# FIFTH, load all the .bed files from the following two directories:
# (the .idx files are index files to make igv faster)
# ~/HALL_2014/tracks/ourTracks
# ~/HALL_2014/tracks/publicTracks

# SIXTH, customize your IGV session.
# I prefer to put all of my custom BED tracks on the top,
# the bam files in the middle, and all the public tracks on the bottom.
# But, you can arrange things however you like.
# You can move tracks up and down by dragging them with the mouse.

# SEVENTH, save the IGV session (under the "file" menu bar)
# for easy reload if (when?) IGV crashes due to memory limitations

#----------------------------------------------------------------------------------------------------------
# 3) SPEND SOME TIME LOOKING AT THE DATA
#----------------------------------------------------------------------------------------------------------

# Do the LUMPY calls make sense based on the raw alignments? 
# When you mouse over the LUMPY call it shows a score; this is the number of alignments that were clustered.

# Try the following features in IGV:

# 1) Click on a BED track to select it, zoom into a ~20kb window, then press "control-f" to flip to the next feature. 
# This is good way to view a lot of calls in a short period of time without IGV freezing up.

# 2) When you find an variant that is between different chromosomes, right click on discordant read
# and select "View mate region in split screen". Right-click near the top to switch back to "standard view".

# Can you find any obvious false positives?
# Can you find any obvious false negatives?

# Do you see any large clusters of LUMPY calls in the "naive" callset that are not present in the "strict" callset?
# Do you any deletion calls that are actually tranposon insertions in the reference genome?
# Do you see any complex variants composed of multiple LUMPY calls that are right next to each other, or overlapping?

#----------------------------------------------------------------------------------------------------------
# 4) LOAD SOME MORE INFORMATIVE TRACKS
#----------------------------------------------------------------------------------------------------------

# ~/HALL_2014/tracks/supplementalTracks contains some additional data including the following:
# a) NA12878.lumpy.bed: these are the LUMPY calls for the entire genome, in BED12 format
# b) NA12878.lumpy.annotated.vcf: these are the LUMPY calls for the entire genome, in VCF format
# c) NA12878.cnvnator.cnvs.bed: these are CNV calls obtained by read-depth analysis using CNVNator. 
# (These are unfiltered so there are likely to be many false positives; e.g., the smallest ones are lowest confidence)
# d) NA12878.cnvnator.readDepth.tdf: this is normalized read-depth data in 100 bp windows. 
# This is very useful for visualizing the raw read-depth data. Your eyes are better than most algorithms.  
    # when loading this last file, you'll need to change the view options:
    # right click on the track name at the very left; change windowing function to "none"; 
    # right click on the track name and change data range min to 0, Mid to 2, and Max to 10;
    # right click on the track name and change Track Height to 100 for better visualization
