# =========================================================
# Full set of commands to generate annotated VCF file of
# structural variants in NA12878 and load into GEMINI.
#
# This annotated VCF file includes:
#     a) LUMPY structural variant calls
#     b) Genotypes calculated by SVTyper
#     c) Read-depths calculated by CNVnator
#     d) Predicted variant impacts from VEP
#
# Note that some of these commands require software that is
# not installed on your computer.
# =========================================================

# -----------------------------------
# 1) Run LUMPY
time ~/HALL_2014/bin/lumpy \
    -mw 7 \
    -tt 0 \
    -x /mnt/thor_pool1/user_data/cc2qe/code/speedseq/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed  \
    -pe bam_file:NA12878_S1.discordants.bam,histo_file:NA12878_S1.bam.x4.histo,mean:319.558535286,stdev:73.8491946995,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,id:10,min_mapping_threshold:20  \
    -sr bam_file:NA12878_S1.splitters.bam,back_distance:10,min_mapping_threshold:20,weight:1,id:11,min_clip:20 \
    > NA12878.lumpy.out

# -----------------------------------
# 2) Convert LUMPY output to VCF
cat NA12878.lumpy.out \
    | ~/HALL_2014/bin/bedpeToVcf -c ../bin/config.txt -f /shared/genomes/b37/full/human_g1k_v37.fasta \
    > NA12878.lumpy.raw.vcf

# -----------------------------------
# 3) Genotype the VCF
~/HALL_2014/bin/svtyper -v NA12878.raw.vcf \
    -B NA12878_S1.bam \
    -S NA12878_S1.splitters.bam \
    > NA12878.lumpy.vcf

# -----------------------------------
# 4) Run CNVnator
# a) create the root file for the genome
cnvnator -root NA12878.root -genome GRCh37 -tree NA12878_S1.bam -unique

# b) generate histogram.
cnvnator -genome GRCh37 -root NA12878.root -his 100 -d /shared/genomes/b37/full/chroms

# c) calculate statistics
echo "3. Calculate statistics"
cnvnator -root NA12878.root -stat 100

# d) read-depth signal processing (most time consuming)
echo "4. Partition"
cnvnator -root NA12878.root -partition 100

# e) call CNVs
cnvnator -root NA12878.root -call 100 > NA12878.cnvs.txt

# -----------------------------------
# 5) Annotate the LUMPY VCF with read-depth from CNVNator
~/HALL_2014/bin/annotate_rd.py \
    --cnvnator $PATH_TO_CNVNATOR \
    -s NA12878 \
    -w 100 \
    -r NA12878_S1.bam.hist.root \
    -v NA12878.lumpy.vcf \
    > NA12878.lumpy.rd.vcf

# -----------------------------------
# 6) Annotate variant impacts with VEP
cat  NA12878.lumpy.rd.vcf \
    | variant_effect_predictor.pl \
    --fork 32 \
    -o /dev/stdout \
    --force_overwrite \
    --format vcf \
    --offline \
    --no_stats \
    --cache \
    --dir_cache /shared/external_bin/ensembl-tools-release-76/cache \
    --species homo_sapiens \
    --sift b \
    --polyphen b \
    --symbol \
    --numbers \
    --biotype \
    --total_length \
    --vcf \
    --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \
    > NA12878.lumpy.annotated.vcf

# -----------------------------------
# 7) Load into GEMINI
# a) compress and index the VCF file
bgzip NA12878.lumpy.annotated.vcf
tabix -p vcf NA12878.lumpy.annotated.vcf.gz

# b) create a GEMINI database
# (adjust --cores to the number of threads available on your system)
gemini load \
    -v NA12878.lumpy.annotated.vcf.gz \
    -t VEP \
    --cores 8 \
    NA12878.lumpy.annotated.vcf.gz.db \

