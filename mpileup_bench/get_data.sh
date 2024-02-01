# ----------------------------------------------------------------------
# Human reference
# I have this locally:
# /nfs/srpipe_references/references/Human/GRCh38_full_analysis_set_plus_decoy_hla/all/fasta/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa

# This has 3366 sequences.
# The BAMs appear to be aligned against 2580 sequences.
# The GIAB reference directory appears to be a few hundred max.
# I'm unsure where their reference used was, but if we stick to the
# already aligned files and stick to the primary chromosomes then
# frankly any GRCh38 should be fine for evaluation purposes.
# (NB: unsure if true for chr21 due to changes, but that may be in
# patch form only)

# ----------------------------------------------------------------------
# HG005 for final evaluation only
# chr20 only

# ---- Truth set
echo "Fetching truth set"
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.bed
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

# ---- Data files
reg=chr20:20000000-21000000
echo "Getting Illumina $reg"
samtools view -o illumina_300x.bam --write-index -@8 https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/HG005.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam $reg

echo "Getting PacBio $reg"
samtools view -o pacbio_50x.bam --write-index -@8 https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/NIST_BGIseq_2x150bp_100x/GRCh38/HG005_GRCh38_BGIseq-2x150-100x_NIST_20211126.bam $reg

echo "Getting BGI $reg"
samtools view -o bgi_100x.bam --write-index -@8 https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/PacBio_CCS_15kb_20kb_chemistry2/GRCh38/GIAB_5mC_CpG/HG005.GRCh38.deepvariant.haplotagged.bam $reg

# ----------------------------------------------------------------------
# HG002 for code modification, tuning, tweaking and round-trips.
# chr1 only

# ---- Truth set
echo "Fetching truth set"
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.noinconsistent.bed

# ---- Data files
# Data is same locations as above, but more data available.  Have a browse
# basically and decide what to test with
reg=chr1
echo "Browse https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/"
# eg https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/Element_AVITI_20231018/HG002_GRCh37_Element-StdInsert_2X150_81x_20231018.bam
