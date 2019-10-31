source ~/.bashrc

#directory holding fastq and mapping files
dir=
#base filename, e.g. for SAM1-62_S1_L001_R1_001.fastq base filename = SAM1-62_S1_L001
bf=

mkdir ${dir}/demultiplexed
cd ${dir}/demultiplexed

# The mapping file I started with was produced by MR DNA, but the only columns 
# we really need for demultiplexing are the first two, which must provide
# the sample name and the barcode. The sabre formatting, as
# laid out on their github, wants 3 tab-delimited columns: 1) barcode; 2)
# name for forward read file; 3) name for reverse read file. So to format file
# to look like this with awk:
awk -v OFS="\t" ' NR > 1 {print $2, $1"_R1.fq", $1"_R2.fq"} ' ${dir}/071918EM515F-mapping2.txt > ${dir}/mapping_file_sabre.txt
# Change the file names to be whatever you want. 

# sabre paired end command. From their github: sabre pe takes two paired-end 
# files and a barcode data file as input and outputs the reads demultiplexed 
# into separate paired-end files using the file names from the data file. The 
# barcodes will be stripped from the reads and the quality values of the barcode 
# bases will also be removed. -f specifies the forward reads, -r is the location 
# of the reverse reads. -b is the location of the barcode file. -u will create a 
# file with reads from the forward reads file with unknown barcodes, -w does the
# same for the reverse reads. -m can be passed to specify a maximum number of 
# mismatches in a barcode allowed (default 0).
/raid1/home/micro/klingesj/grace/local_programs/sabre pe -f ${dir}/${bf}_R1_001.fastq -r ${dir}/${bf}_R2_001.fastq -b ${dir}/mapping_file_sabre.txt -u ${dir}/no_bc_match_R1.fq -w ${dir}/no_bc_match_R2.fq
