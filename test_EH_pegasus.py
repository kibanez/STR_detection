from EH import *

# Generating some tests to run EH to put into the IP
l_bams = ['/genomes/by_date/2017-01-25/HX01879334/LP3000113-DNA_D05/Assembly/LP3000113-DNA_D05.bam','/genomes/by_date/2017-02-16/HX01847202/LP3000111-DNA_C01/Assembly/LP3000111-DNA_C01.bam','/genomes/by_date/2017-01-20/HX01869483/LP3000112-DNA_F03/Assembly/LP3000112-DNA_F03.bam']

output_folder = '/genomes/scratch/kgarikano/GEL_STR/testing_EH_IP'
ref_fasta = '/genomes/resources/genomeref/Illumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa'
specs_folder = '/genomes/scratch/kgarikano/GEL_STR/ExpansionHunter-offtarget-regions/specs_b38'
eh_binary = '/genomes/scratch/kgarikano/GEL_STR/sw/latest_EH_release'

eh_command = EH(l_bams, 'family_test', output_folder, ref_fasta, specs_folder, eh_binary)


# run EH on them
#export LC_ALL=C; unset LANGUAGE

eh_command.run_EH_offtarget_reads()

# create plots for all loci on them
eh_command.visualize_STRs()

# merge vcf files into a multisample vcf file
eh_command.generate_multisample_VCF()