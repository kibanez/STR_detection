from EH import *

# Generating some tests to run EH to put into the IP

l_bams = ['/Users/kibanez/Documents/GEL_STR/testing/LP3000111-DNA_C01.bam',
          '/Users/kibanez/Documents/GEL_STR/testing/LP3000112-DNA_F03.bam',
          '/Users/kibanez/Documents/GEL_STR/testing/LP3000113-DNA_D05.bam']

l_gender = ['female', 'male', 'male']
output_folder = '/Users/kibanez/Documents/GEL_STR/testing'
ref_fasta = '/Users/kibanez//git/fasta/b38/genome.fa'
specs_folder = '/Users/kibanez/Documents/GEL_STR/specs_b38'
eh_binary = '/Users/kibanez/Documents/GEL_STR/EH_releases/ExpansionHunter-v2.5.5-macOS/bin/ExpansionHunter'
eh_command = EH(l_bams, l_gender, 'family_test', output_folder, ref_fasta, specs_folder, eh_binary)

# run EH on them
eh_command.run_EH_offtarget_reads()

# create plots for all loci on them
eh_command.visualize_STRs()

# merge vcf files into a multisample vcf file
#eh_command.generate_multisample_VCF()

