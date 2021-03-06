'''
Created on 19.01.2017

@author: kristina ibanez-garikano

This python script analyses the *vcf.gz Illumina delivers to us in GRCh38
These files are generated by running an old version of EH they have, for 2 loci: FMRI and C9orf72
It creates a *tsv file for each locus, with the AF for each repeat in the cohort analysed

'''

# !/usr/bin/python

import sys, os, vcf

import pandas as pd

import optparse

import logging


# Function to read the configuration file
class OptionParser(optparse.OptionParser):

    def check_required(self, opt):

        option = self.get_option(opt)

        atrib = getattr(self.values, option.dest)

        if atrib is None:
            # self.error("%s option not supplied" % option)
            return False
        else:
            return True


# Function that reads the TSV file containing the information about the individuals retrieved from Catalog
# @input: TSV file with 2 columns: <lp_id>\t<path_to_the_structure>
# It returns a hash table containing the 2 columns - associates the LP id with the BAM file
def read_tsv_file(input):

    data = pd.read_csv(input, sep="\t", header=1)

    l_lp = data.iloc[:, 0].tolist()
    l_paths = data.iloc[:, 1].tolist()

    hash_table = dict(zip(l_lp, l_paths))

    return hash_table




# Function that prints two TSV files: one containing all the allelotype-lobSTR VCF, already enriched somehow
# and the other table, containing only the STR loci that are spotted by having more repetitions that theoretically
# and in practice have been seen
def print_tables(hash_table, f_output):

    l_fields = ['chr', 'start', 'ref', 'alt', 'gene', 'filter', 'repeat1', 'repeat2', 'num_samples', 'AF',
                'list_samples']

    # chr_cyto_sorted = map(lambda i: "chr%i" % (i), range(1, 23)) + ['chrX', 'chrY', 'chrM']
    # list of chr to be printed
    l_chr = set([item[0] for item in hash_table.keys()])

    chr_specials = []
    if 'X' in l_chr:
        l_chr.remove('X')
        chr_specials.append('X')
    elif 'Y' in l_chr:
        l_chr.remove('Y')
        chr_specials.append('Y')

    # l_chr = [int(i) for i in l_chr]
    # l_chr = sorted(l_chr)
    # l_chr = [str(i) for i in l_chr]

    for i in chr_specials:
        l_chr.append(i)

    fo = open(f_output, 'w')
    fo.write('\t'.join(l_fields) + '\n')
    for key in sorted(hash_table.keys()):
        fo.write('\t'.join(map(lambda field: hash_table[key].get(field, '.'), l_fields)) + '\n')
    fo.close()

    return f_output

# Function that reads each VCF that were retrieved from the upload_list, rare diseases and V4
# It analyses the INFO fields, looking for pathogenic expansions
def analyse_expansions_vcf(l_vcf, l_samples, logger):
    # January'2017 : Illumina delivers to us the expansions detected for 2 genes: ALS and FMR1

    total_samples = len(l_vcf)

    hash_table_als = {}
    hash_table_fmr1 = {}

    for j, vcf_input in enumerate(l_vcf):

        if not os.path.exists(vcf_input):
            raise IOError("The VCF file generated by EH does not exist %s" % (vcf_input))

        vcf_reader = vcf.Reader(filename=vcf_input)

        for i, r in enumerate(vcf_reader):

            try:
                hash_fields = dict(r.INFO)

                hash_variant = {}

                chr = str(r.CHROM)
                start = str(r.POS)
                ref = str(r.REF)
                gene = str(r.ID)

                if r.FILTER == []:

                    filter_result = 'PASS'

                else:

                    filter_result = str(r.FILTER[0])

                repeat_count1 = hash_fields.get('REPEAT_COUNT1')
                repeat_count2 = hash_fields.get('REPEAT_COUNT2')

                if repeat_count1 is None:
                    hash_variant['repeat1'] = '.'
                    repeat1_int = 0
                else:
                    # we retrieved the max repeat number
                    hash_variant['repeat1'] = '-'.join(repeat_count1)
                    repeat1_int = max(map(int, repeat_count1))

                if repeat_count2 is None:
                    hash_variant['repeat2'] = '.'
                    repeat2_int = 0
                else:
                    # we retrieved the max repeat number
                    hash_variant['repeat2'] = '-'.join(repeat_count2)
                    repeat2_int = max(map(int, repeat_count2))

                # as ALT we consider the maximum alternate alleles from both alleles detected
                alt = str(max(repeat1_int, repeat2_int))

                # list with the information regarding the 2 alleles
                l_alt = []
                if repeat1_int != 0 :
                    l_alt.append(str(repeat1_int))
                if repeat2_int != 0:
                    l_alt.append(str(repeat2_int))


                hash_variant['chr'] = chr
                hash_variant['start'] = start
                hash_variant['ref'] = ref
                hash_variant['alt'] = alt
                hash_variant['filter'] = filter_result
                hash_variant['gene'] = gene
                hash_variant['list_samples'] = vcf_input

                # We need to record and save the number of the STR on each allele - not only the maximum alternate allele
                if gene == 'ALS':

                    for alt in l_alt:

                        if not hash_table_als.has_key((chr, start, ref, alt)):

                            hash_table_als[(chr, start, ref, alt)] = hash_variant
                            hash_table_als[(chr, start, ref, alt)]['num_samples'] = '1'

                        else:

                            hash_table_als[(chr, start, ref, alt)]['num_samples'] = str(
                                int(hash_table_als.get((chr, start, ref, alt))['num_samples']) + 1)
                            hash_table_als[(chr, start, ref, alt)]['list_samples'] = \
                            hash_table_als.get((chr, start, ref, alt))['list_samples'] + ';' + vcf_input

                elif gene == 'FMR1':

                    for alt in l_alt:

                        if not hash_table_fmr1.has_key((chr, start, ref, alt)):

                            hash_table_fmr1[(chr, start, ref, alt)] = hash_variant
                            hash_table_fmr1[(chr, start, ref, alt)]['num_samples'] = '1'


                        else:

                            hash_table_fmr1[(chr, start, ref, alt)]['num_samples'] = str(
                                int(hash_table_fmr1.get((chr, start, ref, alt))['num_samples']) + 1)
                            hash_table_fmr1[(chr, start, ref, alt)]['list_samples'] = \
                            hash_table_fmr1.get((chr, start, ref, alt))['list_samples'] + ';' + vcf_input

            except:

                raise RuntimeError('analyse_expansions_GRCh38_Illumina: Some error has occurred in variant line')

    # Once all the VCF files are ingested, now we compute the statistics: how many samples have the same alternate allele
    # and how many times is seen, against the total number of samples studying (total_samples)

    for key, value in hash_table_als.iteritems():
        # we calculate the AF for each STR site (for each alternate allele)
        af = float(hash_table_als.get(key)['num_samples']) / float(total_samples)
        hash_table_als[key]['AF'] = str(af)

    for key, value in hash_table_fmr1.iteritems():
        # we calculate the AF for each STR site (for each alternate allele)
        af = float(hash_table_fmr1.get(key)['num_samples']) / float(total_samples)
        hash_table_fmr1[key]['AF'] = str(af)

    return hash_table_als, hash_table_fmr1


def run(argv=None):

    if argv is None: argv = sys.argv

    parser = OptionParser(add_help_option=True, description="")

    parser.add_option("--vcf", default=None,
                      help="A 2 column file including the LP-ID and the corresponding paths to each VCF file - ExpansionHunter for ALS and FMR1 delivered by Illumina",
                      dest="f_vcf")
    parser.add_option("--o", default=None,
                      help="Folder or directory in which the TSV files regading each locus will be saved or created",
                      dest="f_output")

    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:
        sys.exit(0)

    if not parser.check_required("--vcf"):
        raise IOError('The text file including the path to the VCF files does not exist')

    try:

        if options.f_vcf <> None:

            vcf_list_file = options.f_vcf

            if not os.path.isfile(vcf_list_file):
                raise IOError('The text file including the path to the VCF files %s does not exist' % (vcf_list_file))

            f_output = options.f_output

            if not os.path.exists(f_output):
                os.mkdir(f_output)

            # Configure logger
            formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
            console = logging.StreamHandler()
            console.setFormatter(formatter)
            console.setLevel(logging.INFO)
            logger = logging.getLogger("preprocess")
            logger.setLevel(logging.INFO)
            logger.addHandler(console)

            logger.info(
                "Analysing the FMR1 and ALS expansions on V4 VCF files generated by running ExpansionHunter - delivered by Illumina")

            logger.info("1 - Reading the information of the individuals to be analysed ...")

            # hash_table contains the LP id for each individual and the path in which is saved (the BAM file)
            hash_table = read_tsv_file(vcf_list_file)

            logger.info("2 - Analysing the expansions on ALS and FMR1 for each sample ...")

            l_vcf = []
            l_samples = []

            for lp_id, path in hash_table.iteritems():
                l_vcf.append(path)
                l_samples.append(lp_id)
            hash_table_als, hash_table_fmr1 = analyse_expansions_vcf(l_vcf, l_samples, logger)

            f_output_als = os.path.join(f_output, "frequency_ALS_GRCh38.tsv")
            print_tables(hash_table_als, f_output_als)

            f_output_fmr1 = os.path.join(f_output, "frequency_FMR1_GRCh38.tsv")
            print_tables(hash_table_fmr1, f_output_fmr1)

            logger.info("Finished analysing the expansions on ALS and FMR1 in V4 - rare diseases")


    except:

        print >> sys.stderr, '\n%s\t%s' % (sys.exc_info()[0], sys.exc_info()[1])
        sys.exit(2)


if __name__ == '__main__':
    run()

