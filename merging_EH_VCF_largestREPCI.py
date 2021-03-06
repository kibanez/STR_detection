'''
Created on 04/01/2017

@author: kristina ibanez-garikano

Function that takes all the *vcf.gz generated by running EH-off-target-region, merges them, and computes the AF
for each locus
'''

# !/usr/bin/python

import sys
import os
import re
import copy
import ConfigParser
import optparse
from os import path as osp
import logging
import vcf
from operator import itemgetter

localModulesBase = osp.dirname(osp.realpath(__file__))

modulesRelDirs = ["../modules/"]

for moduleRelDir in modulesRelDirs:
    sys.path.insert(0, osp.join(localModulesBase, moduleRelDir))


class OptionParser(optparse.OptionParser):
    def check_required(self, opt):

        option = self.get_option(opt)

        atrib = getattr(self.values, option.dest)

        if atrib is None:
            #            self.error("%s option not supplied" % option)
            return False
        else:
            return True


def print_tables(hash_table, f_output):
    """
    Function that prints two tsv files: one containing all the EH VCF, already enriched somehow
    and the other table, containing only the STR loci that are spotted by having more repetitions that theoretically
    and in practice have been seen

    :param hash_table: it contains the frequencies for each locus
    :param f_output: file in where the function will write the info in hash_table
    :return:
    """

    l_fields = ['chr', 'start', 'end', 'allele', 'gene', 'ref', 'alt', 'Repeat_Motif',
                'num_samples', 'AF', 'list_samples']

    l_chr = set([item[0] for item in hash_table.keys()])

    chr_specials = []
    if 'X' in l_chr:
        l_chr.remove('X')
        chr_specials.append('X')
    elif 'Y' in l_chr:
        l_chr.remove('Y')
        chr_specials.append('Y')

    l_chr = sorted(l_chr)
    l_chr = [str(i) for i in l_chr]

    for i in chr_specials:
        l_chr.append(i)

    fo = open(f_output, 'w')
    fo.write('\t'.join(l_fields) + '\n')
    for key in sorted(sorted(hash_table.keys(), key=itemgetter(1)), key=lambda x: l_chr.index(x[0])):
        fo.write('\t'.join(map(lambda field: hash_table[key].get(field, '.'), l_fields)) + '\n')
    fo.close()

    return f_output


def merging_vcf(l_vcf, path_vcf, logger):
    """
    Function that receives a list of individual VCF files to merge them all

    The main objective here is to extract the GT frequencies for each STR loci (being trinucleotide regions)
    It returns a VCF file containing the AF for each STR detected in the whole cohort
    chr   pos ref alt frequency

    :param l_vcf: list of VCF files generated by EH
    :param path_vcf: folder in where l_vcf are located
    :param logger:
    :return:
    """

    hash_table = {}

    total_samples = len(l_vcf)

    for vcf_input in l_vcf:

        name_vcf = vcf_input

        vcf_input = os.path.join(path_vcf, vcf_input)

        if not os.path.isfile(vcf_input):
            raise IOError("The VCF file does not exist %s " % vcf_input)

        logger.info("Analysing the STR regions within %s" % vcf_input)

        # Sometimes the VCF generated by EH is empty. We need to check whether the VCF is empty or not
        if os.path.getsize(vcf_input) == 0:
            print vcf_input + '\n'
            continue

        vcf_reader = vcf.Reader(filename=vcf_input)

        for i, r in enumerate(vcf_reader):

            try:

                hash_fields = dict(r.INFO)
                hash_fields.update(dict(zip(r.samples[0].data._fields, r.samples[0].data)))

                # include in the hash_fields those STR with GT different than 0/0 ==>
                # here as EH is not consistent is better to check the ALT value

                # if (hash_fields.get('GT') != '0/0'):
                if r.ALT != [None]:

                    # [<STR23>]   or  [<STR23>,<STR26>]
                    # Check whether alt does have more than 1 ALT allele
                    alt = str(r.ALT[0])
                    alt = re.sub('^<', '', alt)
                    alt = re.sub('>$', '', alt)
                    alt = str(re.sub('^STR', '', alt))

                    if len(r.ALT) > 1:
                        alt2 = str(r.ALT[1])
                        alt2 = re.sub('^<', '', alt2)
                        alt2 = re.sub('>$', '', alt2)
                        alt2 = str(re.sub('^STR', '', alt2))


                    # We only analyse STR markers >=3 length
                    if len(str(hash_fields.get('RU', 0))) >= 3:

                        pos = str(r.INFO['END'] - 1)
                        gene = str(hash_fields.get('REPID', ""))

                        # we retrieve all the info contained in the INFO fields
                        hash_variant = {}

                        hash_variant['Repeat_Motif'] = str(hash_fields.get('RU', 0))
                        hash_variant['Reference_length_bp'] = str(hash_fields.get('RL', 0))
                        hash_variant['gene'] = gene
                        hash_variant['gt'] = str(hash_fields.get('GT', ""))
                        hash_variant['chr'] = r.CHROM
                        hash_variant['start'] = str(r.POS)
                        hash_variant['end'] = str(r.INFO['END'])
                        hash_variant['ref'] = r.REF
                        hash_variant['alt'] = alt

                        total_length_STR = re.sub('^STR', '', alt)
                        total_length_REF = str(hash_fields.get('REF', ""))

                        hash_variant['alt_size'] = str(total_length_STR)
                        hash_variant['ref_size'] = str(total_length_REF)
                        hash_variant['num_samples'] = '1'
                        hash_variant['list_samples'] = name_vcf

                        largest_ci_1 = hash_fields['REPCI'].split('/')[0].split('-')[1]
                        largest_ci_2 = hash_fields['REPCI'].split('/')[1].split('-')[1]



                        # Before: each variant consists of CHROM, POS, with the REF and its ALT alleles,
                        # in which alt is not the alternate allele, but the STR repeats
                        # Now: depending on the genotype (GT) we will check which alleles are frequent in our cohort

                        # We will call `allele` to each allele size detected in each genome
                        # For loci in autosomal chromosomes, we expect to have 2 repeat sizes corresponding to both alleles
                        # For loci in sexual chromosomes, we expect to have 2 repeat size alleles if the genome is female,
                        # or 1 repeat size allele if the genome is male (i.e. FMR1, AR). We take the max value of both
                        # alleles estimation



                        if hash_variant.get('gt') == '0/0':
                            allele = largest_ci_1
                            hash_variant['allele'] = allele

                            if hash_table.has_key((r.CHROM, pos, gene, allele)):
                                # we add info to an existing key - repeat size allele
                                hash_table[(r.CHROM, pos, gene, allele)]['num_samples'] = str(int(hash_table.get((r.CHROM, pos, gene, allele))['num_samples']) + 2)

                                # we add the LP id in which has been detected the expansion (for extract later on this info)
                                hash_table[(r.CHROM, pos, gene, allele)]['list_samples'] = \
                                    hash_table.get((r.CHROM, pos, gene, allele))['list_samples'] + ';' + name_vcf + '_x2'

                            else:
                                # we add the new alternative allele info
                                hash_variant['num_samples'] = '2'

                                # we update the number of the genome + '_x2)
                                hash_variant['list_samples'] = name_vcf + '`_2x'

                                hash_table[(r.CHROM, pos, gene, allele)] = hash_variant



                        elif hash_variant.get('gt') == '1/1':
                            allele = largest_ci_1
                            hash_variant['allele'] = allele

                            if hash_table.has_key((r.CHROM, pos, gene, allele)):
                                # we add info to an existing key - repeat size allele
                                hash_table[(r.CHROM, pos, gene, allele)]['num_samples'] = str(int(hash_table.get((r.CHROM, pos, gene, allele))['num_samples']) + 2)

                                # we add the LP id in which has been detected the expansion (for extract later on this info)
                                hash_table[(r.CHROM, pos, gene, allele)]['list_samples'] = \
                                    hash_table.get((r.CHROM, pos, gene, allele))['list_samples'] + ';' + name_vcf + '_x2'

                            else:
                                # we add the new alternative allele info
                                hash_variant['num_samples'] = '2'
                                # we update the number of the genome + '_x2)
                                hash_variant['list_samples'] = name_vcf + '`_2x'
                                hash_table[(r.CHROM, pos, gene, allele)] = hash_variant


                        elif hash_variant.get('gt') == '0/1':
                            allele_ref = largest_ci_1
                            allele_alt = largest_ci_2

                            if hash_table.has_key((r.CHROM, pos, gene, allele_ref)):
                                # we add info to an existing key - repeat size allele
                                hash_table[(r.CHROM, pos, gene, allele_ref)]['num_samples'] = str(
                                    int(hash_table.get((r.CHROM, pos, gene, allele_ref))['num_samples']) + 1)

                                # we add the LP id in which has been detected the expansion (for extract later on this info)
                                hash_table[(r.CHROM, pos, gene, allele_ref)]['list_samples'] = \
                                    hash_table.get((r.CHROM, pos, gene, allele_ref))['list_samples'] + ';' + name_vcf

                            else:
                                # we specify the allele
                                hash_variant['allele'] = allele_ref

                                # we add the new alternative allele info
                                hash_table[(r.CHROM, pos, gene, allele_ref)] = hash_variant

                            if hash_table.has_key((r.CHROM, pos, gene, allele_alt)):
                                # we add info to an existing key - repeat size allele
                                hash_table[(r.CHROM, pos, gene, allele_alt)]['num_samples'] = str(
                                    int(hash_table.get((r.CHROM, pos, gene, allele_alt))['num_samples']) + 1)

                                # we add the LP id in which has been detected the expansion (for extract later on this info)
                                hash_table[(r.CHROM, pos, gene, allele_alt)]['list_samples'] = \
                                    hash_table.get((r.CHROM, pos, gene, allele_alt))['list_samples'] + ';' + name_vcf

                            else:
                                # we specify the allele
                                # we create a new hash_variant_alt to have different info (!!)
                                hash_variant_alt = copy.deepcopy(hash_variant)
                                hash_variant_alt['allele'] = allele_alt

                                # we add the new alternative allele info
                                hash_table[(r.CHROM, pos, gene, allele_alt)] = hash_variant_alt

                        elif hash_variant.get('gt') == '1/2':
                            allele_alt1 = largest_ci_1
                            allele_alt2 = largest_ci_2

                            if hash_table.has_key((r.CHROM, pos, gene, allele_alt1)):
                                # we add info to an existing key - repeat size allele
                                hash_table[(r.CHROM, pos, gene, allele_alt1)]['num_samples'] = str(
                                    int(hash_table.get((r.CHROM, pos, gene, allele_alt1))['num_samples']) + 1)

                                # we add the LP id in which has been detected the expansion (for extract later on this info)
                                hash_table[(r.CHROM, pos, gene, allele_alt1)]['list_samples'] = \
                                    hash_table.get((r.CHROM, pos, gene, allele_alt1))['list_samples'] + ';' + name_vcf

                            else:
                                # we specify the allele
                                hash_variant['allele'] = allele_alt1

                                # we add the new alternative allele info
                                hash_table[(r.CHROM, pos, gene, allele_alt1)] = hash_variant

                            if hash_table.has_key((r.CHROM, pos, gene, allele_alt2)):
                                # we add info to an existing key - repeat size allele
                                hash_table[(r.CHROM, pos, gene, allele_alt2)]['num_samples'] = str(
                                    int(hash_table.get((r.CHROM, pos, gene, allele_alt2))['num_samples']) + 1)

                                # we add the LP id in which has been detected the expansion (for extract later on this info)
                                hash_table[(r.CHROM, pos, gene, allele_alt2)]['list_samples'] = \
                                    hash_table.get((r.CHROM, pos, gene, allele_alt2))['list_samples'] + ';' + name_vcf

                            else:
                                # we specify the allele
                                hash_variant_alt = copy.deepcopy(hash_variant)
                                hash_variant_alt['allele'] = allele_alt2

                                # we add the new alternative allele info
                                hash_table[(r.CHROM, pos, gene, allele_alt2)] = hash_variant_alt

            except:

                raise RuntimeError(
                    'run_merging_EH_VCF.merging_vcf: Some error has occurred in variant line')

    # Once all the VCF files are ingested, now we compute the statistics: how many samples have the same alternate allele
    # and how many times is seen, against the total number of samples studying (total_samples)

    for key, value in hash_table.iteritems():
        # we calculate the AF for each STR site (for each alternate allele)
        af = float(hash_table.get(key)['num_samples']) / float(total_samples)
        hash_table[key]['AF'] = str(af)

    logger.info("Computing the allele frequencies after having digested all the individual VCF files")

    return hash_table


def read_cfg_file(cfg_filename):
    '''
    Function that reads the configuration file which includes the paths to all the fundamental input

    :param cfg_filename: configuration file
    :return: hash_table containing info in the config file
    '''
    fi = open(cfg_filename, 'r')

    config = ConfigParser.ConfigParser()
    config.readfp(fi)

    hash_cfg = {}

    for field in config.options('GENERAL'):
        hash_cfg[field] = config.get('GENERAL', field)

    for field in config.options('OUTPUT'):
        hash_cfg[field] = config.get('OUTPUT', field)

    for field in config.options('REFERENCE'):
        hash_cfg[field] = config.get('REFERENCE', field)

    for field in config.options('SOFTWARE'):
        hash_cfg[field] = config.get('SOFTWARE', field)

    fi.close()

    return hash_cfg


def run(argv=None):
    if argv is None: argv = sys.argv

    parser = OptionParser(add_help_option=True, description="")
    parser.add_option("--s", default=None, help="The path in which the resulting EH VCF files are", dest="f_samples")
    parser.add_option("--o", default=None, help="The output VCF name in which the merged VCF will be write",
                      dest="f_output")
    parser.add_option("--O", default=None, help="The output directory in which the merged VCF will be write",
                      dest="d_output")
    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:
        sys.exit(0)

    if not parser.check_required("--s"):
        raise IOError('The path to the VCF files generated by running EH is missing')

    if not parser.check_required("--o"):
        raise IOError('The merged VCF file is missing')

    if not parser.check_required("--O"):
        raise IOError('The output directory in which the merged VCF files will be write is missing')

    try:

        if options.f_samples != None:

            path_samples = options.f_samples

            if not os.path.exists(path_samples):
                raise IOError('The path to the VCF files generated by running EH %s does not exist') % path_samples

            # Output folder in which the output of EH will be saved
            output_folder = options.d_output

            if not os.path.exists(output_folder):
                os.mkdir(output_folder)

            merged_vcf = options.f_output

            if merged_vcf is None:
                raise IOError('The name for the merged VCF is missing %s') % merged_vcf

            output_file = os.path.join(output_folder, merged_vcf)


            # Configure logger
            formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
            console = logging.StreamHandler()
            console.setFormatter(formatter)
            console.setLevel(logging.INFO)
            logger = logging.getLogger("preprocess")
            logger.setLevel(logging.INFO)
            logger.addHandler(console)

            logger.info("The process of merging all the individual VCF files has started")

            # Here we retrieve all the VCF files generated by running EH algorithm in a cohort
            # EH for each sample generates a *.json, *.vcf and *.log file
            # we are interested in the VCF files

            l_vcf = [f for f in os.listdir(path_samples) if f.endswith('.vcf')]

            l_samples = []

            for vcf in l_vcf:
                sample_name = re.sub('^EH_', '', vcf)
                sample_name = re.sub('.vcf$', '', sample_name)
                l_samples.append(sample_name)

            hash_table = merging_vcf(l_vcf, path_samples, logger)

            tsv_file = print_tables(hash_table, output_file)

            logger.info(
                "The merged VCF files has been annotated and enriched. Take a look to the TSV file %s") % tsv_file

    except:
        print >> sys.stderr, '\n%s\t%s' % (sys.exc_info()[0], sys.exc_info()[1])
        sys.exit(2)


########################################################################################################################

if __name__ == '__main__':
    run()

