"""
Created on 16.01.2017

@author: kristina ibanez-garikano
"""

# !/usr/bin/python

import sys
import os
import ConfigParser
import optparse
import logging
from subprocess import Popen, PIPE
import pandas as pd



class OptionParser(optparse.OptionParser):
    """
    Parser function to digest arguments from the user
    """

    def check_required(self, opt):
        """

        :param opt: option of the attribute passed by the user
        :return: the value of that attribute
        """

        option = self.get_option(opt)

        atrib = getattr(self.values, option.dest)

        if atrib is None:
            return False
        else:
            return True


def read_tsv_file(input_tsv):

    """
    It reads a TSV file containing genomes retrieved from Catalog
    :param input_tsv: tsv file with 2 columns: <lp_id>\t<path_to_the_structure>
    :return: hash table containing the 2 columns - associates the LP id with the BAM file
    """

    data = pd.read_csv(input_tsv, sep="\t", header=None)

    l_lp = data.iloc[:, 0].tolist()
    l_paths = data.iloc[:, 1].tolist()
    l_gender = data.iloc[:, 2].tolist()

    hash_path = dict(zip(l_lp, l_paths))
    hash_gender = dict(zip(l_lp, l_gender))

    return hash_path, hash_gender


def read_cfg_file(cfg_filename):

    """
    It reads the config file, containing required input data to run EH
    :param cfg_filename: config file with the GENERAL (the tsv), OUTPUT paths, REFERENCE
    files and their paths, as well as the SOFTWARE (EH binary file)
    :return: hash table with all this information
    """

    file_input = open(cfg_filename, 'r')

    config = ConfigParser.ConfigParser()
    config.readfp(file_input)

    hash_cfg = {}

    for field in config.options('GENERAL'):
        hash_cfg[field] = config.get('GENERAL', field)

    for field in config.options('OUTPUT'):
        hash_cfg[field] = config.get('OUTPUT', field)

    for field in config.options('REFERENCE'):
        hash_cfg[field] = config.get('REFERENCE', field)

    for field in config.options('SOFTWARE'):
        hash_cfg[field] = config.get('SOFTWARE', field)

    file_input.close()

    return hash_cfg


def compute_expansion_hunter_offtarget_reads(eh_path, lp_id, gender, bam_path, fasta_file,
                                             specs_path, output_path, logger):
    """
    Calling EH algorithm, with the off-target functionality
    :param eh_path: path to the EH binary file
    :param lp_id: LP id
    :param bam_path: path to the BAM file
    :param fasta_file: path to the human reference file (GRCh37 or GRCh38)
    :param specs_path: path to the folder containing genomic coordinates for each STR locus defined
    :param output_path: path to the output folder in which EH will write the results
    :param logger: logger object
    :return:
    """

    # We need to take the BAM file from the structure - from bam_path/Assembly/<lp_id>.bam
    bam_path = os.path.join(bam_path, 'Assembly')
    bam_file = os.path.join(bam_path, lp_id + '.bam')

    # Definition of the output files: VCF, JSON and the output LOG files
    vcf_output = os.path.join(output_path, 'EH_' + lp_id + '.vcf')
    json_output = os.path.join(output_path, 'EH_' + lp_id + '.json')
    log_output = os.path.join(output_path, 'EH_' + lp_id + '_alignments_relevant_reads.log')

    if os.path.exists(bam_file):

        args = [eh_path, '--bam', bam_file, 'sex', gender, '--ref-fasta', fasta_file, '--repeat-specs',
                specs_path, '--vcf', vcf_output, '--json', json_output, '--log', log_output]

        eh_output = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True, bufsize=1)
        (out_info, log_data) = eh_output.communicate()
        eh_output.wait()

        logger.info('EH output %s \n%s\n' % (lp_id, out_info))
        logger.info('EH log output:\n%s\n' % log_data)

        if 'error' in log_data:
            if log_data.lower().find('error') != -1:
                raise RuntimeError('run_EH.compute_expansion_hunter: Error in running EH:\n%s'
                                   % log_data)
    else:

        raise IOError("The BAM file %s does not exist" % bam_file)


def run(argv=None):

    if argv is None:
        argv = sys.argv

    parser = OptionParser(add_help_option=True, description="")

    parser.add_option("--cfg", default=None,
                      help="Config file with the complete information of the target regions"
                           "and paths of the files needed for the calling", dest="f_cfg")

    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:
        sys.exit(0)

    if not parser.check_required("--cfg"):
        raise IOError('The cfg file does not exist')

    try:

        if options.f_cfg is not None:

            cfg_file = options.f_cfg

            if not os.path.exists(cfg_file):
                raise IOError('The cfg file %s does not exist' % cfg_file)

            hash_cfg = read_cfg_file(cfg_file)

            f_catalog = hash_cfg.get('catalog_file', '')

            if not os.path.isfile(f_catalog):
                raise IOError('The tsv file including the individuals to be analyse '
                              'does not exist. %s' % f_catalog)

            # Output folder in which the output of popSTR will be saved
            output_path = hash_cfg.get('output_path', '')

            if not os.path.exists(output_path):
                os.mkdir(output_path)

            # Path to the EH algorithm
            eh_path = hash_cfg.get('eh_path', '')

            if not os.path.exists(eh_path):
                raise IOError('The executable binary file with the EH algorithm '
                              'does not exist %s' % eh_path)

            # FASTA file
            fasta_file = hash_cfg.get('fasta_file', '')

            if not os.path.isfile(fasta_file):
                raise IOError('The human genome reference fasta file does not exist %s'
                              % fasta_file)

            specs_path = hash_cfg.get('specs_path', '')

            if not os.path.exists(specs_path):
                raise IOError(
                    'The folder including all the JSON files including the genomic '
                    'positions of each STR loci or markers does not exist %s'
                    % specs_path)

            formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
            console = logging.StreamHandler()
            console.setFormatter(formatter)
            console.setLevel(logging.INFO)
            logger = logging.getLogger("preprocess")
            logger.setLevel(logging.INFO)
            logger.addHandler(console)

            logger.info("The detection of STRs by using EH started ...")

            logger.info("1 - Reading the information of the individuals to be analysed ...")

            hash_path, hash_gender = read_tsv_file(f_catalog)

            logger.info("2 - Detecting the STR for each sample ...")

            for lp_id, path in hash_path.iteritems():

                gender = hash_gender.get(lp_id)

                logger.info("Running EH in %s, path is %s and gender is %s" % (lp_id, path, gender))

                # for every chr in each individual we pass it

                if path == '.':
                    continue

                compute_expansion_hunter_offtarget_reads(eh_path, lp_id, gender, path, fasta_file,
                                                         specs_path, output_path, logger)

                logger.info("... finished")

            logger.info('Finished running EH algorithm')

    except:

        print >> sys.stderr, '\n%s\t%s' % (sys.exc_info()[0], sys.exc_info()[1])
        sys.exit(2)


if __name__ == '__main__':
    run()
    