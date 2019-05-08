"""
Created on 24.05.2018

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
    It reads a TSV file containing in 2 columns the platekeys regarding each duplicate pair
    :param input_tsv: tsv file with 2 columns: <lp_id>\t<lp_id>
    :return: hash table containing the 2 columns - associates the LP id with the BAM file
    """

    data = pd.read_csv(input_tsv, sep="\t")

    iid1 = data.iloc[:, 0].tolist()
    iid2 = data.iloc[:, 1].tolist()

    hash_table = dict(zip(iid1, iid2))

    return hash_table



def run(argv=None):

    if argv is None:
        argv = sys.argv

    parser = OptionParser(add_help_option=True, description="")

    parser.add_option("--folder", default=None,
                      help="Folder in where the genomes are", dest="f_folder")

    parser.add_option("--list_pair_duplicates", default=None,
                      help="File with pair duplicates in 2 columns", dest="f_input")

    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:
        sys.exit(0)

    if not parser.check_required("--folder"):
        raise IOError('The folder or directory in which we need to do the check is required\n')

    if not parser.check_required("--list_pair_duplicates"):
        raise IOError('The input TSV file with genomes to investigate splitted into 2 columns is required\n')

    try:

        f_folder = options.f_folder
        f_input = options.f_input

        if not os.path.exists(f_folder):
            raise IOError('The directory %s does not exist') % f_folder

        if not os.path.isfile(f_input):
            raise IOError('The input TSV file %s does not exist') % f_input


        formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
        console = logging.StreamHandler()
        console.setFormatter(formatter)
        console.setLevel(logging.INFO)
        logger = logging.getLogger("preprocess")
        logger.setLevel(logging.INFO)
        logger.addHandler(console)

        logger.info("Keeping only one of the pair duplicated genomes")

        logger.info("1 - Reading the input file containing duplicate pair genomes")

        hash_table = read_tsv_file(f_input)

        logger.info("2 - Keeping only one of the genomes if both do exist")

        for iid1, iid2 in hash_table.iteritems():
            # If both iid1 and iid2 do exist ==> we only keep 1 (by default we keep the first one)
            # Otherwise ==> we do not do anything
            iid1_path = os.path.join(f_folder, 'EH_' + iid1 + '.vcf')
            iid2_path = os.path.join(f_folder, 'EH_' + iid2 + '.vcf')

            if os.path.isfile(iid1_path) and os.path.isfile(iid2_path):
                print("Removing %s") % str(iid2)
                os.remove(iid2_path)


        logger.info('Finished keeping only one copy of duplicate genomes')

    except:

        print >> sys.stderr, '\n%s\t%s' % (sys.exc_info()[0], sys.exc_info()[1])
        sys.exit(2)


if __name__ == '__main__':
    run()
    