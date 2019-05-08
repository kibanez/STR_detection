"""
Created on 16.01.2017

@author: kristina ibanez-garikano
"""

# !/usr/bin/python

import os
import re
import STR_pileup
import logging
from subprocess import Popen, PIPE

class EH:

    def __init__(self, l_bams, l_gender, familyID, output_folder, ref_fasta, specs_folder, eh_binary):

        self.l_bams = l_bams
        self.l_gender = l_gender
        self.familyID = familyID
        self.output_folder = output_folder
        self.ref_fasta = ref_fasta
        self.specs_folder = specs_folder
        self.eh_binary = eh_binary


    def get_EH_output(self):
        """

        :return: VCF output files after running EH
        """
        l_vcf = []
        for bam_file in self.l_bams:
            lp_id = os.path.basename(bam_file)
            lp_id = re.sub('\.bam$', '', lp_id)

            output_vcf = os.path.join(self.output_folder, 'EH_' + lp_id)
            l_vcf.append(output_vcf)

        return l_vcf


    def run_EH_offtarget_reads(self):
        """

        :return: it runs EH on the list of BAMs (l_bams) generating output VCF, JSON and LOG output files
        l_bams correspond to the BAM files of a complete family to be interpreted
        """

        for bam_file, gender in zip(self.l_bams, self.l_gender):

            lp_id = os.path.basename(bam_file)
            lp_id = re.sub('\.bam$', '', lp_id)

            # Definition of the output files: VCF, JSON and the output LOG files
            vcf_output = os.path.join(self.output_folder, 'EH_' + lp_id + '.vcf')
            json_output = os.path.join(self.output_folder, 'EH_' + lp_id + '.json')
            log_output = os.path.join(self.output_folder, 'EH_' + lp_id + '_alignments_relevant_reads.log')

            if os.path.exists(bam_file):

                args = [self.eh_binary, '--bam', bam_file, '--sex', gender, '--ref-fasta', self.ref_fasta, '--repeat-specs',
                        self.specs_folder, '--vcf', vcf_output, '--json', json_output, '--log', log_output]

                eh_output = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True, bufsize=1)
                (out_info, log_data) = eh_output.communicate()
                eh_output.wait()

                if 'error' in log_data:
                    if log_data.lower().find('error') != -1:
                        raise RuntimeError('EH.run_EH_offtarget_reads: Error in running EH:\n%s'
                                           % log_data)
            else:

                raise IOError("The BAM file %s does not exist" % bam_file)


    def visualize_STRs(self):
        """

        :return: For each genome, it creates a pdf/png for each locus included in the VCF file
        Useful to share with CIPs to manage better false positives
        """

        l_eh_output = self.get_EH_output()
        for output_file in l_eh_output:
            STR_pileup.main(['--json', output_file + '.json', '--read_align', output_file + '_alignments_relevant_reads.log',
                             '--color', '--output_folder', self.output_folder])


    def generate_multisample_VCF(self):
        """

        :return: it generates a multisample VCF file merging EH output files from a family
        """
        l_eh_output = self.get_EH_output()

        # 1 -  bgzip all VCF files
        for file in l_eh_output:
            vcf_file = file + '.vcf'
            try:
                if os.path.exists(vcf_file):
                    args = ['bgzip', vcf_file]

                    bgzip_output = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True, bufsize=1)
                    (out_info, log_data) = bgzip_output.communicate()
                    bgzip_output.wait()

                    if 'error' in log_data:
                        if log_data.lower().find('error') != -1:
                            raise RuntimeError('EH.bgzip: Error running bgzip:\n%s'
                                               % log_data)
            except Exception, e:
                logging.error("The EH vcf file {file} does not exist:{error}".format(file=vcf_file, error=str(e)))

        # 2 - tabix them
        for file in l_eh_output:
            vcfgz_file = file + '.vcf.gz'
            try:
                if os.path.exists(vcfgz_file):
                    args = ['tabix', vcfgz_file]

                    tabix_output = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True, bufsize=1)
                    (out_info, log_data) = tabix_output.communicate()
                    tabix_output.wait()

                    if 'error' in log_data:
                        if log_data.lower().find('error') != -1:
                            raise RuntimeError('EH.tabix: Error running tabix:\n%s'
                                               % log_data)
            except Exception, e:
                logging.error("The EH vcf.gz file {file} does not exist:{error}".format(file=vcfgz_file, error=str(e)))

        # 3 - merge all VCF file into a single multisample VCF file
        l_vcfgz = [s + '.vcf.gz' for s in l_eh_output]

        try:
            l_exist = []
            for vcf in l_vcfgz:
                l_exist.append(os.path.isfile(vcf))

            if all(l_exist):
                output_multisample = os.path.join(self.output_folder, self.familyID + '.vcf.gz')

                args = ['bcftools', 'merge', '-m', 'all', '-Oz', '-o', output_multisample, '--force-samples']

                merge_output = Popen(args + l_vcfgz, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True, bufsize=1)
                (out_info, log_data) = merge_output.communicate()
                merge_output.wait()

                if log_data != '':
                    if log_data.lower().find('error') != -1:
                        raise RuntimeError('EH.bcftools_merge: Error running bcftools merge:\n%s'
                                           % log_data)

            else:
                raise RuntimeError("There are some compressed vcf files that do not exist")

        except Exception, e:
            logging.error("Problems when running bcftools merge:{error}".format(error=str(e)))

        # 4 - Tabix multisample VCF file
        try:
            if os.path.exists(output_multisample):
                args = ['tabix', output_multisample]

                tabix_output = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True, bufsize=1)
                (out_info, log_data) = tabix_output.communicate()
                tabix_output.wait()

                if 'error' in log_data:
                    if log_data.lower().find('error') != -1:
                        raise RuntimeError('EH.tabix: Error running tabix of the multisample vcf file:\n%s'
                                           % log_data)
        except Exception, e:
            logging.error("The multisample EH vcf.gz file {file} does not exist:{error}".format(file=output_multisample, error=str(e)))

