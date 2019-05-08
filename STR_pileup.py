#!/usr/bin/env python 
#
# STR_pileup
# Copyright(c) 2018 Illumina, Inc.
# Code adapted by Kristina Ibanez-Garikano from
#
# Author: Viraj Deshpande < vdeshpande@illumina.com >
# Concept: Egor Dolzhenko < edolzhenko@illumina.com >, Michael Eberle < meberle@illumina.com >
#
# This program is free software: you can redistribute it and / or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see < http: // www.gnu.org / licenses / >.
#

import argparse
import json
import math
import os
import sys
from time import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from functools import total_ordering


TSTART = time()

# Repeat units for each repeat_id
repeat_units = {}

# Basepair properties
basepair_color = {'A': 'g', 'C': 'b', 'G': 'orange', 'T': 'r', 'N': 'k', 'R': 'k'}
complementary_nucleotide = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                            'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                            'N': 'N', 'n': 'n'}
def reverse_complement(seq):
    """Returns reverse complemented sequence of input string
    """
    return ''.join([complementary_nucleotide[c] for c in seq[::-1]])

@total_ordering
class read_alignment():
    """Alignment of a single read in the read_align file.
    Attributes:
        read_name:  Read name.
        read_seq:   Read sequence as aligned to reference +ve strand.
            High/low quality bases in upper/lower case.
        ref_seq:    Reference sequence in read alignment. Repeat marked by R.
            Optional: Default = ''.
        alignment:  Alignment string between read_seq and ref_seq.
            Matching positions marked using '|' else ' '.
        repeat_id:  Repeat ID of STR location.
            Optional: Default = ''.
        repeat_unit: Repeat unit of STR.
            Optional: Default = ''.
        variant:    Alignment type as reported in the read_align file.
            e.g. "SPANNING_23", "FLANKING", "FLANKING_35", "INREPEAT_50".
        repeat_length:  Length of STR repeat alignment in basepairs.
        left_flank_size:    Length of left flank of the alignment.
        right_flank_size:   Length of right flank of the alignment.
        gtype:  Genotype inferred for SPANNING reads only. Default = -1.
        offset: Offset of the first STR base-pair of read within STR sequence.
            Maybe be non-zero for INREPEAT and right FLANKING reads.
    Ordered by:
        read_seq.upper()
    """
    def __init__(self, variant, read_seq, alignment, ref_seq, repeat_id='', read_name=''):
        """Initializes read alignment class object.
        Args:
            variant
            read_seq
            alignment
            ref_seq
            repeat_id (optional)
            read_name (optional)
        """
        self.variant = variant
        self.read_seq = read_seq
        self.alignment = alignment
        self.ref_seq = ref_seq
        self.repeat_id = repeat_id
        self.read_name = read_name
        self.repeat_unit = ''
        self.gtype = -1
        self.left_flank_size = 0
        self.right_flank_size = 0
        if repeat_id != '' and self.repeat_id in repeat_units:
            self.repeat_unit = repeat_units[self.repeat_id].upper()
        if 'FLANKING' in variant:
            if ref_seq[-1] == 'R':
                self.left_flank_size = self.ref_seq.find('R')
            else:
                self.right_flank_size = len(
                    self.ref_seq) - self.ref_seq.rfind('R') - 1
            self.gtype = -1
        elif 'SPANNING' in self.variant:
            self.gtype = int(self.variant.split('_')[1])
            self.left_flank_size = ref_seq.find('R')
            self.right_flank_size = len(ref_seq) - ref_seq.rfind('R') - 1
        else:
            self.gtype = -1
            self.ref_seq = 'R' * len(self.read_seq)
        self.offset = (len(self.repeat_unit) - (self.read_seq.upper().find(self.repeat_unit) -
                                                self.left_flank_size)) % len(self.repeat_unit)
        self.repeat_length = self.ref_seq.count('R')

    def __eq__(self, y):
        """Equality of two read_align objects based on equality of read_seq attribute
        """
        self.read_seq.upper() == y.read_seq.upper()

    def __lt__(self, y):
        """Less than comparison of two read_align objects bases on comparison of read_seq attribute
        """
        self.read_seq.upper() < y.read_seq.upper()

def parse_read_align(file_name, repeat_id_list=[]):
    """Parses YAML-formatted read alignment file (log file) into read_align objects.
    Args:
        file_name: Name of log file (include path)
        repeat_id_list (optional): List of repeat_ids to be parsed (Default=all)
    Returns:
        A dict mapping repeat IDs to corresponding lists of read alignment objects.
    """
    ll =[l.replace('  ', '\t', 4).split('\t') for l in open(file_name) if l.strip() != '']
    read_aligns = {}
    read_align_iter = iter(ll)
    repeat_id = ''
    while True:
        try:
            l = next(read_align_iter)
        except:
            break
        if l[0] != '':
            if repeat_id != '' and (repeat_id_list == [] or repeat_id in repeat_id_list):
                read_aligns[repeat_id] = site_reads
            repeat_id = l[0].strip().strip(':')
            site_reads = {}
            continue
        elif repeat_id_list != [] and repeat_id not in repeat_id_list:
            continue
        elif l[1] != '':
            variant = l[1].strip().strip(':')
            site_reads[variant] = []
            continue
        elif l[2] != '':
            read_seq = ''
        elif l[2] == '' and l[3] != '' and l[3].split()[0] == 'name:':
            read_name = l[3].split()[1]
        elif l[2] == '' and l[3] != '' and l[3].split()[0] == 'irr:':
            read_seq = l[3].split()[1]
            if read_seq.count(repeat_units[repeat_id]) < reverse_complement(read_seq).count(repeat_units[repeat_id]):
                read_seq = reverse_complement(read_seq)
            site_reads[variant].append(read_alignment(variant, read_seq, '','R'*len(read_seq), repeat_id, read_name=''))
        elif l[2] == '' and l[3] == '' and l[4] != '':
            if read_seq == '':
                read_seq = l[4].strip()
                read_align = ''
            elif read_align == '':
                read_align = l[4].strip()
            else:
                ref_seq = l[4].strip()
                site_reads[variant].append(read_alignment(variant, read_seq, read_align, ref_seq, repeat_id, read_name=''))
    if repeat_id_list == [] or repeat_id in repeat_id_list:
        read_aligns[repeat_id] = site_reads
    return read_aligns

def get_flanks(repeat_id, read_align_dict):
    """Determines partition sizes of flanks and repeat unit for reads mapping to an STR from multiple samples.

    Partitions reads into 3 sections: Left, Repeat and Right and returns the largest size for each partition.

    Args:
        repeat_id: Repeat ID for which to determine flank sizes
        read_align_dict: dict mapping a sample to dict of repeat_ids mapping to list of read_align objects
            e.g. D[sample][repeat_id] = [align1, align2, align3, ..]
    Returns:
        3-tuple: (Largest left flank, Largest Repeat String (rounded up by 1 unit), Largest right flank)
    """
    ref_units = []
    for sample in read_align_dict:
        sample_read_aligns = read_align_dict[sample]
        if repeat_id not in sample_read_aligns:
            continue
        left_flank_size = max([ra.left_flank_size for v in sample_read_aligns[repeat_id] for ra in sample_read_aligns[repeat_id][v]] + [0])
        right_flank_size = max([ra.right_flank_size for v in sample_read_aligns[repeat_id]
                                for ra in sample_read_aligns[repeat_id][v]] + [0])
        if left_flank_size > 0:
            left_flank_seq = [ra for v in sample_read_aligns[repeat_id] for ra in sample_read_aligns[repeat_id][v] if ra.left_flank_size == left_flank_size][0].ref_seq[:left_flank_size]
        else:
            left_flank_seq = ''
        max_right_flank_ra_list = [ra for v in sample_read_aligns[repeat_id] for ra in sample_read_aligns[repeat_id][v]
                                   if ra.right_flank_size == right_flank_size]
        if len(max_right_flank_ra_list) > 0:
            max_right_flank_seq = max_right_flank_ra_list[0].ref_seq
            right_flank_seq = max_right_flank_seq[len(max_right_flank_seq) - right_flank_size:]
        else:
            right_flank_seq = ''

        max_repeat_seq = max([ra.repeat_length + ra.offset for v in sample_read_aligns[repeat_id]
                              for ra in sample_read_aligns[repeat_id][v]] + [2])
        rep_len = int(len(repeat_units[repeat_id]) * math.ceil(float(max_repeat_seq) / len(repeat_units[repeat_id])))
        repeat_seq = 'R' * rep_len
        ref_units.append((left_flank_seq, right_flank_seq, repeat_seq))
    if len(ref_units) == 0:
        return ('', '', '')
    max_ref_units = [max([(len(u[i]), u[i]) for u in ref_units]) for i in range(3)]
    return (max_ref_units[0][1], max_ref_units[1][1], max_ref_units[2][1])

def get_ref_units(repeat_id, read_align_dict):
    """Returns units of reference sequence divided into left, repeat and right flanks

    For all read_aligns, determines the largest left flank, repeat string, and right flank.
    Generates ref_units corresponding to reference sequences of each partition.

    Args:
        repeat_id: Repeat ID for which to compute ref_units.
        read_align_dict: dict mapping sample_name to dict of repeat_ids mapping to list of read_align objects
            e.g. D[sample][repeat_id] = [align1, align2, align3, ..]
    
    Returns:
        List of 3 ref_units corresponding to for left, repeat and right sequences of reference.
        ref_unit: 2-tuple: (index of first base-pair of ref_unit, sequence of ref_unit).
        e.g. 'GGGRRRRTTT' will be represented as [(0, 'GGG'), (3, 'RRRR'), (7, 'TTT')]
    """
    site_flanks = get_flanks(repeat_id, read_align_dict)
    # x position of each partition in the read pileup
    ref_units = []
    ref_units.append((0, site_flanks[0]))
    ref_units.append((len(site_flanks[0]) + 2, site_flanks[2]))
    ref_units.append(
        (len(site_flanks[0]) + len(site_flanks[2]) + 4, site_flanks[1]))
    return ref_units

def pileup_coordinates(repeat_id, ref_units, read_align_list):
    """Determines the leftmost coordinate for a read mapping to an STR for given sample

    read_alignments are split by genotype if spanning else left, right, flanking categories
    read_alignments within each category are sorted by leftmost mapping coordinate

    Args:
        repeat_id: Repeat ID
        ref_units: Reference partitions to be used (generated using get_ref_units())
        read_align_list: List of read_align objects.
    
    Returns:
        dict mapping from genotype category to list of 2-tuples of reads.
        Each 2-tuple consists of (leftmost coordinate, read_align)
        Genotype categories are:
            int repeat unit for spanning reads
            (-1, 'irr'), (-1, 'left'), (-1, 'right') for in-repeat, left flanking and right flanking reads
        e.g.
        {
            1: [(4, read_align1), (5, read_align2)],
            (-1, 'left'): [(2, read_align3)]
            (-1, 'right'): [(20, read_align4)]
            (-1, 'irr'): [(9, read_align5), (10, read_align6)]
        }
    """
    left_flank_size = len(ref_units[0][1])
    gtypes = {
        ra.gtype for v in read_align_list for ra in read_align_list[v]}
    coords = {}
    # print gtypes
    for g in gtypes:
        if g == -1:
            continue
        read_aligns = [ra for v in read_align_list
                       for ra in read_align_list[v] if ra.gtype == g]
        read_aligns.sort(key=lambda x: x.left_flank_size - x.right_flank_size)
        g_positions = []
        for ra in read_aligns:
            if ra.left_flank_size > 0:
                position = left_flank_size - ra.left_flank_size
            elif ra.right_flank_size > 0:
                position = left_flank_size + g * \
                    len(ra.repeat_unit) + ra.right_flank_size - len(ra.ref_seq)
            g_positions.append((position, ra))
        coords[g] = g_positions
    if -1 in gtypes:
        read_aligns = [ra for v in read_align_list
                       for ra in read_align_list[v] if ra.gtype == -1]
        read_aligns.sort(key=lambda x: x.left_flank_size - x.right_flank_size)
        left_reads = [ra for ra in read_aligns if ra.gtype == -1 and ra.left_flank_size != 0]
        right_reads = [ra for ra in read_aligns if ra.gtype == -1 and ra.right_flank_size != 0]
        ir_reads = [ra for ra in read_aligns
                    if ra.gtype == -1 and ra.right_flank_size == 0 and ra.left_flank_size == 0]
        coords[(-1, 'left')] = [(left_flank_size - ra.left_flank_size, ra)
                                for ra in left_reads]
        g_positions = []
        max_right = max([0] + [len(ra.read_seq) - ra.right_flank_size for ra in right_reads])
        for ra in right_reads:
            position = max_right - (len(ra.read_seq) - ra.right_flank_size)
            g_positions.append((position, ra))
        coords[(-1, 'right')] = g_positions
        #TODO: decided by largest unassigned in repeat unit + possible offset for largest repeat unit
        coords[(-1, 'irr')] = [(left_flank_size + ra.offset, ra)
                               for ra in ir_reads]
    for g in coords:
        coords[g].sort()
    return coords

def get_pileup(repeat_id, read_align_list, ref_units=None):
    """ Get genotype_pileup_list for given sample

    Args:
        repeat_id: Repeat ID for which to generate pileup
        read_align_list: List of read_alignments to repeat_id (see parse_read_align)
        ref_units: (Optional) Use this option for multiple samples. See ref_units.

    Returns:
        (ref_units, genotype_pileup_list)
        where:
        ref_units: See get_ref_units
        genotype_pileup_list: list of 2-tuples (genotype, pileup_positions)
            genotype: 2 tuple: ((gtype, 0 if SPANNING else 1))
                gtype: largest exact genotype upto the size of the repeat length.
                    e.g. if the repeat unit has length 3 and spanning reads have gtypes 2 and 5            
                    then for a flanking with repeat sequence of length 11, gtype = 2
            pileup_positions: list of 3-tuples for 3-partitions each read. 
                partition corresponds to left, repeat and right.
                [(xposition, sequence), ..] for each unit within read
                sort order of reads:
                    first by read type: (left, spanning, right, in-repeat)
                    next by left-most mapping coordinates
                includes a gap of 2 between consecutive partions
    """
    if ref_units is None:
        ref_units = get_ref_units(repeat_id, {None:{repeat_id:read_align_list}})

    coords = pileup_coordinates(repeat_id, ref_units, read_align_list)

    gtype_pileup_list = []
    gprev = 0
    gtype_list = [g for g in coords if type(g) == int]
    gtype_list.sort()
    for g in gtype_list:
        gtype_pileup = []
        for gg in [(-1, 'left'), (-1, 'right'), (-1, 'irr')]:
            if gg not in coords:
                continue
            for c in coords[gg]:
                if c[1].repeat_length >= g * len(repeat_units[repeat_id]) or c[1].repeat_length < gprev * len(repeat_units[repeat_id]):
                    continue
                read_units = []
                if c[1].left_flank_size == 0 and c[1].right_flank_size > 0:
                    clen = ref_units[2][0] - 2 - c[1].repeat_length
                    read_units.append((0, ''))
                    read_units.append(
                        (clen, c[1].read_seq[:c[1].repeat_length]))
                elif c[1].left_flank_size > 0:
                    read_units.append(
                        (c[0], c[1].read_seq[:c[1].left_flank_size]))
                    read_units.append(
                        (ref_units[1][0], c[1].read_seq[c[1].left_flank_size: c[1].left_flank_size + c[1].repeat_length]))
                if c[1].right_flank_size > 0:
                    read_units.append(
                        (ref_units[2][0], c[1].read_seq[-1 * c[1].right_flank_size:]))
                else:
                    read_units.append((ref_units[2][0], ''))
                gtype_pileup.append(read_units)
        if len(gtype_pileup) > 0:
            gtype_pileup_list.append(((gprev, 1), gtype_pileup))
        gtype_pileup = []
        for c in coords[g]:
            read_units = []
            if c[1].left_flank_size == 0 and c[1].right_flank_size > 0:
                clen = ref_units[2][0] - 2 - c[1].repeat_length
                read_units.append((0, ''))
                read_units.append(
                    (clen, c[1].read_seq[:c[1].repeat_length]))
            else:
                read_units.append(
                    (c[0], c[1].read_seq[:c[1].left_flank_size]))
                read_units.append(
                    (ref_units[1][0], c[1].read_seq[c[1].left_flank_size: c[1].left_flank_size + c[1].repeat_length]))
            if c[1].right_flank_size > 0:
                read_units.append(
                    (ref_units[2][0], c[1].read_seq[-1 * c[1].right_flank_size:]))
            else:
                read_units.append((ref_units[2][0], ''))
            gtype_pileup.append(read_units)
        gprev = g
        gtype_pileup_list.append(((g, 0), gtype_pileup))

    gtype_pileup = []
    for gg in [(-1, 'left'), (-1, 'right'), (-1, 'irr')]:
        if gg not in coords:
            continue
        for c in coords[gg]:
            if c[1].repeat_length < gprev * len(repeat_units[repeat_id]):
                continue
            read_units = []
            if c[1].left_flank_size == 0 and c[1].right_flank_size > 0:
                clen = ref_units[2][0] - 2 - c[1].repeat_length
                read_units.append((0, ''))
                read_units.append(
                    (clen, c[1].read_seq[:c[1].repeat_length]))
            elif c[1].left_flank_size > 0:
                read_units.append(
                    (c[0], c[1].read_seq[:c[1].left_flank_size]))
                read_units.append(
                    (ref_units[1][0], c[1].read_seq[c[1].left_flank_size: c[1].left_flank_size + c[1].repeat_length]))
            else:
                read_units.append((0, ''))
                read_units.append(
                    (ref_units[1][0] + c[1].offset, c[1].read_seq))
            if c[1].right_flank_size > 0:
                read_units.append(
                    (ref_units[2][0], c[1].read_seq[-1 * c[1].right_flank_size:]))
            else:
                read_units.append((ref_units[2][0], ''))
            gtype_pileup.append(read_units)
    if len(gtype_pileup) > 0:
        gtype_pileup_list.append(((gprev, 1), gtype_pileup))

    return (ref_units, gtype_pileup_list)

def get_pileup_list(repeat_id, all_sample_read_aligns, sample_names):
    """ Generate genotype_pileups for all sample

    Args:
        repeat_id: Repeat ID for which to generate pileup
        all_sample_read_aligns: dict: {sample_name : {repeat_id: [list of read aligns]}}
        sample_names: ordered list of sample names in pileup

    
    Returns:
        List of 3-tuples with pileups of each genotype within each sample
            e.g. [(sample_name, ref_units, genotype_pileup_list)]
        For ref_units: see get_ref_units()
        For genotype_pileup_list: see get_pileup()
    """
    ref_units = get_ref_units(repeat_id, all_sample_read_aligns)

    pileup_list = []
    for sample_name in sample_names:
        if sample_name in all_sample_read_aligns and repeat_id in all_sample_read_aligns[sample_name]:
            genotype_pileup_list = get_pileup(repeat_id, all_sample_read_aligns[sample_name][repeat_id], ref_units)
        else:
            genotype_pileup_list = get_pileup(repeat_id, [], ref_units)
        pileup_list.append((sample_name, ref_units, genotype_pileup_list[1]))
    return pileup_list

def plot_pileup(repeat_id, pileup_list, genotypes, color=False, output_prefix = ''):
    """Plot read pileups for all samples for give repeat_id

    Args:
        repeat_id: Repeat ID for which to plot pileup
        pileup_list: See get_pile_list()
        genotypes: Genotypes reported by EH for each sample
        color: (optional): Set True for IGV scheme. Default: black+red(mismatches).
        output_prefix: (optional) Default ''

    Outputs:
        PNG and PDF files: output_prefix + "_" + repeat_id
    """
    xscale = 1.0 / 10
    yscale = 2.0 / 10
    fontsize = 8
    margin = 0.5
    xlen = max([sum([len(ref_unit[1]) for ref_unit in p[1]]) + 4
               for p in pileup_list])
    ylen = 1
    for site_pileup in pileup_list:
        ylen += 4 + 3 * len(site_pileup[2]) + sum([len(k[1]) for k in site_pileup[2]])
    
    fig = plt.figure(figsize=(xlen * xscale + 2 *
                                margin,  ylen * yscale + 2 * margin))
    fig.suptitle(output_prefix + " - Repeat ID:" + repeat_id + ' - Repeat Unit:' +
                 repeat_units[repeat_id], fontsize=4 * fontsize, y=1 - margin / (ylen * yscale + 2 * margin), va='bottom')
    ax = fig.add_axes([margin / (xlen * xscale + 1), margin / (ylen * yscale + 1),
                       xlen * xscale / (xlen * xscale + 2 * margin), ylen * yscale / (ylen * yscale + 2 * margin)])
    ax.set_xlim(0, xlen)
    ax.set_ylim(ylen, 0)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ypos = 1
    ax.axhline(ypos + 1, linestyle='--', c='0.5')
    for site_pileup in pileup_list:
        ref_units = site_pileup[1]
        if site_pileup[0] in genotypes and repeat_id in genotypes[site_pileup[0]]:
            gtype_str = '/'.join(map(str, genotypes[site_pileup[0]][repeat_id]))
        else:
            gtype_str = "None"
        sample_str = site_pileup[0] + ' - Genotypes predicted:' + gtype_str
        ax.text(xlen / 2, ypos + 3, sample_str, ha='center',
                va='bottom', fontsize=3 * fontsize)
        ypos += 3
        for gtype_pileup in site_pileup[2]:
            ax.axhline(ypos + 0.5, linestyle=':', c='0.25')
            if gtype_pileup[0][1] == 0:
                gtype_str = "Genotype = " + str(gtype_pileup[0][0])
            else:
                gtype_str = "Genotype >= " + str(gtype_pileup[0][0])
            ax.text(xlen / 2, ypos + 2, gtype_str, ha='center', va='bottom', fontsize=2 * fontsize)
            ypos += 2
            pxpos = -1
            for units in site_pileup[1]:
                if pxpos != -1:
                    for x in range(pxpos, units[0]):
                        ax.text(x + 0.5, ypos + 1, '-', ha='center', va='center',
                                fontsize=fontsize, color='k', weight='bold')
                for i in enumerate(units[1]):
                    if color:
                        c = basepair_color[i[1].upper()]
                        ax.text(units[0] + i[0] + 0.5, ypos + 1, i[1], ha='center',
                                va='center', fontsize=fontsize, weight='bold', color=c)
                    else:
                        ax.text(units[0] + i[0] + 0.5, ypos + 1, i[1], ha='center',
                                va='center', fontsize=fontsize, weight='bold')
                    pxpos = units[0] + i[0] + 1
            for i in range(int(len(site_pileup[1][1][1]) / len(repeat_units[repeat_id])) + 1):
                if len(repeat_units[repeat_id]) <= 3 and i % 5 != 0:
                    continue
                xpos = site_pileup[1][1][0] + i * len(repeat_units[repeat_id])
                ax.plot((xpos, xpos), (ypos, ypos + 1.5), color='k')
                ax.text(xpos - 0.5, ypos + 0.25, str(i), ha='right', va='center', fontsize=fontsize)
            ypos += 1
            for read_units in gtype_pileup[1]:
                pxpos = -1
                for (ref_unit, unit) in zip(ref_units, read_units):
                    if pxpos != -1 and len(unit[1]) > 0:
                        for x in range(pxpos, unit[0]):
                            ax.text(x + 0.5, ypos + 1, '-', ha='center', va='center', fontsize=fontsize, color='k', weight='bold')
                    for i in enumerate(unit[1]):
                        if color:
                            c = basepair_color[i[1].upper()]
                        elif ref_unit[1][i[0] + unit[0] - ref_unit[0]] != 'R' and ref_unit[1][i[0] + unit[0] - ref_unit[0]] != i[1].upper():
                            c = 'r'
                        elif ref_unit[1][i[0] + unit[0] - ref_unit[0]] == 'R' and repeat_units[repeat_id][(i[0] + unit[0] - ref_unit[0]) % len(repeat_units[repeat_id])] != i[1].upper():
                            c = 'r'
                        elif i[1].islower():
                            c = '0.5'
                        else:
                            c = 'k'
                        ax.text(unit[0] + i[0] + 0.5, ypos + 1, i[1], ha='center',
                                va='center', fontsize=fontsize, color=c, weight='bold')
                        pxpos = unit[0] + i[0] + 1
                ypos += 1
        ax.axhline(ypos + 1, linestyle='--', c='0.5')
    fig.savefig(output_prefix + '_' + repeat_id + '.pdf')
    fig.savefig(output_prefix + '_' + repeat_id + '.png')
    plt.close(fig)

def get_args(argv):
    """ Get arguments -- parser.parse_args(['--foo', 'FOO'])
    """
    parser = argparse.ArgumentParser(
        description="Creates pileup view of reads from ExpansionHunter read alignment log file")
    parser.add_argument("--json", dest='json_file', type=str,
                        help="Json output from EH. Required with argument --read_align")
    read_align_group = parser.add_mutually_exclusive_group(required=True)
    read_align_group.add_argument(
        "--read_align", dest='read_align_file', type=str, help="Read alignment log file from EH output")
    read_align_group.add_argument(
        "--read_align_list", dest='read_align_file_list', type=str, help="3 column file with list of read alignment outputs from EH for multiple samples. Column1: Sample name, Column2: Json file from EH output, Column3: Read alignment log file from EH output")
    parser.add_argument(
        "--output_folder", dest='output_folder', type=str, help="Output folder where PDF and PNG will be created")
    parser.add_argument("--color", action="store_true",
                        help="Flag to color by nucleotide")
    parser.add_argument("--repeat_id", dest='repeat_id',
                        type=str, help="Comma-separated Repeat IDs for which to plot pileup. Default: All repeats", default="")

    return parser.parse_args(argv)

def get_sample_paths(args):
    """Parse arguments to get sample paths and output prefix
    """
    # List of 3-tuples for each samples: (sample_name, read_align_file, json_file)
    sample_list = []
    # output_prefix: prefix of final output images
    # Populate sample_list and output_prefix
    if args.read_align_file is not None:
        if args.json_file is None:
            raise("EH json output not provided")
        sample_name = os.path.splitext(
            os.path.basename(args.read_align_file))[0]
        sample_list = [(sample_name, args.read_align_file, args.json_file)]
        output_prefix = sample_name
    else:
        in_root = os.path.dirname(args.read_align_file_list)
        for l in open(args.read_align_file_list):
            ll = l.strip().split()
            sample_list.append((ll[0], os.path.join(
                in_root, ll[2]), os.path.join(in_root, ll[1])))
        output_prefix = os.path.splitext(
            os.path.basename(args.read_align_file_list))[0]
    return (sample_list, output_prefix)

def populate_repeat_units(sample_list):
    """Parse EH JSON outputs from list of samples to populate repeat unit for each repeat ID
    """
    for sample in sample_list:
        eh_json = json.load(open(sample[2]))
        for repeat_id in eh_json:
            if repeat_id != 'BamStats':
                repeat_units[repeat_id] = eh_json[repeat_id]['RepeatUnit']
    return

def get_EH_genotypes(sample_list):
    """Parse EH JSON outputs from list of samples to obtain reported genotypes for each sample.
    """
    # List of genotypes for each sample for each repeat_id
    # genotype[sample][repeat_id] = "genotype1/genotype2"
    genotypes = {}
    for sample in sample_list:
        eh_json = json.load(open(sample[2]))
        # extract genotypes from json
        sample_genotypes = {repeat_id: eh_json[repeat_id]['Genotype']
                            for repeat_id in eh_json if repeat_id != 'BamStats'}
        # convert genotypes to list of integers
        sample_genotypes = {repeat_id: map(int, [g for g in sample_genotypes[repeat_id].split(
            '/') if g != '']) for repeat_id in sample_genotypes}
        genotypes[sample[0]] = sample_genotypes
    return genotypes

def main(argv=sys.argv[1:]):
    args = get_args(argv)
    (sample_list, output_prefix) = get_sample_paths(args)
    populate_repeat_units(sample_list)
    genotypes = get_EH_genotypes(sample_list)

    # dict mapping from sample to pileup_coordinates for each (selected) STR region in the sample
    all_sample_read_aligns = {s[0]: parse_read_align(s[1], args.repeat_id.split(',') if args.repeat_id != '' else []) for s in sample_list}
    # set of all (selected) sites reported in all samples
    site_list = set([])
    for sample in all_sample_read_aligns:
        for repeat_id in all_sample_read_aligns[sample]:
            site_list.add(repeat_id)

    for repeat_id in site_list:
        pileup_list = get_pileup_list(repeat_id, all_sample_read_aligns, [s[0] for s in sample_list])
        str_plots = os.path.join(args.output_folder, "str_plots")
        if not os.path.exists(str_plots):
            os.mkdir(str_plots)
        output_dir = os.path.join(str_plots, output_prefix)
        plot_pileup(repeat_id, pileup_list, genotypes, args.color, output_dir)



if __name__ == '__main__':
    main()
