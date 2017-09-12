#!/usr/bin/env python3
import argparse
import collections
import copy
import logging
import os
import os.path
import re
import subprocess as sp
import sys

__version__ = 0.1


def parse_params(args):
    """"""
    opts = collections.defaultdict(dict)
    opts.update({
        'pdf_plots': 1,
        'use_sample_names': 0,
        'verbose': 1,
        'make_pdf': 1,
        'make_plots': 1,
        'merge': 0,
        'args': ' '.join(os.path.basename(arg) for arg in sys.argv),
        'img_width': 11 / 2.54,
        'img_height': 10 / 2.54,
        'id2col': ['orange', 'red', 'darkgreen'],
        'tex': {
            'slide3v': {'height1': '7cm', 'height2': '7cm',  'height3': '4.5cm'},
            'slide3h': {'width1': '15cm', 'width2': '10cm', 'width3': '8cm'}},
        # for file version sanity check
        'sections': [
            {
                'id': 'ID',
                'header': 'Definition of sets',
                'exp': '# ID\t[2]id\t[3]tab-separated file names'
            },
            {
                'id': 'SN',
                'header': 'SN, mysummary numbers',
                'exp': '# SN\t[2]id\t[3]key\t[4]value'
            },
            {
                'id': 'TSTV',
                'header': '# TSTV, transition/transversions:',
                'exp': '# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)'
            },
            {
                'id': 'SiS',
                'header': 'Sis, Singleton stats',
                'exp':
                    '# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t'
                    '[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t'
                    '[9]repeat-inconsistent\t[10]not applicable'
            },
            {
                'id': 'AF',
                'header': 'AF, Stats by non-reference allele frequency',
                'exp':
                    '# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t'
                    '[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent'
                    '\t[10]not applicable'
            },
            {
                'id': 'IDD',
                'header': 'IDD, InDel distribution',
                'exp': '# IDD\t[2]id\t[3]length (deletions negative)\t[4]count'
            },
            {
                'id': 'ST',
                'header': 'ST, Substitution types',
                'exp': '# ST\t[2]id\t[3]type\t[4]count'
            },
            {
                'id': 'GCsAF',
                'header': 'GCsAF, Genotype concordance by non-reference allele frequency (SNPs)',
                'exp':
                    '# GCsAF\t[2]id\t[3]allele frequency\t[4]RR Hom matches\t[5]RA Het matches\t[6]AA Hom matches\t'
                    '[7]RR Hom mismatches\t[8]RA Het mismatches\t[9]AA Hom mismatches\t[10]dosage r-squared\t'
                    '[11]number of genotypes'
            },
            {
                'id': 'GCiAF',
                'header': 'GCiAF, Genotype concordance by non-reference allele frequency (indels)',
                'exp':
                    '# GCiAF\t[2]id\t[3]allele frequency\t[4]RR Hom matches\t[5]RA Het matches\t'
                    '[6]AA Hom matches\t[7]RR Hom mismatches\t[8]RA Het mismatches\t[9]AA Hom mismatches\t'
                    '[10]dosage r-squared\t[11]number of genotypes'
            },
            {
                'id': 'NRDs',
                'header': 'Non-Reference Discordance (NRD), SNPs',
                'exp':
                    '# NRDs\t[2]id\t[3]NRD\t[4]Ref/Ref discordance\t[5]Ref/Alt discordance\t[6]Alt/Alt discordance'
            },
            {
                'id': 'NRDi',
                'header': 'Non-Reference Discordance (NRD), indels',
                'exp':
                    '# NRDi\t[2]id\t[3]NRD\t[4]Ref/Ref discordance\t[5]Ref/Alt discordance\t[6]Alt/Alt discordance'
            },
            {
                'id': 'GCsS',
                'header': 'GCsS, Genotype concordance by sample (SNPs)',
                'exp':
                    '# GCsS\t[2]id\t[3]sample\t[4]non-reference discordance rate\t[5]RR Hom matches\t['
                    '6]RA Het matches\t[7]AA Hom matches\t[8]RR Hom mismatches\t[9]RA Het mismatches\t'
                    '[10]AA Hom mismatches\t[11]dosage r-squared'
            },
            {
                'id': 'GCiS',
                'header': 'GCiS, Genotype concordance by sample (indels)',
                'exp':
                    '# GCiS\t[2]id\t[3]sample\t[4]non-reference discordance rate\t[5]RR Hom matches\t'
                    '[6]RA Het matches\t[7]AA Hom matches\t[8]RR Hom mismatches\t[9]RA Het mismatches\t'
                    '[10]AA Hom mismatches\t[11]dosage r-squared'
            },
            {
                'id': 'PSC',
                'header': 'PSC, Per-sample counts',
                'exp':
                    '# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t'
                    '[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons'
            },
            {
                'id': 'PSI',
                'header': 'PSI, Per-sample Indels',
                'exp':
                    '# PSI\t[2]id\t[3]sample\t[4]in-frame\t[5]out-frame\t[6]not applicable\t'
                    '[7]out/(in+out) ratio\t[8]nHets\t[9]nAA'
            },
            {
                'id': 'DP',
                'header': 'DP, Depth distribution',
                'exp':
                    '# DP\t[2]id\t[3]bin\t[4]number of genotypes\t[5]fraction of genotypes (%)\t'
                    '[6]number of sites\t[7]fraction of sites (%)'
            },
            {
                'id': 'FS',
                'header': 'FS, Indel frameshifts',
                'exp':
                    '# FS\t[2]id\t[3]in-frame\t[4]out-frame\t[5]not applicable\t[6]out/(in+out) ratio\t'
                    '[7]in-frame (1st ALT)\t[8]out-frame (1st ALT)\t[9]not applicable (1st ALT)\t'
                    '[10]out/(in+out) ratio (1st ALT)'
            },
            {
                'id': 'ICS',
                'header': 'ICS, Indel context mysummary',
                'exp':
                    '# ICS\t[2]id\t[3]repeat-consistent\t[4]repeat-inconsistent\t[5]not applicable\t'
                    '[6]c/(c+i) ratio'
            },
            {
                'id': 'ICL',
                'header': 'ICL, Indel context by length',
                'exp':
                    '# ICL\t[2]id\t[3]length of repeat element\t[4]repeat-consistent deletions)\t'
                    '[5]repeat-inconsistent deletions\t[6]consistent insertions\t[7]inconsistent insertions\t'
                    '[8]c/(c+i) ratio'
            },
            {
                'id': 'QUAL',
                'header': 'QUAL, Stats by quality',
                'exp':
                    '# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t'
                    '[6]number of transversions (1st ALT)\t[7]number of indels'
            },
            {
                'id': 'HWE',
                'header': 'HWE',
                'exp':
                    '# HWE\t[2]id\t[3]1st ALT allele frequency\t[4]Number of observations\t[5]25th percentile\t'
                    '[6]median\t[7]75th percentile'}],
        'SN_keys': [
            'number of samples:',
            'number of records:',
            'number of no-ALTs:',
            'number of SNPs:',
            'number of MNPs:',
            'number of indels:',
            'number of others:',
            'number of multiallelic sites:',
            'number of multiallelic SNP sites:']})
    for sec in opts['sections']:
        opts['exp'][sec['id']] = sec['exp']
        opts['id2sec'][sec['id']] = sec
    if args.no_PDF:
        opts['make_pdf'] = 0
    if args.rasterize:
        opts['rasterize'] = 1
        opts['pdf_plots'] = 0
    if args.merge:
        opts['make_plots'] = 0
        opts['make_pdf'] = 0
        opts['merge'] = 1
    if args.sample_names:
        opts['use_sample_names'] = 1
    if args.title:
        opts['titles'] = args.title
    if args.main_title:
        opts['main_title'] = args.main_title
    if args.prefix:
        opts['prefix'] = args.prefix
    else:
        opts['prefix'] = '.'
    opts['vcfstats'] = args.stats_files
    opts['dir'] = opts['prefix']
    opts['logfile'] = 'plot-vcfstats.log'
    if not os.path.isdir(opts['dir']):
        os.mkdir(opts['dir'])
    return opts


def parse_vcfstats(opts):
    """"""
    for i in range(len(opts['vcfstats'])):
        parse_vcfstats1(opts, i)

    # Check sanity
    if 0 not in opts['dat']['ID']:
        error("Sanity check failed: no stats found by vcfstats??")

    # Set titles
    file2title = {}
    title2file = {}
    if 'titles' in opts:
        for i in range(len(opts['titles'])):
            if i not in opts['dat']['ID']:
                continue
            file2title[opts['dat']['ID'][i][0][0]] = opts['titles'][i]
            title2file[opts['titles'][i]] = opts['dat']['ID'][i][0][0]
    for i_d in file_ids(opts):
        if len(opts['dat']['ID'][i_d][0]) > 1:
            continue
        file = opts['dat']['ID'][i_d][0][0]
        if file not in file2title:  # create short title
            bname = file
            bname = re.sub('^.*/', '', bname)
            bname = re.sub('\.vcf\.gz$', '', bname, flags=re.IGNORECASE)
            if len(bname) > 5:
                bname = bname[0:5]
            i = 0
            title = bname
            while title in title2file:
                title = bname + chr(66 + i)
                i += 1
            file2title[file] = title
            title2file[title] = file
    for i_d in file_ids(opts):
        titles = []
        for file in opts['dat']['ID'][i_d][0]:
            if file in file2title:
                titles.append(file2title[file])
        opts['title'][i_d] = ' + '.join(titles)

    # mapping from file names to list of IDs
    for i_d in file_ids(opts):
        for file in opts['dat']['ID'][i_d][0]:
            if 'file2ids' not in opts:
                opts['file2ids'] = {}
            if file not in opts['file2ids']:
                opts['file2ids'][file] = []
            opts['file2ids'][file].append(i_d)

    # check sanity of the file merge: were the correct files merged?
    if 'coalesced_files' in opts and opts['verbose']:
        logging.info('The vcfstats outputs have been merged as follows:')
        printed = {}
        for i_d in opts['coalesced_files']:
            for i in range(len(opts['coalesced_files'][i_d])):
                for j in range(len(opts['coalesced_files'][i_d][i])):
                    if opts['dat']['ID'][i_d][i][j]:
                        continue
                    logging.info('\t' + opts['dat']['ID'][i_d][i][j])
                    for file in opts['coalesced_files'][i_d][i][j]:
                        n = opts['coalesced_files'][i_d][i][j][file]
                        logging.info('\t\t' + file + ('\t..\t' + n + 'x') if n > 1 else '')
                    printed[opts['dat']['ID'][i_d][i][j]] = 1


def parse_vcfstats1(opts, i):
    """"""
    file = opts['vcfstats'][i]
    logging.info('Parsing bcftools stats output: {}'.format(file.name))
    header = next(file)
    if not re.search(r'^# This file was produced by \S*', header):
        error('Sanity check failed: was this file generated by bcftools stats?')
    dat = collections.defaultdict(dict)
    def_line_regex = re.compile('^#\s+(\S+)\t')
    for line in file:
        line = line.strip()
        def_line = re.search(def_line_regex, line)
        if def_line:
            opts['def_line'][def_line.group(1)] = line
            continue
        if line.startswith('#'):
            continue
        items = line.split('\t')
        for idx, item in enumerate(items):
            try:
                item = int(item)
            except ValueError:
                try:
                    item = float(item)
                except ValueError:
                    pass
            items[idx] = item
        if items[0] == 'SN':
            dat[items[1]][items[2]] = int(items[3])
            continue
        if items[1] not in dat[items[0]]:
            dat[items[0]][items[1]] = []
        dat[items[0]][items[1]].append(items[2:])
    for a in dat:
        if a not in opts['dat']:  # First vcfstats file
            opts['dat'][a] = dat[a]
            continue
        for b in dat[a]:  # Merging multiple vcfstats files. Honestly, this is quite hacky.
            if b not in opts['dat'][a]:  # copy all, first occurrence
                opts['dat'][a][b] = dat[a][b]
                continue
            if a == 'ID':
                merge_id(opts, opts['dat'][a], dat[a], b)
            elif not isinstance(dat[a][b], list):  # SN, mysummary numbers, do not mysum sample counts
                if b != 'number of samples:':
                    opts['dat'][a][b] += dat[a][b]
            elif a == 'NRDs':
                add_to_avg(opts['dat'][a][b], dat[a][b], i)
            elif a == 'NRDi':
                add_to_avg(opts['dat'][a][b], dat[a][b], i)
            elif a == 'DP':
                merge_dp(opts['dat'][a][b], dat[a][b])
            elif a == 'GCsS':
                merge_GCsS(opts['dat'][a][b], dat[a][b], i)
            elif a == 'GCiS':
                merge_GCsS(opts['dat'][a][b], dat[a][b], i)
            elif a == 'GCsAF':
                merge_GCsAF(opts['dat'][a][b], dat[a][b])
            elif a == 'GCiAF':
                merge_GCsAF(opts['dat'][a][b], dat[a][b])
            elif a == 'ST':
                add_to_values(opts['dat'][a][b], dat[a][b], cmp_str)
            elif a == 'PSC':
                merge_PSC(opts['dat'][a][b], dat[a][b], i)
            elif a == 'PSI':
                merge_PSI(opts['dat'][a][b], dat[a][b])
            elif a == 'IDD':
                add_to_values(opts['dat'][a][b], dat[a][b], cmp_num)
            elif a == 'FS':
                merge_FS(opts['dat'][a][b], dat[a][b])
            elif a == 'ICS':
                merge_ICS(opts['dat'][a][b], dat[a][b])
            elif a == 'ICL':
                merge_ICL(opts['dat'][a][b], dat[a][b])
            elif a == 'TSTV':
                merge_TSTV(opts['dat'][a][b], dat[a][b])
            elif a == 'DBG':
                continue
            else:
                add_to_values(opts['dat'][a][b], dat[a][b], cmp_num_op)


def merge_id(opts, dst, src, i_d):
    """"""
    for i in range(len(src[i_d])):
        for j in range(len(src[i_d][i])):
            gname = rglob(dst[i_d][i][j], src[i_d][i][j])
            dst[i_d][i][j] = gname
            if i_d not in opts['coalesced_files']:
                opts['coalesced_files'][i_d] = []
            if i + 1 > len(opts['coalesced_files'][i_d]):
                opts['coalesced_files'][i_d].append([])
            if j + 1 > len(opts['coalesced_files'][i_d][i]):
                opts['coalesced_files'][i_d][i].append({src[i_d][i][j]: 0})
            opts['coalesced_files'][i_d][i][j][src[i_d][i][j]] += 1


def merge_dp(a, b):
    """"""
    add_to_values(a, b, cmp_num_op)
    # recalculate fraction of GTs and fraction of sites, cannot be simply mysummed
    gsum = 0  # genotype mysum
    ssum = 0  # site mysum
    for i in range(len(a)):
        gsum += a[i][1]
        if len(a[i]) > 3:
            ssum += a[i][3]
        else:
            #  older stats files will not have last 2 columns for (number of sites, fraction of sites), so fill in as zero
            a[i].extend((0, 0))
    for i in range(len(a)):
        a[i][2] = a[i][1] * 100 / gsum if gsum else 0
        a[i][4] = a[i][3] * 100 / ssum if ssum else 0


def merge_GCsS(a, b, n):
    """average the non-ref discordance rate"""
    for i in range(len(a)):
        a[i][1] *= n
    add_to_sample_values(a, b)
    for i in range(len(a)):
        a[i][1] /= n + 1


def merge_GCsAF(a, b):
    """recalculate r2"""
    for i in range(len(a)):
        a[i][7] *= a[i][8]
    for i in range(len(b)):
        a[i][7] *= a[i][8]
    add_to_values(a, b, cmp_num_op)
    for i in range(len(a)):
        a[i][7] /= a[i][8]


def merge_PSC(a, b, n):
    """"""
    for i in range(len(a)):
        a[i][7] *= n
    add_to_sample_values(a, b)
    for i in range(len(a)):
        a[i][7] /= n + 1


def merge_PSI(a, b):
    """"""
    add_to_sample_values(a, b)
    for i in range(len(b)):
        a[i][4] = float('{:.2f}'.format(a[i][2] / (a[i][1] + a[i][2]) if a[i][1] + a[i][2] else 0))


def merge_FS(a, b):
    """"""
    for i in range(len(a)):
        for j in range(3):
            a[i][j] += b[i][j]
        a[i][3] = float('{:.2f}'.format(a[i][1] / (a[i][0] + a[i][1]) if a[i][0] + a[i][1] else 0))
        for j in range(4, 7):
            a[i][j] += b[i][j]
        a[i][7] = float('{:.2f}'.format(a[i][5] / (a[i][4] + a[i][5]) if a[i][4] + a[i][5] else 0))


def merge_ICS(a, b):
    """"""
    for i in range(len(a)):
        for j in range(3):
            a[i][j] += b[i][j]
        a[i][3] = '{:.4f}'.format(a[i][0] / (a[i][0] + a[i][1]) if a[i][0] + a[i][1] else 0)


def merge_ICL(a, b):
    """"""
    for i in range(len(a)):
        for j in range(1, 5):
            a[i][j] += b[i][j]
        a[i][5] = '{:.4f}'.format((a[i][1] + a[i][3]) / (a[i][1] + a[i][2] + a[i][3] + a[i][4]) if a[i][2] + a[i][4] else 0)


def merge_TSTV(a, b):
    """"""
    for i in range(len(a)):
        for j in range(2):
            a[i][j] += b[i][j]
        a[i][2] = float('{:.2f}'.format(a[i][0] / a[i][1] if a[i][1] else 0))
        for j in range(3, 5):
            a[i][j] += b[i][j]
        a[i][5] = float('{:.2f}'.format(a[i][3] / a[i][4] if a[i][4] else 0))


def add_to_avg(dst, src, n):
    """"""
    for i in range(len(src)):
        if isinstance(dst[i], list):
            for j in range(len(dst[i])):
                dst[i][j] = (n * dst[i][j] + src[i][j]) / (n + 1)
        else:
            dst[i] = (n * dst[i] + src[i]) / (n + 1)


def rglob(a, b):
    """"""
    if a == b:
        return a
    a = re.sub(r'\\*', '', a)
    la = len(a)
    lb = len(b)
    i = 0
    while (i < la) & (i < lb) & (a[i] == b[i]):
        i += 1
    la -= 1
    lb -= 1
    while (la > i) & (lb > i) & (a[la] == b[lb]):
        la -= 1
        lb -= 1
    la = 1 if la == i and lb == i else la - i
    a = a[:i] + '*' + a[i + la:]
    return a


def add_to_values(dst, src, cmp):
    """"""
    i_d = 0
    i_s = 0
    while i_s < len(src):
        while i_d < len(dst) and cmp(src[i_s][0], dst[i_d][0]) > 0:
            i_d += 1
        if i_d < len(dst) and not cmp(src[i_s][0], dst[i_d][0]):
            for j in range(1, len(src[i_s])):
                dst[i_d][j] += src[i_s][j]
        else:
            dst.insert(i_d, src[i_s])
        i_s += 1


def cmp_num_op(a, b):
    """numeric compare with operators
    Cases like <3, >500 make it complicated
    """
    xa = '='
    xb = '='
    a = str(a)
    b = str(b)
    char_start = re.compile('^(\D+)(.*)')
    a_match = re.search(char_start, a)
    if a_match:
        xa = a_match.group(1)
        a = a_match.group(2)
    b_match = re.search(char_start, b)
    if b_match:
        xb = b_match.group(1)
        b = b_match.group(2)
    if a == b:
        return (xa > xb) - (xa < xb)
    try:
        a = int(a)
    except ValueError:
        try:
            a = float(a)
        except ValueError:
            pass
    try:
        b = int(b)
    except ValueError:
        try:
            b = float(b)
        except ValueError:
            pass
    return (a > b) - (a < b)


def cmp_str(a, b):
    """"""
    return (a > b) - (a < b)


def cmp_num(a, b):
    """"""
    return (a > b) - (a < b)


def add_to_sample_values(dst, src):
    """"""
    id2i = {}
    for i in range(len(dst)):
        id2i[dst[i][0]] = i
    for i in range(len(src)):
        if src[i][0] not in id2i:
            error('Whoops, no such dst sample: {}'.format(src[i][0]))
        di = id2i[src[i][0]]
        for j in range(len(src[i])):
            dst[di][j] += src[i][j]


def file_ids(opts):
    i_d = 0
    out = []
    while 'ID' in opts['dat'] and i_d in opts['dat']['ID']:
        out.append(i_d)
        i_d += 1
    return out


def merge_vcfstats(opts):
    """"""
    fh = open('merge.chk', 'w') if not opts['merge'] else sys.stdout
    fh.write('# This file was produced by plot-vcfstats, the command line was:\n'
             '#   ' + opts['args'] + '\n#\n')
    for sec in opts['sections']:
        sid = sec['id']
        if sid not in opts['dat']:
            continue
        fh.write('# ' + sec['header'] + '\n' + sec['exp'] + '\n')
        for i_d in sorted(opts['dat'][sid]):
            for rec in opts['dat'][sid][i_d]:
                fh.write('{}\t{}\t{}\n'.format(sid, i_d, '\t'.join(str(r) for r in rec)))
        if sid == 'ID':
            fh.write('# ' + opts['id2sec']['SN']['header'] + '\n' + opts['id2sec']['SN']['exp'] + '\n')
            # output mysummary numbers here
            for i_d in opts['dat']:
                if i_d in opts['dat']:
                    continue
                for key in opts['SN_keys']:
                    if key not in opts['dat'][i_d]:
                        continue
                    fh.write('SN\t' + 'id' + '\t' + key + '\t' + opts['dat'][i_d][key] + '\n')
    fh.close()


def init_plots(opts):
    """"""
    opts['plt_file'] = 'plot.py'
    titles = '# Title abbreviations:\n'
    for i_d in file_ids(opts):
        titles += '# \t {} .. {} .. {}\n'.format(i_d, opts['title'][i_d], opts['dat']['ID'][i_d][0][0])
    titles += '#'
    fh = open(opts['plt_file'], 'w')
    tprint(fh, '''
        # This file was produced by plot-vcfstats, the command line was:
        #   {args}
        #
        # Edit as necessary and recreate the plots by running
        #   python {plt_file}
        #
        {titles}

        # Set to 1 to plot in PDF instead of PNG
        pdf_plots = {pdf_plots}

        # Use logarithimic X axis for allele frequency plots
        af_xlog = 0

        # Plots to generate, set to 0 to disable
        plot_venn_snps = 1
        plot_venn_indels = 1
        plot_tstv_by_sample = 1
        plot_hethom_by_sample = 1
        plot_snps_by_sample = 1
        plot_indels_by_sample = 1
        plot_singletons_by_sample = 1
        plot_depth_by_sample = 1
        plot_SNP_count_by_af = 1
        plot_Indel_count_by_af = 1
        plot_SNP_overlap_by_af = 1
        plot_Indel_overlap_by_af = 1
        plot_dp_dist = 1
        plot_hwe = 1
        plot_concordance_by_af = 1
        plot_r2_by_af = 1
        plot_discordance_by_sample = 1
        plot_tstv_by_af = 1
        plot_indel_dist = 1
        plot_tstv_by_qual = 1
        plot_substitutions = 1


        # Set to 1 to use sample names for xticks instead of numeric sequential IDs
        #   and adjust margins and font properties if necessary
        sample_names   = {use_sample_names}
        sample_margins = {{'right':0.98, 'left':0.07, 'bottom':0.2}}
        sample_font    = {{'rotation':45, 'ha':'right', 'fontsize':8}}

        if sample_names==0: sample_margins=(); sample_font=();


        #-------------------------------------------------


        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt

        import csv
        csv.register_dialect('tab', delimiter='\\t', quoting=csv.QUOTE_NONE)

        import numpy
        def smooth(x,window_len=11,window='hanning'):
        \\tif x.ndim != 1: raise ValueError("The function 'smooth' only accepts 1 dimension arrays.")
        \\tif x.size < window_len: return x
        \\tif window_len<3: return x
        \\tif not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']: raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
        \\ts = numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
        \\tif window == 'flat': # moving average
        \\t\\tw = numpy.ones(window_len,'d')
        \\telse:
        \\t\\tw = eval('numpy.'+window+'(window_len)')
        \\ty = numpy.convolve(w/w.sum(),s,mode='valid')
        \\treturn y[(window_len/2-1):-(window_len/2)]

    '''.format(
        args=opts['args'], plt_file=opts['plt_file'], titles=titles, pdf_plots=opts['pdf_plots'],
        use_sample_names=opts['use_sample_names']))
    opts['plt_fh'] = fh


def plot_venn_bars(opts):
    """"""
    ids = file_ids(opts)
    if len(ids) != 3:
        return

    snps = []
    indels = []
    tstv = []
    snp_titles = []
    indel_titles = []
    for i_d in range(3):
        snps.append(get_value(opts, i_d, 'number of SNPs:'))
        indels.append(get_value(opts, i_d, 'number of indels:'))
        tstv.append(float('{:.2f}'.format(get_values(opts, i_d, 'TSTV', 0, 5))))
        snp_titles.append('{}\\nts/tv {}\\n'.format(opts['title'][i_d], tstv[i_d]) + bignum(snps[i_d]))
        fs = get_values(opts, i_d, 'FS')
        fss = 'frm{}\\n'.format(fs[0][3]) if len(fs) else ''
        indel_titles.append('{}\\n{}'.format(opts['title'][i_d], fss) + bignum(indels[i_d]))

    fh = opts['plt_fh']
    tprint(fh, '''
    
                if plot_venn_snps:
            \\tfig = plt.figure(figsize=({img_width},{img_height}))
            \\tax1 = fig.add_subplot(111)
            \\tax1.bar([1,2,3],[{snps_0},{snps_2},{snps_1}],align='center',color='{id2col_0}',width=0.3)
            \\tax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            \\tax1.set_xlim(0.5,3.5)
            \\tplt.xticks([1,2,3],('{snp_titles_0}','{snp_titles_2}','{snp_titles_1}'))
            \\tplt.title('Number of SNPs')
            \\tplt.subplots_adjust(right=0.95,bottom=0.15)
            \\tplt.savefig('venn_bars.snps.png')
            \\tif pdf_plots: plt.savefig('venn_bars.snps.pdf')
            \\tplt.close()


            if plot_venn_indels:
            \\tfig = plt.figure(figsize=({img_width},{img_height}))
            \\tax1 = fig.add_subplot(111)
            \\tax1.bar([1,2,3],[{indels_0},{indels_2},{indels_1}],align='center',color='{id2col_1}',width=0.3)
            \\tax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            \\tax1.set_xlim(0.5,3.5)
            \\tplt.xticks([1,2,3],('{indel_titles_0}','{indel_titles_2}','{indel_titles_1}'))
            \\tplt.title('Number of indels')
            \\tplt.subplots_adjust(right=0.95,bottom=0.15)
            \\tplt.savefig('venn_bars.indels.png')
            \\tif pdf_plots: plt.savefig('venn_bars.indels.pdf')
            \\tplt.close()
    
    '''.format(img_width=opts['img_width'], img_height=opts['img_height'], snps_0=snps[0], snps_1=snps[1],
               snps_2=snps[2], id2col_0=opts['id2col'][0], snp_titles_0=snp_titles[0], snp_titles_1=snp_titles[1],
               snp_titles_2=snp_titles[2], indels_0=indels[0], indels_1=indels[1], indels_2=indels[2],
               id2col_1=opts['id2col'][1], indel_titles_0=indel_titles[0], indel_titles_1=indel_titles[1],
               indel_titles_2=indel_titles[2]))


def plot_counts_by_AF(opts):
    """"""
    plot_counts_by_AF_col(opts, 1, 'SNP')
    plot_counts_by_AF_col(opts, 4, 'Indel')


def plot_counts_by_AF_col(opts, col, title):
    """"""
    fh = opts['plt_fh']
    img = 'counts_by_af.' + title.lower() + 's'
    tfh = open('{}.dat'.format(img), 'w')
    tfh.write('# [1]id\t[2]Nonref Allele Frequency\t[3]Number of sites\n')
    for i_d in file_ids(opts):
        tmp = get_values(opts, i_d, 'AF')
        vals = rebin_values(tmp, 1, 0)
        for val in vals:
            if not val[col]:
                continue
            tfh.write('{}\t{}\t{}\n'.format(i_d, val[0], val[col]))
    tfh.close()

    tprint(fh, '''

            dat = {{}}
            with open('{img}.dat', 'r') as f:
            \\treader = csv.reader(f, 'tab')
            \\tfor row in reader:
            \\t\\tif row[0][0] == '#': continue
            \\t\\tid = int(row[0])
            \\t\\tif id not in dat: dat[id] = []
            \\t\\tdat[id].append([float(row[1]),float(row[2])])

            if plot_{title}_count_by_af:
            \\tfig = plt.figure(figsize=(2*{img_width},{img_height}*0.7))
            \\tax1 = fig.add_subplot(111)
            \\tax1.set_ylabel('Number of sites')
            \\tax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            \\tax1.set_yscale('log')
            \\tif af_xlog: ax1.set_xscale('log')
            \\tax1.set_xlabel('Non-reference allele frequency')
            \\tax1.set_xlim(-0.05,1.05)
            \\thas_data = 0
    '''.format(img=img, title=title, img_width=opts['img_width'], img_height=opts['img_height']))
    for i_d in file_ids(opts):
        tprint(fh, '''
            \\tif {id} in dat and len(dat[{id}])>2:
            \\t\\tax1.plot([row[0] for row in dat[{id}]], [row[1] for row in dat[{id}]], '-o',markersize=3, color='{id2col}',mec='{id2col}',label='{title}')
            \\t\\thas_data = 1
        '''.format(id=i_d, id2col=opts['id2col'][i_d], title=opts['title'][i_d]))
    tprint(fh, '''
            \\tif has_data:
            \\t\\tax1.legend(numpoints=1,markerscale=1,loc='best',prop={{'size':10}},frameon=False)
            \\t\\tplt.title('{title} count by AF')
            \\t\\tplt.subplots_adjust(bottom=0.2,left=0.1,right=0.95)
            \\t\\tplt.savefig('{img}.png')
            \\t\\tif pdf_plots: plt.savefig('{img}.pdf')
            \\t\\tplt.close()


    '''.format(title=title, img=img))


def plot_overlap_by_AF(opts):
    """"""
    plot_overlap_by_AF_col(opts, 1, 'SNP')
    plot_overlap_by_AF_col(opts, 4, 'Indel')


def plot_overlap_by_AF_col(opts, col, title):
    """"""
    ids = file_ids(opts)
    if len(ids) != 3:
        return
    ia = ib = iab = 0
    for i in range(len(ids)):
        if len(opts['dat']['ID'][ids[i]][0]) > 1:
            iab = i
            continue
        if not ia:
            ia = i
            continue
        ib = i
    fh = opts['plt_fh']
    img = 'overlap_by_af' + title.lower() + 's'
    vals_a = get_values(opts, ia, 'AF')
    vals_ab = get_values(opts, iab, 'AF')
    afs = {}
    af_a = {}
    af_ab = {}
    for val in vals_a:
        afs[val[0]] = val[col]
        af_a[val[0]] = val[col]
    for val in vals_ab:
        afs[val[0]] = val[col]
        af_ab[val[0]] = val[col]
    tfh = open('{}.dat'.format(img), 'w')
    tfh.write('# [1]Allele frequency\t[2]Fraction of sites from {title_1} also in {title_2}\t[3]Number of sites\n'.format(
        title_1=opts['title'][ids[ia]], title_2=opts['title'][ids[ib]]))
    for af in sorted(afs):
        a = af_a[af] if af_a[af] else 0
        ab = af_ab[af] if af_ab[af] else 0
        yval = ab * 100 / (a + ab) if (a + ab) else 0
        tfh.write('{}\t{}\t{}\n'.format(af, yval, a + ab))
    tfh.close()
    tprint(fh, '''

        dat = []

        with open('{img}.dat', 'r') as f:
        \\treader = csv.reader(f, 'tab')
        \\tfor row in reader:
        \\t\\tif row[0][0] != '#': dat.append(row)

        if plot_{title}_overlap_by_af and len(dat)>1:
        \\tfig = plt.figure(figsize=(2*{img_width},{img_height}*0.7))
        \\tax1 = fig.add_subplot(111)
        \\tax1.plot([row[0] for row in dat], [row[1] for row in dat],'-o',markersize=3, color='{id2col_1}',mec='{id2col_1}')
        \\tax1.set_ylabel('Fraction found in {title_ib} [%]')
        \\tax1.set_xscale('log')
        \\tax1.set_xlabel('Non-reference allele frequency in {title_ia}')
        \\tax1.set_xlim(0,1.01)
        \\tplt.title('{title} overlap by AF')
        \\tplt.subplots_adjust(bottom=0.2,left=0.1,right=0.95)
        \\tplt.savefig('{img}.png')
        \\tif pdf_plots: plt.savefig('{img}.pdf')
        \\tplt.close()
    '''.format(img=img, title=title, img_width=opts['img_width'], img_height=opts['img_height'],
               id2col_1=opts['id2col'][1], title_ia=opts['title'][ia], title_ib=opts['title'][ib]))


def plot_concordance_by_AF(opts):
    """"""
    vals = get_values(opts, 2, 'GCsAF')
    if not vals:
        return

    # create a local copy and prepare r2 for rebinning
    vals = copy.deepcopy(vals)
    for i in range(len(vals)):
        vals[i][7] *= vals[i][8]
    vals = rebin_values(vals, 0.01, 0)
    fh = opts['plt_fh']
    img = 'gts_by_af'
    img2 = 'r2_by_af'
    tfh = open('{}.dat'.format(img), 'w')
    tfh.write('# [1]Allele Frequency\t[2]RR concordance\t[3]RA concordance\t[4]AA concordance\t[5]nRR\t[6]nRA\t[7]nAA\t[8]R^2\t[9]Number of genotypes\n')
    for i in range(len(vals)):
        tfh.write('{:f}\t{:f}\t{:f}\t{:f}\t{:d}\t{:d}\t{:d}\t{:f}\t{:d}\n'.format(
            vals[i][0],
            vals[i][1] / (vals[i][1] + vals[i][4]) if vals[i][1] + vals[i][4] else 1,
            vals[i][2] / (vals[i][2] + vals[i][5]) if vals[i][2] + vals[i][5] else 1,
            vals[i][3] / (vals[i][3] + vals[i][6]) if vals[i][3] + vals[i][6] else 1,
            vals[i][1] + vals[i][4],
            vals[i][2] + vals[i][5],
            vals[i][3] + vals[i][6],
            vals[i][7] / vals[i][8] if vals[i][8] else 1,
            vals[i][8]))
    tfh.close()

    tprint(fh, '''

            dat = []
            with open('{img}.dat', 'r') as f:
            \\treader = csv.reader(f, 'tab')
            \\tfor row in reader:
            \\t\\tif row[0][0] != '#': dat.append(row)

            if plot_concordance_by_af and len(dat)>1:
            \\tfig = plt.figure(figsize=({img_width}*1.2,{img_height}))
            \\tax1 = fig.add_subplot(111)
            \\tax1.plot([row[0] for row in dat], [row[1] for row in dat],'.',color='{id2col_1}',label='Hom RR')
            \\tax1.plot([row[0] for row in dat], [row[2] for row in dat],'.',color='{id2col_0}',label='Het RA')
            \\tax1.plot([row[0] for row in dat], [row[3] for row in dat],'.',color='k',label='Hom AA')
            \\tax1.set_xlabel('Non-ref allele frequency')
            \\tax1.set_ylabel('Concordance')
            \\tleg = ax1.legend(title='Concordance:',numpoints=1,markerscale=2,loc='best',prop={{'size':9}})
            \\tleg.draw_frame(False)
            \\tplt.setp(leg.get_title(),fontsize=9)
            \\tax2 = ax1.twinx()
            \\tax2.plot([row[0] for row in dat], [row[4] for row in dat],color='{id2col_1}')
            \\tax2.plot([row[0] for row in dat], [row[5] for row in dat],color='{id2col_0}')
            \\tax2.plot([row[0] for row in dat], [row[6] for row in dat],color='k')
            \\tax2.set_ylabel('Number of genotypes')
            \\tax2.set_yscale('log')
            \\tif af_xlog: ax1.set_xscale('log')
            \\tif af_xlog: ax2.set_xscale('log')
            \\tplt.subplots_adjust(left=0.15,right=0.83,bottom=0.11)
            \\tplt.savefig('{img}.png')
            \\tif pdf_plots: plt.savefig('{img}.pdf')
            \\tplt.close()

            if plot_r2_by_af and len(dat)>1:
            \\tfig = plt.figure(figsize=({img_width}*1.3,{img_height}))
            \\tax1 = fig.add_subplot(111)
            \\tax2 = ax1.twinx()
            \\tax1.set_zorder(ax2.get_zorder()+1)
            \\tax1.patch.set_visible(False)
            \\tax2.plot([row[0] for row in dat], [row[8] for row in dat], '-o', color='r',mec='r',markersize=3)
            \\tax1.plot([row[0] for row in dat], [row[7] for row in dat], '-^', color='k',markersize=3)
            \\tfor tl in ax2.get_yticklabels(): tl.set_color('r')
            \\tax2.set_ylabel('Number of genotypes', color='r')
            \\tax2.set_yscale('log')
            \\tif af_xlog: ax1.set_xscale('log')
            \\tif af_xlog: ax2.set_xscale('log')
            \\tax1.set_ylabel('Aggregate allelic R\$^2\$', color='k')
            \\tax1.set_xlabel('Non-ref allele frequency')
            \\tplt.subplots_adjust(left=0.19,right=0.83,bottom=0.11)
            \\tplt.savefig('{img2}.png')
            \\tif pdf_plots: plt.savefig('{img2}.pdf')
            \\tplt.close()

    '''.format(img=img, img2=img2, img_width=opts['img_width'], img_height=opts['img_height'],
               id2col_0=opts['id2col'][0], id2col_1=opts['id2col'][1]))


def plot_concordance_by_sample(opts):
    """"""
    vals = get_values(opts, 2, 'GCsS')
    if not vals:
        return
    fh = opts['plot_fh']
    img = 'gts_by_sample'
    tfh = open('{}.dat'.format(img), 'w')
    tfh.write('# [1]Sample ID\t[2]Discordance\t[3]Sample Name\n')
    for i in range(len(vals)):
        tfh.write('{:d}\t{:f}\t{:s}\n'.format(i, vals[i][1], vals[i][0]))
    tfh.close()

    tprint(fh, '''

            dat = []
            with open('{img}.dat', 'r') as f:
            \\treader = csv.reader(f, 'tab')
            \\tfor row in reader:
            \\t\\tif row[0][0] != '#': dat.append(row)

            if plot_discordance_by_sample:
            \\tfig = plt.figure(figsize=(2*{img_width},{img_height}*0.7))
            \\tax1 = fig.add_subplot(111)
            \\tax1.plot([row[0] for row in dat], [row[1] for row in dat],'.',color='orange')
            \\tax1.set_ylabel('Non-ref discordance')
            \\tax1.set_ylim(0,)
            \\tif sample_names:
            \\t\\t     plt.xticks([int(row[0]) for row in dat],[row[2] for row in dat],**sample_font)
            \\t\\t     plt.subplots_adjust(**sample_margins)
            \\telse:
            \\t\\t     plt.subplots_adjust(right=0.98,left=0.07,bottom=0.17)
            \\t\\t     ax1.set_xlabel('Sample ID')
            \\tplt.savefig('{img}.png')
            \\tif pdf_plots: plt.savefig('{img}.pdf')
            \\tplt.close()


    '''.format(img=img, img_width=opts['img_width'], img_height=opts['img_height']))


def plot_tstv_by_AF(opts, i_d):
    """"""
    vals = get_values(opts, i_d, 'AF')
    if not vals:
        return
    fh = opts['plt_fh']
    img = 'tstv_by_af.{}'.format(i_d)
    vals = rebin_values(vals, 8, 0)

    tfh = open('{}.dat'.format(img), 'w')
    tfh.write('# [1]Allele frequency\t[2]Number of sites\t[3]ts/tv\n')
    for i in range(len(vals)):
        if vals[i][2] + vals[i][3] == 0:
            continue
        tfh.write('{:f}\t{:d}\t{:f}\n'.format(
            vals[i][0],
            vals[i][2] + vals[i][3],
            vals[i][2] / vals[i][3] if vals[i][3] else 0))
    tfh.close()

    tprint(fh, '''

            dat = []
            with open('{img}.dat', 'r') as f:
            \\treader = csv.reader(f, 'tab')
            \\tfor row in reader:
            \\t\\tif row[0][0] != '#': dat.append([float(x) for x in row])


            if plot_tstv_by_af and len(dat)>2:
            \\tfig = plt.figure(figsize=({img_width},{img_height}))
            \\tax1 = fig.add_subplot(111)
            \\tax1.plot([row[0] for row in dat], [row[1] for row in dat], '-o',color='k',mec='k',markersize=3)
            \\tax1.set_ylabel('Number of sites',color='k')
            \\tax1.set_yscale('log')
            \\t#ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            \\tfor tl in ax1.get_yticklabels(): tl.set_color('k')
            \\tax1.set_xlabel('Non-ref allele frequency')
            \\tax2 = ax1.twinx()
            \\tax2.plot([row[0] for row in dat], [row[2] for row in dat], '-o',color='{id2col_id}',mec='{id2col_id}',markersize=3)
            \\tax2.set_ylabel('Ts/Tv',color='{id2col_id}')
            \\tax2.set_ylim(0,0.5+max(3,max(row[2] for row in dat)))
            \\tax1.set_xlim(0,1)
            \\tfor tl in ax2.get_yticklabels(): tl.set_color('{id2col_id}')
            \\tplt.subplots_adjust(right=0.88,left=0.15,bottom=0.11)
            \\tplt.title('{title_id}')
            \\tplt.savefig('{img}.png')
            \\tif pdf_plots: plt.savefig('{img}.pdf')
            \\tplt.close()

    '''.format(img=img, img_width=opts['img_width'], img_height=opts['img_height'],
               id2col_id=opts['id2col'][i_d], title_id=opts['title'][i_d]))


def plot_tstv_by_QUAL(opts, i_d):
    """"""
    vals = get_values(opts, i_d, 'QUAL')
    if not vals:
        return
    fh = opts['plt_fh']
    img = 'tstv_by_qual.{}'.format(i_d)
    tfh = open('{}.dat'.format(img), 'w')
    tfh.write('# [1]Quality\t[2]Number of sites\t[3]Marginal Ts/Tv\n')

    dat = []
    ntot = 0
    for val in vals:
        dat.append([val[0], val[2], val[3]])  # qual, nts, ntv
        ntot += val[2] + val[3]
    sdat = sorted(dat, reverse=True)
    sdat.append([-1])
    dn = ntot * 0.05
    qprev = sdat[0][0]
    nts = 0
    ntv = 0
    nout = 0
    for rec in sdat:
        if rec[0] == -1 or nts + ntv > dn:
            if ntv:
                tfh.write('{}\t{:d}\t{:f}\n'.format(qprev, nts + ntv + nout, nts / ntv))
            if rec[0] == -1:
                break
            nout += nts + ntv
            nts = 0
            ntv = 0
            qprev = rec[0]
        nts += rec[1]
        ntv += rec[2]
    tfh.close()

    tprint(fh, '''

            dat = []
            with open('{img}.dat', 'r') as f:
            \\treader = csv.reader(f, 'tab')
            \\tfor row in reader:
            \\t\\tif row[0][0] != '#': dat.append([float(x) for x in row])

            if plot_tstv_by_qual and len(dat)>2:
            \\tfig = plt.figure(figsize=({img_width},{img_height}))
            \\tax1 = fig.add_subplot(111)
            \\tax1.plot([row[1] for row in dat], [row[2] for row in dat], '^-', ms=3, mec='{id2col_id}', color='{id2col_id}')
            \\tax1.set_ylabel('Ts/Tv',fontsize=10)
            \\tax1.set_xlabel('Number of sites\\n(sorted by QUAL, descending)',fontsize=10)
            \\tax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='x')
            \\tax1.set_ylim(min(2,min(row[2] for row in dat))-0.3,0.3+max(2.2,max(row[2] for row in dat)))

            \\tplt.subplots_adjust(right=0.88,left=0.15,bottom=0.15)
            \\tplt.title('{title_id}')
            \\tplt.savefig('{img}.png')
            \\tif pdf_plots: plt.savefig('{img}.pdf')
            \\tplt.close()

    '''.format(img=img, img_width=opts['img_width'], img_height=opts['img_height'],
               id2col_id=opts['id2col'][i_d], title_id=opts['title'][i_d]))


def plot_indel_distribution(opts, i_d):
    """"""
    vals = get_values(opts, i_d, 'IDD')
    if not vals:
        return
    # Set xlim to show 99 of indels but ignore outliers
    tmp = []
    for i_d in file_ids(opts):
        v = get_values(opts, i_d, 'IDD')
        for sv in v:
            while len(tmp) < abs(sv[0]) + 1:
                tmp.append(0)
            tmp[abs(sv[0])] += sv[1]
    n = 0
    for t in tmp:
        n += t if t else 0
    mysum = xlim = 0
    for xlim in range(len(tmp)):
        mysum += tmp[xlim] if tmp[xlim] else 0
        if mysum / n >= 0.99:
            break
    if xlim < 20:
        xlim = 20
    fh = opts['plt_fh']
    img = 'indels.{}'.format(i_d)

    tfh = open('{}.dat'.format(img), 'w')
    tfh.write('# [1]Indel length\t[2]Count\n')
    for val in vals:
        tfh.write('{}\t{}\n'.format(val[0], val[1]))
    tfh.close()

    tprint(fh, '''

            dat = []
            with open('{img}.dat', 'r') as f:
            \\treader = csv.reader(f, 'tab')
            \\tfor row in reader:
            \\t\\tif row[0][0] != '#': dat.append([float(x) for x in row])

            if plot_indel_dist and len(dat)>0:
            \\tfig = plt.figure(figsize=({img_width},{img_height}))
            \\tax1 = fig.add_subplot(111)
            \\tax1.bar([row[0]-0.5 for row in dat], [row[1] for row in dat], color='{id2col_0}')# , edgecolor='{id2col_0}')
            \\tax1.set_xlabel('InDel Length')
            \\tax1.set_ylabel('Count')
            \\tax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            \\tax1.set_xlim(-{xlim},{xlim})
            \\tplt.subplots_adjust(bottom=0.17)
            \\tplt.title('{title_id}')
            \\tplt.savefig('{img}.png')
            \\tif pdf_plots: plt.savefig('{img}.pdf')
            \\tplt.close()
    '''.format(img=img, img_width=opts['img_width'], img_height=opts['img_height'], id2col_0=opts['id2col'][0],
               xlim=xlim, title_id=opts['title'][i_d]))


def plot_substitutions(opts, i_d):
    """"""
    vals = get_values(opts, i_d, 'ST')
    if not vals:
        return
    fh = opts['plt_fh']
    img = 'substitutions.{}'.format(i_d)
    tprint(fh, '''
            dat = [
    ''')
    for i in range(len(vals)):
        val = vals[i]
        tprint(fh, '\\t[{i},\'{val_0}\',{val_1}],\n'.format(i=i, val_0=val[0], val_1=val[1]))
    tprint(fh, ''']

            if plot_substitutions:
            \\tfig = plt.figure(figsize=({img_width},{img_height}))
            \\tcm  = mpl.cm.get_cmap('autumn')
            \\tn = 12
            \\tcol = []
            \\tfor i in range(n): col.append(cm(1.*i/n))
            \\tax1 = fig.add_subplot(111)
            \\tax1.bar([row[0] for row in dat], [row[2] for row in dat], color=col)
            \\tax1.set_ylabel('Count')
            \\tax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            \\tax1.set_xlim(-0.5,n+0.5)
            \\tplt.xticks([row[0] for row in dat],[row[1] for row in dat],rotation=45)
            \\tplt.title('{title_id}')
            \\tplt.savefig('{img}.png')
            \\tif pdf_plots: plt.savefig('{img}.pdf')
            \\tplt.close()

    '''.format(img_width=opts['img_width'], img_height=opts['img_height'], title_id=opts['title'][i_d], img=img))


def plot_per_sample_stats(opts, i_d):
    vals = get_values(opts, i_d, 'PSC')
    if not vals:
        return
    fh = opts['plt_fh']
    img = 'tstv_by_sample.'.format(i_d)
    img2 = 'hets_by_sample.'.format(i_d)
    img3 = 'snps_by_sample.'.format(i_d)
    img4 = 'indels_by_sample.'.format(i_d)
    img5 = 'singletons_by_sample.'.format(i_d)
    img6 = 'dp_by_sample.'.format(i_d)

    tfh = open('{}.dat'.format(img), 'w')
    tfh.write('# [1]Sample ID\t[2]ts/tv\t[3]het/hom\t[4]nSNPs\t[5]nIndels\t[6]Average depth\t[7]nSingletons\t[8]Sample name\n')
    for i in range(len(vals)):
        tstv = vals[i][4] / vals[i][5] if vals[i][5] else 0
        hethom = vals[i][3] / vals[i][2] if vals[i][2] else 0
        tfh.write('{:d}\t{:f}\t{:f}\t{:d}\t{:d}\t{:f}\t{:d}\t{:s}\n'.format(
            i, tstv, hethom, vals[i][4] + vals[i][5], vals[i][6], vals[i][6], vals[i][7], vals[i][8], vals[i][0]))
    tfh.close()

    tprint(fh, '''

            dat = []
            with open('{img}.dat', 'r') as f:
            \\treader = csv.reader(f, 'tab')
            \\tfor row in reader:
            \\t\\tif row[0][0] != '#': dat.append(row)

            if plot_tstv_by_sample:
            \\tfig = plt.figure(figsize=(2*{img_width},{img_height}*0.7))
            \\tax1 = fig.add_subplot(111)
            \\tax1.plot([row[0] for row in dat], [row[1] for row in dat], 'o', color='{id2col_id}',mec='$$opts{id2col_id}')
            \\tax1.set_ylabel('Ts/Tv')
            \\tax1.set_ylim(min(float(row[1]) for row in dat)-0.1,max(float(row[1]) for row in dat)+0.1)
            \\tif sample_names:
            \\t\\t     plt.xticks([int(row[0]) for row in dat],[row[7] for row in dat],**sample_font)
            \\t\\t     plt.subplots_adjust(**sample_margins)
            \\telse:
            \\t\\t     plt.subplots_adjust(right=0.98,left=0.07,bottom=0.17)
            \\t\\t     ax1.set_xlabel('Sample ID')
            \\tplt.title('{title_id}')
            \\tplt.savefig('{img}.png')
            \\tif pdf_plots: plt.savefig('{img}.pdf')
            \\tplt.close()


            if plot_hethom_by_sample:
            \\tfig = plt.figure(figsize=(2*{img_width},{img_height}*0.7))
            \\tax1 = fig.add_subplot(111)
            \\tax1.plot([row[0] for row in dat], [row[2] for row in dat], 'o', color='{id2col_id}',mec='{id2col_id}')
            \\tax1.set_ylabel('nHet(RA) / nHom(AA)')
            \\tax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            \\tif sample_names:
            \\t\\t     plt.xticks([int(row[0]) for row in dat],[row[7] for row in dat],**sample_font)
            \\t\\t     plt.subplots_adjust(**sample_margins)
            \\telse:
            \\t\\t     plt.subplots_adjust(right=0.98,left=0.07,bottom=0.17)
            \\t\\t     ax1.set_xlabel('Sample ID')
            \\tplt.title('{title_id}')
            \\tplt.savefig('{img2}.png')
            \\tif pdf_plots: plt.savefig('{img2}.pdf')
            \\tplt.close()


            if plot_snps_by_sample:
            \\tfig = plt.figure(figsize=(2*{img_width},{img_height}*0.7))
            \\tax1 = fig.add_subplot(111)
            \\tax1.plot([row[0] for row in dat], [row[3] for row in dat], 'o', color='{id2col_id}',mec='{id2col_id}')
            \\tax1.set_ylabel('Number of SNPs')
            \\tax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            \\tif sample_names:
            \\t\\t     plt.xticks([int(row[0]) for row in dat],[row[7] for row in dat],**sample_font)
            \\t\\t     plt.subplots_adjust(**sample_margins)
            \\telse:
            \\t\\t     plt.subplots_adjust(right=0.98,left=0.07,bottom=0.17)
            \\t\\t     ax1.set_xlabel('Sample ID')
            \\tplt.title('{title_id}')
            \\tplt.savefig('{img3}.png')
            \\tif pdf_plots: plt.savefig('{img3}.pdf')
            \\tplt.close()


            if plot_indels_by_sample:
            \\tfig = plt.figure(figsize=(2*{img_width},{img_height}*0.7))
            \\tax1 = fig.add_subplot(111)
            \\tax1.plot([row[0] for row in dat], [row[4] for row in dat], 'o', color='{id2col_id}',mec='{id2col_id}')
            \\tax1.set_ylabel('Number of indels')
            \\tax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            \\tif sample_names:
            \\t\\t     plt.xticks([int(row[0]) for row in dat],[row[7] for row in dat],**sample_font)
            \\t\\t     plt.subplots_adjust(**sample_margins)
            \\telse:
            \\t\\t     plt.subplots_adjust(right=0.98,left=0.07,bottom=0.17)
            \\t\\t     ax1.set_xlabel('Sample ID')
            \\tplt.title('{title_id}')
            \\tplt.savefig('{img4}.png')
            \\tif pdf_plots: plt.savefig('{img4}.pdf')
            \\tplt.close()


            if plot_singletons_by_sample:
            \\tfig = plt.figure(figsize=(2*{img_width},{img_height}*0.7))
            \\tax1 = fig.add_subplot(111)
            \\tax1.plot([row[0] for row in dat], [row[6] for row in dat], 'o', color='{id2col_id}',mec='{id2col_id}')
            \\tax1.set_ylabel('Number of singletons')
            \\tax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            \\tif sample_names:
            \\t\\t     plt.xticks([int(row[0]) for row in dat],[row[7] for row in dat],**sample_font)
            \\t\\t     plt.subplots_adjust(**sample_margins)
            \\telse:
            \\t\\t     plt.subplots_adjust(right=0.98,left=0.07,bottom=0.17)
            \\t\\t     ax1.set_xlabel('Sample ID')
            \\tplt.title('{title_id}')
            \\tplt.savefig('{img5}.png')
            \\tif pdf_plots: plt.savefig('{img5}.pdf')
            \\tplt.close()


            if plot_depth_by_sample:
            \\tfig = plt.figure(figsize=(2*{img_width},{img_height}*0.7))
            \\tax1 = fig.add_subplot(111)
            \\tax1.plot([row[0] for row in dat], [row[5] for row in dat], 'o', color='{id2col_id}',mec='{id2col_id}')
            \\tax1.set_ylabel('Average depth')
            \\tax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            \\tif sample_names:
            \\t\\t     plt.xticks([int(row[0]) for row in dat],[row[7] for row in dat],**sample_font)
            \\t\\t     plt.subplots_adjust(**sample_margins)
            \\telse:
            \\t\\t     plt.subplots_adjust(right=0.98,left=0.07,bottom=0.17)
            \\t\\t     ax1.set_xlabel('Sample ID')
            \\tplt.title('{title_id}')
            \\tplt.savefig('{img6}.png')
            \\tif pdf_plots: plt.savefig('{img6}.pdf')
            \\tplt.close()

    '''.format(img=img, img2=img2, img3=img3, img4=img4, img5=img5, img6=img6, img_width=opts['img_width'],
               img_height=opts['img_height'], id2col_id=opts['id2col'][i_d], title_id=opts['title'][i_d]))


def plot_DP(opts, i_d):
    """"""
    vals = get_values(opts, i_d, 'DP')
    if not vals:
        return
    fh = opts['plt_fh']
    img = 'depth.{}'.format(i_d)

    tfh = open('{}.dat'.format(img), 'w')
    tfh.write('# [1]Depth\t[2]Cumulative number of genotypes\t[3]Number of genotypes\n')
    mysum = 0
    dp_regex = re.compile('^\d+$')
    for i in range(len(vals)):
        if mysum > 99:
            break
        if not re.search(dp_regex, str(vals[i][0])):  # DP ">500" case
            continue
        mysum += vals[i][2]
        tfh.write('{:d}\t{:f}\t{:f}\n'.format(vals[i][0], mysum, vals[i][2]))
    tfh.close()

    tprint(fh, '''

            dat = []
            with open('{img}.dat', 'r') as f:
            \\treader = csv.reader(f, 'tab')
            \\tfor row in reader:
            \\t\\tif row[0][0] != '#': dat.append(row)

            if plot_dp_dist:
            \\tfig = plt.figure(figsize=({img_width},{img_height}))
            \\tax1 = fig.add_subplot(111)
            \\tax1.plot([row[0] for row in dat], [row[2] for row in dat], '-^', color='k')
            \\tax1.set_ylabel('Number of genotypes [%]',color='k')
            \\tax1.set_xlabel('Depth')
            \\tax2 = ax1.twinx()
            \\tax2.plot([row[0] for row in dat], [row[1] for row in dat], '-o', color='{id2col_id}')
            \\tax2.set_ylabel('Cumulative number of genotypes [%]',color='{id2col_id}')
            \\tfor tl in ax2.get_yticklabels(): tl.set_color('{id2col_id}')
            \\tplt.subplots_adjust(left=0.15,bottom=0.15,right=0.87)
            \\tplt.title('{title_id}')
            \\tplt.savefig('{img}.png')
            \\tif pdf_plots: plt.savefig('{img}.pdf')
            \\tplt.close()

    '''.format(img=img, img_width=opts['img_width'],  img_height=opts['img_height'], id2col_id=opts['id2col'][i_d],
               title_id=opts['title'][i_d]))


def plot_hwe(opts, i_d):
    """"""
    vals = get_values(opts, i_d, 'HWE')
    if not vals:
        return
    fh = opts['plt_fh']
    img = 'hwe.{}'.format(i_d)

    tfh = open('{}.dat'.format(img), 'w')
    tfh.write('# [1]Allele Frequency\t[2]Depth\t[3]Number of hets (median)\t[4]Number of hets (25-75th percentile)\n')
    for i in range(len(vals)):
        if not vals[i][1]:
            continue
        tfh.write('\t'.join(vals[i]) + '\n')
    tfh.close()

    tprint(fh, '''


            dat = []
            with open('{img}.dat', 'r') as f:
            \\treader = csv.reader(f, 'tab')
            \\tfor row in reader:
            \\t\\tif row[0][0] != '#': dat.append(row)

            if plot_hwe and len(dat)>1:
            \\tx  = [float(row[0]) for row in dat]
            \\ty1 = smooth(numpy.array([float(row[2]) for row in dat]),40,'hanning')
            \\ty2 = smooth(numpy.array([float(row[3]) for row in dat]),40,'hanning')
            \\ty3 = smooth(numpy.array([float(row[4]) for row in dat]),40,'hanning')
            \\tdp = smooth(numpy.array([float(row[1]) for row in dat]),40,'hanning')
            \\thwe = []
            \\tfor af in x: hwe.append(2*af*(1-af))

            \\tfig = plt.figure(figsize=({img_width},{img_height}))
            \\tax1 = fig.add_subplot(111)
            \\tplots  = ax1.plot(x,hwe,'--',color='#ff9900',label='HWE')
            \\tplots += ax1.plot(x,y2,color='#ff9900',label='Median')
            \\tplots += ax1.plot(x,y3,color='#ffe0b2',label='25-75th percentile')
            \\tax1.fill_between(x,y1,y3, facecolor='#ffeacc',edgecolor='#ffe0b2')
            \\tax1.set_ylabel('Fraction of hets',color='#ff9900')
            \\tax1.set_xlabel('Allele frequency')
            \\tfor tl in ax1.get_yticklabels(): tl.set_color('#ff9900')
            \\tax2 = ax1.twinx()
            \\tplots += ax2.plot(x,dp, 'k', label='Number of sites')
            \\tax2.set_ylabel('Number of sites')
            \\tax2.set_yscale('log')
            \\tif af_xlog: ax1.set_xscale('log')
            \\tif af_xlog: ax2.set_xscale('log')
            \\tlabels = [l.get_label() for l in plots]
            \\tplt.legend(plots,labels,numpoints=1,markerscale=2,loc='center',prop={{'size':9}},frameon=False)
            \\tplt.subplots_adjust(left=0.15,bottom=0.15,right=0.86)
            \\tplt.title('{title_id}')
            \\tplt.savefig('{img}.png')
            \\tif pdf_plots: plt.savefig('{img}.pdf')
            \\tplt.close()

    '''.format(img=img, img_width=opts['img_width'], img_height=opts['img_height'], title_id=opts['title'][i_d]))


def plot(opts):
    """"""
    if not opts['plt_fh']:
        return
    opts['plt_fh'].close()
    cmd = 'python {}'.format(opts['plt_file'])
    logging.info('Plotting graphs: {}'.format(cmd))
    retcode = os.system(cmd)
    if retcode != 0:
        error('The command exited with non-zero status {}:\n\t{}\n\n'.format(retcode, cmd))


def create_pdf(opts):
    """"""
    ids = file_ids(opts)
    width = '25.4cm'  # todo: move all this to $$opts{tex}
    height = '19cm'
    height1 = '13cm'
    width1 = '23cm'
    width2 = '10.5cm' if ids == 3 else '10.5cm'
    fmt = 'png' if opts['rasterize'] else 'pdf'
    ext = 'type={fmt},ext=.{fmt},read=.{fmt}'.format(fmt=fmt)
    opts['fmt'] = fmt
    opts['ext'] = ext

    # Check that xcolor is available
    cmd = 'kpsewhich xcolor.sty'
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    out, err = proc.communicate()
    has_color = out.decode("utf-8").strip()
    if not out:
        logging.warning('Note: The xcolor.sty package not available, black and white tables only...\n\n')

    tex_file = 'summary.tex'
    pdf_file = 'summary.pdf'
    tex = open(tex_file, 'w')
    tprint(tex, '''
            % This file was produced by plot-vcfstats, the command line was:
            %   {args}
            %
            % Edit as necessary and recreate the PDF by running
            %   pdflatex {tex_file}
            %

            % Slides style and dimensions
            %
            \\nonstopmode
            \\documentclass[17pt]{{memoir}}
            \\setstocksize{{{height}}}{{{width}}}
            \\settrimmedsize{{\\stockheight}}{{\\stockwidth}}{{*}}
            \\settrims{{0pt}}{{0pt}}
            \\setlrmarginsandblock{{1cm}}{{*}}{{*}}
            \\setulmarginsandblock{{1.5cm}}{{*}}{{*}}
            \\setheadfoot{{1mm}}{{1cm}}
            \\setlength{{\\parskip}}{{0pt}}
            \\setheaderspaces{{*}}{{1mm}}{{*}}
            \\setmarginnotes{{1mm}}{{1mm}}{{1mm}}
            \\checkandfixthelayout[fixed]
            \\usepackage{{charter}}   % font
            \\pagestyle{{plain}}
            \\makeevenfoot{{plain}}{{}}{{}}{{\\thepage}}
            \\makeoddfoot{{plain}}{{}}{{}}{{\\thepage}}
            \\usepackage{{graphicx}}

            % For colored tables. If xcolor.sty is not available on your system,
            % download xcolor.sty.gz LaTeX class style from
            %   http://www.ukern.de/tex/xcolor.html
            % Unpack and install system-wide or place elsewhere and make available by
            % setting the TEXINPUTS environment variable (note the colon)
            %   export TEXINPUTS=some/dir:
            % The list of the recognised path can be obtained by running `kpsepath tex`
            %
            \\usepackage{{multirow}}
            \\setlength{{\\tabcolsep}}{{0.6em}}
            \\renewcommand{{\\arraystretch}}{{1.2}}
    '''.format(args=opts['args'], tex_file=tex_file, height=height, width=width))
    if has_color:
        tprint(tex, '\\usepackage[table]{xcolor}')
    else:
        tprint(tex, '''
            \\newcommand{\\definecolor}[3]{}
            \\newcommand{\\columncolor}[1]{}
            \\newcommand{\\rowcolors}[4]{}
            \\newcommand{\\arrayrulecolor}[1]{}
        ''')
    tprint(tex, '''
            \\definecolor{{hcol1}}{{rgb}}{{1,0.6,0}}
            \\definecolor{{hcol2}}{{rgb}}{{1,0.68,0.2}}
            \\definecolor{{row1}}{{rgb}}{{1,0.88,0.7}}
            \\definecolor{{row2}}{{rgb}}{{1,0.92,0.8}}    % #FFEBCC
            \\setlength{{\\arrayrulewidth}}{{1.5pt}}

            % Slide headings
            \\newcommand*{{\\head}}[1]{{{{\\Large\\centerline{{#1}}\\vskip0.5em}}}}

            % Slide definition
            \\newcommand*{{\\hslide}}[2]{{%
                    \\head{{#1}}%
                    \\begin{{vplace}}[0.5]\\centerline{{#2}}\\end{{vplace}}\\newpage}}
            \\newcommand{{\\pdf}}[2]{{\\IfFileExists{{#2.{fmt}}}{{\\includegraphics[#1]{{#2}}}}{{}}}}


            % The actual slides
            \\begin{{document}}
    '''.format(fmt=fmt))
    slide = '''
            \\begin{minipage}{\\textwidth}\centering
            \small \\rowcolors*{3}{row2}{row1} \\arrayrulecolor{black}
            \\begin{tabular}{l | r r r | r r | r | r}
            \multicolumn{1}{>{\columncolor{hcol1}}l|}{}
            & \multicolumn{3}{>{\columncolor{hcol1}}c|}{SNPs}
            & \multicolumn{2}{>{\columncolor{hcol1}}c|}{indels}
            & \multicolumn{1}{>{\columncolor{hcol1}}c|}{MNPs}
            & \multicolumn{1}{>{\columncolor{hcol1}}c}{others}  \\\\
            %
            \multicolumn{1}{>{\columncolor{hcol2}}l|}{Callset}
            & \multicolumn{1}{>{\columncolor{hcol2}}c}{n}
            & \multicolumn{1}{>{\columncolor{hcol2}}c }{ts/tv}
            & \multicolumn{1}{>{\columncolor{hcol2}}c|}{\\footnotesize(1st ALT)}
            & \multicolumn{1}{>{\columncolor{hcol2}}c}{n}
            & \multicolumn{1}{>{\columncolor{hcol2}}c}{frm$^*$}
            & \multicolumn{1}{>{\columncolor{hcol2}}c|}{}
            & \multicolumn{1}{>{\columncolor{hcol2}}c}{} \\\\ \hline
            '''
    tex_titles = {}
    for i_d in ids:
        snps = get_value(opts, i_d, 'number of SNPs:')
        indels = get_value(opts, i_d, 'number of indels:')
        mnps = get_value(opts, i_d, 'number of MNPs:')
        others = get_value(opts, i_d, 'number of others:')
        tstv = float('{:.2f}'.format(get_values(opts, i_d, 'TSTV', 0, 2)))
        tstv1 = float('{:.2f}'.format(get_values(opts, i_d, 'TSTV', 0, 5)))
        frsh = get_values(opts, i_d, 'FS')
        frshlen = frsh[0][3] if frsh else '--'
        title = opts['title'][i_d]
        title = re.sub('_', '\\_', title)
        title = re.sub('^\s*\*', '\$*\$', title)
        tex_titles[i_d] = title
        slide += ' {} & '.format(title) + bignum(snps) + ' & {} & {} & '.format(tstv, tstv1) + bignum(indels) +\
                 ' & {} & '.format(frshlen) + bignum(mnps) + ' & ' + bignum(others) + ' \\\\ \n'
    slide += '''%
    \multicolumn{8}{r}{$^*$ frameshift ratio: out/(out+in)} \\\\
    \end{tabular}
    \\\\ \\vspace{1em}
    \\begin{tabular}{l | r r r | r r}
    \multicolumn{1}{>{\columncolor{hcol1}}l|}{}
    & \multicolumn{3}{>{\columncolor{hcol1}}c|}{singletons {\\footnotesize(AC=1)}}
    & \multicolumn{2}{>{\columncolor{hcol1}}c}{multiallelic}  \\\\
    %
    \multicolumn{1}{>{\columncolor{hcol2}}l|}{Callset}
    & \multicolumn{1}{>{\columncolor{hcol2}}c}{SNPs}
    & \multicolumn{1}{>{\columncolor{hcol2}}c}{ts/tv}
    & \multicolumn{1}{>{\columncolor{hcol2}}c}{indels}
    & \multicolumn{1}{>{\columncolor{hcol2}}c}{sites}
    & \multicolumn{1}{>{\columncolor{hcol2}}c}{SNPs} \\\\ \hline
    '''
    for i_d in ids:
        s = singletons(opts, i_d)
        mals = get_value(opts, i_d, 'number of multiallelic sites:')
        msnps = get_value(opts, i_d, 'number of multiallelic SNP sites:')
        title = tex_titles[i_d]
        slide += ' {} & {}\\% & {} & {}\\% & '.format(title, s['snps'], s['tstv'], s['indels']) + bignum(mals) +\
                 ' & ' + bignum(msnps) + ' \\\\ \n'
    slide += ''' \\end{tabular}
        \\vspace{2em}
        \\begin{itemize}[-]
        \\setlength{\\itemsep}{0pt}
    '''
    for i_d in ids:
        fname = opts['dat']['ID'][i_d][0][0]
        if re.search(' \+ ', opts['title'][i_d]):
            continue
        fname = re.sub('.{80}', '$&\\\\\\hskip2em ', fname)
        fname = re.sub('_', '\\_', fname)
        slide += '\\item {} .. \\texttt{{\\footnotesize {}}}\n'.format(tex_titles[i_d], fname)
    slide += '\\end{itemize}\\end{minipage}'
    title = opts['main_title'] if 'main_title' in opts else 'Summary Numbers'
    tprint(tex, '''
        % Table with summary numbers
        %
        \\hslide{{{title}}}{{{slide}}}

    '''.format(title=title, slide=slide))

    # Venn bars
    if len(ids) == 3:
        tprint(tex, '''

            % Venn numbers
            %
            \\hslide{{Total counts}}{{%
                \\includegraphics[{ext},width={width2}]{{venn_bars.snps}}%
                \\includegraphics[{ext},width={width2}]{{venn_bars.indels}}
            }}
        '''.format(ext=ext, width2=width2))

    tprint(tex, fmt_slide3v(opts, 'tstv_by_sample', 'Ts/Tv by sample'))
    tprint(tex, fmt_slide3v(opts, 'hets_by_sample', 'Hets vs non-ref Homs by sample'))
    tprint(tex, fmt_slide3v(opts, 'singletons_by_sample', 'Singletons by sample {\normalsize(hets and homs)}'))
    tprint(tex, fmt_slide3v(opts, 'dp_by_sample', 'Average depth by sample'))
    tprint(tex, fmt_slide3v(opts, 'snps_by_sample', 'Number of SNPs by sample'))
    tprint(tex, fmt_slide3v(opts, 'indels_by_sample', 'Number of indels by sample'))
    if get_values(opts, 2, 'GCsS'):
        tprint(tex, '''
            % Genotype discordance by sample
            %
            \\hslide{{Genotype discordance by sample}}{{\\pdf{{{ext},width={width1}}}{{gts_by_sample}}}}

        '''.format(ext=ext, width1=width1))
    if get_values(opts, 2, 'GCsAF'):
        vals = get_values(opts, 2, 'NRDs')
        nrd = float('{:.2f}'.format(vals[0][0]))
        rr = float('{:.2f}'.format(vals[0][1]))
        ra = float('{:.2f}'.format(vals[0][2]))
        aa = float('{:.2f}'.format(vals[0][3]))
        table = '''%
            {{\\small
            \\rowcolors*{{1}}{{row2}}{{row1}}\\arrayrulecolor{{black}}
            \\begin{{tabular}}{{c | c | c | c }}
            \\multicolumn{{1}}{{>{{\\columncolor{{hcol1}}}}c|}}{{REF/REF}} &
            \\multicolumn{{1}}{{>{{\\columncolor{{hcol1}}}}c|}}{{REF/ALT}} &
            \\multicolumn{{1}}{{>{{\\columncolor{{hcol1}}}}c|}}{{ALT/ALT}} &
            \\multicolumn{{1}}{{>{{\\columncolor{{hcol1}}}}c}}{{NRDs}} \\\\ \\hline
            {rr}\\% & {ra}\\% & {aa}\\% & {nrd}\\% \\\\
            \\end{{tabular}}}}]'''.format(rr=rr, ra=ra, aa=aa, nrd=nrd)
        tprint(tex, '''
                % Genotype discordance by AF
                %
                \\head{{Genotype discordance by AF}}\\begin{{vplace}}[0.7]\\centerline{{{table}}}%
                \\centerline{{\\pdf{{{ext},height={height1}}}{{gts_by_af}}}}\\end{{vplace}}
                \\newpage

                % dosage r2 by AF
                %
                \\hslide{{Allelic R\$^2\$ by AF}}{{\\pdf{{{ext},height={height1}}}{{r2_by_af}}}}
        '''.format(table=table, ext=ext, height1=height1))
    if os.path.isfile('counts_by_af.snps.{}'.format(fmt)) and os.path.isfile('counts_by_af.indels.{}'.format(fmt)):
        tprint(tex, '''
            % SNP and indel counts by AF
            %
            \\hslide{{}}{{\\vbox{{\\noindent\\includegraphics[{ext},width={width1}]{{counts_by_af.snps}}\\\\%
                \\noindent\\includegraphics[{ext},width={width1}]{{counts_by_af.indels}}}}
            '''.format(ext=ext, width1=width1))
    if os.path.isfile('overlaps_by_af.snps.{}'.format(fmt)) and os.path.isfile('overlaps_by_af.indels.{}'.format(fmt)):
        tprint(tex, '''
            % SNP and indel overlap by AF
            %
            \\hslide{{}}{{\\vbox{{\\noindent\\includegraphics[{ext},width={width1}]{{overlap_by_af.snps}}\\\\%
                \\noindent\\includegraphics[{ext},width={width1}]{{overlap_by_af.indels}}}}}}
            '''.format(ext=ext, width1=width1))
    tprint(tex, fmt_slide3h(opts, 'tstv_by_af', 'Ts/Tv by AF'))
    tprint(tex, fmt_slide3h(opts, 'tstv_by_qual', 'Ts/Tv stratified by QUAL'))
    tprint(tex, fmt_slide3h(opts, 'indels', 'Indel distribution'))
    tprint(tex, fmt_slide3h(opts, 'depth', 'Depth distribution'))
    tprint(tex, fmt_slide3h(opts, 'hwe', 'Number of HETs by AF'))
    tprint(tex, fmt_slide3h(opts, 'substitutions', 'Substitution types'))
    # tprint(tex, fmt_slide3h(opts, 'irc_by_af', 'Indel Repeat Consistency by AF'))
    # tprint(tex, fmt_slide3h(opts, 'irc_by_rlen', 'Indel Consistency by Repeat Type'))

    tprint(tex, '\n\n\\end{document}\n')
    tex.close()

    tex_file = re.sub('^.+/', '', tex_file)
    cmd = 'pdflatex {tex_file} >{opts_logfile} 2>&1'.format(tex_file=tex_file, opts_logfile=opts['logfile'])
    logging.info('Creating PDF: {}\n'.format(cmd))
    retcode = os.system(cmd)
    if retcode != 0:
        error('The command exited with non-zero status, please consult the output of pdflatex: {}{}\n\n'.format(
            opts['dir'], opts['logfile']))
    logging.info('Finished: {}/{}'.format(opts['dir'], pdf_file))


def singletons(opts, i_d):
    """"""
    si_vals = get_values(opts, i_d, 'SiS')
    si_snps = si_vals[0][1]
    si_indels = si_vals[0][4]
    si_irc = float('{:.3f}'.format(si_vals[0][5] / (si_vals[0][5] + si_vals[0][6]) if si_vals[0][6] else 0))
    si_tstv = float('{:.2f}'.format(si_vals[0][2] / si_vals[0][3] if si_vals[0][3] else 0))
    all_vals = get_values(opts, i_d, 'AF')
    nsnps = 0
    nindels = 0
    for val in all_vals:
        nsnps += val[1]
        nindels += val[4]
    si_snps = '{:.1f}'.format(si_snps * 100 / nsnps if nsnps else 0)
    si_indels = '{:.1f}'.format(si_indels * 100 / nindels if nindels else 0)
    return {'snps': si_snps, 'indels': si_indels, 'tstv': si_tstv, 'irc': si_irc}


def fmt_slide3v(opts, image, title):
    """"""
    n = 0
    for i_d in range(2):
        if os.path.isfile('{}.{}.{}'.format(image, i_d, opts['fmt'])):
            n += 1
    if not n:
        return ''
    h = opts['tex']['slide3v']['height{}'.format(n)]
    slide = '\\vbox{'
    for i_d in range(2):
        if not os.path.isfile('{}.{}.{}'.format(image, i_d, opts['fmt'])):
            continue
        slide += '\\centerline{{\\includegraphics[{opts_ext},height={h}]{{{image}.{i_d}}}}'.format(
            opts_ext=opts['ext'], h=h, image=image, i_d=i_d)
    return '''
            % {title}
            %
            \\hslide{{{title}}}{{{slide}}}
        '''.format(title=title, slide=slide)


def fmt_slide3h(opts, image, title):
    """"""
    n = 0
    for i_d in range(2):
        if os.path.isfile('{}.{}.{}'.format(image, i_d, opts['fmt'])):
            n += 1
    if not n:
        return ''
    w = opts['tex']['slide3h']['width{}'.format(n)]
    slide = ''
    for i_d in range(2):
        if not os.path.isfile('{}.{}.{}'.format(image, i_d, opts['fmt'])):
            continue
        slide += '\\includegraphics[{opts_ext},width={w}]{{{image}.{i_d}}}'.format(
            opts_ext=opts['ext'], w=w, image=image, i_d=i_d)
    return '''
            % {title}
            %
            \\hslide{{{title}}}{{{slide}}}
        '''.format(title=title, slide=slide)


def tabreplace(matchobj):
    """"""
    return matchobj.group(1) + '\t'


def tprint(fh, txt):
    """"""
    txt = re.sub('^[ \t]+', '', txt, flags=re.MULTILINE)  # eat leading tabs
    while re.search('^\t*\\\\t', txt, flags=re.MULTILINE):
        txt = re.sub('^(\t*)\\\\t', '\g<1>\t', txt, flags=re.MULTILINE)  # replace escaped tabs (\\t) with tabs
    fh.write(txt)


def get_value(opts, i_d, key):
    """"""
    if i_d not in opts['dat']:
        return
    if key not in opts['dat'][i_d]:
        return
    return opts['dat'][i_d][key]


def get_values(opts, i_d, key, i=None, j=None):
    """"""
    if key not in opts['dat']:
        return
    if i_d not in opts['dat'][key]:
        return
    if key not in opts['dat']:
        error('todo: sanity check for {}'.format(key))
    if key in opts['def_line'] and opts['def_line'][key] != opts['exp'][key] and not opts['def_line_warned'][key]:
        logging.warning(
            'Warning: Possible version mismatch, the definition line differs\n'
            '\texpected: {expected}\n'
            '\tfound:    {found}'.format(expected=opts['exp'][key], found=opts['def_line'][key]))
        opts['def_line_warned'][key] = 1
    if i is not None:
        if j is not None:
            return opts['dat'][key][i_d][i][j]
        return opts['dat'][key][i_d][i]
    return opts['dat'][key][i_d]


def bignum(num):
    """"""
    num = str(num)
    if not num:
        return ''
    num_re = re.compile('^\d+$')
    if not re.search(num_re, num):
        return int(num)
    length = len(num)
    out = ''
    for i in range(length):
        out += num[i]
        if i + 1 < length and not (length - i - 1) % 3:
            out += ','
    return out


def rebin_values(vals, bin_size, col):
    """"""
    avg = {}
    prev = vals[0][col]
    iout = 0
    nsum = 0
    out = []
    dat = []
    for i in range(len(vals)):
        for icol in range(len(vals[i])):
            if icol == col:
                continue
            while len(dat) < icol + 1:
                dat.append(0)
            dat[icol] += vals[i][icol]
        nsum += 1
        if i + 1 < len(vals) and vals[i][col] - prev < bin_size:
            continue
        dat[col] = prev
        for icol in range(len(vals[i])):
            while len(out) < iout + 1:
                out.append(None)
            if not out[iout]:
                out[iout] = [None] * len(vals[i])
            out[iout][icol] = dat[icol] if dat[icol] else 0
            if icol in avg and nsum:
                out[iout][icol] /= nsum
        nsum = 0
        dat = []
        iout += 1
        prev = vals[i][col]
    return out


def error(msg):
    logging.error(msg)
    raise Exception(msg)


def main():
    """"""
    parser = argparse.ArgumentParser(description='Work with stats generated by bcftools stats')
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='plot-vcfstats v{}'.format(__version__))
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '-m', '--merge',
        help='Merge outputs from a single VCF that was broken up',
        action='store_true')
    group.add_argument(
        '-p', '--prefix',
        metavar='<dir>',
        help='Output directory')
    parser.add_argument(
        '-T', '--main-title',
        help='Main title for the PDF',
        metavar='<title>')
    parser.add_argument(
        '-t', '--title',
        help='Identify files by these titles in plots (can be given multiple times)',
        metavar='<title>',
        action='append')
    parser.add_argument(
        '-s', '--sample-names',
        help=' Use sample names for xticks rather than numeric IDs',
        action='store_true')
    parser.add_argument(
        '-P', '--no-PDF',
        help='Skip the PDF creation step',
        action='store_true')
    parser.add_argument(
        '-r', '--rasterize',
        help='Rasterize PDF images for fast rendering',
        action='store_true')
    parser.add_argument(
        'stats_files',
        help='File(s) containing output from bcftools stats',
        metavar='<file>',
        type=argparse.FileType('r'),
        nargs='+')
    args = parser.parse_args()
    opts = parse_params(args)
    parse_vcfstats(opts)
    if len(opts['vcfstats']) > 1:
        merge_vcfstats(opts)
    os.chdir(opts['dir'])
    if opts['make_plots']:
        init_plots(opts)
        plot_venn_bars(opts)
        plot_counts_by_AF(opts)
        plot_overlap_by_AF(opts)
        plot_concordance_by_AF(opts)
        plot_concordance_by_sample(opts)
        for i_d in file_ids(opts):
            plot_tstv_by_AF(opts, i_d)
            plot_tstv_by_QUAL(opts, i_d)
            plot_indel_distribution(opts, i_d)
            plot_substitutions(opts, i_d)
            plot_per_sample_stats(opts, i_d)
            plot_DP(opts, i_d)
            plot_hwe(opts, i_d)
        plot(opts)
    if opts['make_pdf']:
        create_pdf(opts)


if __name__ == '__main__':
    main()
