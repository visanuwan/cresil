#!/usr/bin/env python
import mappy
import datetime
import pathlib
import subprocess
import intervaltree as tree
import argparse, re, sys
import pandas as pd
import numpy as np
from os import path
from tqdm import tqdm
from itertools import chain
from collections import Counter, namedtuple
from multiprocessing import cpu_count, Pool
from Bio import SeqIO
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def cal_feature_region(value, ini_start):
    str_ = (value['acc_length'] + ini_start) - value['length']
    end_ = value['acc_length'] + ini_start
    return_str = "{}-{}".format(str_, end_)
    return return_str

def getNonOverlap(ref, name, seq, mapq, list_ex_chr):
    refDict = {}
    mergedRead = {}
    i=0
    nonOvl = []
    fst = True
    for hit in ref.map(seq): # alignments
        if hit.is_primary and hit.mapq >= mapq and not hit.ctg in list_ex_chr:
            if fst:
                q_st_fst = hit.q_st
                q_en_fst = hit.q_en
                r_st_fst = hit.r_st
                r_en_fst = hit.r_en
                fst = False
                refDict.setdefault(hit.ctg, tree.IntervalTree()).add(tree.Interval(hit.r_st, hit.r_en))
                nonOvl.append([name,hit.q_st,hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq, hit.strand])
            else:
                if hit.ctg in refDict:
                    checkOvl = refDict[hit.ctg].overlaps(hit.r_st, hit.r_en) 
                else:
                    checkOvl = False
                if checkOvl:
                    for region in refDict[hit.ctg][hit.r_st:hit.r_en]:
                        st = region.begin
                        en = region.end
                        sdata = sorted([hit.r_st, hit.r_en, st, en])
                        ovl = sdata[2]-sdata[1]
                        lenR = en-st
                        lenQ = hit.r_en-hit.r_st
                        covR = ovl/float(lenR)
                        covQ = ovl/float(lenQ)
                        if covQ < 0.05:
                            refDict.setdefault(hit.ctg, tree.IntervalTree()).add(tree.Interval(hit.r_st, hit.r_en, i))
                            nonOvl.append([name,hit.q_st,hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq, hit.strand])
                else:
                    nonOvl.append([name,hit.q_st,hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq, hit.strand])
    return nonOvl

def getMerge(nonOvl_o, gap_len=10):

    min_start = min([x[1] for x in nonOvl_o])
    max_end = max([x[2] for x in nonOvl_o])

    this_nonOvl_o = []
    if len(nonOvl_o) > 1:
        for row in nonOvl_o:

            mod_row = list(row)
            mod_row[1] = max(min_start, mod_row[1] - gap_len)
            mod_row[2] = min(max_end, mod_row[2] + gap_len)

            this_nonOvl_o.append(mod_row)
    else:
        this_nonOvl_o = nonOvl_o
            
    newNonOvl = []
    mergedRead = tree.IntervalTree()
    for row in this_nonOvl_o: 
        mergedRead.add(tree.Interval(row[1], row[2]))

    mergedRead.merge_overlaps()
    
    read_mapped_pos = mergedRead[nonOvl_o[0][1]:nonOvl_o[0][2]]
    read_mapped_pos_begin = list(read_mapped_pos)[0].begin
    read_mapped_pos_end = list(read_mapped_pos)[0].end

    return_start = read_mapped_pos_begin + gap_len if read_mapped_pos_begin != min_start else read_mapped_pos_begin
    return_end = read_mapped_pos_end - gap_len if read_mapped_pos_end != max_end else read_mapped_pos_end

    return_read_mapped_pos = tree.Interval(return_start, return_end)
    
    newNonOvl = []
    for row in this_nonOvl_o:
        if row[1] + gap_len == return_read_mapped_pos.begin:
            row[1] += gap_len
        if row[2] - gap_len == return_read_mapped_pos.end:
            row[2] -= gap_len
        if return_read_mapped_pos.overlaps(row[1], row[2]):
            newNonOvl.append(row)

    return(return_read_mapped_pos, newNonOvl)

def readFaidx(fastxFile):
    if path.isfile("{}.fai".format(fastxFile)):
        count=0
        with open("{}.fai".format(fastxFile)) as inf:
            for line in inf: count += 1
        return count
    else:
        exit("please index fasta/fastq using samtools faidx/fqidx")

class pdRegion(object):
    def __init__(self, pdObject):
        self.readid = pdObject.readid
        self.q_len = pdObject.q_len
        self.q_start = pdObject.q_start
        self.q_end = pdObject.q_end
        self.ref = pdObject.ref
        self.r_start = pdObject.r_start
        self.r_end = pdObject.r_end
        self.matchLen = pdObject.matchLen
        self.blockLen = pdObject.blockLen
        self.mapq = pdObject.mapq
        self.strand = pdObject.strand
    def cal_q_len(self):
        return self.q_end - self.q_start
    def cal_r_len(self):
        return self.r_end - self.r_start
    def get_list_attr(self):
        return [self.readid, self.q_len, self.q_start, self.q_end, self.ref,
               self.r_start, self.r_end, self.matchLen, self.blockLen, self.mapq, self.strand]

class Region(object):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
    def get_str_region(self):
        return "{}_{}_{}".format(self.chrom, self.start, self.end)

def combine_region_strand(value):
    return "{}_{}_{}_{}".format(value['ref'], value['r_start'], value['r_end'], value['strand'])

def make_mappedRegion(value):
    mappedRegion = "{}_{}_{}".format(value['ref'], value['r_start'], value['r_end'])
    return mappedRegion
    
def is_overlapping(obj, offset, chrom, start, end):
    list_bool = [False, False]
    if obj.chrom == chrom:
        if abs(obj.start - start) <= offset:
            list_bool[0] = True
        if abs(obj.end - end) <= offset:
            list_bool[1] = True
    return all(list_bool)

def any_overlapping_range(start1, end1, start2, end2):
    if start1 <= end2 and start2 <= end1:
        return True
    else:
        return False

def split_list_step_window(list_check, step, window):
    return [list_check[i:i + window] for i in range(0, len(list_check), window - step)]


def make_group_order(list_test, offset):
    
    list_group_order = []
    list_group = []

    for region_strand in list_test:

        chrom, start, end, strand = region_strand.rsplit('_', 3)

        if len(list_group) == 0:

            obj_region = Region(chrom, int(start), int(end))
            list_group.append(obj_region)

            mod_region = "{}_{}_{}_{}".format(chrom, start, end, strand)
            list_group_order.append(mod_region)

        else:

            found_flag = 0
            for obj in list_group:
                if is_overlapping(obj, offset, chrom, int(start), int(end)):

                    found_flag = 1

                    mod_region = "{}_{}".format(obj.get_str_region(), strand)
                    list_group_order.append(mod_region)
                    break

            if found_flag == 0:
                obj_region = Region(chrom, int(start), int(end))
                list_group.append(obj_region)

                mod_region = "{}_{}_{}_{}".format(chrom, start, end, strand)
                list_group_order.append(mod_region)
                
    return list_group_order


def get_sum_pattern(list_group_order):
    
    list_sum_pattern = []

    list_most_common = Counter(list_group_order).most_common()
    if list_most_common[0][1] > 1:
        if len(list_most_common) == 1:
            list_sum_pattern = [True, list_most_common[0][0], list_most_common[0][1]]
        else:
            idx_start = list_group_order.index(list_most_common[0][0])
            list_trim = list_group_order[idx_start:]

            list_trim_most_common = Counter(list_trim).most_common()
            if len(list_trim_most_common) == 1:
                list_sum_pattern = [True, list_trim_most_common[0][0], list_trim_most_common[0][1]]
            else:
                list_sum_pattern = [False, '', 0]

                for i in range(1, len(list_trim)):
                    list_split = split_list_step_window(list_trim, 0, i)
                    list_region = [','.join(x) for x in list_split]

                    list_rev_split_most_common = Counter(list_region).most_common()[::-1]

                    for tup_ in list_rev_split_most_common:
                        if tup_[1] > 1:
                            list_sum_pattern = [True, tup_[0], tup_[1]]
    else:
        list_sum_pattern = [False, '', 0]
    
    return list_sum_pattern

def get_idx_longest_pattern(list_group_order, list_pattern_region):

    len_list_group_order = len(list_group_order)
    len_list_pattern_region = len(list_pattern_region)

    list_idx_pattern = []
    idx = 0
    while True:
        if idx + len_list_pattern_region > len_list_group_order:
            break
        idx_end = idx + len_list_pattern_region
        if list_group_order[idx:idx_end] == list_pattern_region:
            list_idx_pattern.append(list(range(idx,idx_end)))
            idx = idx_end
        else:
            idx += 1

    list_merge_idx_pattern = []

    for list_ in list_idx_pattern:

        if len(list_merge_idx_pattern) > 0:

            if list_merge_idx_pattern[-1][-1] + 1 == list_[0]:
                list_merge_idx_pattern[-1] += list_
            else:
                list_merge_idx_pattern.append(list_)
        else:
            list_merge_idx_pattern.append(list_)

    return max(list_merge_idx_pattern, key=len)


def check_read_pattern(ref, name, seq, mapq, list_ex_chr):

    list_result = []
    len_seq = len(seq)

    list_hit = []
    for hit in ref.map(seq):
        if hit.is_primary and hit.mapq >= mapq and not hit.ctg in list_ex_chr:
            list_local_hit = [name, hit.is_primary, hit.mapq, len_seq, hit.q_st, hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.strand]
            list_hit.append(list_local_hit)

    if len(list_hit) > 1:
        header = ['readid', 'is_primary', 'mapq', 'q_len', 'q_start', 'q_end', 'ref', 'r_start', 'r_end', 'matchLen', 'blockLen','strand']
        df_mapping = pd.DataFrame(list_hit, columns=header)

        df_mapping['mappedRegion'] = df_mapping.apply(combine_region_strand, axis=1)
        df_mapping = df_mapping.sort_values(by='q_start')
        df_mapping.reset_index(drop=True, inplace=True)

        list_check_pattern = df_mapping['mappedRegion'].tolist()
        list_group_order = make_group_order(list_check_pattern, 50)
        list_sum_pattern = get_sum_pattern(list_group_order)

        list_idx_pattern = []

        if list_sum_pattern[0]:
            
            list_pattern_region = list_sum_pattern[1].split(',')
            list_idx_pattern = get_idx_longest_pattern(list_group_order, list_pattern_region)

            if len(list_idx_pattern) > len(list_pattern_region):
                
                selected_cols = ['readid', 'q_len', 'q_start', 'q_end', 'ref', 'r_start', 'r_end', 'matchLen', 'blockLen', 'mapq', 'strand']
                df_mod_mapping = df_mapping[df_mapping.index.isin(list_idx_pattern)][selected_cols].copy()
                df_mod_mapping['qlenTrimmed'] = df_mod_mapping['q_end'].tolist()[-1] - df_mod_mapping['q_start'].tolist()[0]

                list_result = df_mod_mapping.values.tolist()
            
    return list_result


def getMergeAll(ref, name, seq, mapq, list_ex_chr, allow_gap, allow_overlap):

    list_merged_result = []
    len_seq = len(seq)

    list_hit = []
    for hit in ref.map(seq):
        if hit.is_primary and hit.mapq >= mapq and not hit.ctg in list_ex_chr:
            list_local_hit = [name, len_seq, hit.q_st, hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq, hit.strand]
            list_hit.append(list_local_hit)

    if len(list_hit) > 0:
        header = ['readid', 'q_len', 'q_start', 'q_end', 'ref', 'r_start', 'r_end', 'matchLen', 'blockLen', 'mapq','strand']
        df_mapping = pd.DataFrame(list_hit, columns=header)

        df_mapping = df_mapping.sort_values(by='q_start')
        df_mapping.reset_index(drop=True, inplace=True)

        if len(df_mapping) > 1:

            #### make initial hit objects

            list_init_obj = []
            for index, value in df_mapping.iterrows():
                list_init_obj.append(pdRegion(value))

            #### merge overlapping regions

            list_result_obj = [list_init_obj[0]]

            for idx, obj_ in enumerate(list_init_obj[1:], start=1):
                
                #### deletion compensation

                if any_overlapping_range(list_result_obj[-1].q_start, list_result_obj[-1].q_end, obj_.q_start, obj_.q_end):
                    
                    if abs(list_result_obj[-1].q_end - obj_.q_start) >= allow_overlap:

                        if (list_result_obj[-1].ref == obj_.ref) and (list_result_obj[-1].strand == obj_.strand):

                            if abs(list_result_obj[-1].r_end - obj_.r_start) <= 1000:

                                series_new = pd.Series([obj_.readid, obj_.q_len, min(list_result_obj[-1].q_start, obj_.q_start), max(list_result_obj[-1].q_end, obj_.q_end),
                                                        obj_.ref, min(list_result_obj[-1].r_start, obj_.r_start), max(list_result_obj[-1].r_end, obj_.r_end),
                                                        max(list_result_obj[-1].matchLen, obj_.matchLen), max(list_result_obj[-1].blockLen, obj_.blockLen),
                                                        max(list_result_obj[-1].mapq, obj_.mapq), obj_.strand], index = df_mapping.columns)
                                obj_new = pdRegion(series_new)
                                list_result_obj[-1] = obj_new
                            else:
                                list_result_obj.append(obj_)
                        else:
                            if obj_.cal_q_len() > list_result_obj[-1].cal_q_len():
                                list_result_obj[-1] = obj_
                    else:
                        list_result_obj.append(obj_)
                else:
                    
                    #### insertion compensation   
                    
                    if any_overlapping_range(list_result_obj[-1].q_start, list_result_obj[-1].q_end + allow_gap, obj_.q_start, obj_.q_end):
                        
                        if (list_result_obj[-1].ref == obj_.ref) and (list_result_obj[-1].strand == obj_.strand):

                            if abs(list_result_obj[-1].r_end - obj_.r_start) <= 1000:
                            
                                series_new = pd.Series([obj_.readid, obj_.q_len, min(list_result_obj[-1].q_start, obj_.q_start), max(list_result_obj[-1].q_end, obj_.q_end),
                                                        obj_.ref, min(list_result_obj[-1].r_start, obj_.r_start), max(list_result_obj[-1].r_end, obj_.r_end),
                                                        max(list_result_obj[-1].matchLen, obj_.matchLen), max(list_result_obj[-1].blockLen, obj_.blockLen),
                                                        max(list_result_obj[-1].mapq, obj_.mapq), obj_.strand], index = df_mapping.columns)
                                obj_new = pdRegion(series_new)
                                list_result_obj[-1] = obj_new
                            else:
                                list_result_obj.append(obj_)
                        else:
                            list_result_obj.append(obj_)

            list_merged_result = [obj.get_list_attr() for obj in list_result_obj]
        else:
            list_merged_result = df_mapping.values.tolist()
        
    return list_merged_result


def cal_pattern_trim(record):
    global ref
    global mapq
    global list_ex_chr
    global allow_gap
    global allow_overlap
    
    list_result = []
    
    name = record.id
    seq = str(record.seq)
    
    list_chk_pattern = check_read_pattern(ref, name, seq, mapq, list_ex_chr)
    if len(list_chk_pattern) > 0:
        for idx, rTab in enumerate(list_chk_pattern, start=0):
            oTab = rTab+["{:.2f}".format((rTab[3]-rTab[2])/rTab[11]), idx, True]
            list_result.append(oTab)
    else:
        list_merged_result = getMergeAll(ref, name, seq, mapq, list_ex_chr, allow_gap, allow_overlap)
        if len(list_merged_result) > 0:
            
            len_trim_read = list_merged_result[-1][3] - list_merged_result[0][2]

            for idx, rTab in enumerate(list_merged_result, start=0):
                oTab = rTab+[len_trim_read, "{:.2f}".format( (rTab[3] - rTab[2]) / len_trim_read ), idx, False]
                list_result.append(oTab)
                    
    return list_result


def argparser():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='Find and trim potential eccDNA regions from ONT reads')
    general = parser.add_argument_group(title='General options')
    general.add_argument('-fq', "--fq-input",  dest='fqinput',
                            help="input fasta/fastq",
                            type=str, default=None)
    general.add_argument('-mapq', "--map-quality", dest='mapq',  
                            help="mapping quality [30]",
                            type=int, default=30)
    general.add_argument('-t', "--threads", dest='threads',
                            help="Number of threads [all CPU cores]",
                            type=int, default=0)
    general.add_argument('-g', "--gap", dest='allow_gap', 
                            help="allowing gap length [10]",
                            type=int, default=10)
    general.add_argument('-l', "--overlap", dest='allow_overlap', 
                            help="allowing overlapping length [10]",
                            type=int, default=10)
    general.add_argument('-e', "--exclude", dest='ex_chr',  
                            help="exclude chromosome(s) separating with comma [None]",
                            type=str, default='')
    general.add_argument('-r', "--ref",  
                            help="reference Minimap2 index .mmi",
                            type=str, default=None)
    general.add_argument('-o', "--output",  
                            help="output directory",
                            type=str, default="cresil_result")
    return parser

def main(args):

    if not args.fqinput:
        sys.stderr.write("[ABORT] Fastq reads are needed\n")
        exit(1)

    global ref
    global mapq
    global list_ex_chr
    global allow_gap
    global allow_overlap

    threads = cpu_count() if args.threads == 0 else args.threads
    fname = args.fqinput
    fref = args.ref
    mapq = args.mapq
    allow_gap = args.allow_gap
    allow_overlap = args.allow_overlap
    bname, ext = path.splitext(fname)

    list_ex_chr = []
    if args.ex_chr != '':
        list_ex_chr = [x.strip() for x in args.ex_chr.split(',')]

    if ext in ['.fasta', '.fna', '.fa']:
        filetype = 'fasta'
    elif ext in ['.fastq', '.fq']:
        filetype = 'fastq'
    else:
        exit("please input type (fasta/fastq)")

    pathlib.Path(args.output).mkdir(parents=True, exist_ok=True)

    oname_map = "{}/trim.txt".format(args.output)

    list_record = []
    iter_ = 10000
    flag = 0

    print("\n######### CReSIL : start trimming process (thread : {})".format(threads), flush=True)

    ## check a fastq index
    fname_index = "{}.fai".format(fname)
    if not path.isfile(fname_index):

        ct = datetime.datetime.now()
        print("[{}] indexing an input file".format(ct), flush=True)

        cmd = "samtools fqidx {}".format(fname)
        process = subprocess.call("{}".format(cmd), shell=True)

        ct = datetime.datetime.now()
        print("[{}] finised indexing an input file".format(ct), flush=True)
    
    numSeq = readFaidx(fname)

    ct = datetime.datetime.now()
    print("[{}] total read : {}".format(ct, numSeq), flush=True)

    ## load Aligner ############################################################

    MM_F_NO_LJOIN = 0x400
    ref = mappy.Aligner(fref, preset="map-ont", extra_flags=MM_F_NO_LJOIN)  # load or build index
    if not ref:
        raise Exception("ERROR: failed to load/build index")

    ## trimming process ########################################################
    
    for counter, record in enumerate(SeqIO.parse(fname, filetype), start=1):
        
        if counter % iter_ == 0:
            ct = datetime.datetime.now()
            print("[{}] {}".format(ct, counter), flush=True)
        
        if counter % iter_ == 0:
            
            list_record.append(record)
            
            pool = Pool(processes = threads)
            list_pool_result = pool.map(cal_pattern_trim, list_record)
            pool.close()
            pool.join()
            
            list_result = list(chain.from_iterable(list_pool_result))
            
            list_cols = ['readid', 'q_len', 'q_start', 'q_end', 'ref', 'r_start', 'r_end',
                        'match', 'mapBlock', 'mapq', 'strand', 'qlenTrimmed', 'freqCov', 'order', 'has_repeat_pattern']
            df_trim = pd.DataFrame(list_result, columns=list_cols)
            
            if flag == 0:
                df_trim.to_csv(oname_map, sep='\t', mode='w', header=True, index=None)
                flag = 1
            else:
                df_trim.to_csv(oname_map, sep='\t', mode='a', header=False, index=None)
                
            list_record = []
        else:
            list_record.append(record)
                
    if len(list_record) > 0:

        pool = Pool(processes = threads)
        list_pool_result = pool.map(cal_pattern_trim, list_record)
        pool.close()
        pool.join()
        
        list_result = list(chain.from_iterable(list_pool_result))
        
        list_cols = ['readid', 'q_len', 'q_start', 'q_end', 'ref', 'r_start', 'r_end',
                    'match', 'mapBlock', 'mapq', 'strand', 'qlenTrimmed', 'freqCov', 'order', 'has_repeat_pattern']
        df_trim = pd.DataFrame(list_result, columns=list_cols)

        if path.exists(oname_map):
            df_trim.to_csv(oname_map, sep='\t', mode='a', header=False, index=None)
        else:
            df_trim.to_csv(oname_map, sep='\t', mode='w', header=True, index=None)

    ct = datetime.datetime.now()
    print("[{}] finished trimming process\n".format(ct), flush=True)

