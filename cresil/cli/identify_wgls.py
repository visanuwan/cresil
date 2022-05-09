#!/usr/bin/env python
import os, sys, signal
import math
import datetime
import pathlib
import subprocess
import operator
from turtle import right
from xmlrpc.client import FastParser
import pandas as pd
import numpy as np
import pysam, shutil
import pybedtools as bt
import graphviz as gv
import networkx as nx
import matplotlib.pyplot as plt
from itertools import chain
from multiprocessing import cpu_count, Pool
from Bio.Seq import Seq
from Bio import SeqIO
from numpy.lib.arraysetops import _unique_dispatcher
from networkx.drawing.nx_agraph import to_agraph
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

#####################################################################
## utility fuctions #################################################

def lenLoci(loci):
    length = 0
    for region in loci.split(','):
        chrom, start, end = region.split('_')[:3]
        length += int(end)-int(start)
    return length

def format_merge_region(merge_region):
    list_formatted = []
    for region in merge_region.split(','):
        chrom, start, end, strand = region.split('_')
        list_formatted.append("{}:{}-{}_{}".format(chrom, int(start) + 1, end, strand))
    return ','.join(list_formatted)

def draw_graph(graph, filename="test"):
    import graphviz as gv
    d = gv.Digraph(name=filename)
    for k in graph:
        node = eval(k)
        d.edge(node[0], node[1], label=str(graph[k]))
    d.render()
    return filename

def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c)

def majority_strand(df_merged):
    strand = '+'
    list_str_strand = df_merged['strand'].astype(str).tolist()
    if list_str_strand.count('-1') > list_str_strand.count('1'):
        strand = '-'
    return strand

def chk_ovl_nodes(df_final_read_mergeid, list_nodes):
    list_uniq_nodes = list(set(list_nodes))
    l1 = df_final_read_mergeid['mergeid'].isin(list_uniq_nodes)
    df_x = df_final_read_mergeid[l1].groupby(by='readid')[['mergeid']].nunique().sort_values(by='mergeid', ascending=False)
    list_readid = list(set(df_x[df_x['mergeid'] == len(list_uniq_nodes)].index))
    
    chk_return = False
    if len(list_readid) > 0:
        df_y = df_final_read_mergeid[df_final_read_mergeid['readid'].isin(list_readid)].groupby(by='mergeid').agg({'ovl_5end':['sum'], 'ovl_3end':['sum']})
        list_boo = []
        for node in list_uniq_nodes:
            list_boo.append(df_y.loc[node,'ovl_3end']['sum'] > 0)
            list_boo.append(df_y.loc[node,'ovl_5end']['sum'] > 0)
        if all(list_boo):
            chk_return = True
        
    return chk_return

def chk_circular_subgraph(G, graph, dict_pair_strand):
    list_nodes = list(graph.nodes)
    num_nodes = len(list_nodes)
    print_nodes = ",".join(list_nodes)
    sl = True if len(list(nx.selfloop_edges(graph))) > 0 else False
    solved = True
    cyclic = True
    if num_nodes == 2:
        rm_sl_subgraph = graph.copy()
        rm_sl_subgraph.remove_edges_from(nx.selfloop_edges(rm_sl_subgraph))
        solved = all(x <= 2 for x in [rm_sl_subgraph.degree[node] for node in list_nodes])
        key_forward = repr(list_nodes)
        key_reverse = repr(list_nodes[::-1])
        cyclic = False
        if all(x in dict_pair_strand.keys() for x in [key_forward, key_reverse]):
            for fstrands in dict_pair_strand[key_forward]:
                fl, fr  = fstrands.split('_')
                key_check = "{}{}{}".format(fr, key_reverse ,fl)
                for rstrands in dict_pair_strand[key_reverse]:
                    rl, rr  = rstrands.split('_')
                    key_reverse_check = "{}{}{}".format(rl, key_reverse ,rr)
                    if key_check == key_reverse_check:
                        cyclic = True
                        break
    elif num_nodes > 2:
        rm_sl_subgraph = graph.copy()
        rm_sl_subgraph.remove_edges_from(nx.selfloop_edges(rm_sl_subgraph))
        solved = all(x <= 2 for x in [rm_sl_subgraph.degree[node] for node in list_nodes])
        cyclic = True if len(nx.cycle_basis(nx.DiGraph(rm_sl_subgraph).to_undirected())) > 0 else False
    else:
        cyclic = True if len(nx.cycle_basis(nx.DiGraph(graph).to_undirected())) > 0 else False
        
    if not solved:
        cyclic = False
    
    test_graph = graph.copy()
    test_graph.remove_edges_from(nx.selfloop_edges(test_graph))
    list_traversal = nx.cycle_basis(nx.DiGraph(test_graph).to_undirected())
    if len(list_traversal) > 0:
        list_traversal_ini = list_traversal[0]
        if len(list_traversal_ini) == len(list_nodes):
            print_nodes = ",".join(list_traversal_ini)

    return (print_nodes, num_nodes, solved, sl, cyclic)

def reverse_strand(strand):
    return ''.join(['+' if x == '-' else '-' if x == '+' else x for x in list(strand)])

def get_longest_node(chk_subgraph):
    list_nodes = list(chk_subgraph.nodes)
    longest_node = list_nodes[0]
    for node in list_nodes[1:]:
        if lenLoci(node) > lenLoci(longest_node):
            longest_node = node
    return longest_node
    
def get_dict_weight_subgraph(directed_graph, chk_subgraph):
    dict_return = {}
    for node_name in list(chk_subgraph.nodes):
        for k, v in directed_graph.adj[node_name].items():
            key = "{},{}".format(node_name, k)
            dict_return[key] = v[0]['weight']
    return dict_return

def nodes_to_merge_regions(list_nodes):
    list_merge_regions = []
    for node in list_nodes:
        chr_, start_, end_ = node.split('_')
        region_ = "{}:{}-{}".format(chr_, int(start_) + 1, end_)
        list_merge_regions.append(region_)
    return ','.join(list_merge_regions)

def get_fasta_length(fasta_path):
    total_len = 0
    for record in SeqIO.parse(fasta_path, 'fasta'):
        total_len += len(str(record.seq))
    return total_len

def get_consensus_status(value, dict_medaka_seq_len):
    status = '0|Unknown'
    id_ = value['id']
    merge_len = value['merge_len']
    consensus_len = dict_medaka_seq_len.get(id_, None)
    if consensus_len:
        if merge_len > consensus_len:
            status = "{}|Deletion".format(consensus_len)
        elif merge_len < consensus_len:
            status = "{}|Insertion".format(consensus_len)
        else:
            status = "{}|No_indel".format(consensus_len)
    return status

def any_overlapping_range(start1, end1, start2, end2):
    if start1 <= end2 and start2 <= end1:
        return True
    else:
        return False

def check_breakpoint_direction(df_check):
    
    list_pairs = []

    df_check.reset_index(drop=True, inplace=True)

    list_order_check = [(i , i+1) for i in df_check.index.tolist() if i + 1 < len(df_check)]

    for tup_check in list_order_check:

        l1 = df_check.index.isin(tup_check)

        q_starts = df_check[l1]['q_start'].tolist()
        q_ends = df_check[l1]['q_end'].tolist()

        if any_overlapping_range(q_starts[0], q_ends[0] + 50, q_starts[1], q_ends[1]):

            list_mergeid = df_check[l1]['mergeid'].tolist()
            list_check = df_check[l1]['strand'].tolist()

            if list_check[0] == 1 and list_check[1] == 1:
                list_5 = df_check[l1]['ovl_5end'].tolist()
                list_3 = df_check[l1]['ovl_3end'].tolist()
                if list_3[0] == 1 and list_5[1] == 1:
                    pair = (list_mergeid[0], list_mergeid[1], '+_+' , True)
                    list_pairs.append(pair)

            if list_check[0] == -1 and list_check[1] == -1:
                list_5 = df_check[l1]['ovl_5end'].tolist()
                list_3 = df_check[l1]['ovl_3end'].tolist()
                if list_3[1] == 1 and list_5[0] == 1:
                    pair = (list_mergeid[1], list_mergeid[0], '-_-', True)
                    list_pairs.append(pair)


            if list_check[0] == -1 and list_check[1] == 1:
                list_5 = df_check[l1]['ovl_5end'].tolist()
                list_3 = df_check[l1]['ovl_3end'].tolist()
                if list_5[0] == 1 and list_5[1] == 1:
                    pair = (list_mergeid[0], list_mergeid[1], '-_+', True)
                    list_pairs.append(pair)


            if list_check[0] == 1 and list_check[1] == -1:
                list_5 = df_check[l1]['ovl_5end'].tolist()
                list_3 = df_check[l1]['ovl_3end'].tolist()
                if list_3[0] == 1 and list_3[1] == 1:
                    pair = (list_mergeid[0], list_mergeid[1], '+_-', True)
                    list_pairs.append(pair)
        else:
            break

    return list_pairs

def check_abs_ovl(value):
    
    cut_off=1
    
    r_start = value['r_start']
    r_end = value['r_end']
    mergeid = value['mergeid']
    merge_5end = int(mergeid.split('_')[1])
    merge_3end = int(mergeid.split('_')[2])
    
    ovl_5end_value = abs(r_start - merge_5end)
    ovl_3end_value = abs(r_end - merge_3end)
    
    ovl_5end = 1 if ovl_5end_value <= cut_off else 0
    ovl_3end = 1 if ovl_3end_value <= cut_off else 0
    sum_ends = ovl_5end + ovl_3end
    
    return ovl_5end, ovl_3end, sum_ends

def prepare_identify_seq(tup_value):

    gname, merge_region, num_region, is_cyclic, assemGraph, fastaRef, fastaName = tup_value

    fa_ref = pysam.Fastafile(fastaRef)
    fa_reads = pysam.FastaFile(fastaName)

    lengthLoci = lenLoci(merge_region)
    eccdna_status = 'cyclic' if is_cyclic == True else 'non-cyclic'

    ## create a subgraph folder
    assemFol = "{}/{}".format(assemGraph, gname)
    pathlib.Path(assemFol).mkdir(parents=True, exist_ok=True)

    ## write a reference sequence from regions
    ref_fasta_path = "{}/reference_regions.fa".format(assemFol)
    seq_ = ""
    with open(ref_fasta_path, 'w') as w_f:
        w_f.write(">{}\n".format(gname))
        for node in merge_region.split(','):
            chr_, start_, end_, strand_ = node.split('_')
            chr_ = str(chr_)
            if strand_ == '-':
                seq_ += str(Seq(fa_ref.fetch(chr_, int(start_), int(end_))).reverse_complement()).upper()
            else:
                seq_ += str(fa_ref.fetch(chr_, int(start_), int(end_))).upper()
        w_f.write(seq_ + "\n")

    ## write all regions in reads associating eccDNA regions
    total_base = 0
    list_readid = []

    out_filename_final = "{}/final_reads.fa".format(assemFol)
    with open(out_filename_final, mode='w') as foutF:
        for node in merge_region.split(','):
            l1 = read_merged_ins_df['mergeid'] == node[:-2]
            readO = read_merged_ins_df[l1].copy()
            for readO_idx, readO_row in readO.iterrows():
                readO_readid = str(readO_row['readid'])
                readO_start = readO_row['q_start']
                readO_end = readO_row['q_end']
                readO_name = "{}_{}_{}_{}".format(gname, readO_readid, readO_start, readO_end)
                readO_seq = str(fa_reads.fetch(readO_readid, readO_start, readO_end)).upper()
                foutF.write(">{}\n{}\n".format(readO_name, readO_seq))

                total_base += len(readO_seq)
                list_readid.append(readO_readid)

    expectCov = "{:.2f}".format((float(total_base) / lengthLoci))

    list_tup_result = [(gname, merge_region, lengthLoci, num_region, eccdna_status, len(set(list_readid)), total_base, expectCov)]

    return list_tup_result

def dist_work(cmd):
    process = subprocess.call("{}".format(cmd), shell=True)

def argparser():

    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='Identify and verify eccDNA from whole-genome long-read (WGLS) sequencing data')
    general = parser.add_argument_group(title='General options')
    general.add_argument('-t', "--threads",
                            help="Number of threads [all CPU cores]",
                            type=int, default=0)
    general.add_argument('-m', "--mode", dest='mode',
                            help="mode to identify potential regions for eccDNAs : linkage or depth [depth]",
                            type=str, default='depth')
    general.add_argument('-r', "--ref",  
                            help="reference Minimap2 index .mmi",
                            type=str, default=None)
    general.add_argument('-fa', "--fa-ref", dest='faref',
                            help="reference sequence .fasta",
                            type=str, default=None)
    general.add_argument('-fai', "--fa-index", dest='fai',
                            help="genome size file e.g. hg19.chrom.sizes or file.fasta.fai",
                            type=str, default=None)
    general.add_argument('-fq', "--fq-input",  dest='fqinput',
                            help="input fasta/fastq",
                            type=str, default=None)
    general.add_argument('-trim', "--trim-input", dest='triminput',
                            help="CReSIL trim table",
                            type=str, default=None)
    general.add_argument('-minsize', "--minimum-eccdna-size", dest='minsize',
                            help="minimum size of eccDNA [200]",
                            type=int, default=200)
    general.add_argument('-ovl', "--ovl-size", dest='ovlsize',
                            help="size of 5' and 3' regions on reads for a breakpiont overlapping check [50]",
                            type=int, default=50)
    general.add_argument('-break', "--brakpoint-depth", dest='breakpointdepth',
                            help="lowest number of supported breakpoints of eccDNA [3]",
                            type=int, default=3)
    general.add_argument('-depth', "--average-depth", dest='averagedepth',
                            help="lowest average depth of eccDNA [5]",
                            type=int, default=5)
    general.add_argument('-mwindow', "--window-merge", dest='windowmerge',
                            help="maximum size of poteintial eccDNA regions (effective with -m linkage) [100000]",
                            type=int, default=100000)
    general.add_argument('-mbreak', "--brakpoint-merge", dest='breakpointmerge',
                            help="lowest number of breakpoints for estimating poteintial eccDNA regions (effective with -m depth) [2]",
                            type=int, default=2)
    general.add_argument('-moff', "--offset-merge", dest='offsetmerge',
                            help="offset below average depth for estimating poteintial eccDNA regions (effective with -m depth) [2]",
                            type=int, default=2)
    return parser

def main(args):

    threads = cpu_count() if args.threads == 0 else args.threads
    adj_threads = max(4, int(threads / 4))
    bed_threads = max(4, int(threads / 5))
    mode = 'linkage' if args.mode == 'linkage' else 'depth'
    fref = args.ref
    fastaRef = args.faref
    chromSizes = args.fai
    fastaName = args.fqinput
    minsize = args.minsize
    check_ovl_size = args.ovlsize
    breakpointdepth = args.breakpointdepth
    regiondepth = args.averagedepth
    windowmerge = args.windowmerge
    breakpointmerge = args.breakpointmerge
    offsetmerge = args.offsetmerge

    if not args.triminput:
        sys.stderr.write("[ABORT] a table from trimming step (trim.txt) is needed\n")
        exit(1)

    if args.threads == 0:
        threads = os.cpu_count()

    bname = str(pathlib.Path(args.triminput).parent)

    readTrim = pd.read_csv(args.triminput, sep="\t", header=0)

    outDir="{}/cresil_run".format(bname)
    tmpDir="{}/tmp".format(outDir)
    assemGraph = "{}/assemGraph".format(outDir)

    ## input ############################################################

    #####################################################################
    ## output ###########################################################

    pathlib.Path(outDir).mkdir(parents=True, exist_ok=True)
    pathlib.Path(tmpDir).mkdir(parents=True, exist_ok=True)
    pathlib.Path(assemGraph).mkdir(parents=True, exist_ok=True)

    ## output ###########################################################
    
    ct = datetime.datetime.now()
    print("\n######### CReSIL : Start identifying process : WGLS - {} (thread : {})\n[{}] total trimmed region : {}".format(mode, threads, ct, len(readTrim)), flush=True)
    if len(readTrim) == 0:
        sys.exit("[ABORT] zero trimmed region\n")

    #####################################################################
    ## Identify potential eccDNA regions ################################

    ## Calculate read coverage and filter out restuion with depth <= 5x
    ord_header = ['ref', 'r_start', 'r_end', 'readid', 'q_start', 'q_end', 'match',
                  'mapBlock', 'mapq', 'strand', 'qlenTrimmed', 'freqCov', 'order']
    readTrim = readTrim.loc[:,ord_header]

    ## prepare aligned reads
    aln_reads = bt.BedTool.from_dataframe(readTrim).sort()

    ## calculate reads coverage by bedtools genomecov
    genome_cov = aln_reads.genome_coverage(bg=True, g=chromSizes)

    if mode == 'linkage':
        ## identify potential eccDNA regions by linkage
        ct = datetime.datetime.now()
        print("[{}] identifying potential ecccDNA regions by linkage".format(ct), flush=True)

        trim_sup = readTrim.query('order>0')

        list_dfs = []
        l1 = trim_sup['strand'] == 1
        l2 = trim_sup['strand'] == -1
        trim_sup_plus = trim_sup[l1][['ref','r_start','r_start','readid']]
        trim_sup_plus.columns = [0,1,2,3]
        trim_sup_plus[2] = trim_sup_plus[2].apply(lambda x: x + 5)
        list_dfs.append(trim_sup_plus)

        trim_sup_minus = trim_sup[l2][['ref','r_end','r_end','readid']]
        trim_sup_minus.columns = [0,1,2,3]
        trim_sup_minus[1] = trim_sup_minus[1].apply(lambda x: x - 5)
        list_dfs.append(trim_sup_minus)

        df_breaks_bed = pd.concat(list_dfs, ignore_index=False).sort_values(by=[0,1])
        breaks_bed = bt.BedTool.from_dataframe(df_breaks_bed)

        df_read_break_bed = trim_sup[['ref','r_start','r_end','readid']]
        read_break_bed = bt.BedTool.from_dataframe(df_read_break_bed)

        breaks_bdg = breaks_bed.genome_coverage(bg=True, g=chromSizes)

        df_filtered_breaks_bdg = breaks_bdg.to_dataframe()
        df_filtered_breaks_bdg = df_filtered_breaks_bdg[df_filtered_breaks_bdg['name']>0]
        filtered_breaks_bdg = bt.BedTool.from_dataframe(df_filtered_breaks_bdg)

        df_filtered_breaks_bdg = breaks_bdg.to_dataframe()
        df_filtered_breaks_bdg = df_filtered_breaks_bdg[df_filtered_breaks_bdg['name']>0]
        filtered_breaks_bdg = bt.BedTool.from_dataframe(df_filtered_breaks_bdg)

        break_min1x_merge_bed = filtered_breaks_bdg.merge()

        df_ref_chromSizes = pd.read_csv(chromSizes, sep='\t', header=None)
        df_ref_chrom_bed = pd.DataFrame(list(zip(df_ref_chromSizes[0],[0 for x in range(len(df_ref_chromSizes))],df_ref_chromSizes[1])))
        ref_chrom_bed = bt.BedTool.from_dataframe(df_ref_chrom_bed)

        ref_chrom_subtract_filtered_breaks_bed = ref_chrom_bed.subtract(filtered_breaks_bdg)
        df_ref_chrom_subtract_filtered_breaks_bed = ref_chrom_subtract_filtered_breaks_bed.to_dataframe()
        df_ref_chrom_subtract_filtered_breaks_bed['len_region'] = df_ref_chrom_subtract_filtered_breaks_bed['end'] - df_ref_chrom_subtract_filtered_breaks_bed['start']
        df_ref_chrom_subtract_filtered_breaks_bed = df_ref_chrom_subtract_filtered_breaks_bed[df_ref_chrom_subtract_filtered_breaks_bed['len_region'] <= windowmerge]
        filtered_subtract_regions = bt.BedTool.from_dataframe(df_ref_chrom_subtract_filtered_breaks_bed)

        df_rem_bg_aln_reads_merge = filtered_subtract_regions.intersect(read_break_bed, wa=True, f=0.9, u=True).to_dataframe()
        df_rem_bg_aln_reads_merge = df_rem_bg_aln_reads_merge[['chrom','start','end']]

        ct = datetime.datetime.now()
        print("[{}] finished identifying potential regions : {}".format(ct, len(df_rem_bg_aln_reads_merge)), flush=True)
    else:
        ## identify potential eccDNA regions by depth
        ct = datetime.datetime.now()
        print("[{}] identifying potential eccDNA regions by depth\n[{}] creating a BAM file".format(ct, ct), flush=True)

        cmd = "minimap2 -t {} --no-long-join -a {} {} | samtools sort -o {}/ref_aln.bam; samtools index {}/ref_aln.bam".format(threads, fref, fastaName, tmpDir, tmpDir)
        process = subprocess.call(cmd, shell=True)

        ct = datetime.datetime.now()
        print("[{}] finished creating a BAM file\n[{}] estimating read depth".format(ct, ct), flush=True)

        cmd = "mosdepth -n -t 4 --fast-mode --by 10 ref_aln ref_aln.bam"
        process = subprocess.call(cmd, shell=True, cwd=tmpDir)

        ct = datetime.datetime.now()
        print("[{}] finished estimating read depth\n[{}] identifying potential regions".format(ct, ct), flush=True)

        mosdepth_summary_path = "{}/ref_aln.mosdepth.summary.txt".format(tmpDir)
        mosdepth_region_bed_gz_path = "{}/ref_aln.regions.bed.gz".format(tmpDir)
        filtered_region_depth_path = "{}/ref_aln.filtered_region_depth.bdg".format(tmpDir)

        df_mosdepth_summary = pd.read_csv(mosdepth_summary_path, sep='\t')
        df_mosdepth_summary = df_mosdepth_summary[~df_mosdepth_summary['chrom'].str.contains('region')]
        df_mosdepth_summary['mean'] = df_mosdepth_summary['mean'].apply(lambda x: max(0, float(int(x) - offsetmerge)))
        dict_chrom_avg_depth = dict(zip(df_mosdepth_summary['chrom'], df_mosdepth_summary['mean']))

        chunk_size = 10**7
        for df in pd.read_csv(mosdepth_region_bed_gz_path, compression='gzip', sep='\t', header=None, chunksize=chunk_size):
            df[4] = df[0].apply(lambda x: dict_chrom_avg_depth.get(x, 0.0))
            l1 = df[3] >= df[4]
            l2 = df[0].isin(dict_chrom_avg_depth.keys())
            filtered_df = df[l1 & l2][[0,1,2,3]]
            
            if len(filtered_df) > 0:
                filtered_df.to_csv(filtered_region_depth_path, mode='a', sep='\t', header=None, index=None)

        filtered_region_depth_bed = bt.BedTool(filtered_region_depth_path)
        merged_filtered_region_depth = filtered_region_depth_bed.merge()

        trim_sup = readTrim.query('order>0')

        list_dfs = []
        l1 = trim_sup['strand'] == 1
        l2 = trim_sup['strand'] == -1
        trim_sup_plus = trim_sup[l1][['ref','r_start','r_start','readid']]
        trim_sup_plus.columns = [0,1,2,3]
        trim_sup_plus[2] = trim_sup_plus[2].apply(lambda x: x + 5)
        list_dfs.append(trim_sup_plus)

        trim_sup_minus = trim_sup[l2][['ref','r_end','r_end','readid']]
        trim_sup_minus.columns = [0,1,2,3]
        trim_sup_minus[1] = trim_sup_minus[1].apply(lambda x: x - 5)
        list_dfs.append(trim_sup_minus)

        df_breaks_bed = pd.concat(list_dfs, ignore_index=False).sort_values(by=[0,1])
        breaks_bed = bt.BedTool.from_dataframe(df_breaks_bed)

        df_read_break_bed = trim_sup[['ref','r_start','r_end','readid']]
        read_break_bed = bt.BedTool.from_dataframe(df_read_break_bed)

        breaks_bdg = breaks_bed.genome_coverage(bg=True, g=chromSizes)

        df_filtered_breaks_bdg = breaks_bdg.to_dataframe()
        df_filtered_breaks_bdg = df_filtered_breaks_bdg[df_filtered_breaks_bdg['name'] >= breakpointmerge]
        filtered_breaks_bdg = bt.BedTool.from_dataframe(df_filtered_breaks_bdg)

        break_min2x_merge_bed = filtered_breaks_bdg.merge()

        sub_breaks2x_bed = merged_filtered_region_depth.subtract(break_min2x_merge_bed)

        merged_window = sub_breaks2x_bed.window(break_min2x_merge_bed, w=1)

        df_count_merged_window = merged_window.groupby(g=[1, 2, 3], c=4, o=['count']).to_dataframe()
        df_count_merged_window = df_count_merged_window[df_count_merged_window['name'] > 1]
        count_merged_window = bt.BedTool.from_dataframe(df_count_merged_window)

        df_rem_bg_aln_reads_merge = count_merged_window.slop(b=5, g=chromSizes).to_dataframe()
        df_rem_bg_aln_reads_merge = df_rem_bg_aln_reads_merge[['chrom','start','end']]

        ct = datetime.datetime.now()
        print("[{}] finished identifying potential regions : {}".format(ct, len(df_rem_bg_aln_reads_merge)), flush=True)

    if len(df_rem_bg_aln_reads_merge) == 0:
        sys.exit("[ABORT] no potential merge region was detected\n")

    aln_reads_merge = bt.BedTool.from_dataframe(df_rem_bg_aln_reads_merge)

    ## join merged regions with read coverage and perform filtering
    merge_genomecov = aln_reads_merge.intersect(genome_cov, output='{}/merge_genomecov.txt'.format(tmpDir), wo=True)
    merge_genomecov_df = pd.read_csv("{}/merge_genomecov.txt".format(tmpDir),sep="\t",header=None,dtype=str)
    merge_genomecov_df.columns = ['m_chrom','m_start','m_end','bg_chrom','bg_start','bg_end','depth','d_ovl']
    merge_genomecov_df[['depth','bg_start','d_ovl']] = merge_genomecov_df[['depth','bg_start','d_ovl']].apply(pd.to_numeric, errors='coerce')

    ## generate mergeID
    merge_genomecov_df['mergeid'] = merge_genomecov_df['m_chrom'] + "_" + merge_genomecov_df['m_start'] + "_" + merge_genomecov_df['m_end']

    ## trim low-coverage regions
    merge_genomecov_df_filt = merge_genomecov_df[merge_genomecov_df['depth'] >= regiondepth]
    merge_genomecov_df = merge_genomecov_df.sort_values(['mergeid','bg_start'])
    ## group reads coverage of each mergeid and calculate mean coverage
    merge_bg_filt = merge_genomecov_df_filt.groupby(['mergeid']).agg({'bg_chrom': np.max,'bg_start' : np.min,
                                                           'bg_end' : np.max,'depth' : np.mean}).reset_index()
    ## create final potential eccDNA regions and 
    ## re-create mergeid based on trimming low-coverage regions
    merge_bg_filt[['bg_start','bg_end']] = merge_bg_filt[['bg_start','bg_end']].apply(pd.to_numeric, errors='coerce')
    merge_bg_filt['length'] = merge_bg_filt.bg_end - merge_bg_filt.bg_start
    merge_bg_filt["bg_start"]= merge_bg_filt["bg_start"].astype(str) 
    merge_bg_filt["bg_end"]= merge_bg_filt["bg_end"].astype(str) 
    merge_bg_filt['mergeid'] = merge_bg_filt.bg_chrom+"_"+merge_bg_filt.bg_start+"_"+merge_bg_filt.bg_end
    merge_bg_filt = merge_bg_filt.loc[:,['bg_chrom','bg_start','bg_end','depth','length','mergeid']]

    # filter regions with the length >= 200 bp
    merge_bg_filt = merge_bg_filt[merge_bg_filt['length'] >= minsize]

    ## Identify potential eccDNA regions ################################

    #####################################################################
    ## Assign readid to each potential eccDNA regions (mergeid) #########

    ct = datetime.datetime.now()
    print("[{}] calculating breakpoints and merging regions".format(ct), flush=True)

    ## create merge region
    aln_reads_merge = bt.BedTool.from_dataframe(merge_bg_filt)
    read_merged_intersect = aln_reads_merge.intersect(aln_reads, output='{}/reads_merge_intersect.bed'.format(tmpDir), wo=True)
    ## read merge region and create mergeid

    if os.stat("{}/reads_merge_intersect.bed".format(tmpDir)).st_size == 0:
        sys.exit("[ABORT] no merge region was detected (lower -breakpoint, -depth; or region too short)\n")


    read_merged_ins_df_org = pd.read_csv("{}/reads_merge_intersect.bed".format(tmpDir),sep="\t",header=None,dtype=str)
    header = list(merge_bg_filt.columns)+list(readTrim.columns)+["ovl"]
    read_merged_ins_df_org.columns = header
    header_select = ['ref','r_start','r_end','mergeid','depth','length','readid', 
                     'q_start','q_end','match','mapBlock','mapq','strand','qlenTrimmed', 
                     'freqCov','order','ovl']
    read_merged_ins_df_org = read_merged_ins_df_org[header_select]
    ## Assign readid to each potential eccDNA regions (mergeid) #########

    #####################################################################
    ## Annotate 5'end and 3'end region of each readid ###################

    ## collect 200bp from 5'end and 3'end of each potential eccDNA
    end5 = open("{}/end5_merge_region.bed".format(tmpDir),'w')
    end3 = open("{}/end3_merge_region.bed".format(tmpDir),'w')
    for idx_, row in merge_bg_filt.iterrows():
        chrom = row.bg_chrom
        start = row.bg_start
        end = row.bg_end
        end_size = check_ovl_size if minsize >= 200 else int(round(minsize * 0.3, 0))
        end5.write("{}\n".format("\t".join(map(str, [chrom, start, int(start) + end_size]))))
        end3.write("{}\n".format("\t".join(map(str, [chrom, int(end) - end_size, end]))))
    end5.close()
    end3.close()
    bt_5end = bt.BedTool("{}/end5_merge_region.bed".format(tmpDir))
    bt_3end = bt.BedTool("{}/end3_merge_region.bed".format(tmpDir))
    ## Annotate each read that has 5'end and 3'end overlap
    ## add column ovl_5end and ovl_3end (0=no overlap, 1=presence of overlap) 
    read_merged_ins_df_org_bt = bt.BedTool.from_dataframe(read_merged_ins_df_org)
    read_merged_final_bt = read_merged_ins_df_org_bt.annotate(files=["{}/end5_merge_region.bed".format(tmpDir),"{}/end3_merge_region.bed".format(tmpDir)], counts=True)

    global read_merged_ins_df
    read_merged_ins_df = read_merged_final_bt.to_dataframe(names=header_select + ['ovl_5end','ovl_3end'])
    ## sum 5'end and 3'end overlap (0=no overlap, 1=either 5'end or 5'end overlap, 2=both ends overlap)
    read_merged_ins_df['sum_ends'] = read_merged_ins_df.ovl_5end + read_merged_ins_df.ovl_3end

    ## create dictionary of nodes and their temporary strands
    dict_majority_strand = {}
    for group_name, group in read_merged_ins_df.groupby(by='mergeid'):
        dict_majority_strand[group_name] = majority_strand(group)

    ## create dictionary of node pairs and all graphs
    ct = datetime.datetime.now()
    print("[{}] analyzing graphs".format(ct), flush=True)

    dict_pair_strand = {}
    graph = {}

    for readid, group in read_merged_ins_df.groupby(by='readid'):
        
        df_check = group.sort_values(by='order').copy()
        df_check.reset_index(drop=True, inplace=True)
        
        if len(df_check) > 1:
            for tup_result in check_breakpoint_direction(df_check):
                if tup_result[3] == True:

                    repr_mergeid = repr([tup_result[0], tup_result[1]])
                    
                    if repr_mergeid in dict_pair_strand.keys():
                        dict_pair_strand[repr_mergeid].add(tup_result[2])
                    else:
                        dict_pair_strand[repr_mergeid] = set()
                        dict_pair_strand[repr_mergeid].add(tup_result[2])

                    if repr_mergeid in graph.keys():
                        graph[repr_mergeid] += 1
                    else:
                        graph[repr_mergeid] = 1

    ## filter graphs with low breakpoints
    graphFilt = {k: v for k, v in graph.items() if v >= breakpointdepth}

    ## create all graph objects
    G = nx.MultiDiGraph()
    for k in graphFilt:
        nodes = eval(k)
        G.add_edge(nodes[0], nodes[1], weight = graphFilt[k])

    ## create list of subgraphs
    subgraphs = list(connected_component_subgraphs(G.to_undirected()))

    ct = datetime.datetime.now()
    print("[{}] initial subgraphs : {}".format(ct, len(subgraphs)), flush=True)

    ## correct eccDNA strands and create a graph summary file
    list_graph_summary = []

    for idx, graph in enumerate(subgraphs, start=0):
        nodes = list(graph.nodes())
        gname = "ec{}".format(idx + 1)

        regions, num_nodes, can_be_solved, contain_selfloop, is_cyclic = chk_circular_subgraph(G, graph, dict_pair_strand)
        
        if len(nodes) == 1:
            select_repr = repr([nodes[0], nodes[0]])
            regions = "{}_{}".format(eval(select_repr)[0], list(dict_pair_strand[select_repr])[0][0])

        elif len(nodes) == 2:
            select_repr = repr(nodes)
            list_region = eval(select_repr)
            
            if not select_repr in dict_pair_strand.keys():
                select_repr = repr(nodes[::-1])
                
            list_strand = list(dict_pair_strand[select_repr])[0].split('_')
            regions = ','.join("{}_{}".format(x[0], x[1]) for x in list(zip(list_region, list_strand)))
            
        else:
            test_graph = graph.copy()
            test_graph.remove_edges_from(nx.selfloop_edges(test_graph))
            
            list_traversal = nx.cycle_basis(nx.DiGraph(test_graph).to_undirected())
            if len(list_traversal) == 0:
                regions = ','.join(["{}_{}".format(node, dict_majority_strand[node]) for node in nodes])

            else:
                list_traversal = list_traversal[0]

                list_order_check = [(i , i+1) for i in list(range(len(list_traversal))) if i + 1 < len(list_traversal)]

                list_temp_all = []
                for tup_check in list_order_check:
                    l_region = list_traversal[tup_check[0]]
                    r_region = list_traversal[tup_check[1]]
                    repr_mergeid = repr([l_region, r_region])

                    list_this_level = []

                    if repr_mergeid in dict_pair_strand.keys():

                        for str_strand in list(dict_pair_strand[repr_mergeid]):
                            list_str_strand = ['_'.join(x) for x in list(zip([l_region, r_region], str_strand.split('_')))]
                            list_this_level.append(list_str_strand)
                    else:
                        rev_repr_mergeid = repr([r_region, l_region])

                        for str_strand in [ reverse_strand(x) for x  in list(dict_pair_strand[rev_repr_mergeid]) ]:
                            list_str_strand = ['_'.join(x) for x in list(zip([l_region, r_region], str_strand.split('_')))]
                            list_this_level.append(list_str_strand)

                    list_temp_all.append(list_this_level)


                list_traverse = list_temp_all[0]
                for list_level in list_temp_all[1:]:
                    for list_ in list_level:
                        for index, list_traverse_level in enumerate(list_traverse, start=0):
                            if list_[0] == list_traverse_level[-1]:
                                list_traverse[index].append(list_[1])

                regions = ','.join(max(list_traverse, key=len))
        
        tup_graph_summary = (gname, regions, num_nodes, can_be_solved, contain_selfloop, is_cyclic)
        
        list_graph_summary.append(tup_graph_summary)
        
    if len(list_graph_summary) > 0:
        header = ['id', 'regions', 'num_nodes', 'can_be_solved', 'contain_selfloop', 'is_cyclic']
        df_graph_summary = pd.DataFrame(list_graph_summary, columns=header)
    else:
        sys.exit("[ABORT] no eccDNA was detected (lower -breakpoint, -depth; or region too short)\n")

    ## Run consensus sequence and variants
    ct = datetime.datetime.now()
    print("[{}] preparing data from subGraphs : {}".format(ct, len(subgraphs)), flush=True)

    list_identify_results = []
    list_tup_value = []
    iter_ = 100

    for counter, (idx, value) in enumerate(df_graph_summary.iterrows(), start=1):

        if counter % iter_ == 0:
            ct = datetime.datetime.now()
            print("[{}] {}".format(ct, counter), flush=True)

        if counter % iter_ == 0:
            gname = value['id']
            merge_region = value['regions']
            num_region = value['num_nodes']
            is_cyclic = value['is_cyclic']

            tup_value = (gname, merge_region, num_region, is_cyclic, assemGraph, fastaRef, fastaName)
            list_tup_value.append(tup_value)

            pool = Pool(processes = bed_threads)
            list_pool_result = pool.map(prepare_identify_seq, list_tup_value)
            pool.close()
            pool.join()

            list_identify_results += list(chain.from_iterable(list_pool_result))

            list_tup_value = []
        else:
            gname = value['id']
            merge_region = value['regions']
            num_region = value['num_nodes']
            is_cyclic = value['is_cyclic']

            tup_value = (gname, merge_region, num_region, is_cyclic, assemGraph, fastaRef, fastaName)
            list_tup_value.append(tup_value)

    if len(list_tup_value) > 0:
        pool = Pool(processes = bed_threads)
        list_pool_result = pool.map(prepare_identify_seq, list_tup_value)
        pool.close()
        pool.join()

        list_identify_results += list(chain.from_iterable(list_pool_result))

    ## write a subGraph summary file
    header = ['id', 'merge_region', 'merge_len', 'num_region', 'eccdna_status', 'numreads', 'totalbase', 'coverage']
    selected_header = ['id', 'merge_region', 'merge_len', 'num_region', 'can_be_solved', 'contain_selfloop', 'eccdna_status', 'numreads', 'totalbase', 'coverage']
    df_temp_graph_summary = pd.DataFrame(list_identify_results, columns=header)
    df_final_graph_summary = pd.merge(left=df_graph_summary, right=df_temp_graph_summary, left_on='id', right_on='id', how='inner')
    df_final_graph_summary = df_final_graph_summary[selected_header]
    df_final_graph_summary['merge_region'] = df_final_graph_summary['merge_region'].apply(lambda x: format_merge_region(x))

    write_path = "{}/subGraphs.summary.txt".format(outDir)
    df_final_graph_summary.to_csv(write_path, sep='\t', index=None)

    ## write identified eccDNA summary
    header = ['id', 'merge_region', 'merge_len', 'num_region', 'eccdna_status', 'numreads', 'totalbase', 'coverage']
    df_identify_summary = df_final_graph_summary[header].copy()
    eccdna_final_path = "{}/eccDNA_final.txt".format(bname)
    df_identify_summary.to_csv(eccdna_final_path, sep='\t', index=None)

    ct = datetime.datetime.now()
    print("[{}] finished identifying process : WGLS - {}\n".format(ct, mode), flush=True)
    sys.exit()