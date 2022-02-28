#!/usr/bin/env python
import os
import shutil
import pathlib
import datetime
import fileinput
import subprocess
import pandas as pd
import numpy as np
import re, math
import pysam
import pybedtools
import intervaltree as tree
from Bio import SeqIO
from itertools import chain
from multiprocessing import cpu_count
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def format_merge_based1_region(merge_based1_region):
    list_formatted = []
    for region in merge_based1_region.split(','):
        chrom = region.split(':')[0]
        start, end = region.split(':')[1].split('_')[0].split('-')
        strand = region.split(':')[1].split('_')[1]
        list_formatted.append('_'.join([chrom, start, end, strand]))
    return ','.join(list_formatted)

def lenLoci(loci):
    length = 0
    for region in loci.split(','):
        chrom, start, end = region.split('_')[:3]
        length += int(end)-int(start)
    return length

def any_overlapping_range(start1, end1, start2, end2):
    if start1 <= end2 and start2 <= end1:
        return True
    else:
        return False

class modRegion(object):
    def __init__(self, tup_mod_region):
        self.chrom, self.start, self.end, self.strand, self.len_region, self.mod_start, self.mod_end = tup_mod_region

class mergedRegion(object):
    def __init__(self, merge_based1_regions):
        self.list_obj_region = []
        mod_start, mod_end = (0,0)
        for region in merge_based1_regions.split(','):
            mod_region = format_merge_based1_region(region)
            len_region = lenLoci(mod_region) + 1
            mod_end += len_region

            chrom, start, end, strand = mod_region.split('_')
            tup_mod_region = (chrom, int(start), int(end), strand, len_region, mod_start, mod_end)
            obj_ = modRegion(tup_mod_region)
            self.list_obj_region.append(obj_)
            
            mod_start += len_region
            
    def get_list_rel_region(self):
        return [(obj_.mod_start, obj_.mod_end) for obj_ in self.list_obj_region]

    def cal_overlap_rel_region(self, based1_region):
        list_overlap_rel_region = []
        q_chrom, q_start, q_end = based1_region.split('_')
        q_start = int(q_start)
        q_end = int(q_end)
        for obj_ in self.list_obj_region:
            if obj_.chrom == q_chrom:
                if any_overlapping_range(obj_.start, obj_.end, q_start, q_end):
                    rel_range_start = 0
                    rel_range_end = 0
                    
                    # cal start
                    if q_start <= obj_.start:
                        rel_range_start = obj_.mod_start
                    else:
                        diff_start = q_start - obj_.start
                        rel_range_start = (obj_.mod_start + 1) + diff_start
                    # cal end
                    if q_end >= obj_.end:
                        rel_range_end = obj_.mod_end
                    else:
                        diff_end = obj_.end - q_end
                        rel_range_end = obj_.mod_end - diff_end
                    
                    tup_overlap_rel_range = (rel_range_start, rel_range_end)
                    list_overlap_rel_region.append(tup_overlap_rel_range)
        
        return list_overlap_rel_region

def recalculate_windowing(value=None, type='max', window=100, step=1):
    if type=='real':
        return value
    global mem
    if value == None:
        return 0
    tmp = [[value]]*(window/step)
    if len(mem)==0:
        out = value
        mem = tmp
        mem.pop(0)
        return out
    else:
        mem.append([])
        for j in range(len(mem)):
            mem[j] = mem[j] + tmp[j]
        out = mem.pop(0)
        if type=='min':
            return min(out)
        elif type=='max':
            return max(out)
        elif type=='mean':
            if len(out) > 0:
                return sum(out)/float(len(out))
            else:
                return 0
        else:
            exit("only min, max, mean, real type available")

def wigToBedGraph(wigFile, genomeSizeFile, log2=False):
    genomeSize = {}
    for l in open(genomeSizeFile).readlines():
        L = l.strip().split("\t")
        genomeSize[L[0]] = int(L[1])

    list_bdg_result = []
    for l in open(wigFile).readlines():
        if re.search("variableStep", l):
            m = re.search(r"variableStep chrom=(\S+) span=(\d+)", l)
            if m:
                chrom = m.group(1)
                span = int(m.group(2))
                continue
        if re.search(r"^\d+", l):
            L = l.strip().split("\t")
            start = int(L[0])-1
            end = int(L[0])-1+span

            if genomeSize[chrom] < end:
                end = genomeSize[chrom]
            if log2:
                list_bdg_result.append(list(map(str,[chrom, start, end, math.log(float(L[1]),2)])))
            else:
                value = float(L[1])
                valueOut = recalculate_windowing(value,type='real',window=50,step=1)
                list_bdg_result.append(list(map(str,[chrom, start, end, valueOut])))

    return list_bdg_result

def bamfilterByReadid(infile_, outfile_):
    infile = pysam.AlignmentFile(infile_, "rb")
    outfile = pysam.AlignmentFile("{}.bam".format(outfile_), "wb", template=infile)

    query = set(open(outfile_).read().strip().split("\n"))
    
    for read in infile.fetch():
        if read.qname in query:
            outfile.write(read)

    outfile.close()
    infile.close()

def argparser():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='Identified eccDNA visualization')
    general = parser.add_argument_group(title='General options')
    general.add_argument('-t', "--threads",
                            help="Number of threads [all CPU cores]",
                            type=int, default=0)
    general.add_argument('-uc', "--unit-circos", dest='unitcircos',
                            help="size of unit in bp to display in Circos [1000]",
                            type=int, default=1000)
    general.add_argument('-mc', "--mode-circos", dest='modecircos',
                            help="mode to display ctc/non-ctc reads in Circos : ctc or proportion [proportion]",
                            type=str, default='proportion')
    general.add_argument('-c', "--eccdna-id",  dest='eccdna_id',
                            help="ID to visualize",
                            type=str, default=None)
    general.add_argument('-identify', "--identify-input", dest='identifyinput',
                            help="identified eccDNA table",
                            type=str, default=None)
    return parser


def main(args):

    threads = cpu_count() if args.threads == 0 else args.threads
    unitcircos = args.unitcircos
    modecircos = 'ctc' if args.modecircos == 'ctc' else 'proportion'
    ec_id = args.eccdna_id
    identify_path = args.identifyinput
    cresil_fol = str(pathlib.Path(identify_path).parent)

    annotate_fol = "{}/cresil_gAnnotation".format(cresil_fol)
    ec_fol = "{}/cresil_run/assemGraph/{}".format(cresil_fol, ec_id)

    circos_fol = "{}/for_Circos/{}".format(cresil_fol, ec_id)
    circos_data_fol = "{}/data".format(circos_fol)
    circos_temp_fol = "{}/tmp".format(circos_fol)

    reference_region_path = "{}/reference_regions.fa".format(ec_fol)
    final_read_path = "{}/final_reads.fa".format(ec_fol)

    ct = datetime.datetime.now()
    print("\n######### CReSIL : Start visualizing process (thread : {})\n[{}] initiating data and folder structure : {}".format(threads, ct, ec_id), flush=True)

    ## make a folder structure for visualizing
    pathlib.Path(circos_fol).mkdir(parents=True, exist_ok=True)
    pathlib.Path(circos_data_fol).mkdir(parents=True, exist_ok=True)
    pathlib.Path(circos_temp_fol).mkdir(parents=True, exist_ok=True)

    circos_conf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'template/circos.conf')
    ideogram_conf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'template/ideogram.conf')
    ticks_conf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'template/ticks.conf')

    ## copy Circos config files
    shutil.copy(circos_conf, circos_fol)
    shutil.copy(ideogram_conf, circos_fol)
    shutil.copy(ticks_conf, circos_fol)
    
    if unitcircos != 1000:
        circos_conf_path = "{}/circos.conf".format(circos_fol)
        with fileinput.FileInput(circos_conf_path, inplace=True) as mf:
            for line in mf:
                print(line.replace("chromosomes_units = 1000", "chromosomes_units = {}".format(unitcircos)), end='')

    ## prepare regions to plot
    df_identify = pd.read_csv(identify_path, sep='\t')
    l1 = df_identify['id'] == ec_id
    merge_region = df_identify[l1]['merge_region'].values[0]
    obj_ec = mergedRegion(merge_region)

    ## karyotype.conf
    write_path = "{}/karyotype.conf".format(circos_fol)
    list_range = []
    for obj_ in obj_ec.list_obj_region:
        list_range.append((obj_.mod_start, obj_.mod_end))
    str_karyotype = "chr - contig_1 contig_1 0 {} grey".format(list_range[-1][1])
    with open(write_path, 'w') as w_f:
        write_line = "{}\n".format(str_karyotype)
        w_f.write(write_line)

    ## chrom.tiles
    write_path = "{}/chrom.tiles".format(circos_data_fol)
    list_str_result = []
    for index, obj_ in enumerate(obj_ec.list_obj_region, start=1):
        if index > 1:
            mod_start = obj_.mod_start + 1
        else:
            mod_start = obj_.mod_start
        str_rel_region = "contig_1 {} {} color={}".format(mod_start, obj_.mod_end, obj_.chrom)
        list_str_result.append(str_rel_region)
    with open(write_path, 'w') as w_f:
        for str_result in list_str_result:
            write_line = "{}\n".format(str_result)
            w_f.write(write_line)

    ## text.labels.txt
    write_path = "{}/text.labels.txt".format(circos_data_fol)
    list_str_result = []
    for index, obj_ in enumerate(obj_ec.list_obj_region, start=1):
        if index > 1:
            mod_start = obj_.mod_start + 1
        else:
            mod_start = obj_.mod_start
        str_rel_region = "contig_1 {} {} {}:{}-{}:({})".format(mod_start, obj_.mod_end, obj_.chrom, obj_.start, obj_.end, obj_.strand)
        list_str_result.append(str_rel_region)
    with open(write_path, 'w') as w_f:
        for str_result in list_str_result:
            write_line = "{}\n".format(str_result)
            w_f.write(write_line)

    ## repeat.eccDNA.tiles
    write_path = "{}/repeat.eccDNA.tiles".format(circos_data_fol)
    repeat_annotate_path = "{}/repeat.annotate.txt".format(annotate_fol)
    df_repeat_annotate = pd.read_csv(repeat_annotate_path, sep='\t')
    l1 = df_repeat_annotate['id'] == ec_id
    l2 = df_repeat_annotate['ovl_cov'] > 0
    df_repeat_annotate = df_repeat_annotate[l1 & l2].copy()

    list_str_result = []
    if len(df_repeat_annotate) > 0:
        df_repeat_annotate['str_region'] = df_repeat_annotate.agg(lambda x: "{}_{}_{}".format(x['rep_chrom'], x['rep_start'], x['rep_end']), axis=1)

        list_uniq_region = list(set(df_repeat_annotate['str_region'].tolist()))

        merged_inverval = tree.IntervalTree()
        for region in list_uniq_region:
            list_ = obj_ec.cal_overlap_rel_region(region)
            for rel_region in list_:
                merged_inverval.add(tree.Interval(rel_region[0], rel_region[1]))

        merged_inverval.merge_overlaps()
        
        prev_end = None
        for interval in sorted(merged_inverval):
            if prev_end:
                if prev_end != interval[0]:
                    str_rel_region = "contig_1 {} {}".format(interval[0], interval[1])
                    list_str_result.append(str_rel_region)
                    prev_end = interval[1]
                else:
                    str_rel_region = "contig_1 {} {}".format(interval[0] + 1, interval[1])
                    list_str_result.append(str_rel_region)
                    prev_end = interval[1]
            else:
                str_rel_region = "contig_1 {} {}".format(interval[0], interval[1])
                list_str_result.append(str_rel_region)
                prev_end = interval[1]

    with open(write_path, 'w') as w_f:
        for str_result in list_str_result:
            write_line = "{}\n".format(str_result)
            w_f.write(write_line)

    ## cpg.eccDNA.tiles
    write_path = "{}/cpg.eccDNA.tiles".format(circos_data_fol)
    cpg_annotate_path = "{}/CpG.annotate.txt".format(annotate_fol)
    df_cpg_annotate = pd.read_csv(cpg_annotate_path, sep='\t')
    l1 = df_cpg_annotate['id'] == ec_id
    l2 = df_cpg_annotate['ovl_cov'] > 0
    df_cpg_annotate = df_cpg_annotate[l1 & l2].copy()

    list_str_result = []
    if len(df_cpg_annotate) > 0:
        df_cpg_annotate['str_region'] = df_cpg_annotate.agg(lambda x: "{}_{}_{}".format(x['cpg_chrom'], x['cpg_start'], x['cpg_end']), axis=1)

        list_uniq_region = list(set(df_cpg_annotate['str_region'].tolist()))

        merged_inverval = tree.IntervalTree()
        for region in list_uniq_region:
            list_ = obj_ec.cal_overlap_rel_region(region)
            for rel_region in list_:
                merged_inverval.add(tree.Interval(rel_region[0], rel_region[1]))

        merged_inverval.merge_overlaps()
            
        prev_end = None
        for interval in sorted(merged_inverval):
            if prev_end:
                if prev_end != interval[0]:
                    str_rel_region = "contig_1 {} {}".format(interval[0], interval[1])
                    list_str_result.append(str_rel_region)
                    prev_end = interval[1]
                else:
                    str_rel_region = "contig_1 {} {}".format(interval[0] + 1, interval[1])
                    list_str_result.append(str_rel_region)
                    prev_end = interval[1]
            else:
                str_rel_region = "contig_1 {} {}".format(interval[0], interval[1])
                list_str_result.append(str_rel_region)
                prev_end = interval[1]
            
    with open(write_path, 'w') as w_f:
        for str_result in list_str_result:
            write_line = "{}\n".format(str_result)
            w_f.write(write_line)

    ## snp_indel.eccDNA.tiles
    var_qual_cutoff = 20

    write_path = "{}/snp_indel.eccDNA.tiles".format(circos_data_fol)
    snp_annotate_path = "{}/variant.annotate.txt".format(annotate_fol)
    df_snp_annotate = pd.read_csv(snp_annotate_path, sep='\t')
    df_snp_annotate['quality'] = df_snp_annotate['quality'].apply(lambda x: x if x != '.' else '0.0')
    df_snp_annotate['quality'] = df_snp_annotate['quality'].astype(float)
    l1 = df_snp_annotate['id'] == ec_id
    l2 = df_snp_annotate['quality'] >= var_qual_cutoff
    df_snp_annotate = df_snp_annotate[l1 & l2].copy()

    list_str_result = []
    if len(df_snp_annotate) > 0:
        for index, value in df_snp_annotate.iterrows():
            str_rel_region = "contig_1 {} {} type={}".format(value['position'], min(int(value['position']) + 1, obj_ec.list_obj_region[-1].mod_end), 'snp' if (len(value['ref_base']) == 1 and len(value['alt_base']) == 1) else 'indels')
            list_str_result.append(str_rel_region)
            
    with open(write_path, 'w') as w_f:
        for str_result in list_str_result:
            write_line = "{}\n".format(str_result)
            w_f.write(write_line)

    ## gene.eccDNA.tile
    ## gene.eccDNA.exon.tiles
    ## gene.eccDNA.intron.tiles
    ## gene.eccDNA.intron.strand.tiles

    gene_annotate_path = "{}/gene.annotate.txt".format(annotate_fol)
    df_gene_annotate = pd.read_csv(gene_annotate_path, sep='\t')
    df_gene_annotate.fillna('unknown', inplace=True)
    l1 = df_gene_annotate['id'] == ec_id
    l2 = df_gene_annotate['strand'] != 'unknown'
    df_gene_annotate = df_gene_annotate[l1 & l2].copy()

    list_gene_result = []
    list_exon_result = []
    list_intron_result = []
    list_intron_strand_result = []
    for group_name, group in df_gene_annotate.groupby(by='name'):
        sorted_group = group.sort_values(by='match_start', ascending=True)
        chrom = sorted_group.iloc[0]['match_chrom']
        strand = sorted_group.iloc[0]['strand']
        str_region  = "{}_{}_{}".format(chrom, sorted_group.iloc[0]['match_start'], sorted_group.iloc[-1]['match_end'])
        list_rel_region = obj_ec.cal_overlap_rel_region(str_region)
        
        ## gene
        str_rel_region = "contig_1 {} {} {}".format(list_rel_region[0][0], list_rel_region[-1][1], group_name.replace('_','-'))
        list_gene_result.append(str_rel_region)
        
        ## exon
        l1 = group['region_type'].str.startswith('Intron')
        df_exon = group[~l1].copy()
        if len(df_exon) > 0:
            for exon_index, exon_value in df_exon.iterrows():
                exon_check_region = "{}_{}_{}".format(chrom, exon_value['match_start'], exon_value['match_end'])
                list_exon_region = obj_ec.cal_overlap_rel_region(exon_check_region)
                str_exon_region = "contig_1 {} {} type={}".format(list_exon_region[0][0], list_exon_region[0][1], 'exon' if exon_value['region_type'].startswith('Exon') else 'utr')
                list_exon_result.append(str_exon_region)
        
        ## intron
        str_rel_region = "contig_1 {} {}".format(list_rel_region[0][0], list_rel_region[-1][1])
        list_intron_result.append(str_rel_region)
        
        ## intron_strand
        rel_region_start = list_rel_region[0][0]
        rel_region_end = list_rel_region[-1][1]
        arrow_pos = rel_region_start
        arrow_shape = '>' if strand == '+' else '<'
        while ( arrow_pos + 300 ) < rel_region_end:
            arrow_pos += 300
            str_intron_strand = "contig_1 {} {} {}".format(arrow_pos, arrow_pos + 1, arrow_shape)
            list_intron_strand_result.append(str_intron_strand)

    write_path = "{}/gene.eccDNA.tiles".format(circos_data_fol)
    with open(write_path, 'w') as w_f:
        for str_result in list_gene_result:
            write_line = "{}\n".format(str_result)
            w_f.write(write_line)
            
    write_path = "{}/gene.eccDNA.exon.tiles".format(circos_data_fol)
    with open(write_path, 'w') as w_f:
        for str_result in list_exon_result:
            write_line = "{}\n".format(str_result)
            w_f.write(write_line)
            
    write_path = "{}/gene.eccDNA.intron.tiles".format(circos_data_fol)
    with open(write_path, 'w') as w_f:
        for str_result in list_intron_result:
            write_line = "{}\n".format(str_result)
            w_f.write(write_line)
            
    write_path = "{}/gene.eccDNA.intron.strand.tiles".format(circos_data_fol)
    with open(write_path, 'w') as w_f:
        for str_result in list_intron_strand_result:
            write_line = "{}\n".format(str_result)
            w_f.write(write_line)
    
    ct = datetime.datetime.now()
    print("[{}] preparing ctc/non-ctc reads".format(ct), flush=True)

    ## reads
    if all([os.path.isfile(reference_region_path), os.path.isfile(final_read_path)]):
        circos_reference_region_path = "{}/reference.fa".format(circos_temp_fol)
        circos_reference_size_path = "{}/reference.chrom.sizes".format(circos_temp_fol)
        
        ## reference
        chrom_limit = 0
        with open(circos_reference_region_path, 'w') as w_f:
            for record in SeqIO.parse(reference_region_path, 'fasta'):
                seq = str(record.seq)
                chrom_limit = len(seq)
                w_f.write('>contig_1\n{}\n'.format(str(record.seq)))
                
        with open(circos_reference_size_path, 'w') as w_f:
            w_f.write('contig_1\t{}\n'.format(chrom_limit))
        
        ## reads
        circos_final_read_path = "{}/reads.fa".format(circos_temp_fol)
        
        list_result = []
        for record in SeqIO.parse(final_read_path, 'fasta'):
            readid = record.id
            list_split_id = readid.split('_')[1:]
            trim_start = int(list_split_id[-2])
            trim_end = int(list_split_id[-1])
            comb_readid = '_'.join(list_split_id[:-2])
            seq = str(record.seq)
            tup_result = (readid, trim_start, trim_end, comb_readid, seq)
            list_result.append(tup_result)
        cols = ['readid','trim_start','trim_end','comb_readid','seq']
        df_readid_trim_info = pd.DataFrame(list_result, columns=cols)
        
        with open(circos_final_read_path, 'w') as w_f:
            for comb_readid, group in df_readid_trim_info.groupby(by='comb_readid'):
                seq = ''
                for index, value in group.sort_values(by='trim_start', ascending=True).iterrows():
                    seq += value['seq']
                write_line = ">{}\n{}\n".format(comb_readid, seq)
                w_f.write(write_line)

        ## minimap2 & igv & samtools index
        cmd = "minimap2 -t {} -ax map-ont reference.fa reads.fa | samtools view -bq 30 | samtools sort -o read.aligned.bam; igvtools count -w 50 -f max read.aligned.bam read.aligned.bam.wig reference.chrom.sizes; samtools index read.aligned.bam".format(threads)
        process = subprocess.call(cmd, shell=True, cwd=circos_temp_fol)

        wigFile = "{}/read.aligned.bam.wig".format(circos_temp_fol)
        genomeSizeFile = "{}/reference.chrom.sizes".format(circos_temp_fol)
        out_file = "{}/read.aligned.bam.wig.fill_ends.log10.bdg".format(circos_temp_fol)

        bdg_data = wigToBedGraph(wigFile, genomeSizeFile)

        if len(bdg_data) > 0:
            bdg_data = [[bdg_data[0][0], bdg_data[0][1], bdg_data[0][1], bdg_data[0][3]]] + bdg_data + [[bdg_data[-1][0], bdg_data[-1][2], bdg_data[-1][2], bdg_data[-1][3]]]
            
        df_bdg_data = pd.DataFrame(bdg_data)
        df_bdg_data[1] = df_bdg_data[1].astype(int)
        df_bdg_data[2] = df_bdg_data[2].astype(int)
        df_bdg_data[3] = df_bdg_data[3].astype(float)
        df_bdg_data[4] = df_bdg_data[3].apply(lambda x: math.log10(x if x >= 1.0 else 1.0))
        df_bdg_data = df_bdg_data[df_bdg_data[1] <= chrom_limit].copy()
        df_bdg_data[2] = df_bdg_data[2].apply(lambda x: x if x <= chrom_limit else chrom_limit)
        df_bdg_data.drop(columns=[3], axis=1, inplace=True)
        df_bdg_data.to_csv(out_file, sep='\t', header=None, index=None)

        ## generate 20-depth bed file
        cmd = "minimap2 -t {} -ax map-ont reference.fa reads.fa | samtools view -bq 30 | samtools sort -O BAM | bedtools bamtobed -i stdin > depth20.bed".format(threads)
        process = subprocess.call(cmd, shell=True, cwd=circos_temp_fol)

        depth20bed_path = "{}/depth20.bed".format(circos_temp_fol)
        out_file = "{}/reads.depth20.tiles".format(circos_temp_fol)
        out_closegap = "{}/depth20.close_gaps.txt".format(circos_temp_fol)

        cols = ['chrom', 'start', 'end', 'readid', 'mapq', 'strand']
        df_depth20bed = pd.read_csv(depth20bed_path, sep='\t', header=None, names=cols)
        df_depth20bed['end'] = df_depth20bed['end'].apply(lambda x: x if x <= chrom_limit else chrom_limit)
        df_group_depth20bed = df_depth20bed.groupby(by='readid').size().reset_index(name='count')
        l1 = df_group_depth20bed['count'] > 1
        df_breakpoint_group_depth20bed = df_group_depth20bed[l1]

        max_viz_reads = 1000
        if modecircos == 'proportion':

            ct = datetime.datetime.now()
            print("[{}] sampling reads : proportion".format(ct), flush=True)

            list_breakpoint_readid = []
            if len(df_breakpoint_group_depth20bed) > 0:
                set_breakpoint_readid = set(df_breakpoint_group_depth20bed['readid'].tolist())

                list_breakpoint_readid = list(set_breakpoint_readid)

                len_break = len(set_breakpoint_readid)
                len_non_break = len(df_group_depth20bed) - len_break
                if len_non_break > 0:
                    if len_break < len_non_break:
                        prop_break = 1
                        prop_non_break = round(len_non_break / len_break)
                    else:
                        prop_non_break = 1
                        prop_break = round(len_break / len_non_break)

                    mult = round(max_viz_reads / ( prop_break + prop_non_break ))

                    # gap close
                    l1 = df_depth20bed['readid'].isin(set_breakpoint_readid)
                    df_reads_depth20_tiles_close = df_depth20bed[l1].copy()
                    df_reads_depth20_tiles_close['type'] = df_reads_depth20_tiles_close['readid'].apply(lambda x: "gap=close,id={}".format(x))
                    inv_num = max(1, round(len(df_reads_depth20_tiles_close) / (prop_break * mult)))

                    df_out = df_reads_depth20_tiles_close.iloc[::inv_num]
                    selected_cols = ['chrom','start','end','type']
                    df_out[selected_cols].to_csv(out_file, sep=' ', mode='a', header=None, index=None)

                    # gap none
                    l1 = df_depth20bed['readid'].isin(set_breakpoint_readid)
                    df_reads_depth20_tiles_none = df_depth20bed[~l1].copy()
                    df_reads_depth20_tiles_none['type'] = df_reads_depth20_tiles_none['readid'].apply(lambda x: "gap=none,id={}".format(x))
                    inv_num = max(1, round(len(df_reads_depth20_tiles_none) / (prop_non_break * mult)))
                    
                    df_out = df_reads_depth20_tiles_none.iloc[::inv_num]
                    selected_cols = ['chrom','start','end','type']
                    df_out[selected_cols].to_csv(out_file, sep=' ', mode='a', header=None, index=None)

                else:
                    df_reads_depth20_tiles_close = df_depth20bed.copy()
                    df_reads_depth20_tiles_close['type'] = df_reads_depth20_tiles_close['readid'].apply(lambda x: "gap=close,id={}".format(x))
                    inv_num = max(1, round(len(df_reads_depth20_tiles_close) / max_viz_reads))
                    df_out = df_reads_depth20_tiles_close.iloc[::inv_num]
                    selected_cols = ['chrom','start','end','type']
                    df_out[selected_cols].to_csv(out_file, sep=' ', header=None, index=None)

            else:
                # all gap none
                df_reads_depth20_tiles_none = df_depth20bed.copy()
                df_reads_depth20_tiles_none['type'] = df_reads_depth20_tiles_none['readid'].apply(lambda x: "gap=none,id={}".format(x))
                inv_num = max(1, round(len(df_reads_depth20_tiles_none) / max_viz_reads))
                df_out = df_reads_depth20_tiles_none.iloc[::inv_num]
                selected_cols = ['chrom','start','end','type']
                df_out[selected_cols].to_csv(out_file, sep=' ', header=None, index=None)

            with open(out_closegap, 'w') as w_f:
                for readid in list_breakpoint_readid:
                    w_f.write(str(readid) + '\n')
        else:

            ct = datetime.datetime.now()
            print("[{}] sampling reads : ctc".format(ct), flush=True)

            list_breakpoint_readid = []
            if len(df_breakpoint_group_depth20bed) > 0:
                set_breakpoint_readid = set(df_breakpoint_group_depth20bed['readid'].tolist())
                list_breakpoint_readid = list(set_breakpoint_readid)
                list_dfs = []

                # gap close
                l1 = df_depth20bed['readid'].isin(set_breakpoint_readid)
                df_reads_depth20_tiles_close = df_depth20bed[l1].copy()
                df_reads_depth20_tiles_close['type'] = df_reads_depth20_tiles_close['readid'].apply(lambda x: "gap=close,id={}".format(x))
                if len(df_reads_depth20_tiles_close) > max_viz_reads:
                    inv_num = max(1, round(len(df_reads_depth20_tiles_close) / max_viz_reads))
                    df_out = df_reads_depth20_tiles_close.iloc[::inv_num]
                    selected_cols = ['chrom','start','end','type']
                    df_out[selected_cols].to_csv(out_file, sep=' ', header=None, index=None)

                else:
                    list_dfs.append(df_reads_depth20_tiles_close)

                    # gap none
                    l1 = df_depth20bed['readid'].isin(set_breakpoint_readid)
                    df_reads_depth20_tiles_none = df_depth20bed[~l1].copy()
                    df_reads_depth20_tiles_none['type'] = df_reads_depth20_tiles_none['readid'].apply(lambda x: "gap=none,id={}".format(x))
                    list_dfs.append(df_reads_depth20_tiles_none)

                    df_concat = pd.concat(list_dfs, ignore_index=True)
                    selected_cols = ['chrom','start','end','type']
                    df_concat[selected_cols].head(max_viz_reads).to_csv(out_file, sep=' ', header=None, index=None)

            else:
                # all gap none
                df_reads_depth20_tiles_none = df_depth20bed.copy()
                df_reads_depth20_tiles_none['type'] = df_reads_depth20_tiles_none['readid'].apply(lambda x: "gap=none,id={}".format(x))
                inv_num = max(1, round(len(df_reads_depth20_tiles_none) / max_viz_reads))
                df_out = df_reads_depth20_tiles_none.iloc[::inv_num]
                selected_cols = ['chrom','start','end','type']
                df_out[selected_cols].to_csv(out_file, sep=' ', header=None, index=None)

            with open(out_closegap, 'w') as w_f:
                for readid in list_breakpoint_readid:
                    write_line = "{}\n".format(readid)
                    w_f.write(write_line)

        ct = datetime.datetime.now()
        print("[{}] filtering reads for visualizing".format(ct), flush=True)

        ## filter bam file
        fbam = "{}/read.aligned.bam".format(circos_temp_fol)
        fname = "{}/depth20.close_gaps.txt".format(circos_temp_fol)
        bamfilterByReadid(fbam, fname)

        ## make wig file
        cmd = "igvtools count -w 50 -f max depth20.close_gaps.txt.bam read.close_gaps.aligned.bam.wig reference.chrom.sizes"
        process = subprocess.call(cmd, shell=True, cwd=circos_temp_fol)

        wigFile = "{}/read.close_gaps.aligned.bam.wig".format(circos_temp_fol)
        genomeSizeFile = "{}/reference.chrom.sizes".format(circos_temp_fol)
        out_file = "{}/read.close_gaps.aligned.bam.wig.fill_ends.log10.bdg".format(circos_temp_fol)

        bdg_data = wigToBedGraph(wigFile, genomeSizeFile)

        if len(bdg_data) > 0:
            bdg_data = [[bdg_data[0][0], bdg_data[0][1], bdg_data[0][1], bdg_data[0][3]]] + bdg_data + [[bdg_data[-1][0], bdg_data[-1][2], bdg_data[-1][2], bdg_data[-1][3]]]
            
        df_bdg_data = pd.DataFrame(bdg_data)
        df_bdg_data[1] = df_bdg_data[1].astype(int)
        df_bdg_data[2] = df_bdg_data[2].astype(int)
        df_bdg_data[3] = df_bdg_data[3].astype(float)
        df_bdg_data[4] = df_bdg_data[3].apply(lambda x: math.log10(x if x >= 1.0 else 1.0))
        df_bdg_data = df_bdg_data[df_bdg_data[1] <= chrom_limit].copy()
        df_bdg_data[2] = df_bdg_data[2].apply(lambda x: x if x <= chrom_limit else chrom_limit)
        df_bdg_data.drop(columns=[3], axis=1, inplace=True)
        df_bdg_data.to_csv(out_file, sep='\t', header=None, index=None)

        depth20_tiles_path = "{}/reads.depth20.tiles".format(circos_temp_fol)
        read_bdg_path = "{}/read.aligned.bam.wig.fill_ends.log10.bdg".format(circos_temp_fol)
        read_bdg_close_gaps_path = "{}/read.close_gaps.aligned.bam.wig.fill_ends.log10.bdg".format(circos_temp_fol)

        shutil.copy(depth20_tiles_path, circos_data_fol)
        shutil.copy(read_bdg_path, circos_data_fol)
        shutil.copy(read_bdg_close_gaps_path, circos_data_fol)
    else:
        depth20_tiles_path = "{}/reads.depth20.tiles".format(circos_data_fol)
        read_bdg_path = "{}/read.aligned.bam.wig.fill_ends.log10.bdg".format(circos_data_fol)
        read_bdg_close_gaps_path = "{}/read.close_gaps.aligned.bam.wig.fill_ends.log10.bdg".format(circos_data_fol)

        pathlib.Path(depth20_tiles_path).touch()
        pathlib.Path(read_bdg_path).touch()
        pathlib.Path(read_bdg_close_gaps_path).touch()

    ct = datetime.datetime.now()
    print("[{}] finished visualizing process : {}\n".format(ct, ec_id), flush=True)