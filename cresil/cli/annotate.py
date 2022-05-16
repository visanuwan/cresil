#!/usr/bin/env python

import gzip
import pathlib
import datetime
import pybedtools
import os, sys
import pandas as pd
from itertools import chain
from multiprocessing import cpu_count, Pool
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def loadRmsk(fname, extension):
    if extension == 'gzip':
        rmsk = pd.read_csv(fname, compression='gzip', header=None, sep='\t', quotechar='"', on_bad_lines='skip')
    else:
        rmsk = pd.read_csv(fname, header=None, sep='\t', quotechar='"', on_bad_lines='skip')
    rmskbed = rmsk.loc[:,[5,6,7,10,11,12,]]
    rmskbed.loc[:,6] = rmskbed.loc[:,6]
    rmskbed.loc[:,13] = rmskbed.loc[:,5]+"_"+rmskbed.loc[:,6].astype(str)+"_" \
                        +(rmskbed.loc[:,7]-rmskbed.loc[:,6]).astype(str)
    rmskbed.loc[:,14] = rmskbed.loc[:,7]-rmskbed.loc[:,6]
    rmskbed.columns = ['chrom','start','end','repName','repClass','repFamily','repid','size']

    return rmskbed

def loadCpG(fname, extension):
    if extension == 'gzip':
        cpgbed = pd.read_csv(fname, compression='gzip', header=None, sep='\t', quotechar='"', on_bad_lines='skip')
    else:
        cpgbed = pd.read_csv(fname, header=None, sep='\t', quotechar='"', on_bad_lines='skip')
    cpgbed.columns = ['chrom', 'start', 'end', 'cpg_id']
    cpgbed['cpg_size'] = cpgbed['end'] - cpgbed['start']

    return cpgbed
        
def readEccData(fname):
    ecc = pd.read_csv(fname, sep='\t')
    if len(ecc) == 0:
        return []
    
    list_result = []
    for idx, value in ecc.iterrows():
        merge_region = value['merge_region']
        merge_len = value['merge_len']
        eccdna_status = value['eccdna_status']
        id_ = value['id']

        for mrow in merge_region.split(","):
            chrom = mrow.split(':')[0]
            start = int(mrow.split(':')[1].split('-')[0]) - 1
            end = int(mrow.split(':')[1].split('-')[1].split('_')[0])
            size = end - start
            freq = "%.2f"%(float(size)/float(merge_len))
            
            
            tup_result = (chrom, start, end, id_, size, merge_len, eccdna_status, freq)
            list_result.append(tup_result)

    header = ['chrom','start','end','id','region_size','merge_len','eccdna_status','freq']
    eccbed = pd.DataFrame(list_result, columns=header)

    return eccbed

def trim_exon_intron_regions(df_):
    r_start = df_['r_start']
    r_end = df_['r_end']
    q_start = df_['q_start']
    q_end = df_['q_end']
    
    b_start = r_start
    if r_start < q_start:
        b_start = q_start
    b_end = q_end
    if r_end < q_end:
        b_end = r_end
        
    df_['b_start'] = b_start
    df_['b_end'] = b_end
    
    return df_

def split_detail(df_):
    r_name = df_['r_name']
    name, region_type, strand = r_name.split('|')
    df_['name'] = name
    df_['region_type'] = region_type
    df_['strand'] = strand
    
    return df_

def get_match_intron_exon(df_overlap_temp, query_start, query_end, query):

    list_dfs = []

    for idx, value in df_overlap_temp.iterrows():

        chrom = value['chrom']
        chromStart = value['chromStart']
        chromEnd = value['chromEnd']
        name = value['name']
        strand = value['strand']
        thickStart = value['thickStart']
        thickEnd = value['thickEnd']
        blockCount = value['blockCount']
        blockSize = value['blockSize']
        blockStart = value['blockStart']

        if blockSize[-1] == ',':
            blockSize = blockSize[:-1]
        if blockStart[-1] == ',':
            blockStart = blockStart[:-1]

        ## EXON coordinates

        list_exon_block = list(zip(list(map(int, blockSize.split(','))), list(map(int, blockStart.split(',')))))
        len_exon_block = len(list_exon_block)

        list_exon_intron = []

        for counter, tup_ in enumerate(list_exon_block, start=1):
            start_ = chromStart + tup_[1]
            end_ = start_ + tup_[0]

            if strand == '+':
                name_this = "{}|Exon_{}|{}".format(name, counter, strand)
            else:
                name_this = "{}|Exon_{}|{}".format(name, len_exon_block - (counter - 1), strand)
            

            tup_result = (chrom, start_, end_, name_this)
            list_exon_intron.append(tup_result)

        list_temp_start = []
        list_temp_end = []
        for tup_ in list_exon_intron:
            list_temp_start.append(tup_[1])
            list_temp_end.append(tup_[2])

        ## INTRON coordinates
            
        list_intron_block = list(zip(list_temp_end[:-1], list_temp_start[1:]))
        len_intron_block = len(list_intron_block)

        for counter, tup_ in enumerate(list_intron_block, start=1):
            start_ = tup_[0] + 1
            end_ = tup_[1] - 1
            
            if start_ > end_:
                start_ = start_ - 1
                end_ = start_

            if strand == '+':
                name_this = "{}|Intron_{}|{}".format(name, counter, strand)
            else:
                name_this = "{}|Intron_{}|{}".format(name, len_intron_block - (counter - 1), strand)

            tup_result = (chrom, start_, end_, name_this)
            list_exon_intron.append(tup_result)

        ## UTR cooridinates

        if strand == '+':
            tup_utr5 = (chrom, chromStart, thickStart, "{}|UTR5|{}".format(name, strand))
            tup_utr3 = (chrom, thickEnd, chromEnd, "{}|UTR3|{}".format(name, strand))
        else:
            tup_utr5 = (chrom, thickEnd, chromEnd, "{}|UTR5|{}".format(name, strand))
            tup_utr3 = (chrom, chromStart, thickStart, "{}|UTR3|{}".format(name, strand))

        if thickStart != thickEnd:
            if tup_utr5[1] != tup_utr5[2]:
                list_exon_intron.append(tup_utr5)
            if tup_utr3[1] != tup_utr3[2]:
                list_exon_intron.append(tup_utr3)

        ## combine gene coordinates

        list_cols = ['chrom', 'chromStart', 'chromEnd', 'name_this']
        df_this_target = pd.DataFrame(list_exon_intron, columns=list_cols)
        
        this_target  = pybedtools.BedTool.from_dataframe(df_this_target)

        list_cols = ['r_chrom','r_start','r_end','r_name','q_chrom','q_start','q_end','match_size']
        df_overlap_this_temp = pybedtools.BedTool.to_dataframe(this_target.intersect(query, wo=True ), disable_auto_names=True, header=None, names=list_cols)
        
        if len(df_overlap_this_temp) > 0:
            df_overlap_this_temp = df_overlap_this_temp.sort_values(by=['r_start', 'r_end'])
            df_overlap_this_temp = df_overlap_this_temp.apply(trim_exon_intron_regions, axis=1)
            list_dfs.append(df_overlap_this_temp)

    ## combine all overlapping gene/exons/introns
    df_return = []
    if len(list_dfs) > 0:
        df_return  = pd.concat(list_dfs, ignore_index=True)
        df_return  = df_return.apply(split_detail, axis=1)
        df_return.drop(columns = ['r_name'], inplace=True)
        return df_return
    else:
        return df_return

def get_gene_annotate(tup_record):

    list_result = []
        
    chrom, chromStart, chromEnd, id_, seqSize, eccSize, eccdna_status, freq = tup_record

    try:
        df_test_query = pd.DataFrame([(chrom, chromStart, chromEnd)], columns=['chrom', 'chromStart', 'chromEnd'])

        query  = pybedtools.BedTool.from_dataframe(df_test_query)

        df_overlap_temp = pd.read_table(target.intersect(query, u=False, wa=True).fn, header=None)
        df_overlap_temp = df_overlap_temp[~df_overlap_temp.isnull().any(axis=1)]
        df_overlap_temp.reset_index(drop=True, inplace=True)
        
        if len(df_overlap_temp) > 0:
            list_cols = ['chrom','chromStart','chromEnd','name','score','strand','thickStart', 'thickEnd','itemRgb','blockCount', 'blockSize','blockStart']
            df_overlap_temp.columns = list_cols

            df_ovl_exon_intron = get_match_intron_exon(df_overlap_temp, chromStart, chromEnd, query)

            if len(df_ovl_exon_intron) > 0:
                for idx, value in df_ovl_exon_intron.iterrows():
                    
                    match_size = value['match_size']
                    r_chrom = value['r_chrom']
                    b_start = value['b_start']
                    b_end = value['b_end']
                    name = value['name']
                    region_type = value['region_type']
                    strand = value['strand']
                    
                    tup_result = (chrom, chromStart, chromEnd, id_, seqSize, eccSize, eccdna_status, freq, match_size, r_chrom, b_start, b_end, name, region_type, strand)
                    list_result.append(tup_result)
            else:
                tup_result = (chrom, chromStart, chromEnd, id_, seqSize, eccSize, eccdna_status, freq, 0, '', 0, 0, '', '', '')
                list_result.append(tup_result)
        else:
            tup_result = (chrom, chromStart, chromEnd, id_, seqSize, eccSize, eccdna_status, freq, 0, '', 0, 0, '', '', '')
            list_result.append(tup_result)

    except pd.errors.EmptyDataError:
        pass

    return list_result

def argparser():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='Annotate identified eccDNA')
    general = parser.add_argument_group(title='General options')
    general.add_argument('-t', "--threads",
                            help="Number of threads [all CPU cores]",
                            type=int, default=0)
    general.add_argument('-rp', "--repeat-bed", dest= 'rmsk_bed',
                            help="repeat masker bed file",
                            type=str, default=None)
    general.add_argument('-cg', "--cpg-bed", dest='cpg_bed',
                            help="CpG islands bed file",
                            type=str, default=None)
    general.add_argument('-gb', "--gene-bed", dest='gene_bed',
                            help="gene annotation bed file",
                            type=str, default=None)
    general.add_argument('-identify', "--identify-input", dest='identifyinput',
                            help="identified eccDNA table",
                            type=str, default=None)
    return parser

def main(args):

    if not args.rmsk_bed:
        sys.stderr.write("[WARNING] a repeat masker bed file cannot be found (continue without it)\n")

    if not args.cpg_bed:
        sys.stderr.write("[WARNING] a CpG islands bed file cannot be found (continue without it)\n")

    if not args.gene_bed:
        sys.stderr.write("[WARNING] a gene annotation bed file cannot be found (continue without it)\n")

    if not args.identifyinput:
        sys.stderr.write("[ABORT] an identified eccDNA table is needed\n")
        exit(1)

    threads = cpu_count() if args.threads == 0 else args.threads
    rmsk_bed = args.rmsk_bed
    cpg_bed = args.cpg_bed
    gene_bed = args.gene_bed
    fname = args.identifyinput
    bname = str(pathlib.Path(fname).parent)
    assemGraph = "{}/cresil_run/assemGraph".format(bname)
    outDir="{}/cresil_gAnnotation".format(bname)

    #####################################################################
    ## output ###########################################################

    pathlib.Path(outDir).mkdir(parents=True, exist_ok=True)

    ## output ###########################################################

    ## load eccDNA table
    eccbed = readEccData(fname)

    ct = datetime.datetime.now()
    print("\n######### CReSIL : Start annotating process (thread : {})\n[{}] preparing data for eccDNA regions".format(threads, ct), flush=True)

    if rmsk_bed:
        ## prepare data
        extension = rmsk_bed.split('.')[-1]
        rmskbed = loadRmsk(rmsk_bed, extension)
        rmskbed_bt = pybedtools.BedTool.from_dataframe(rmskbed)
        eccbed_bt = pybedtools.BedTool.from_dataframe(eccbed)

        if len(eccbed) == 0:
            sys.exit("[ABORT] zero identified eccDNA\n")

        ct = datetime.datetime.now()
        print("[{}] finished preparing data : {}".format(ct, len(eccbed)), flush=True)

        ## repeat ###########################################################

        ct = datetime.datetime.now()
        print("[{}] calculating repeats".format(ct), flush=True)

        try:
            eccrmsk_intersect = pybedtools.BedTool.to_dataframe(eccbed_bt.intersect(rmskbed_bt, wao=True ), disable_auto_names=True, header=None, dtype={15: str}) 
            eccrmsk_intersect[9] = eccrmsk_intersect[9].apply(lambda x: 0 if x < 0 else x)
            eccrmsk_intersect[10] = eccrmsk_intersect[10].apply(lambda x: 0 if x < 0 else x)
            eccrmsk_intersect[15] = eccrmsk_intersect[15].apply(lambda x: 0 if x == '.' else int(x))
            eccrmsk_intersect[17] = eccrmsk_intersect[16].divide(eccrmsk_intersect[15].astype(int))
            eccrmsk_intersect[17] = eccrmsk_intersect[17].apply(lambda x: round(x, 2))
            eccrmsk_intersect[17] = eccrmsk_intersect[17].fillna(value = 0)
            eccrmsk_intersect.replace('.','', inplace=True)
            eccrmsk_intersect.columns = ['chrom','start','end','id','region_size','merge_len','eccdna_status','freq',
                                        'rep_chrom','rep_start','rep_end','rep_name','rep_class','rep_family','rep_id',
                                        'rep_size','rep_ovl_ecc','ovl_cov']
            
            write_path = "{}/repeat.annotate.txt".format(outDir)
            eccrmsk_intersect.to_csv(write_path, sep='\t', index=None)

        except pd.errors.EmptyDataError:
            ct = datetime.datetime.now()
            print("[{}] no repeat overlapping eccDNA".format(ct), flush=True)
            pass
    else:
        cols = ['chrom','start','end','id','region_size','merge_len','eccdna_status','freq',
                'rep_chrom','rep_start','rep_end','rep_name','rep_class','rep_family',
                'rep_id','rep_size','rep_ovl_ecc','ovl_cov']
        eccrmsk_intersect = pd.DataFrame([], columns=cols)
        write_path = "{}/repeat.annotate.txt".format(outDir)
        eccrmsk_intersect.to_csv(write_path, sep='\t', index=None)


    ## CpG island ###########################################################

    ct = datetime.datetime.now()
    print("[{}] calculating CpG islands".format(ct), flush=True)

    if cpg_bed:
        ## prepare data
        extension = cpg_bed.split('.')[-1]
        cpgbed = loadCpG(cpg_bed, extension)
        cpgbed_bt = pybedtools.BedTool.from_dataframe(cpgbed)

        try:
            ecccpg_intersect = pybedtools.BedTool.to_dataframe(eccbed_bt.intersect(cpgbed_bt, wao=True ), disable_auto_names=True, header=None)
            ecccpg_intersect[9] = ecccpg_intersect[9].apply(lambda x: 0 if x < 0 else x)
            ecccpg_intersect[10] = ecccpg_intersect[10].apply(lambda x: 0 if x < 0 else x)
            ecccpg_intersect[12] = ecccpg_intersect[12].apply(lambda x: 0 if x < 0 else x)
            ecccpg_intersect[14] = ecccpg_intersect[13].divide(ecccpg_intersect[12].astype(int))
            ecccpg_intersect[14] = ecccpg_intersect[14].apply(lambda x: round(x, 2))
            ecccpg_intersect[14] = ecccpg_intersect[14].fillna(value = 0)
            ecccpg_intersect.replace('.','', inplace=True)
            ecccpg_intersect.columns = ['chrom', 'start', 'end', 'id', 'region_size', 'merge_len', 'eccdna_status', 'freq',
                                        'cpg_chrom', 'cpg_start', 'cpg_end', 'cpg_id', 'cpg_size',
                                        'cpg_ovl_ecc', 'ovl_cov']
            
            write_path = "{}/CpG.annotate.txt".format(outDir)
            ecccpg_intersect.to_csv(write_path, sep='\t', index=None)

        except pd.errors.EmptyDataError:
            ct = datetime.datetime.now()
            print("[{}] no CpG islands overlapping eccDNA".format(ct), flush=True)
            pass
    else:
        cols = ['chrom', 'start', 'end', 'id', 'region_size', 'merge_len', 'eccdna_status', 'freq',
                'cpg_chrom', 'cpg_start', 'cpg_end', 'cpg_id', 'cpg_size',
                'cpg_ovl_ecc', 'ovl_cov']
        ecccpg_intersect = pd.DataFrame([], columns=cols)
        write_path = "{}/CpG.annotate.txt".format(outDir)
        ecccpg_intersect.to_csv(write_path, sep='\t', index=None)

    ## gene/exons/introns ###########################################################

    ct = datetime.datetime.now()
    print("[{}] calculating gene/exons/introns".format(ct), flush=True)

    if gene_bed:
        ## prepare data
        global target
        target = pybedtools.BedTool(gene_bed)

        write_path = "{}/gene.annotate.txt".format(outDir)

        list_record = []
        iter_ = 1000
        flag = 0

        for counter, (index, value) in enumerate(eccbed.iterrows(), start=1):
                
            if counter % iter_ == 0:
                ct = datetime.datetime.now()
                print("[{}] {}".format(ct, counter), flush=True)

            chrom = value['chrom']
            chromStart = value['start']
            chromEnd = value['end']
            id_ = value['id']
            seqSize = value['region_size']
            eccSize = value['merge_len']
            eccdna_status = value['eccdna_status']
            freq = value['freq']
            
            if counter % iter_ == 0:
                
                tup_record = (chrom, chromStart, chromEnd, id_, seqSize, eccSize, eccdna_status, freq)
                list_record.append(tup_record)
                
                pool = Pool(processes = threads)
                list_pool_result = pool.map(get_gene_annotate, list_record)
                pool.close()
                pool.join()
                
                list_result = list(chain.from_iterable(list_pool_result))
                
                list_cols = ['chrom', 'start', 'end', 'id', 'region_size', 'merge_len', 'eccdna_status', 'freq', 'match_size', 'match_chrom',
                            'match_start', 'match_end', 'name', 'region_type', 'strand']
                
                df_eccdna_exon_intron = pd.DataFrame(list_result, columns=list_cols)
                
                if flag == 0:
                    df_eccdna_exon_intron.to_csv(write_path, sep='\t', mode='w', header=True, index=None)
                    flag = 1
                else:
                    df_eccdna_exon_intron.to_csv(write_path, sep='\t', mode='a', header=False, index=None)
                    
                list_record = []
            else:
                tup_record = (chrom, chromStart, chromEnd, id_, seqSize, eccSize, eccdna_status, freq)
                list_record.append(tup_record)

        if len(list_record) > 0:

            pool = Pool(processes = threads)
            list_pool_result = pool.map(get_gene_annotate, list_record)
            pool.close()
            pool.join()
            
            list_result = list(chain.from_iterable(list_pool_result))
            
            cols = ['chrom', 'start', 'end', 'id', 'region_size', 'merge_len', 'eccdna_status', 'freq', 'match_size', 'match_chrom',
                    'match_start', 'match_end', 'name', 'region_type', 'strand']

            df_eccdna_exon_intron = pd.DataFrame(list_result, columns=cols)

            if os.path.exists(write_path):
                df_eccdna_exon_intron.to_csv(write_path, sep='\t', mode='a', header=False, index=None)
            else:
                df_eccdna_exon_intron.to_csv(write_path, sep='\t', mode='w', header=True, index=None)
    else:
        cols = ['chrom', 'start', 'end', 'id', 'region_size', 'merge_len', 'eccdna_status', 'freq', 'match_size', 'match_chrom',
                'match_start', 'match_end', 'name', 'region_type', 'strand']
        df_eccdna_exon_intron = pd.DataFrame([], columns=cols)
        write_path = "{}/gene.annotate.txt".format(outDir)
        df_eccdna_exon_intron.to_csv(write_path, sep='\t', index=None)

    ## variant ###########################################################

    ct = datetime.datetime.now()
    print("[{}] combining variants".format(ct), flush=True)

    ## prepare data
    df_identify = pd.read_csv(fname, sep='\t')
    list_id = df_identify['id'].tolist()

    write_path = "{}/variant.annotate.txt".format(outDir)
    with open(write_path, 'w') as w_f:
        w_f.write('id\tmerge_len\tnum_region\teccdna_status\tposition\tidentifier\tref_base\talt_base\tquality\tfilter\tinfo\tformat\tsample\n')

    for id_ in list_id:

        merge_len = df_identify[df_identify['id'] == id_]['merge_len'].values[0]
        num_region = df_identify[df_identify['id'] == id_]['num_region'].values[0]
        eccdna_status = df_identify[df_identify['id'] == id_]['eccdna_status'].values[0]

        variant_path = "{}/{}/{}_variant.vcf".format(assemGraph, id_, id_)
        if os.path.isfile(variant_path):
            list_to_write = []
            with open(variant_path, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        temp_line = [id_, str(merge_len), str(num_region), eccdna_status] + line.strip().split('\t')[1:]
                        write_line = "\t".join(temp_line)
                        list_to_write.append(write_line)
            if len(list_to_write) == 0:
                list_to_write = ["{}\t{}\t{}\t{}\t.\t.\t.\t.\t.\t.\t.\t.\t.".format(id_, merge_len, num_region, eccdna_status)]

            with open(write_path, 'a') as w_f:
                for line_ in list_to_write:
                    w_f.write("{}\n".format(line_))
        else:
            print("[WARNING] file not found : {}".format(variant_path), flush=True)

    ct = datetime.datetime.now()
    print("[{}] finished annotating process\n".format(ct), flush=True)

