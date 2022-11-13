import appscript as appscript
import pandas as pd
import xlrd
import yaml
from phenoapt import PhenoApt
import numpy as np
import subprocess
import re
import json
from suds import client

import urllib.request as urlre
from functools import reduce
import requests, sys
from openpyxl import load_workbook
from tqdm import tqdm
import appscript
import time
import allel
from pandas.testing import assert_frame_equal
import networkx
import obonet

from script.file_prepare_script import dir_of_Filtered_TSV_MUltiple_source
from script.process_file_script import get_case_id_file

import sys,os
custom_module_path = os.path.dirname(os.path.dirname (os.path.abspath ('__file__')))
print(f'Adding extra dir to sys.path: {custom_module_path}')
sys.path.append(custom_module_path)
print(sys.path)

from phrank.phrank import Phrank
# import phrank

def output_columns(hpo,tools,intersect,REVEL_thresh,CADD_thresh,collect_variants_tofile,collect_organ):
    columns_fill = ['CaseID', 'hpo', 'hpo_id_input', 'Symbol']
    for tool in tools:
        columns_fill = columns_fill + [f'{tool}_rank']
        if intersect:
            columns_fill = columns_fill + [f'{tool}_rank_intersect']
        if len(REVEL_thresh) != 0:
            for thresh in REVEL_thresh:
                columns_fill = columns_fill + [f'{tool}_rank_revel_{thresh}']
        if len(CADD_thresh) != 0:
            for thresh in CADD_thresh:
                columns_fill = columns_fill + [f'{tool}_rank_cadd_{thresh}']
        if tool == 'phenoapt' and hpo == 'Weight':
            for k in ['all_hpo_plus_weight', 'only_weight_hpo_weight', 'only_weight_hpo_no_weight']:
                columns_fill = columns_fill + [f'{k}_phenoapt_rank']
                if intersect: columns_fill = columns_fill + [f'{k}_phenoapt_rank_intersect']
                if len(REVEL_thresh) != 0:
                    for thresh in REVEL_thresh:
                        columns_fill = columns_fill + [f'{k}_phenoapt_rank_revel_{thresh}']
                if len(CADD_thresh) != 0:
                    for thresh in CADD_thresh:
                        columns_fill = columns_fill + [f'{k}_phenoapt_rank_cadd_{thresh}']

    if collect_variants_tofile:
        columns_fill = columns_fill + ['Frequency','Variant_class','Inheritance_ADAR']
    if collect_organ:
        columns_fill = columns_fill+[f'{hpo}_organ_system',f'{hpo}_organ_system_number']
    return columns_fill


def get_phenomizer_rank(correct_gene, searchlist):
    correct_rank = {}
    time = 0
    for k in range(len(searchlist)):
        if f'{correct_gene} (' in str(searchlist[k]):
            correct_rank[time] = k
            time = time + 1
    if time == 0:
        rank = 'NA'
        rr = 0
        rp = 'NA'
    else:
        rank = correct_rank[0]
        rr = 1 / (rank + 1)
        rp = (rank+1)/len(searchlist)
    return rank, rr , rp

def load_txt_list(filedir):
    with open(filedir) as f:
        list = (f.read().replace('\n','')).split(',')
    return list

def process_rank_result_table(tool,result_mode,result_list,case_id,i,pathogenic_gene_dict,pwd,filename,intersect,REVEL_thresh,CADD_thresh,part_table,part_rr_table,part_rp_table):
    txt_dir_dict = {'Symbol':'candidategene','Ensembl Gene ID':'candidategene_ensemblid','HGNC':'candidategene_hgncid'}
    tsv_candidate_gene_list = load_txt_list(f'{pwd}/{filename}/{txt_dir_dict[result_mode]}/{case_id}.txt')
    Filter_txt_dir_dict = {'Symbol': '', 'Ensembl Gene ID': '_ensemblid','HGNC': '_hgncid'}
    if len(REVEL_thresh) != 0:
        tsv_candidate_gene_revel_dict = {}
        for thresh in REVEL_thresh:
            tsv_candidate_gene_revel_dict[thresh] = load_txt_list(f'{pwd}/{filename}/revel_{thresh}{Filter_txt_dir_dict[result_mode]}/{case_id}.txt')

    if len(CADD_thresh) != 0:
        tsv_candidate_gene_cadd_dict = {}
        for thresh in CADD_thresh:
            tsv_candidate_gene_cadd_dict[thresh] = load_txt_list(f'{pwd}/{filename}/cadd_{thresh}{Filter_txt_dir_dict[result_mode]}/{case_id}.txt')

    if tool == 'phenomizer':
        part_table.loc[i, f'{tool}_rank'], part_rr_table.loc[i, f'{tool}_rank'], part_rp_table.loc[
            i, f'{tool}_rank'] = get_phenomizer_rank(pathogenic_gene_dict[result_mode], result_list)
    else:
        part_table.loc[i, f'{tool}_rank'], part_rr_table.loc[i, f'{tool}_rank'], part_rp_table.loc[
            i, f'{tool}_rank'] = get_rank(pathogenic_gene_dict[result_mode], result_list)
    if intersect:
        if tool =='phenomizer':
            tool_gene_intersect_list =[]
            for y in result_list:
                for x in tsv_candidate_gene_list:
                    if f'{x} (' in str(y):
                        tool_gene_intersect_list = tool_gene_intersect_list+[str(y)]
                        break
            part_table.loc[i, f'{tool}_rank_intersect'], part_rr_table.loc[i, f'{tool}_rank_intersect'], \
            part_rp_table.loc[i, f'{tool}_rank_intersect'] = get_phenomizer_rank(pathogenic_gene_dict[result_mode],
                                                                      tool_gene_intersect_list)
        else:
            tool_gene_intersect_list = [x for x in result_list if
                                        x in tsv_candidate_gene_list]
            part_table.loc[i, f'{tool}_rank_intersect'], part_rr_table.loc[i, f'{tool}_rank_intersect'], \
            part_rp_table.loc[i, f'{tool}_rank_intersect'] = get_rank(pathogenic_gene_dict[result_mode], tool_gene_intersect_list)
    if len(REVEL_thresh) != 0:
        for thresh in REVEL_thresh:
            if tool =='phenomizer':
                tool_gene_revel_list = []
                for y in result_list:
                    for x in tsv_candidate_gene_revel_dict[thresh]:
                        if f'{x} (' in str(y):
                            tool_gene_revel_list = tool_gene_revel_list + [str(y)]
                            break
                part_table.loc[i, f'{tool}_rank_revel_{thresh}'], part_rr_table.loc[
                    i, f'{tool}_rank_revel_{thresh}'], part_rp_table.loc[
                    i, f'{tool}_rank_revel_{thresh}'] = get_phenomizer_rank(pathogenic_gene_dict[result_mode], tool_gene_revel_list)
            else:
                tool_gene_revel_list = [x for x in result_list if
                                        x in tsv_candidate_gene_revel_dict[thresh]]
                part_table.loc[i, f'{tool}_rank_revel_{thresh}'], part_rr_table.loc[
                    i, f'{tool}_rank_revel_{thresh}'], part_rp_table.loc[
                    i, f'{tool}_rank_revel_{thresh}'] = get_rank(pathogenic_gene_dict[result_mode], tool_gene_revel_list)
    if len(CADD_thresh) != 0:
        for thresh in CADD_thresh:
            if tool == 'phenomizer':
                tool_gene_cadd_list = []
                for y in result_list:
                    for x in tsv_candidate_gene_cadd_dict[thresh]:
                        if f'{x} (' in str(y):
                            tool_gene_cadd_list = tool_gene_cadd_list + [str(y)]
                            break
                part_table.loc[i, f'{tool}_rank_cadd_{thresh}'], part_rr_table.loc[
                    i, f'{tool}_rank_cadd_{thresh}'], part_rp_table.loc[
                    i, f'{tool}_rank_cadd_{thresh}'] = get_phenomizer_rank(pathogenic_gene_dict[result_mode],tool_gene_cadd_list)
            else:
                tool_gene_cadd_list = [x for x in result_list if
                                       x in tsv_candidate_gene_cadd_dict[thresh]]
                part_table.loc[i, f'{tool}_rank_cadd_{thresh}'], part_rr_table.loc[
                    i, f'{tool}_rank_cadd_{thresh}'], part_rp_table.loc[
                    i, f'{tool}_rank_cadd_{thresh}'] = get_rank(pathogenic_gene_dict[result_mode],
                                                                     tool_gene_cadd_list)
    return part_table, part_rr_table, part_rp_table

def processdict(file_dir_dict, mode,hpo,case_id,inheritance):
    file_dir_dict['Exomiser'][
        mode] = f'{hpo}_Exomiseroutput{mode}/{case_id}_{inheritance}.genes.tsv'
    file_dir_dict['Exomiser13.1.0'][
        mode] = f'{hpo}_Exomiser1310output{mode}/{case_id}.genes.tsv'
    file_dir_dict['LIRICAL'][
        mode] = f'{hpo}_LIRICALoutput{mode}/{case_id}.tsv'


def statistic_to_table(df,pwd,filename,hpo,tools,intersect,REVEL_thresh,CADD_thresh,collect_variants_tofile,collect_organ):
    columns_fill =output_columns(hpo,tools,intersect,REVEL_thresh,CADD_thresh,collect_variants_tofile,collect_organ)
    part_table = pd.DataFrame(columns=columns_fill)
    part_rr_table = pd.DataFrame(columns=columns_fill)
    part_rp_table = pd.DataFrame(columns=columns_fill)
    part_table_for_R = pd.DataFrame(columns=columns_fill)
    part_rr_table_for_R = pd.DataFrame(columns=columns_fill)
    part_rp_table_for_R = pd.DataFrame(columns=columns_fill)
    ##for R是后面去掉NA的。
    log = pd.DataFrame(columns=['case_id','hpo_id_input','Symbol']+tools)
    svcf_vcf_adress_df = pd.read_csv(f'{pwd}/{filename}/log/standard_vcf_log_file.csv').set_index('case_id')

    for i in tqdm(range(len(df))):
        case_id = df.loc[i, 'Blood ID']
        log.loc[i,'case_id'] = case_id
        try:
            if df.loc[i, hpo][0] == '[':
                hpo_id = df.loc[i, hpo][2:-2]
                hpo_id_input = [k for k in hpo_id.split("', '")]
            else:
                hpo_id = df.loc[i, hpo]
                hpo_id_input = [k for k in hpo_id.split(";")]
        except Exception as e:
            log.loc[i,'hpo_id_input'] = e
            continue


        symbol = df.loc[i, "Symbol"]
        inheritance = df.loc[i, 'Inheritance_ADAR']
        pathogenic_gene_dict = (df[['Symbol', 'entrezGeneId', 'entrezId', 'Ensembl Gene ID', 'HGNC']].loc[i]).to_dict()
        pathogenic_gene_dict['HGNC'] = str(int(pathogenic_gene_dict['HGNC']))
        pathogenic_gene_dict['Symbol_inheritance_ADAR'] = symbol + '_' + inheritance



        filedir,tsv_source_mode = dir_of_Filtered_TSV_MUltiple_source(pwd,filename,case_id)
        if tsv_source_mode=='no Filtered TSV anywhere':
            print(f'{case_id} no tsv file, make sure you have manual diagnosis')
            log.loc[i,'Symbol'] = f'no tsv file'
            continue

        part_table.loc[i, 'CaseID'],part_rr_table.loc[i, 'CaseID'],part_rp_table.loc[i, 'CaseID'] = case_id,case_id,case_id
        part_table.loc[i, 'hpo'],part_rr_table.loc[i, 'hpo'],part_rp_table.loc[i, 'hpo'] = hpo,hpo,hpo
        part_table.loc[i, 'hpo_id_input'],part_rr_table.loc[i, 'hpo_id_input'],part_rp_table.loc[i, 'hpo_id_input'] = str(hpo_id_input),str(hpo_id_input),str(hpo_id_input)
        part_table.loc[i, 'Symbol'],part_rr_table.loc[i, 'Symbol'],part_rp_table.loc[i, 'Symbol'] = symbol,symbol,symbol

        if collect_variants_tofile:
            part_rr_table.loc[i, 'Frequency'] = df.loc[i,'Frequency']
            part_rr_table.loc[i, 'Variant_class'] = df.loc[i, 'Variant_class']
            part_rr_table.loc[i, 'Inheritance_ADAR'] = df.loc[i, 'Inheritance_ADAR']
        if collect_organ:
            part_rr_table.loc[i, f'{hpo}_organ_system'] =df.loc[i,f'{hpo}_organ_system']
            #part_rr_table.loc[i, f'{hpo}_organ_system_number'] = df.loc[i, f'{hpo}_organ_system_number']

        mode_dict = {'result': {'PhenoGenius': 'Symbol', 'phrank': 'Ensembl Gene ID', 'phenolyzer': 'Symbol',
                                'GADO': 'Ensembl Gene ID','phenoapt':'Symbol','HANRD':'HGNC','PhenoRank':'Ensembl Gene ID',
                                'phenomizer':'Symbol'},
                     'dir': {'PhenoGenius': f'{hpo}_PhenoGenius_output/{case_id}.tsv',
                             'phrank': f'phrankoutput/{case_id}_{hpo}_nogenelist_phrank_rank.tsv',
                             'phenolyzer': f'phenolyzeroutput/{case_id}_{hpo}.final_gene_list',
                             'GADO': f'GADOoutput/{case_id}.txt','phenoapt':f'phenoaptoutput/{case_id}_{hpo}_phenoapt_rank.tsv',
                             'HANRD':f'HANRDoutput/output/{case_id}_{hpo}.tsv',
                             'PhenoRank':f'PhenoRankoutput/{case_id}_{hpo}.tsv',
                             'phenomizer':f'phenomizeroutput/{hpo}/{case_id}.tsv'},
                     'column': {'PhenoGenius': 'gene_symbol', 'phrank': '1', 'phenolyzer': 'Gene',
                                'GADO': 'Ensg','phenoapt':'gene_symbol',
                                'HANRD':'HGNCid',
                                'PhenoRank':'GENE',
                                'phenomizer':'Disease-Name'}
                     }

        if hpo == 'Weight':
            def get_phenoapt_result_dir(k):
                if k == 'all_hpo_plus_weight':
                    pheno_result_dir = f'phenoaptoutput/{case_id}_{hpo}_phenoapt_rank.tsv'
                    return pheno_result_dir
                if k == 'only_weight_hpo_no_weight':
                    pheno_result_dir = f'phenoaptoutput/{case_id}_only_{hpo}_hpo_no_weight_phenoapt_rank.tsv'
                    return pheno_result_dir
                if k == 'only_weight_hpo_weight':
                    pheno_result_dir = f'phenoaptoutput/{case_id}_only_{hpo}_hpo_weight_phenoapt_rank.tsv'
                    return pheno_result_dir
            for k in ['all_hpo_plus_weight', 'only_weight_hpo_weight', 'only_weight_hpo_no_weight']:
                tools = tools + [f'{k}_phenoapt']
                mode_dict['result'][f'{k}_phenoapt'] = 'Symbol'
                mode_dict['dir'][f'{k}_phenoapt'] = get_phenoapt_result_dir(k)
                mode_dict['column'][f'{k}_phenoapt'] = 'gene_symbol'


        ## pheno-only tools table
        for tool in tools:
            if tool in mode_dict['result'].keys():
                try:
                    tool_rank_tsv = pd.read_csv(f"{pwd}/{filename}/{mode_dict['dir'][tool]}", sep='\t')
                    tool_gene_list = list(tool_rank_tsv[mode_dict['column'][tool]])
                    if tool=='HANRD':
                        tool_gene_list = [id.split(':')[1] for id in tool_gene_list]
                    part_table, part_rr_table, part_rp_table = process_rank_result_table(tool, mode_dict['result'][tool], tool_gene_list, case_id, i, pathogenic_gene_dict, pwd, filename,
                                              intersect, REVEL_thresh, CADD_thresh, part_table, part_rr_table, part_rp_table)
                    print(f'{tool} ok')
                except Exception as e:
                    log.loc[i,tool] = e
        ## 输入结果直接包含筛选，多种output，不需要和txt intersect
        file_dir_dict = {
            'Exomiser':{'original':f'{hpo}_Exomiseroutput/{case_id}_{inheritance}.genes.tsv',
                        'column':'ENTREZ_GENE_ID',
                        'result':'entrezId'},
            'phen2gene':{'original':f'phen2geneoutput_nolist/{case_id}-phen2gene-nolist-output.txt',
                         'column':'ID',
                         'result':'entrezId'},
            'Exomiser13.1.0': {'original': f'{hpo}_Exomiser1310output/{case_id}.genes.tsv',
                             'column':'ID',
                             'result':'Symbol_inheritance_ADAR'},
            'LIRICAL': {'original': f'{hpo}_LIRICALoutput/{case_id}.tsv',
                        'column':'entrezGeneId',
                        'result':'entrezGeneId'}
        }
        if intersect:
            file_dir_dict['phen2gene']['intersect'] = f'phen2geneoutput/{case_id}-phen2gene-candidategene-output.txt'
            processdict(file_dir_dict, 'intersect',hpo,case_id,inheritance)
        if len(REVEL_thresh)!=0:
            for thresh in REVEL_thresh:
                file_dir_dict['phen2gene'][f'revel_{thresh}'] = f'phen2geneoutput/{case_id}-phen2gene-revel_{thresh}-output.txt'
                processdict(file_dir_dict, f'revel_{thresh}',hpo,case_id,inheritance)
        if  len(CADD_thresh)!=0:
            for thresh in CADD_thresh:
                file_dir_dict['phen2gene'][
                    f'cadd_{thresh}'] = f'phen2geneoutput/{case_id}-phen2gene-cadd_{thresh}-output.txt'
                processdict(file_dir_dict, f'cadd_{thresh}',hpo,case_id,inheritance)

            ##准备好了integ需要的dict,处理table
        for tool in tools:
            if tool in file_dir_dict.keys():
                if tool in log.columns:
                    log.drop(tool, inplace=True, axis=1)
                for mode in file_dir_dict[list(file_dir_dict.keys())[0]].keys():
                    if mode not in ['column','result']:
                        try:
                            tool_rank_tsv = pd.read_csv(f"{pwd}/{filename}/{file_dir_dict[tool][mode]}", sep='\t')
                            tool_gene_list = list(tool_rank_tsv[file_dir_dict[tool]['column']])
                            if mode =='original':
                                part_table.loc[i, f'{tool}_rank'], part_rr_table.loc[i, f'{tool}_rank'], part_rp_table.loc[
                                    i, f'{tool}_rank'] = get_rank(pathogenic_gene_dict[file_dir_dict[tool]['result']], tool_gene_list)
                            else:
                                part_table.loc[i, f'{tool}_rank_{mode}'], part_rr_table.loc[i, f'{tool}_rank_{mode}'], part_rp_table.loc[
                                    i, f'{tool}_rank_{mode}'] = get_rank(pathogenic_gene_dict[file_dir_dict[tool]['result']], tool_gene_list)
                            print(f'{tool} ok')
                        except Exception as e:
                            if f'{tool}_{mode}' not in log.columns:
                                loginfo = ['' for k in range(len(log)-1)] + [str(e)]
                                log[f'{tool}_{mode}'] = loginfo
                            else:
                                log.loc[i,f'{tool}_{mode}'] = str(e)

        if 'Phen_gen' in tools:
            try:
                part_table.loc[i, 'Phen_gen_rank'], part_rr_table.loc[i, 'Phen_gen_rank'] = phen_gen(case_id, pwd, hpo)
                print('phen_gen ok')
            except Exception as e:
                log.loc[i, 'Phen_gen'] = e
            ## 转化成作图数据

        if svcf_vcf_adress_df.loc[case_id,'vcf_file_address'] != 0: ##'NA'是没有结果，nan是找结果的过程中有bug
            part_table_for_R.loc[i] = part_table.loc[i]
            part_rr_table_for_R.loc[i] = part_rr_table.loc[i]
            part_rp_table_for_R.loc[i] = part_rp_table.loc[i]
            ## 替换排名
            for column in part_table_for_R.columns:
                if 'rank' in column:
                    rank_value =part_table_for_R.loc[i, column]
                    if str(rank_value) != 'NA' and str(rank_value) != 'nan':
                        rank_value = int(rank_value)
                        if rank_value == 0:
                            part_table_for_R.loc[i, column] = 'TOP1'
                        if rank_value >= 1 and rank_value <= 4:
                            part_table_for_R.loc[i, column] = 'TOP5'
                        if rank_value > 4 and rank_value <= 9:
                            part_table_for_R.loc[i, column] = 'TOP10'
                        if rank_value > 9 and rank_value <= 19:
                            part_table_for_R.loc[i, column] = 'TOP20'
                        if rank_value > 19 and rank_value <= 49:
                            part_table_for_R.loc[i, column] = 'TOP50'
                        if rank_value > 49 and rank_value <= 99:
                            part_table_for_R.loc[i, column] = 'TOP100'
                        if rank_value > 99:
                            part_table_for_R.loc[i, column] = '>100'
                    print(f'{column} rank count ok')
            ## 替换NA
            if 'phenoapt' in tools and hpo == 'Weight':
                tools = tools+[f'{k}_phenoapt' for k in ['all_hpo_plus_weight', 'only_weight_hpo_weight', 'only_weight_hpo_no_weight']]
            for tool in tools:
                if part_table_for_R.loc[i, f'{tool}_rank'] == 'NA':
                    for column in part_table_for_R.columns:
                        if tool in str(column):
                            part_table_for_R.loc[i, column] = 'not ranked'
                else:
                    if intersect:
                        if part_table_for_R.loc[i, f'{tool}_rank_intersect'] == 'NA':
                            for column in part_table_for_R.columns:
                                if tool in str(column):
                                    if column != f'{tool}_rank':
                                        part_table_for_R.loc[i, column] = 'not_in_filtered_tsv'
                        else:
                            if len(REVEL_thresh) != 0:
                                thresh_na = []
                                for thresh in REVEL_thresh:
                                    # part_table_for_R.loc[i, f'{tool}_rank_REVEL_{thresh}']
                                    if part_table_for_R.loc[i, f'{tool}_rank_revel_{thresh}'] == 'NA':
                                        thresh_na = thresh_na + [thresh]
                                if len(thresh_na) != 0:
                                    filter_thresh = min(thresh_na)
                                    for thresh in thresh_na:
                                        part_table_for_R.loc[
                                            i, f'{tool}_rank_REVEL_{thresh}'] = f'filtered_by_revel_{filter_thresh}'
                            if len(CADD_thresh) != 0:
                                thresh_na = []
                                for thresh in CADD_thresh:
                                    if part_table_for_R.loc[i, f'{tool}_rank_cadd_{thresh}'] == 'NA':
                                        thresh_na = thresh_na + [thresh]
                                if len(thresh_na) != 0:
                                    filter_thresh = min(thresh_na)
                                    for thresh in thresh_na:
                                        part_table_for_R.loc[
                                            i, f'{tool}_rank_cadd_{thresh}'] = f'filtered_by_cadd_{filter_thresh}'
                print(f'{tool} NA replace ok')
                if tool == 'Phen_gen':
                    if part_table_for_R.loc[i, f'{tool}_rank'] == 'NA':
                        part_table_for_R.loc[i, f'{tool}_rank'] = f'not ranked'
                        print(f'Phen_gen NA replace ok')
                ##if TOPpercent:
            print(f'{case_id} statistic ok')
        log.to_csv(f'{pwd}/{filename}/log/statistic.csv')
        part_table.to_csv(f'{hpo}_{len(tools)}_rank_len_{len(df)}.tsv', sep='\t')
        part_rr_table.to_csv(f'{hpo}_{len(tools)}_rr_len_{len(df)}.tsv', sep='\t')
        part_rp_table.to_csv(f'{hpo}_{len(tools)}_rp_len_{len(df)}.tsv', sep='\t')
        part_table_for_R.to_csv(f'{hpo}_{len(tools)}_rank_with_vcf_len_{len(df)}.tsv', sep='\t')
        part_rr_table_for_R.to_csv(f'{hpo}_{len(tools)}_rr_with_vcf_len_{len(df)}.tsv', sep='\t')
        part_rp_table_for_R.to_csv(f'{hpo}_{len(tools)}_rp_with_vcf_len_{len(df)}.tsv', sep='\t')


def get_rank(pathogene,rank_list):
    if pathogene in rank_list:
        rank = rank_list.index(pathogene)
        rr = 1 / (rank + 1)
        rp = (rank+1)/(len(rank_list))
    else:
        rank = 'NA'
        rr = 0
        rp= 'NA'
    return rank, rr,rp


def get_case_id_integ_file_map(pwd, filename, integ, case_ids, hpo):
    case_id_tsv_file_map = {}
    if hpo == 'hpo_id':
        tsv_dir = f'{pwd}/{filename}/{integ}output'
    else:
        tsv_dir = f'{pwd}/{filename}/{hpo}_{integ}output'
    files = os.listdir(tsv_dir)
    for case_id in case_ids:
        for file_name in files:
            if file_name.endswith(f'.tsv'):
                if case_id in file_name:
                    case_id_tsv_file_map[case_id] = os.path.join(tsv_dir, file_name)
    return case_id_tsv_file_map


def read_xlsx(path, sheet):
    wb = load_workbook(path)
    ws = wb[sheet]
    sh = pd.DataFrame(ws.values)
    df = sh.rename(columns=sh.loc[0])[1:].reset_index(drop=True)
    return df


def containtype(biao, type):
    data_mut = biao[biao["Mutation_type"].str.contains(type)]
    return data_mut
    # 包含字符串type=“”的rows，Na的单元格自动去掉

def get_HANRD(df,pwd,filename,hpo,refresh=False,refresh_seperate_file=False):
    #print('start')
    HANRD_input=pd.DataFrame(columns=['ID','HPO','Orphanet ID (Diagnosed Disorder)','HGNC (Casual Gene)'])
    for i in range(len(df)):
        if df.loc[i, hpo][0] == '[':
            hpo_id = df.loc[i, hpo][2:-2]
            hpo_id_input = [str(k) for k in hpo_id.split("', '")]
        else:
            hpo_id = df.loc[i, hpo]
            hpo_id_input = [k for k in hpo_id.split(";")]
        #print(hpo_id_input)
        HPO=''
        for k in hpo_id_input:
            HPO = HPO+str(k)+', '
        print(HPO[:-2])
        HANRD_input.loc[i,'HPO'] = HPO[:-2]
        HANRD_input.loc[i,'ID'] = df.loc[i, 'Blood ID']
        HANRD_input.loc[i,'HGNC (Casual Gene)'] = '-'
        HANRD_input.loc[i,'Orphanet ID (Diagnosed Disorder)'] = '-'
        # HANRD_input.set_index('ID')
    if 'HANRDoutput' not in os.listdir(f'{pwd}/{filename}'):
        os.system('mkdir HANRDoutput')
    if f'{hpo}_gVCF-diagnosed-cohort-2022-HANRD.tsv' in os.listdir(f'{pwd}/{filename}/HANRDoutput/'):
        ##已经有输入file了
        old_input_file = pd.read_csv(f'{pwd}/{filename}/HANRDoutput/{hpo}_gVCF-diagnosed-cohort-2022-HANRD.tsv', sep='\t')
        if old_input_file.equals(HANRD_input) and refresh==False:
            print('using HANRD result in cache')
        else:
            os.system(f'cp {pwd}/{filename}/HANRDoutput/{hpo}_gVCF-diagnosed-cohort-2022-HANRD.tsv {pwd}/{filename}/HANRDoutput/{hpo}_gVCF-diagnosed-cohort-2022-HANRD-last-input.tsv')
            HANRD_input.to_csv(f'{pwd}/{filename}/HANRDoutput/{hpo}_gVCF-diagnosed-cohort-2022-HANRD.tsv', sep='\t',index=False)
            print('new HANRD file ready')
            os.system(f'cd ~/Software/HANRD/gcas && java -Xmx8096M -jar gcas.jar {pwd}/{filename}/HANRDoutput/{hpo}_gVCF-diagnosed-cohort-2022-HANRD.tsv')
            print('HANRD done')
    else:
        HANRD_input.to_csv(f'{pwd}/{filename}/HANRDoutput/{hpo}_gVCF-diagnosed-cohort-2022-HANRD.tsv', sep='\t',index=False)
        print('new HANRD file ready')
        os.system(f'cd ~/Software/HANRD/gcas && java -Xmx8096M -jar gcas.jar {pwd}/{filename}/HANRDoutput/{hpo}_gVCF-diagnosed-cohort-2022-HANRD.tsv && mv ./data/output/{hpo}_gVCF-diagnosed-cohort-2022-HANRD_genes.out {pwd}/{filename}/HANRDoutput/{hpo}_gVCF-diagnosed-cohort-2022-HANRD_genes.csv')
        print('HANRD done')
    rankdf = pd.read_csv(f'{pwd}/{filename}/HANRDoutput/{hpo}_gVCF-diagnosed-cohort-2022-HANRD_genes.out', header=None,
                         sep='\t')
    rankdf.columns=['id', 'rank', 'HGNCid', 'score', '-', 'gene_note', 'chr', 'note']
    ## 拆分成单个文件
    for i in range(len(df)):
        case_id = df.loc[i,'Blood ID']
        HANRD_seperate_df = pd.DataFrame(columns=['id', 'rank', 'HGNCid', 'score', '-', 'gene_note', 'chr', 'note'])
        if f'{case_id}_{hpo}.tsv' not in os.listdir(f'{pwd}/{filename}/HANRDoutput/output/') or refresh_seperate_file:
            for k in range(len(rankdf)):
               if rankdf.iloc[k, 0] == case_id:
                    ##iloc可以用数字索引
                    HANRD_seperate_df.loc[k] = list(rankdf.loc[k])
            HANRD_seperate_df.to_csv(f'{pwd}/{filename}/HANRDoutput/output/{case_id}_{hpo}.tsv', sep='\t',index=False)
            print(f'{case_id} result saved')
        else:
            print(f'{case_id} result exist')

def pubcasefinder(df,pwd,filename,hpo,refresh=False):
    ##日本时间晚上11点到凌晨五点之间
    if 'pubcasefinderoutput' not in os.listdir(f'{pwd}/{filename}'):
        os.system(f'mkdir {pwd}/{filename}/pubcasefinderoutput')
    for i in range(len(df)):
        if df.loc[i, hpo][0] == '[':
            hpo_id = df.loc[i, hpo][2:-2]
            hpo_id_input = [str(k) for k in hpo_id.split("', '")]
        else:
            hpo_id = df.loc[i, hpo]
            hpo_id_input = [k for k in hpo_id.split(";")]
        #print(hpo_id_input)
        HPO=''
        for k in hpo_id_input:
            HPO = HPO+str(k)+','
        url = f'https://pubcasefinder.dbcls.jp/api/get_ranked_list?target=gene&format=json&hpo_id={HPO[:-1]}'
        print(url)
        case_id = df.loc[i,'Blood ID']
        if f'{case_id}_{hpo}.json' not in os.listdir(f'{pwd}/{filename}/pubcasefinderoutput') or refresh:
            os.system(f'cd {pwd}/{filename}/pubcasefinderoutput && wget -O {case_id}_{hpo}.json {url}')
            time.sleep(10)

def get_PhenoRank(df,pwd,filename,hpo,refresh=False):
    for i in range(len(df)):
        if len(df.loc[i, f'{hpo}_PhenoRank'])!=0:
            hpo_reset = f'{hpo}_PhenoRank'
            if df.loc[i, hpo_reset][0] == '[':
                hpo_id = df.loc[i, hpo_reset][2:-2]
                hpo_id_input = [str(k) for k in hpo_id.split("', '")]
            else:
                hpo_id = df.loc[i, hpo_reset]
                hpo_id_input = [k for k in hpo_id.split(";")]
        else:

                if df.loc[i, hpo][0] == '[':
                    hpo_id = df.loc[i, hpo][2:-2]
                    hpo_id_input = [str(k) for k in hpo_id.split("', '")]
                else:
                    hpo_id = df.loc[i, hpo]
                    hpo_id_input = [k for k in hpo_id.split(";")]
        #需要手动删掉搜不到的HPO写在hpo_id_PhenoRank里
        HPO=''
        for k in hpo_id_input:
            HPO = HPO+str(k)+';'
        case_id = df.loc[i,'Blood ID']
        if 'PhenoRankoutput' not in os.listdir(f'{pwd}/{filename}'):
            os.system(f'mkdir {pwd}/{filename}/PhenoRankoutput')
        if f'{case_id}_{hpo}.tsv' not in os.listdir(f'{pwd}/{filename}/PhenoRankoutput') or refresh:
            cmd = f"cd ~/Software/PhenoRank && zsh && source ~/opt/anaconda3/etc/profile.d/conda.sh && conda activate phenorank && python run_PhenoRank.py -o {pwd}/{filename}/PhenoRankoutput/{case_id}_{hpo}.tsv -p '{HPO[:-1]}' && exit"
            appscript.app('Terminal').do_script(cmd)
            time.sleep(15)
            print('new phenorank ok')
        else:
            print('phenorank exist')

def GADO(df,pwd,filename,hpo,refresh=False):
    hpotxt = {}
    for i in range(len(df)):
        case_id = df['Blood ID'][i]
        if df.loc[i, hpo][0] == '[':
            hpo_id = df.loc[i, hpo][2:-2]
            hpo_id_input = [k for k in hpo_id.split("', '")]
        else:
            hpo_id = df.loc[i, hpo]
            hpo_id_input = [k for k in hpo_id.split(";")]
        ##所有case记录在hpotxt{}里
        hpotxt[case_id] = hpo_id_input

    file = open(f"./{filename}/hpotxt/totalhpo_{hpo}.txt", "w")
    for key, value in hpotxt.items():
        hpoline = ''
        for k in value:
            hpoline = hpoline + '\t' + k
        file.write('%s%s\n' % (key, hpoline))
        print(file)
    file.close()

    # present_in_put_hpo = pd.read_csv(f'./{filename}/hpotxt/totalhpo_{hpo}.txt',sep='\t')
    # if present_in_put_hpo.equals(previous_input_hpo) and (refresh == False):
    #     print('GADO already finish')
    # else:
    cmd1 = f'cd ~/Software/GadoCommandline-1.0.1 && java -jar GADO.jar \
     --mode PROCESS \
     --output {pwd}/{filename}/GADOoutput/hpoProcessed.txt \
     --caseHpo {pwd}/{filename}/hpotxt/totalhpo_{hpo}.txt \
     --hpoOntology data/hp.obo \
     --hpoPredictionsInfo data/predictions_auc_bonf.txt'
    cmd2 = f'cd ~/Software/GadoCommandline-1.0.1 && java -jar GADO.jar \
     --mode PRIORITIZE \
     --output {pwd}/{filename}/{hpo}_GADOoutput/ \
     --caseHpoProcessed {pwd}/{filename}/GADOoutput/hpoProcessed.txt \
     --genes data/hpo_prediction_genes.txt \
     --hpoPredictions data/genenetwork_bonf_spiked/genenetwork_bonf_spiked'
    os.system(cmd1)
    os.system(cmd2)
    ##按照case_id输出；可以总是刷新，hpo_id结果也没有
    ##print('GADO finish with new input')

def phen2gene_rank(df, pwd, filename, hpo='hpo_id',profile='candidategene',refresh=False):
    # #phen2gene intersect
    # df = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx', 'Sheet1')
    # df = df[:1].reset_index(drop=True)
    for i in range(len(df)):
        case_id = df['Blood ID'][i]
        if df.loc[i, hpo][0] == '[':
            hpo_id = df.loc[i, hpo][2:-2]
            hpo_id_input = [k for k in hpo_id.split("', '")]
        else:
            hpo_id = df.loc[i, hpo]
            hpo_id_input = [k for k in hpo_id.split(";")]
        ## 没有有基因列表
        temp = ''
        for the_hpo_id in hpo_id_input:
            temp = temp + ' ' + the_hpo_id
        HPO_for_Phen2Gene = temp
        if profile == 'no_gene_list':
            if refresh or (f'{case_id}-phen2gene-nolist-output.txt' not in os.listdir(
                    f'{pwd}/{filename}/phen2geneoutput_nolist/')):
                command_phen2gene = f'cd ~/Phen2Gene && python3 phen2gene.py -m{HPO_for_Phen2Gene} -out {pwd}/{filename}/phen2geneoutput_nolist && mv {pwd}/{filename}/phen2geneoutput_nolist/output_file.associated_gene_list {pwd}/{filename}/phen2geneoutput_nolist/{case_id}-phen2gene-nolist-output.txt'
                os.system(command_phen2gene)
                print('run phen2gene')
        else:
            if refresh or (f'{case_id}-phen2gene-{profile}-output.txt' not in os.listdir(f'{pwd}/{filename}/phen2geneoutput/')):
                command_phen2gene_intersect = f'cd ~/Phen2Gene && python3 phen2gene.py -m{HPO_for_Phen2Gene} -out {pwd}/{filename}/phen2geneoutput -l {pwd}/{filename}/{profile}/{case_id}.txt && mv {pwd}/{filename}/phen2geneoutput/output_file.associated_gene_list {pwd}/{filename}/phen2geneoutput/{case_id}-phen2gene-{profile}-output.txt'
                os.system(command_phen2gene_intersect)
                print('run phen2gene')
        ## 有基因列表

def Phenogenius_rank(df,pwd,filename,hpo,refresh):
    print('----------Phenogenius running-------------\n')
    for i in range(len(df)):
        case_id = df['Blood ID'][i]
        outputdir = f'{hpo}_PhenoGenius_output'
        if outputdir not in os.listdir(f'{pwd}/{filename}'):
            os.system(f'mkdir {pwd}/{filename}/{outputdir}')
        if f'{case_id}.tsv' not in os.listdir(f'{pwd}/{filename}/{outputdir}') or refresh:
            if df.loc[i, hpo][0] == '[':
                hpo_id = df.loc[i, hpo][2:-2]
                hpo_id_input = [k for k in hpo_id.split("', '")]
            else:
                hpo_id = df.loc[i, hpo]
                hpo_id_input = [k for k in hpo_id.split(";")]
            hpo_id_list_str =''
            for hpo_id_str in hpo_id_input:
                hpo_id_list_str = hpo_id_list_str + hpo_id_str+','
            hpo_id_list_str = hpo_id_list_str[:-1]
            cmd = f"/bin/zsh -c 'source ~/opt/anaconda3/etc/profile.d/conda.sh && conda activate PhenoGenius && cd  /Users/liyaqi/Software/PhenoGenius && python phenogenius_cli.py --hpo_list {hpo_id_list_str} --result_file {pwd}/{filename}/{outputdir}/{case_id}.tsv'"
            print(cmd)
            os.system(cmd)
        print(f"No. {i} case id {case_id}: PhenoGenius done")
    print('----------Phenogenius has finished-------------\n')


def hpotxt_phenolyzer(df,pwd,filename,hpo,refresh_hpo_txt):
    #df = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx', 'Sheet1')
    # df = df[:1].reset_index(drop=True)
    if 'hpotxt' not in os.listdir(f'{pwd}/{filename}'):
        os.system(f'mkdir {pwd}/{filename}/hpotxt')
    for i in range(len(df)):
        case_id = df['Blood ID'][i]
        if (f'{case_id}_{hpo}.txt' not in os.listdir(f'{pwd}/{filename}/hpotxt/')) or refresh_hpo_txt:
            if df.loc[i, hpo][0] == '[':
                hpo_id = df.loc[i, hpo][2:-2]
                hpo_id_input = [k for k in hpo_id.split("', '")]
            else:
                hpo_id = df.loc[i, hpo]
                hpo_id_input = [k for k in hpo_id.split(";")]
            ##每个case单独一个txt-phenolyzer
            hpo_id_input_df = pd.DataFrame(hpo_id_input)
            hpo_id_input_df.to_csv(
                f'{pwd}/{filename}/hpotxt/{case_id}_{hpo}.txt',
                sep='\n',
                index=False, header=False)
            print(f"{i}txt done")


def phenolyzer(df,pwd,filename,hpo,refresh=False):
    ##需要新修改既往hpotxt for phenolyzer时打开refresh，单纯新增case不需要refresh
    hpotxt_phenolyzer(df,pwd,filename,hpo,refresh)
    ##运行phenolyzer
    #df = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx', 'Sheet1')
    # df = df[:1].reset_index(drop=True)
    if 'phenolyzeroutput' not in os.listdir(f'{pwd}/{filename}'):
        os.system(f'mkdir {pwd}/{filename}/phenolyzeroutput')
    for i in tqdm(range(len(df))):
        case_id = df['Blood ID'][i]
        if hpo !='hpo_id':
            hpo_txt_dir = f'{pwd}/{filename}/hpotxt/{case_id}_{hpo}.txt'
            if (f'{case_id}_{hpo}' not in ''.join(os.listdir(f'{pwd}/{filename}/phenolyzeroutput/'))) or refresh:
                cmd_phenolyzer = f'cd ~/Software/phenolyzer && perl disease_annotation.pl {hpo_txt_dir} -file -prediction -phenotype -logistic -out {pwd}/{filename}/phenolyzeroutput/{case_id}_{hpo} -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 && exit'
                print(cmd_phenolyzer)
                appscript.app('Terminal').do_script(cmd_phenolyzer)
                time.sleep(30)
        else:
            hpo_txt_dir = f'{pwd}/{filename}/hpotxt/{case_id}_{hpo}.txt'
            ##如果有新case，还是用新产生的——hpo_id后缀的hpo.txt
            if (f'{case_id}.' not in ''.join(os.listdir(f'{pwd}/{filename}/phenolyzeroutput/'))) or refresh:
                cmd_phenolyzer = f'cd ~/Software/phenolyzer && perl disease_annotation.pl {hpo_txt_dir} -file -prediction -phenotype -logistic -out {pwd}/{filename}/phenolyzeroutput/{case_id} -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 && exit'
                print(cmd_phenolyzer)
                appscript.app('Terminal').do_script(cmd_phenolyzer)
                time.sleep(30)
    print('phenolyzer ready')
    return True


def phrank_rank(df, pwd, filename, mode, hpo='hpo_id', filter='candidategene_ensemblid', refresh=False):
    DAG = "./phrank/demo/data/hpodag.txt"
    DISEASE_TO_PHENO = "./phrank/demo/data/disease_to_pheno.build127.txt"
    DISEASE_TO_GENE = "./phrank/demo/data/gene_to_disease.build127.txt"
    GENE_TO_PHENO = "./phrank/demo/data/gene_to_pheno.amelie.txt"
    p_hpo = Phrank(DAG, diseaseannotationsfile=DISEASE_TO_PHENO, diseasegenefile=DISEASE_TO_GENE)
    if 'phrankoutput' not in os.listdir(f'{pwd}/{filename}'):
        os.system(f'mkdir {pwd}/{filename}/phrankoutput')
    if mode == 'intersect':
        for i in range(len(df)):
            case_id = df['Blood ID'][i]
            if (f'{case_id}_{hpo}_{filter}_phrank_rank.tsv' not in os.listdir(
                    f'{pwd}/{filename}/phrankoutput')) or refresh:
                if df.loc[i, hpo][0] == '[':
                    hpo_id = df.loc[i, hpo][2:-2]
                    hpo_id_input = [k for k in hpo_id.split("', '")]
                else:
                    hpo_id = df.loc[i, hpo]
                    hpo_id_input = [k for k in hpo_id.split(";")]
                patient_genes = pd.read_csv(f'./{filename}/{filter}/{case_id}.txt', header=None)
                phrank_gene_ranking = p_hpo.rank_genes_directly(patient_genes[0], hpo_id_input)
                phrank_gene_ranking = pd.DataFrame(phrank_gene_ranking)
                phrank_gene_ranking.to_csv(f'{pwd}/{filename}/phrankoutput/{case_id}_{hpo}_{filter}_phrank_rank.tsv',sep='\t')
                #print(f'phrank {hpo} {filter} finish')

    else:
        if mode == 'OMIM':
            for i in range(len(df)):
                case_id = df['Blood ID'][i]
                if (f'{case_id}_{hpo}_nogenelist_phrank_rank.tsv' not in os.listdir(
                        f'{pwd}/{filename}/phrankoutput')) or refresh:
                    if df.loc[i, hpo][0] == '[':
                        hpo_id = df.loc[i, hpo][2:-2]
                        hpo_id_input = [k for k in hpo_id.split("', '")]
                    else:
                        hpo_id = df.loc[i, hpo]
                        hpo_id_input = [k for k in hpo_id.split(";")]
                    patient_genes_df = pd.read_csv(f'./{filename}/omim-non-empty.txt', header=None)
                    patient_genes = patient_genes_df[0]
                    phrank_gene_ranking = p_hpo.rank_genes_directly(patient_genes, hpo_id_input)
                    phrank_gene_ranking = pd.DataFrame(phrank_gene_ranking)
                    phrank_gene_ranking.to_csv(
                        f'{pwd}/{filename}/phrankoutput/{case_id}_{hpo}_nogenelist_phrank_rank.tsv', sep='\t')
                    #print(f'phrank {hpo} no gene list finish')
        else:
            print('wrong mode input')

def phenoapt_rank(df, pwd, filename, hpo, refresh=False, all_hpo_plus_weight=False, only_weight_hpo_weight=False):
    ## PhenoApt排序
    if 'phenoaptoutput' not in os.listdir(f'{pwd}/{filename}'):
        os.system(f'mkdir {pwd}/{filename}/phenoaptoutput')
    for i in range(len(df)):
        case_id = df['Blood ID'][i]
        if (f'{case_id}_{hpo}_phenoapt_rank.tsv' not in os.listdir(f'{pwd}/{filename}/phenoaptoutput')) or refresh:
            print('generating new phenoapt searching')
            if hpo != 'Weight':
                if df.loc[i, hpo][0] == '[':
                    hpo_id = df.loc[i, hpo][2:-2]
                    hpo_id_input = [k for k in hpo_id.split("', '")]
                else:
                    hpo_id = df.loc[i, hpo]
                    hpo_id_input = [k for k in hpo_id.split(";")]
                weight_1 = [1 for k in range(len(hpo_id_input))]
                client = PhenoApt(token='H0pVk00CX07VzkZbdnvHI$24XiU$u9q')
                pheno_result = (client.rank_gene(phenotype=hpo_id_input, weight=weight_1, n=5000)).rank_frame
                pheno_result = pd.DataFrame(pheno_result)
                pheno_result.to_csv(f'{pwd}/{filename}/phenoaptoutput/{case_id}_{hpo}_phenoapt_rank.tsv', sep='\t')
            else:
                intrisic_weight_df= pd.read_csv(f'{pwd}/weight.csv')
                intrisic_weight_df = intrisic_weight_df.set_index('hpo_id')
                if df.loc[i, hpo][0] == '[':
                    hpo_id = df.loc[i, hpo][2:-2]
                    hpo_id_weight = [k for k in hpo_id.split("', '")]
                else:
                    hpo_id = df.loc[i, hpo]
                    hpo_id_weight = [k for k in hpo_id.split(";")]

                if all_hpo_plus_weight:
                    ##有所有HPO，重点HPO加上非1的intrisic weight，如果没有找到对应weight，就加上5
                    if df.loc[i, 'hpo_id'][0] == '[':
                        hpo_id = df.loc[i, 'hpo_id'][2:-2]
                        hpo_id_input = [k for k in hpo_id.split("', '")]
                    else:
                        hpo_id = df.loc[i, 'hpo_id']
                        hpo_id_input = [k for k in hpo_id.split(";")]
                    weight_dict = {}
                    for k in hpo_id_input:
                        if k not in hpo_id_weight:
                            weight_dict[k]=1
                        else:
                            if k in intrisic_weight_df.index:
                                weight_dict[k]=intrisic_weight_df.loc[k,'intrinsic_weight']
                            else:
                                weight_dict[k]=5
                    weight_1 = [weight_dict[k] for k in hpo_id_input]
                    client = PhenoApt(token='H0pVk00CX07VzkZbdnvHI$24XiU$u9q')
                    pheno_result = (client.rank_gene(phenotype=hpo_id_input, weight=weight_1, n=5000)).rank_frame
                    pheno_result = pd.DataFrame(pheno_result)
                    pheno_result.to_csv(f'{pwd}/{filename}/phenoaptoutput/{case_id}_{hpo}_phenoapt_rank.tsv',
                                        sep='\t')
                else:
                    hpo_id_input = hpo_id_weight
                    weight_dict={}
                    if only_weight_hpo_weight:
                        for k in hpo_id_input:
                            if k in intrisic_weight_df.index:
                                weight_dict[k] = intrisic_weight_df.loc[k, 'intrinsic_weight']
                            else:
                                weight_dict[k] = 5
                            ##只有重点HPO，有加权
                        weight_1 = [weight_dict[k] for k in hpo_id_input]
                        client = PhenoApt(token='H0pVk00CX07VzkZbdnvHI$24XiU$u9q')
                        pheno_result = (client.rank_gene(phenotype=hpo_id_input, weight=weight_1, n=5000)).rank_frame
                        pheno_result = pd.DataFrame(pheno_result)
                        pheno_result.to_csv(f'{pwd}/{filename}/phenoaptoutput/{case_id}_only_{hpo}_hpo_weight_phenoapt_rank.tsv',
                                            sep='\t')
                    else:
                        ##只有重点HPO，没有有加权
                        weight_1 = [1 for k in hpo_id_input]
                        client = PhenoApt(token='H0pVk00CX07VzkZbdnvHI$24XiU$u9q')
                        pheno_result = (client.rank_gene(phenotype=hpo_id_input, weight=weight_1, n=5000)).rank_frame
                        pheno_result = pd.DataFrame(pheno_result)
                        pheno_result.to_csv(f'{pwd}/{filename}/phenoaptoutput/{case_id}_only_{hpo}_hpo_no_weight_phenoapt_rank.tsv', sep='\t')
    print(f'phenoapt ready')


def getyml(case_id_file, hpo_id_input, dir):
    with open(f'../yml/test-analysis-exome.yml', 'r') as f:
        config = yaml.safe_load(f)
        config['analysis']['vcf'] = f'{dir}*vcf/{case_id_file}.vcf'
        config['analysis']['hpoIds'] = hpo_id_input  # add the command as a list for the correct yaml
        config['outputOptions']['outputFormats'] = ['TSV_GENE']
        config['outputOptions']['outputPrefix'] = f'{dir}Exomiseroutput/{case_id_file}'
        ##config['inheritanceModes'] = inheri
        print(config['analysis']['hpoIds'])

    with open(f'{dir}*yml/{case_id_file}.yml', "w") as f:  # open the file in append mode
        f.truncate(0)
        yaml.dump(config, f)  ##default_flow_style=Faulse

def phen_gen(case_id,pwd,hpo):
    phen_gen_hpo_id = read_xlsx(f'{pwd}/phen_gen_rank_rr.xlsx',hpo)
    if case_id in list(phen_gen_hpo_id['case_id']):
        rank_index = list(phen_gen_hpo_id['case_id']).index(case_id)
        rank = phen_gen_hpo_id['rank'][rank_index]
        rr = phen_gen_hpo_id['rr'][rank_index]
    else:
        rank = 'No_vcf_ranked_by_Phen-gen'
        rr = 0
    return rank, rr

def Rscript_brief_benchmark(tools=[],intersect=True,REVEL_thresh=[],CADD_thresh=[]):
    # ##在R script中添加strategy分组
    strategy_script = ""
    ##'data_set$strategy<-gsub("nursing 1st","nursing",data_set$strategy)'
    for tool in tools:
        if tool in ['phenoapt','phrank', 'phenolyzer', 'GADO', 'phen2gene','PhenoRank','HANRD','phenomizer']:
            strategy_script = strategy_script + "data_set$strategy<-gsub(" + f"'{tool}_rank$'" + "," + "'Pheno_only_tools',data_set$strategy)\n"
            if intersect:
                strategy_script = strategy_script + "data_set$strategy<-gsub(" + f"'{tool}_intersect_rank'" + "," + "'PUMCHpipeline',data_set$strategy)\n"
                if len(REVEL_thresh)!=0:
                    for thresh in REVEL_thresh:
                        strategy_script = strategy_script + "data_set$strategy<-gsub(" + f"'{tool}_rank_REVEL_{thresh}'" + "," + f"'REVEL_{thresh}',data_set$strategy)\n"
                if len(CADD_thresh)!=0:
                    for thresh in CADD_thresh:
                        strategy_script = strategy_script + "data_set$strategy<-gsub(" + f"'{tool}_rank_CADD_{thresh}'" + "," + f"'CADD_{thresh}',data_set$strategy)\n"
        else:
            strategy_script = strategy_script + "data_set$strategy<-gsub(" + f"'{tool}_rank$'" + "," + "'Integrated_tools',data_set$strategy)\n"

    tool_script = ""
    ##针对hpo_id做benchmark的代码
    for tool in tools:
        if tool in ['phenoapt', 'phrank', 'phenolyzer', 'GADO', 'phen2gene','PhenoRank','HANRD','phenomizer']:
            tool_script = tool_script + "data_set$variable<-gsub(" + f"'{tool}_rank$'" + "," + f"'{tool}',data_set$variable)\n"
            if intersect:
                tool_script = tool_script + "data_set$variable<-gsub(" + f"'{tool}_intersect_rank'" + "," + f"'{tool}',data_set$variable)\n"
                if len(REVEL_thresh) != 0:
                    for thresh in REVEL_thresh:
                        tool_script = tool_script + "data_set$variable<-gsub(" + f"'{tool}_rank_REVEL_{thresh}'" + "," + f"'{tool}',data_set$variable)\n"
                if len(CADD_thresh) != 0:
                    for thresh in CADD_thresh:
                        tool_script = tool_script + "data_set$variable<-gsub(" + f"'{tool}_rank_CADD_{thresh}'" + "," + f"'{tool}',data_set$variable)\n"
        else:
            tool_script = tool_script + "data_set$variable<-gsub(" + f"'{tool}_rank$'" + "," + f"'{tool}',data_set$variable)\n"

    level = "'Pheno_only_tools'"
    if intersect:
        level = level+", 'PUMCHpipeline'"
        if len(REVEL_thresh)!=0:
            for thresh in REVEL_thresh:
                level = level + f", 'REVEL_{thresh}'"
        if len(CADD_thresh)!=0:
            for thresh in CADD_thresh:
                level = level + f", 'CADD_{thresh}'"
    level = level + ", 'Integrated_tools'"
    strategy_script = strategy_script + f"data_set$strategy<-factor(data_set$strategy,levels=c({level}))"

    ## 输出所有出现过的分类名称，用来命名
    count_class = ['TOP1', 'TOP5', 'TOP10', "TOP20", 'TOP50', 'TOP100', '>100', 'not ranked']
    if intersect:
        count_class = count_class + ['not_in_filtered_tsv']
        if len(REVEL_thresh) != 0:
            for thresh in REVEL_thresh:
                count_class = count_class + [f'filtered_by_REVEL_{thresh}']
        if len(CADD_thresh) != 0:
            for thresh in CADD_thresh:
                count_class = count_class + [f'filtered_by_CADD_{thresh}']
    class_names = 'rowname = c(' + "".join(",'" + k + "'" for k in count_class)[1:] + ')'

    print(f'copy this class name to R script:\n{class_names}\n\n{strategy_script}\n\n{tool_script}')


def collect_variants(pwd,df, filename='scoliosis_gVCF_from_zs_updating'):
    if 'collected_variants.tsv' in os.listdir(f'{pwd}/{filename}'):
        os.system(f'rm -f {pwd}/{filename}/collected_variants.tsv')
    for i in tqdm(range(len(df))):
        symbol = df.loc[i, "Symbol"]
        ensemblid = df.loc[i, 'Ensembl Gene ID']
        print(ensemblid)
        case_id = df.loc[i, 'Blood ID']
        if get_case_id_file(case_id, 'sporadic',pwd,filename,filetype='tsv') != 0:
            filedir = get_case_id_file(case_id, 'sporadic',pwd,filename,filetype='tsv')
        else:
            if get_case_id_file(case_id, 'trio',pwd,filename,filetype='tsv') != 0:
                filedir = get_case_id_file(case_id, 'trio',pwd,filename,filetype='tsv')
            else:
                print('no tsv file')
                continue
        variation = pd.read_csv(filedir, sep="\t")
        filtered_result = variation[variation.Gene_name.isin([symbol])].reset_index()
        filtered_result['CaseID'] = [case_id for gene in filtered_result['Gene_name']]
        cl = filtered_result.columns[-1:].tolist()
        cl = cl + filtered_result.columns[:-1].tolist()
        filtered_result = filtered_result[cl]
        file = open(f"./{filename}/collected_variants.tsv", "a")
        for k in range(len(filtered_result)):
            line_first = ''
            if i == 0:
                for column in filtered_result.columns:
                    line_first = line_first + '\t' + column
                line = line_first + '\n'
            else:
                line = ''
            for content in list(filtered_result.loc[k]):
                line = line + '\t' + str(content)
            file.write('%s\n' % line)
            print(file)
        file.close()


def search_organ(hpo_term):
    url = 'http://purl.obolibrary.org/obo/hp.obo'
    graph = obonet.read_obo(url)
    searchgroup = networkx.descendants(graph, hpo_term)
    return searchgroup
import numpy
def readobo(df,hpo,df_organ,save):
    # Read the taxrank ontology
    groupnames = df_organ['Term']
    if f'{hpo}_organ_system' not in df.columns:
        df[f'{hpo}_organ_system']= ""
        df[f'{hpo}_organ_system_number'] = np.nan
    for j in tqdm(range(len(df))):
        if str(df.loc[j, f'{hpo}_organ_system'])=='nan':
            try:
                groupname = []
                if df.loc[j, hpo][0] == '[':
                    hpo_id = df.loc[j, hpo][2:-2]
                    hpo_id_input = [k for k in hpo_id.split("', '")]
                else:
                    hpo_id = df.loc[j, hpo]
                    hpo_id_input = [k for k in hpo_id.split(";")]
                for hpo_term in hpo_id_input:
                    searchgroup = search_organ(hpo_term)
                    groupname = groupname+[k for k in groupnames if k in searchgroup]
                groupname = list(set(groupname))
                df.loc[j,f'{hpo}_organ_system'] = "".join(['+'+name for name in groupname])[1:]
                df.loc[j,f'{hpo}_organ_system_number'] = int(len(groupname))
                print(groupname)
            except Exception as e:
                print(e)
            if save:
                df.to_csv('TF_IDF_weight_organ_2.tsv', sep='\t', index=False)
    return df
    # # Number of nodes
    # print(len(graph))
    # print(networkx.is_directed_acyclic_graph(graph))
    # print(networkx.descendants(graph,'HP:0000118'))

def cnv_FA(pwd,filename,df,refresh):
    # df = df[:1].reset_index(drop=True)
    if 'cnv_function_annotation' not in os.listdir(f'{pwd}/{filename}'):
        os.system(f'mkdir {pwd}/{filename}/cnv_function_annotation')
    for i in tqdm(range(len(df))):
        sequence_id = df['sequence ID'][i]
        cnvnator_dir = f'{pwd}/{filename}/cnv_function_annotation/{sequence_id}.cnvnator'
        if f'{sequence_id}.cnvnator' not in os.listdir(f'{pwd}/{filename}/cnv_function_annotation/'):
            print(f'yaki error: {sequence_id}.cnvnator not in {pwd}/{filename}/cnv_function_annotation/')
        else:
            if (f'{sequence_id}.cnv' not in (os.listdir(f'{pwd}/{filename}/cnv_function_annotation/'))) or refresh:
                cnvnator = pd.read_csv(cnvnator_dir,sep='\t',header=None)
                #cnvnator.set_axis(1,axis=1)
                cnv = pd.DataFrame(columns=['id','chr','start','end','SVType'])
                for k in range(len(cnvnator)):
                    chr = cnvnator.loc[k, 1].split(':')[0][3:]
                    start = (cnvnator.loc[k, 1].split(':')[1]).split('-')[0]
                    end = (cnvnator.loc[k, 1].split(':')[1]).split('-')[1]
                    if cnvnator.loc[k,0]=='deletion':
                        SVtype = 'DEL'
                    else:
                        SVtype='DUP'
                    if start!='1':
                        cnv.loc[k,'id']=f'chr{chr}s{start}e{end}'
                        cnv.loc[k, 'chr'] =chr
                        cnv.loc[k, 'start'] = start
                        cnv.loc[k, 'end'] = end
                        cnv.loc[k, 'SVType'] = SVtype
                cnv.to_csv(f'{pwd}/{filename}/cnv_function_annotation/{sequence_id}.cnv',sep='\t',index=False)
            if f'{sequence_id}.cnv.reg_disruption' not in os.listdir(f'{pwd}/{filename}/cnv_function_annotation/') or refresh:
                ##使用{case_id}.cnv.reg_disruption检查FA是否做了
                cmd = f'zsh && source ~/opt/anaconda3/etc/profile.d/conda.sh && conda activate CNV_FA && cd ~/Software/CNV_FunctionalAnnotation && sh bin/annotate-cnv.sh ~/Software {pwd}/{filename}/cnv_function_annotation/{sequence_id}.cnv'
                print(cmd)
                os.system(cmd)
            if f'{sequence_id}_standard.bed' not in os.listdir(f'{pwd}/{filename}/cnv_function_annotation/') or refresh:
                ##准备标准CNV的bed文件
                cnvnator = pd.read_csv(cnvnator_dir, sep='\t', header=None)
                # cnvnator.set_axis(1,axis=1)
                cnvbed = pd.DataFrame(columns=['#Chrom', 'Start', 'End', 'SV_type'])
                for k in range(len(cnvnator)):
                    chr = cnvnator.loc[k, 1].split(':')[0][3:]
                    start = (cnvnator.loc[k, 1].split(':')[1]).split('-')[0]
                    end = (cnvnator.loc[k, 1].split(':')[1]).split('-')[1]
                    if cnvnator.loc[k, 0] == 'deletion':
                        SVtype = 'DEL'
                    else:
                        SVtype = 'DUP'
                    cnvbed.loc[k, '#Chrom'] = chr
                    cnvbed.loc[k, 'Start'] = start
                    cnvbed.loc[k, 'End'] = end
                    cnvbed.loc[k, 'SV_type'] = SVtype
                cnvbed.to_csv(f'{pwd}/{filename}/cnv_function_annotation/{sequence_id}_standard.bed', sep='\t', index=False)
    print('cnvFA done')
    return True

def subgroup(df,tools,pwd,filename):
    result_file =f'hpo_id_{len(tools)}_rr_len_{len(df)}'
    result_file_weight= f'Weight_{len(tools)}_rr_len_{len(df)}'
    analysis_result_2 = pd.read_csv(f'{result_file}.tsv', sep='\t')
    analysis_result = pd.read_csv(f'{result_file_weight}.tsv', sep='\t')
    result_annot = analysis_result_2
    for k in tqdm(range(len(analysis_result))):
        noi_wei = analysis_result.loc[k, 'all_hpo_plus_weight_phenoapt_rank']
        n_noi_wei = analysis_result.loc[k, 'only_weight_hpo_weight_phenoapt_rank']
        n_noi_n_wei = analysis_result.loc[k, 'only_weight_hpo_no_weight_phenoapt_rank']
        noi_n_wei = analysis_result_2.loc[k, 'phenoapt_rank']
        analysis_result.loc[k,'phenoapt_rank'] = n_noi_n_wei
        noi = None
        wei = None
        cliInd = None
        concl = None
        if df.loc[k,'hpo_id'] == df.loc[k,'Weight']:
            result_annot.loc[k,'set noise'] = 'N'
            if noi_wei>=noi_n_wei and n_noi_wei>=n_noi_n_wei:
                wei = 1
            if noi_wei<=noi_n_wei and n_noi_wei<=n_noi_n_wei:
                wei = -1
            if noi_wei==noi_n_wei and n_noi_wei==n_noi_n_wei:
                wei = 0
        else:
            result_annot.loc[k, 'set noise'] = 'Y'
            if noi_wei>=noi_n_wei and n_noi_wei>=n_noi_n_wei:
                wei = 1
            if noi_wei<=noi_n_wei and n_noi_wei<=n_noi_n_wei:
                wei = -1
            if noi_wei>=n_noi_wei and noi_n_wei>=n_noi_n_wei:
                noi = 1
            if noi_wei<=n_noi_wei and noi_n_wei<=n_noi_n_wei:
                noi = -1
            if noi_wei==noi_n_wei and n_noi_wei==n_noi_n_wei:
                wei=0
            if noi_wei==n_noi_wei and noi_n_wei==n_noi_n_wei:
                noi=0
            if (noi is not None) and (wei is not None):
                if noi == 1: cliInd = wei
                if noi == 0: cliInd = wei
                if noi == -1: cliInd = 2
                ##2=仅cliInd HPO属于该疾病，可能是综合征=0，HPO是综合征里的高频特征表型=1，重点HPO判断可能有误=-1
                cliInd = int(cliInd)
                noi = int(noi)
                wei = int(wei)
                if (cliInd,noi,wei) == (-1,1,-1): concl = 1
                if (cliInd, noi, wei) == (0, 1, 0) or (cliInd, noi, wei) == (0, 1, 1) or (cliInd, noi, wei) == (1, 1, 1): concl = 4
                if (cliInd, noi, wei) == (0,0,0): concl = 3
                if (cliInd, noi, wei) == (2,-1,0) or (cliInd, noi, wei) == (2,-1,1): concl = 2
                if (noi, wei) == (0,-1) or (noi, wei) == (-1,-1):concl = '加权数据有问题'
                #医生clinInd判断错误=1；拓展新表型可能=2；综合征+补充表型可能或者为高度特异表型综合症=3；综合征+所有输入表型已知=4
        if df.loc[k, 'hpo_id'][0] == '[':
            hpo_id = df.loc[k, 'hpo_id'][2:-2]
            hpo_id_input = [hpo for hpo in hpo_id.split("', '")]
        else:
            hpo_id = df.loc[k, 'hpo_id']
            hpo_id_input = [hpo for hpo in hpo_id.split(";")]
        hpo_all = hpo_id_input
        if df.loc[k, 'Weight'][0] == '[':
            hpo_id = df.loc[k, 'Weight'][2:-2]
            hpo_id_input = [hpo for hpo in hpo_id.split("', '")]
        else:
            hpo_id = df.loc[k, 'Weight']
            hpo_id_input = [hpo for hpo in hpo_id.split(";")]
        hpo_weight = hpo_id_input
        noi_hpo_dict =[hpo_term for hpo_term in hpo_all if hpo_term not in list(hpo_weight)]
        noi_hpo = ''
        for hpo_term in noi_hpo_dict:
            noi_hpo = noi_hpo+str(hpo_term)+';'
        result_annot.loc[k,'noise HPO'] = noi_hpo[:-1]
        result_annot.loc[k, 'noi'] = noi
        result_annot.loc[k, 'wei'] = wei
        result_annot.loc[k, 'clinical indication'] = cliInd
        result_annot.loc[k, 'conclusion'] = concl
        enterzid = df.loc[k,'entrezId']
        new_link_hpo_dict =[]
        print(enterzid)
        print(noi_hpo_dict)
        for hpo_term in noi_hpo_dict:
            if search_hpo_annotation(enterzid,hpo_term,pwd,filename)==1:
                new_link_hpo_dict=new_link_hpo_dict+[hpo_term]
        new_link_hpo = ''
        for hpo_term in new_link_hpo_dict:
            new_link_hpo = new_link_hpo+str(hpo_term)+';'
        result_annot.loc[k,'new link hpo'] = new_link_hpo[:-1]
    analysis_result['conclusion'] = result_annot['conclusion']
    analysis_result = analysis_result[analysis_result['conclusion'].notna()]
    result_annot = result_annot[result_annot['conclusion'].notna()]
    result_annot.to_csv(f'{result_file}_subgrouped.tsv',sep='\t',index=False)
    analysis_result.to_csv(f'{result_file_weight}_subgrouped.tsv',sep='\t',index=False)

import subprocess
def search_hpo_annotation(enterzid,hpo,pwd,filename):
    search = str(int(enterzid))+hpo
    print(search)
    call = subprocess.call(f'grep -x "{search}" {pwd}/{filename}/enterzidhpo.txt | echo $?', shell=True)
    return call

def compare_df_for_R(df,tools,pwd,mode,dMRR_heatmap=False):
    if mode=='hpo_id_Weight':
        hpo_id_result = pd.read_csv(f'{pwd}/hpo_id_{len(tools)}_rr_len_{len(df)}_subgrouped.tsv',sep='\t')
        Weight_result= pd.read_csv(f'{pwd}/Weight_{len(tools)}_rr_len_{len(df)}_subgrouped.tsv',sep='\t')
    if mode=='wilcoxon':
        hpo_id_result = pd.read_csv(f'{pwd}/hpo_id_{len(tools)}_rr_with_vcf_len_{len(df)}.tsv', sep='\t')
        Weight_result = pd.read_csv(f'{pwd}/Weight_{len(tools)}_rr_with_vcf_len_{len(df)}.tsv', sep='\t')
    column = ['CaseID','hpo']
    if mode=='':
        column=column+['conclusion']
    for tool in tools:
        column = column + [f'{tool}_rank']
        if mode=='wilcoxon':
            if tool in ['phenoapt','phrank', 'phenolyzer', 'GADO', 'phen2gene','PhenoRank','HANRD','phenomizer']:
                column = column + [f'{tool}_intersect_rank'] +[f'{tool}_rank_CADD_15']+[f'{tool}_rank_REVEL_0.75']
    combine_result = pd.DataFrame(columns=column)
    for k in range(len(hpo_id_result)):
        for j in column:
            combine_result.loc[k,j] = hpo_id_result.loc[k,j]
    for k in range(len(Weight_result)):
        for j in column:
            if j == 'phenoapt_rank':
                combine_result.loc[k + len(hpo_id_result), j] = Weight_result.loc[k, 'only_weight_hpo_no_weight_phenoapt_rank']
            if j == 'phenoapt_intersect_rank':
                combine_result.loc[k + len(hpo_id_result), j] = Weight_result.loc[k, 'only_weight_hpo_no_weight_phenoapt_intersect_rank']
            if j not in ['phenoapt_rank','phenoapt_intersect_rank']:
                if 'CADD' not in j and 'REVEL' not in j:
                    combine_result.loc[k+len(hpo_id_result), j] = Weight_result.loc[k, j]
    if mode=='hpo_id_Weight':
        combine_result.to_csv(f'combine_subgroup_{len(combine_result)}.tsv',sep='\t',index=False)
    if mode=='wilcoxon':
        combine_result.to_csv(f'combine_subgroup_tools_{len(combine_result)}_intersect.tsv', sep='\t', index=False)
    if dMRR_heatmap:
        dMRR_table(tools,combine_result)


def df0_phenoapt_4_strategy_rr(df,tools,pwd):
    df0 = pd.DataFrame(columns=['patients','Weight_strategy','Reciprocal_rank'])
    Weight_result_all = pd.read_csv(f'{pwd}/Weight_{len(tools)}_rr_len_{len(df)}.tsv',sep='\t')
    hpo_id_result = pd.read_csv(f'{pwd}/hpo_id_{len(tools)}_rr_len_{len(df)}.tsv', sep='\t')
    strategy = ['Clinical indication','Clinical indication+weight','All hpo+weight','All hpo']
    result_strategy = ['only_weight_hpo_no_weight_phenoapt_rank','only_weight_hpo_weight_phenoapt_rank','all_hpo_plus_weight_phenoapt_rank','phenoapt_rank']
    for k in range(len(Weight_result_all)):
        m=k + k*len(strategy)
        for n in range(len(strategy)):
            df0.loc[m+n,'Weight_strategy'] = strategy[n]
            df0.loc[m+n, 'patients'] = Weight_result_all.loc[k, 'CaseID']
            if n != (len(strategy) - 1):
                df0.loc[m+n,'Reciprocal_rank'] = Weight_result_all.loc[k, result_strategy[n]]
            else:
                df0.loc[m + n, 'Reciprocal_rank'] = hpo_id_result.loc[k, result_strategy[n]]
    df0.to_csv('df0_phenoapt_4_strategy_rr.tsv',sep='\t',index=False)

def df_for_all_strategy_hpo_id_rr(df,tools,pwd):##小提琴图不错
    df0 = pd.DataFrame(columns=['patients','Weight_strategy','Reciprocal_rank'])
    #Weight_result_all = pd.read_csv(f'{pwd}/Weight_{len(tools)}_rr_len_{len(df)}.tsv',sep='\t')
    hpo_id_result = pd.read_csv(f'{pwd}/hpo_id_{len(tools)}_rr_with_vcf_len_{len(df)}.tsv', sep='\t')
    strategy = ['Clinical indication','Clinical indication+weight','All hpo+weight','All hpo']
    result_strategy = ['only_weight_hpo_no_weight_phenoapt_rank','only_weight_hpo_weight_phenoapt_rank','all_hpo_plus_weight_phenoapt_rank','phenoapt_rank']
    for k in range(len(hpo_id_result)):
        m=k + k*len(strategy)
        for n in range(len(strategy)):
            df0.loc[m+n,'Weight_strategy'] = strategy[n]
            df0.loc[m+n, 'patients'] = hpo_id_result.loc[k, 'CaseID']
            if n != (len(strategy) - 1):
                df0.loc[m+n,'Reciprocal_rank'] = hpo_id_result.loc[k, result_strategy[n]]
            else:
                df0.loc[m + n, 'Reciprocal_rank'] = hpo_id_result.loc[k, result_strategy[n]]
    df0.to_csv('df_for_all_strategy_hpo_id_rr.tsv',sep='\t',index=False)


def dMRR_table(tools,combine_result):
    dMRR_table = pd.DataFrame(columns=['tools','Group1', 'Group2', 'Group3', 'Group4'])
    for group in [1,2,3,4]:
        for n in range(len(tools)):
            subset=combine_result[combine_result['conclusion']==group]
            subset_hpo_id = subset[subset['hpo']=='hpo_id']
            subset_weight = subset[subset['hpo']=='Weight']
            dMRR_table.loc[n,'tools'] = tools[n]
            dMRR_table.loc[n, f'Group{group}'] = subset_hpo_id[f'{tools[n]}_rank'].mean()-subset_weight[f'{tools[n]}_rank'].mean()
    dMRR_table.to_csv(f'dMRR_{len(combine_result)}.tsv',sep='\t',index=False)
    ##做热图用的dMRR数据，列名固定，需要调整的话，得在R中同步修改

def searchTF_IDF_for_R(df,tools,pwd,filename):
    intrisic_weight_df = pd.read_csv(f'{pwd}/weight.csv')
    intrisic_weight_df = intrisic_weight_df.set_index('hpo_id')
    pus_hpo_df = pd.read_csv(f'hpo_id_{len(tools)}_rr_len_{len(df)}_subgrouped.tsv',sep='\t')
    #pus_hpo_df = pus_hpo_df[11:12].reset_index()
    pus_hpo_tf_idf = pd.DataFrame(columns=['hpo','groups','tf-idf','gene_anno_number'])
    m=0
    for k in tqdm(range(len(pus_hpo_df))):
        try:
            hpo_id = pus_hpo_df.loc[k, 'noise HPO']
            hpo_id_input = [hpo for hpo in hpo_id.split(";")]
            for hpo_n in range(len(hpo_id_input)):
                weight = intrisic_weight_df.loc[hpo_id_input[hpo_n],'intrinsic_weight']
                m = m+hpo_n
                pus_hpo_tf_idf.loc[m, 'hpo'] = hpo_id_input[hpo_n]
                pus_hpo_tf_idf.loc[m,'groups'] = 'Group'+str(int(pus_hpo_df.loc[k, 'conclusion']))
                print(pus_hpo_df.loc[k, 'conclusion'])
                pus_hpo_tf_idf.loc[m, 'tf-idf'] = weight
                cmd =f"grep -E '{hpo_id_input[hpo_n]}' {pwd}/{filename}/hpo_count.tsv"
                count = subprocess.getoutput(cmd)
                if count!='':
                    pus_hpo_tf_idf.loc[m,'gene_anno_number'] = int((count).split(' ')[-2])
                else:
                    pus_hpo_tf_idf.loc[m, 'gene_anno_number'] = None
            m =m+1
        except Exception as e:
            print(e)
    pus_hpo_tf_idf.to_csv('pus_TF_IDF.tsv',sep='\t',index=False)
