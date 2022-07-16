# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


import appscript as appscript
import pandas as pd
import xlrd
from phenoapt import PhenoApt
import numpy as np
import subprocess
import re
import json
from suds import client


from phrank import Phrank
import urllib.request as urlre
from functools import reduce
import requests, sys
from openpyxl import load_workbook
from tqdm import tqdm
import appscript
import time
import allel
from functions import *


def main_2(tools=[],CADD_thresh=[],REVEL_thresh=[],hpo='hpo_id',intersect=False,statistic = True,refresh=False,newcase=False): ##无精炼HPO版本

    filename = 'scoliosis_gVCF_from_zs_updating'
    pwd = '/Users/liyaqi/PycharmProjects/PhenoAptAnalysis'

    df_organ = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx',' organ_system')

    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    #df = df[:1].reset_index(drop=True)
    case_ids = df['Blood ID']

    print(f"{len(df)}")
    if newcase:
        ##所有预先需要的文件准备
        file_prepare(df, pwd, filename, REVEL_thresh=REVEL_thresh, CADD_thresh=CADD_thresh, refresh=refresh)
        if 'phen2gene' in tools:
            phen2gene_rank(df,pwd=pwd, filename=filename, profile='no_gene_list',refresh=refresh)
            ## phen2gene直接rank版本
            if intersect:
                phen2gene_rank(df,pwd=pwd, filename=filename, profile=f'candidategene',refresh=refresh)
            if len(REVEL_thresh)!=0:
                for thresh in REVEL_thresh:
                    phen2gene_rank(df,pwd=pwd, filename=filename, profile=f'revel_{thresh}',refresh=refresh)
            ## phen2gene重新跑REVEL
            if len(CADD_thresh)!=0:
                for thresh in CADD_thresh:
                    phen2gene_rank(df, pwd=pwd, filename=filename, profile=f'cadd_{thresh}',refresh=refresh)
            ## phen2gene重新跑CADD
            print('phen2gene ready')
        if 'LIRICAL' in tools:
            generate_phenopaket_zs_diag_zs_gVCF(df,pwd,filename,hpo,refresh=refresh)
            ## 新的hpo，新的case
            print('LIRICAL ready')
        if 'Exomiser' in tools:
            generateyml_zs_diag_zs_gVCF(df,pwd,filename,hpo,refresh=refresh)
            print('Exomiser ready')
            ## 新的hpo，新的case
        if 'GADO' in tools:
            GADO(df,pwd,filename,hpo,refresh=refresh)
            ## 每次增加case都会重新生成hpototal，相当于每次都refresh了hpo；新的hpo
        if 'phenolyzer' in tools:
            phenolyzer(df,pwd,filename,hpo,refresh=refresh)
            # 注意关闭terminal；各个软件结果都有缓存，有新case时才设置True;需要全部cohort重新生成hpotxt的时候再refresh，单纯增加了case会直接加上,单纯增加一种hpo输入，可直接输入在hpo里面; 需要更新hpo的时候，不论完整版还是weight版任何一行做了修改，也是要refresh整个cohort的
        if 'phrank'in tools:
            phrank_rank(df=df,pwd=pwd,filename=filename, mode='OMIM', hpo=hpo, refresh=refresh)
            ## phrank直接rank版本
            if intersect:
                phrank_rank(df,pwd,filename, mode='intersect', hpo=hpo, filter='candidategene_ensemblid', refresh=refresh)
            if len(REVEL_thresh) != 0:
                for thresh in REVEL_thresh:
                    phrank_rank(df,pwd,filename, mode='intersect', hpo=hpo, filter=f'revel_{thresh}_ensemblid', refresh=refresh)
            ## phrank检查是否有漏跑REVEL
            if len(CADD_thresh) != 0:
                for thresh in CADD_thresh:
                    phrank_rank(df,pwd,filename, mode='intersect', hpo=hpo, filter=f'cadd_{thresh}_ensemblid', refresh=refresh)
            ## phrank检查是否有漏跑CADD
            print('phrank ready')

        if 'phenoapt' in tools:
            if hpo != 'Weight':
                phenoapt_rank(df, pwd, filename, hpo, refresh)
            else:
                phenoapt_rank(df, pwd, filename, hpo, refresh=refresh)#只有重点hpo，不加权
                phenoapt_rank(df,pwd,filename,hpo,refresh=refresh, all_hpo_plus_weight=True)#有所有hpo，加权重点hpo
                phenoapt_rank(df, pwd, filename, hpo, refresh=refresh, only_weight_hpo_weight=True)#只有重点hpo，加权重点hpo


    def statistic_to_table(hpo):
        columns_fill = ['CaseID', 'hpo', 'hpo_id_input', 'Symbol']
        for tool in tools:
            if tool in [ 'phrank', 'phenolyzer', 'GADO', 'phen2gene']:
                columns_fill = columns_fill + [f'{tool}_rank']
                if intersect: columns_fill = columns_fill + [f'{tool}_intersect_rank']
                if len(REVEL_thresh) != 0:
                    for thresh in REVEL_thresh:
                        columns_fill = columns_fill + [f'{tool}_rank_REVEL_{thresh}']
                if len(CADD_thresh) != 0:
                    for thresh in CADD_thresh:
                        columns_fill = columns_fill + [f'{tool}_rank_CADD_{thresh}']
            else:
                if tool =='phenoapt':
                    if hpo =='Weight':
                        for k in ['all_hpo_plus_weight','only_weight_hpo_weight','only_weight_hpo_no_weight']:
                            columns_fill = columns_fill + [f'{k}_phenoapt_rank']
                            if intersect: columns_fill = columns_fill + [f'{k}_phenoapt_intersect_rank']
                            if len(REVEL_thresh) != 0:
                                for thresh in REVEL_thresh:
                                    columns_fill = columns_fill + [f'{k}_phenoapt_rank_REVEL_{thresh}']
                            if len(CADD_thresh) != 0:
                                for thresh in CADD_thresh:
                                    columns_fill = columns_fill + [f'{k}_phenoapt_rank_CADD_{thresh}']
                    else:
                        columns_fill = columns_fill + [f'{tool}_rank']
                        if intersect: columns_fill = columns_fill + [f'{tool}_intersect_rank']
                        if len(REVEL_thresh) != 0:
                            for thresh in REVEL_thresh:
                                columns_fill = columns_fill + [f'{tool}_rank_REVEL_{thresh}']
                        if len(CADD_thresh) != 0:
                            for thresh in CADD_thresh:
                                columns_fill = columns_fill + [f'{tool}_rank_CADD_{thresh}']
                else:
                    columns_fill = columns_fill + [f'{tool}_rank']

        part_table = pd.DataFrame(columns=columns_fill)
        part_rr_table = pd.DataFrame(columns=columns_fill)
        part_table_32 = pd.DataFrame(columns=columns_fill)
        part_rr_table_32 = pd.DataFrame(columns=columns_fill)
        case_id_exomiser_tsv_file_dict = get_case_id_integ_file_map(pwd, filename, 'Exomiser', case_ids, hpo)
        case_id_LIRICAL_tsv_file_dict = get_case_id_integ_file_map(pwd, filename, 'LIRICAL', case_ids, hpo)
        for i in tqdm(range(len(df))):
            if df.loc[i, hpo][0] == '[':
                hpo_id = df.loc[i, hpo][2:-2]
                hpo_id_input = [k for k in hpo_id.split("', '")]
            else:
                hpo_id = df.loc[i, hpo]
                hpo_id_input = [k for k in hpo_id.split(";")]
            symbol = df.loc[i, "Symbol"]
            entrezGeneId = df.loc[i, 'entrezGeneId']
            ## NCBIGene:6911 entrezGeneId的格式
            entrezId = df.loc[i, 'entrezId']
            ensemblid = df.loc[i, 'Ensembl Gene ID']
            print(entrezId)
            case_id = df.loc[i, 'Blood ID']
            if get_case_id_file(case_id, 'sporadic') != 0:
                filedir = get_case_id_file(case_id, 'sporadic')
            else:
                if get_case_id_file(case_id, 'trio') != 0:
                    filedir = get_case_id_file(case_id, 'trio')
                else:
                    print('no tsv file')
                    continue
            variation = pd.read_csv(filedir, sep="\t")
            variation_gene_name = variation['Gene_name']
            part_table.loc[i, 'CaseID'] = case_id
            part_table.loc[i, 'hpo'] = hpo
            part_table.loc[i, 'hpo_id_input'] = f'{hpo_id_input}'
            part_table.loc[i, 'Symbol'] = symbol
            part_rr_table.loc[i, 'CaseID'] = case_id
            part_rr_table.loc[i, 'hpo'] = hpo
            part_rr_table.loc[i, 'hpo_id_input'] = f'{hpo_id_input}'
            part_rr_table.loc[i, 'Symbol'] = symbol
            if 'phrank' in tools:
                phrank_gene_ranking = pd.read_csv(
                    f'{pwd}/{filename}/phrankoutput/{case_id}_{hpo}_nogenelist_phrank_rank.tsv', sep='\t')
                phrank_gene_ranking = list(phrank_gene_ranking['1'])
                part_table.loc[i, 'phrank_rank'], part_rr_table.loc[i, 'phrank_rank'] = get_rank(ensemblid,
                                                                                                 phrank_gene_ranking)
                if intersect:
                    phrank_gene_ranking = pd.read_csv(
                        f'{pwd}/{filename}/phrankoutput/{case_id}_{hpo}_candidategene_ensemblid_phrank_rank.tsv',
                        sep='\t')
                    phrank_gene_ranking = list(phrank_gene_ranking['1'])
                    part_table.loc[i, 'phrank_intersect_rank'], part_rr_table.loc[
                        i, 'phrank_intersect_rank'] = get_rank(
                        ensemblid, phrank_gene_ranking)
                if len(REVEL_thresh) != 0:
                    for thresh in REVEL_thresh:
                        phrank_gene_ranking = pd.read_csv(
                            f'{pwd}/{filename}/phrankoutput/{case_id}_{hpo}_revel_{thresh}_ensemblid_phrank_rank.tsv',
                            sep='\t')
                        phrank_gene_ranking = list(phrank_gene_ranking['1'])
                        part_table.loc[i, f'phrank_rank_REVEL_{thresh}'], part_rr_table.loc[
                            i, f'phrank_rank_REVEL_{thresh}'] = get_rank(ensemblid, phrank_gene_ranking)
                if len(CADD_thresh) != 0:
                    for thresh in CADD_thresh:
                        phrank_gene_ranking = pd.read_csv(
                            f'{pwd}/{filename}/phrankoutput/{case_id}_{hpo}_cadd_{thresh}_ensemblid_phrank_rank.tsv',
                            sep='\t')
                        phrank_gene_ranking = list(phrank_gene_ranking['1'])
                        part_table.loc[i, f'phrank_rank_CADD_{thresh}'], part_rr_table.loc[
                            i, f'phrank_rank_CADD_{thresh}'] = get_rank(ensemblid, phrank_gene_ranking)

            if 'phenolyzer' in tools:
                if hpo == 'hpo_id':
                    phenolyzer_df = pd.read_csv(
                        f'{pwd}/{filename}/phenolyzeroutput/{case_id}.final_gene_list',
                        sep='\t')
                    phenolyzer_result = list(phenolyzer_df['ID'])
                else:
                    phenolyzer_df = pd.read_csv(
                        f'{pwd}/{filename}/phenolyzeroutput/{case_id}_{hpo}.final_gene_list',
                        sep='\t')
                    phenolyzer_result = list(phenolyzer_df['ID'])
                part_table.loc[i, 'phenolyzer_rank'], part_rr_table.loc[i, 'phenolyzer_rank'] = get_rank(entrezId,
                                                                                                         phenolyzer_result)
                if intersect:
                    phenolyzer_intersect_gene_name = (
                        phenolyzer_df[phenolyzer_df.Gene.isin(variation_gene_name)]).reset_index()
                    phenolyzer_intersect_gene_name = list(phenolyzer_intersect_gene_name['Gene'])
                    part_table.loc[i, 'phenolyzer_intersect_rank'], part_rr_table.loc[
                        i, 'phenolyzer_intersect_rank'] = get_rank(symbol, phenolyzer_intersect_gene_name)
                if len(REVEL_thresh) != 0:
                    for thresh in REVEL_thresh:
                        tsv_revel_filter_dir = f'./{filename}/revel_{thresh}/{case_id}.txt'
                        tsv_revel_filter = pd.read_csv(tsv_revel_filter_dir, sep='\t', header=None)
                        phenolyzer_intersect_gene_name = (
                            phenolyzer_df[phenolyzer_df.Gene.isin(tsv_revel_filter[0])]).reset_index()
                        phenolyzer_intersect_gene_name = list(phenolyzer_intersect_gene_name['Gene'])
                        part_table.loc[i, f'phenolyzer_rank_REVEL_{thresh}'], part_rr_table.loc[
                            i, f'phenolyzer_rank_REVEL_{thresh}'] = get_rank(symbol,phenolyzer_intersect_gene_name)
                if len(CADD_thresh) != 0:
                    for thresh in CADD_thresh:
                        tsv_cadd_filter_dir = f'./{filename}/cadd_{thresh}/{case_id}.txt'
                        tsv_cadd_filter = pd.read_csv(tsv_cadd_filter_dir, sep='\t', header=None)
                        phenolyzer_intersect_gene_name = (
                            phenolyzer_df[phenolyzer_df.Gene.isin(tsv_cadd_filter[0])]).reset_index()
                        phenolyzer_intersect_gene_name = list(phenolyzer_intersect_gene_name['Gene'])
                        part_table.loc[i, f'phenolyzer_rank_CADD_{thresh}'], part_rr_table.loc[
                            i, f'phenolyzer_rank_CADD_{thresh}'] = get_rank(symbol,phenolyzer_intersect_gene_name)

            if 'GADO' in tools:
                GADO_result = pd.read_csv(f'{pwd}/{filename}/GADOoutput/{case_id}.txt', sep='\t')
                result_ensemblid = list(GADO_result['Ensg'])
                part_table.loc[i, 'GADO_rank'], part_rr_table.loc[i, 'GADO_rank'] = get_rank(ensemblid,
                                                                                             result_ensemblid)
                if intersect:
                    ensemblid_candidate_gene = pd.read_csv(f'./{filename}/candidategene_ensemblid/{case_id}.txt',
                                                           sep='\t', header=None)
                    GADO_intersect_gene = (
                    GADO_result[GADO_result.Ensg.isin(ensemblid_candidate_gene[0])]).reset_index()
                    GADO_intersect_gene = list(GADO_intersect_gene['Ensg'])
                    part_table.loc[i, 'GADO_intersect_rank'], part_rr_table.loc[
                        i, 'GADO_intersect_rank'] = get_rank(ensemblid, GADO_intersect_gene)
                if len(REVEL_thresh) != 0:
                    for thresh in REVEL_thresh:
                        tsv_revel_filer_ensemblid_dir = f'./{filename}/revel_{thresh}_ensemblid/{case_id}.txt'
                        tsv_revel_filer_ensemblid = pd.read_csv(tsv_revel_filer_ensemblid_dir, sep='\t', header=None)
                        GADO_intersect_gene = (
                            GADO_result[GADO_result.Ensg.isin(tsv_revel_filer_ensemblid[0])]).reset_index()
                        GADO_intersect_gene = list(GADO_intersect_gene['Ensg'])
                        part_table.loc[i, f'GADO_rank_REVEL_{thresh}'], part_rr_table.loc[
                            i, f'GADO_rank_REVEL_{thresh}'] = get_rank(ensemblid, GADO_intersect_gene)
                if len(CADD_thresh) != 0:
                    for thresh in CADD_thresh:
                        tsv_CADD_filer_ensemblid_dir = f'./{filename}/cadd_{thresh}_ensemblid/{case_id}.txt'
                        tsv_CADD_filer_ensemblid = pd.read_csv(tsv_CADD_filer_ensemblid_dir, sep='\t', header=None)
                        GADO_intersect_gene = (
                            GADO_result[GADO_result.Ensg.isin(tsv_CADD_filer_ensemblid[0])]).reset_index()
                        GADO_intersect_gene = list(GADO_intersect_gene['Ensg'])
                        part_table.loc[i, f'GADO_rank_CADD_{thresh}'], part_rr_table.loc[
                            i, f'GADO_rank_CADD_{thresh}'] = get_rank(ensemblid, GADO_intersect_gene)

            if 'phen2gene' in tools:
                phen2gene_result = pd.read_csv(
                    f'{pwd}/{filename}/phen2geneoutput_nolist/{case_id}-phen2gene-nolist-output.txt', sep="\t")
                phen2gene_rank = list(phen2gene_result["ID"])
                part_table.loc[i, 'phen2gene_rank'], part_rr_table.loc[i, 'phen2gene_rank'] = get_rank(entrezId,
                                                                                                       phen2gene_rank)
                if intersect:
                    phen2gene_intersect_result = pd.read_csv(
                        f'{pwd}/{filename}/phen2geneoutput/{case_id}-phen2gene-candidategene-output.txt', sep="\t")
                    phen2gene_intersect_result = list(phen2gene_intersect_result['Gene'])
                    part_table.loc[i, 'phen2gene_intersect_rank'], part_rr_table.loc[
                        i, 'phen2gene_intersect_rank'] = get_rank(symbol, phen2gene_intersect_result)
                if len(REVEL_thresh) != 0:
                    for thresh in REVEL_thresh:
                        phen2gene_intersect_result = pd.read_csv(
                            f'{pwd}/{filename}/phen2geneoutput/{case_id}-phen2gene-revel_{thresh}-output.txt', sep="\t")
                        phen2gene_intersect_result = list(phen2gene_intersect_result['Gene'])
                        part_table.loc[i, f'phen2gene_rank_REVEL_{thresh}'], part_rr_table.loc[
                            i, f'phen2gene_rank_REVEL_{thresh}'] = get_rank(symbol, phen2gene_intersect_result)
                        ##输入列表就是基因名称，所以可能不存在排出来的基因名称、数目和Tsv全体candidate gene不一致的情况。中间有做symbol转换吗可以看看文章

                if len(CADD_thresh) != 0:
                    for thresh in CADD_thresh:
                        phen2gene_intersect_result = pd.read_csv(
                            f'{pwd}/{filename}/phen2geneoutput/{case_id}-phen2gene-cadd_{thresh}-output.txt', sep="\t")
                        phen2gene_intersect_result = list(phen2gene_intersect_result['Gene'])
                        part_table.loc[i, f'phen2gene_rank_CADD_{thresh}'], part_rr_table.loc[
                            i, f'phen2gene_rank_CADD_{thresh}'] = get_rank(symbol, phen2gene_intersect_result)

            if 'phenoapt' in tools:
                if hpo!='Weight':
                    pheno_result = pd.read_csv(f'{pwd}/{filename}/phenoaptoutput/{case_id}_{hpo}_phenoapt_rank.tsv',
                                               sep='\t')
                    pheno_rank = list(pheno_result['gene_symbol'])
                    part_table.loc[i, 'phenoapt_rank'], part_rr_table.loc[i, 'phenoapt_rank'] = get_rank(symbol,
                                                                                                         pheno_rank)
                    if intersect:
                        intersect_gene = (
                        pheno_result[pheno_result.gene_symbol.isin(variation_gene_name)]).reset_index()
                        intersect_gene = list(intersect_gene['gene_symbol'])
                        part_table.loc[i, 'phenoapt_intersect_rank'], part_rr_table.loc[
                            i, 'phenoapt_intersect_rank'] = get_rank(symbol, intersect_gene)
                    if len(REVEL_thresh) != 0:
                        for thresh in REVEL_thresh:
                            tsv_revel_filter_dir = f'./{filename}/revel_{thresh}/{case_id}.txt'
                            tsv_revel_filter = pd.read_csv(tsv_revel_filter_dir, sep='\t', header=None)
                            intersect_gene = (
                            pheno_result[pheno_result.gene_symbol.isin(tsv_revel_filter[0])]).reset_index()
                            intersect_gene = list(intersect_gene['gene_symbol'])
                            part_table.loc[i, f'phenoapt_rank_REVEL_{thresh}'], part_rr_table.loc[
                                i, f'phenoapt_rank_REVEL_{thresh}'] = get_rank(symbol, intersect_gene)
                    if len(CADD_thresh) != 0:
                        for thresh in CADD_thresh:
                            tsv_cadd_filter_dir = f'./{filename}/cadd_{thresh}/{case_id}.txt'
                            tsv_cadd_filter = pd.read_csv(tsv_cadd_filter_dir, sep='\t', header=None)
                            intersect_gene = (
                            pheno_result[pheno_result.gene_symbol.isin(tsv_cadd_filter[0])]).reset_index()
                            intersect_gene = list(intersect_gene['gene_symbol'])
                            part_table.loc[i, f'phenoapt_rank_CADD_{thresh}'], part_rr_table.loc[
                                i, f'phenoapt_rank_CADD_{thresh}'] = get_rank(symbol, intersect_gene)
                else:
                    for k in ['all_hpo_plus_weight','only_weight_hpo_weight','only_weight_hpo_no_weight']:
                        def get_phenoapt_result_dir(k):
                            if k == 'all_hpo_plus_weight':
                                pheno_result_dir = f'{pwd}/{filename}/phenoaptoutput/{case_id}_{hpo}_phenoapt_rank.tsv'
                                return pheno_result_dir
                            if k =='only_weight_hpo_no_weight':
                                pheno_result_dir = f'{pwd}/{filename}/phenoaptoutput/{case_id}_only_{hpo}_hpo_no_weight_phenoapt_rank.tsv'
                                return pheno_result_dir
                            if k =='only_weight_hpo_weight':
                                pheno_result_dir = f'{pwd}/{filename}/phenoaptoutput/{case_id}_only_{hpo}_hpo_weight_phenoapt_rank.tsv'
                                return pheno_result_dir
                        pheno_result = pd.read_csv(get_phenoapt_result_dir(k),sep='\t')
                        pheno_rank = list(pheno_result['gene_symbol'])
                        part_table.loc[i, f'{k}_phenoapt_rank'], part_rr_table.loc[i, f'{k}_phenoapt_rank'] = get_rank(symbol,pheno_rank)
                        if intersect:
                            intersect_gene = (
                            pheno_result[pheno_result.gene_symbol.isin(variation_gene_name)]).reset_index()
                            intersect_gene = list(intersect_gene['gene_symbol'])
                            part_table.loc[i, f'{k}_phenoapt_intersect_rank'], part_rr_table.loc[
                                i, f'{k}_phenoapt_intersect_rank'] = get_rank(symbol, intersect_gene)
                        if len(REVEL_thresh) != 0:
                            for thresh in REVEL_thresh:
                                tsv_revel_filter_dir = f'./{filename}/revel_{thresh}/{case_id}.txt'
                                tsv_revel_filter = pd.read_csv(tsv_revel_filter_dir, sep='\t', header=None)
                                intersect_gene = (
                                pheno_result[pheno_result.gene_symbol.isin(tsv_revel_filter[0])]).reset_index()
                                intersect_gene = list(intersect_gene['gene_symbol'])
                                part_table.loc[i, f'{k}_phenoapt_rank_REVEL_{thresh}'], part_rr_table.loc[
                                    i, f'{k}_phenoapt_rank_REVEL_{thresh}'] = get_rank(symbol, intersect_gene)
                        if len(CADD_thresh) != 0:
                            for thresh in CADD_thresh:
                                tsv_cadd_filter_dir = f'./{filename}/cadd_{thresh}/{case_id}.txt'
                                tsv_cadd_filter = pd.read_csv(tsv_cadd_filter_dir, sep='\t', header=None)
                                intersect_gene = (
                                pheno_result[pheno_result.gene_symbol.isin(tsv_cadd_filter[0])]).reset_index()
                                intersect_gene = list(intersect_gene['gene_symbol'])
                                part_table.loc[i, f'{k}_phenoapt_rank_CADD_{thresh}'], part_rr_table.loc[
                                    i, f'{k}_phenoapt_rank_CADD_{thresh}'] = get_rank(symbol, intersect_gene)
                    # filtered_result = variation[variation.Gene_name.isin(intersect_gene['gene_symbol'])]
                    # filtered_result['pheno_rank'] = [pheno_gene_rank[gene] for gene in filtered_result['Gene_name']]
                    # filtered_result['intersect_rank'] = [intersect_gene_rank[gene] for gene in filtered_result['Gene_name']]
                    # filtered_result['CaseID'] = [case_id for gene in filtered_result['Gene_name']]
                    # filtered_result.to_csv(f"output/{case_id}.csv")
            if 'Exomiser' in tools:
                if case_id in case_id_exomiser_tsv_file_dict:
                    rankdf = pd.read_csv(case_id_exomiser_tsv_file_dict[case_id], sep="\t")
                    patho_gene_ID_Exo = list(rankdf['ENTREZ_GENE_ID'])
                    if entrezId in patho_gene_ID_Exo:
                        rank = patho_gene_ID_Exo.index(entrezId)
                        rr = 1 / (1 + rank)
                    else:
                        rank = 'NA'
                        rr = 0
                else:
                    rank = 'No_vcf_ranked_by_Exo'
                    rr = 0
                part_table.loc[i, 'Exomiser_rank'] = rank
                part_rr_table.loc[i, 'Exomiser_rank'] = rr

            if 'LIRICAL' in tools:
                patho_gene_rank_LIR = {}
                if case_id in case_id_LIRICAL_tsv_file_dict:
                    rankdf2 = pd.read_csv(case_id_LIRICAL_tsv_file_dict[case_id], sep="\t")
                    patho_gene_ID_LIR = rankdf2['entrezGeneId']
                    df_count = pd.DataFrame(enumerate(patho_gene_ID_LIR))
                    for s in range(len(df_count)):
                        k = (df_count[df_count[1].isin([df_count.loc[s, 1]])]).reset_index(drop=True)
                        #print(k)
                        patho_gene_rank_LIR[k.loc[0, 1]] = k.loc[0, 0]
                    #print(patho_gene_rank_LIR)

                else:
                    patho_gene_rank_LIR[entrezGeneId] = 'No_vcf_ranked_by_LIR'
                rank = patho_gene_rank_LIR.get(entrezGeneId, 'NA')
                if rank == 'NA' or rank == 'No_vcf_ranked_by_LIR':
                    rr = 0
                else:
                    rr = 1 / (rank + 1)
                part_table.loc[i, 'LIRICAL_rank'] = rank
                part_rr_table.loc[i, 'LIRICAL_rank'] = rr

            ## 转化成作图数据
            if part_table.loc[i,'LIRICAL_rank']!='No_vcf_ranked_by_LIR':
                part_table_32.loc[i] = part_table.loc[i]
                part_rr_table_32.loc[i] = part_rr_table.loc[i]
                for tool in tools:
                    if tool in ['phrank', 'phenolyzer', 'GADO', 'phen2gene']:
                        if part_table_32.loc[i,f'{tool}_rank'] == 'NA':
                            for column in part_table_32.columns:
                                if tool in str(column):
                                    part_table_32.loc[i, column] = 'not ranked'
                        else:
                            if intersect:
                                if part_table_32.loc[i, f'{tool}_intersect_rank'] == 'NA':
                                    for column in part_table_32.columns:
                                        if tool in str(column):
                                            if f'{tool}_rank' not in column:
                                                part_table_32.loc[i, column] = 'not_in_filtered_tsv'
                                else:
                                    if len(REVEL_thresh)!=0:
                                        thresh_na=[]
                                        for thresh in REVEL_thresh:
                                            # part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}']
                                            if part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}'] == 'NA':
                                                thresh_na=thresh_na+[thresh]
                                        if len(thresh_na)!=0:
                                            filter_thresh = min(thresh_na)
                                            for thresh in thresh_na:
                                                part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}'] = f'filtered_by_REVEL_{filter_thresh}'
                                    if len(CADD_thresh)!=0:
                                        thresh_na=[]
                                        for thresh in CADD_thresh:
                                            # part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}']
                                            if part_table_32.loc[i, f'{tool}_rank_CADD_{thresh}'] == 'NA':
                                                thresh_na = thresh_na + [thresh]
                                        if len(thresh_na)!=0:
                                            filter_thresh = min(thresh_na)
                                            for thresh in thresh_na:
                                                part_table_32.loc[i, f'{tool}_rank_CADD_{thresh}'] = f'filtered_by_CADD_{filter_thresh}'

                    if tool == 'phenoapt':
                        if hpo == 'hpo_id':
                            if part_table_32.loc[i, f'{tool}_rank'] == 'NA':
                                for column in part_table_32.columns:
                                    if tool in str(column):
                                        part_table_32.loc[i, column] = 'not ranked'
                            else:
                                if intersect:
                                    if part_table_32.loc[i, f'{tool}_intersect_rank'] == 'NA':
                                        for column in part_table_32.columns:
                                            if tool in str(column):
                                                if f'{tool}_rank' not in column:
                                                    part_table_32.loc[i, column] = 'not_in_filtered_tsv'
                                    else:
                                        if len(REVEL_thresh) != 0:
                                            thresh_na = []
                                            for thresh in REVEL_thresh:
                                                # part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}']
                                                if part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}'] == 'NA':
                                                    thresh_na = thresh_na + [thresh]
                                            if len(thresh_na)!=0:
                                                filter_thresh = min(thresh_na)
                                                for thresh in thresh_na:
                                                    part_table_32.loc[
                                                        i, f'{tool}_rank_REVEL_{thresh}'] = f'filtered_by_REVEL_{filter_thresh}'
                                        if len(CADD_thresh) != 0:
                                            thresh_na = []
                                            for thresh in CADD_thresh:
                                                # part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}']
                                                if part_table_32.loc[i, f'{tool}_rank_CADD_{thresh}'] == 'NA':
                                                    thresh_na = thresh_na + [thresh]
                                            if len(thresh_na)!=0:
                                                filter_thresh = min(thresh_na)
                                                for thresh in thresh_na:
                                                    part_table_32.loc[
                                                        i, f'{tool}_rank_CADD_{thresh}'] = f'filtered_by_CADD_{filter_thresh}'

                        else:
                            for k in ['all_hpo_plus_weight','only_weight_hpo_weight','only_weight_hpo_no_weight']:
                                if part_table_32.loc[i, f'{k}_{tool}_rank'] == 'NA':
                                    for column in part_table_32.columns:
                                        if tool in str(column):
                                            part_table_32.loc[i, column] = 'not ranked'
                                else:
                                    if intersect:
                                        if part_table_32.loc[i, f'{k}_{tool}_intersect_rank'] == 'NA':
                                            for column in part_table_32.columns:
                                                if tool in str(column):
                                                    if f'{k}_{tool}_rank' not in column:
                                                        part_table_32.loc[i, column] = 'not_in_filtered_tsv'
                                        else:
                                            if len(REVEL_thresh) != 0:
                                                thresh_na = []
                                                for thresh in REVEL_thresh:
                                                    # part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}']
                                                    if part_table_32.loc[i, f'{k}_{tool}_rank_REVEL_{thresh}'] == 'NA':
                                                        thresh_na = thresh_na + [thresh]
                                                if len(thresh_na)!=0:
                                                    filter_thresh = min(thresh_na)
                                                    for thresh in thresh_na:
                                                        part_table_32.loc[
                                                            i, f'{k}_{tool}_rank_REVEL_{thresh}'] = f'filtered_by_REVEL_{filter_thresh}'
                                            if len(CADD_thresh) != 0:
                                                thresh_na = []
                                                for thresh in CADD_thresh:
                                                    # part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}']
                                                    if part_table_32.loc[i, f'{k}_{tool}_rank_CADD_{thresh}'] == 'NA':
                                                        thresh_na = thresh_na + [thresh]
                                                if len(thresh_na)!=0:
                                                    filter_thresh = min(thresh_na)
                                                    for thresh in thresh_na:
                                                        part_table_32.loc[
                                                            i, f'{k}_{tool}_rank_CADD_{thresh}'] = f'filtered_by_CADD_{filter_thresh}'
        part_table.to_csv(f'{hpo}_{len(tools)}_rank.tsv', sep='\t')
        part_rr_table.to_csv(f'{hpo}_{len(tools)}_rr.tsv', sep='\t')
        part_table_32.to_csv(f'{hpo}_{len(tools)}_rank_with_vcf.tsv', sep='\t')
        part_rr_table_32.to_csv(f'{hpo}_{len(tools)}_rr_with_vcf.tsv', sep='\t')

    if statistic:
        statistic_to_table(hpo)


# Press the green button in the gutter to run the script.

    # filename = 'scoliosis_gVCF_from_zs_updating'
    # pwd = '/Users/liyaqi/PycharmProjects/PhenoAptAnalysis'
    # df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    # # df = df[:1].reset_index(drop=True)
    main_2(tools=['phrank','phenolyzer','GADO','phen2gene','phenoapt','LIRICAL','Exomiser'],CADD_thresh=[10,15,20],REVEL_thresh=[0.25,0.5,0.75],hpo='Weight',intersect=True,refresh=False,newcase=False,statistic=True)

    # See PyCharm help at https://www.jetbrains.com/help/pycharm/

