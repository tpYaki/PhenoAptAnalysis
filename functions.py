import os

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
from pandas.testing import assert_frame_equal
import networkx
import obonet

def output_columns(hpo,tools,intersect,REVEL_thresh,CADD_thresh):
    columns_fill = ['CaseID', 'hpo', 'hpo_id_input', 'Symbol']
    for tool in tools:
        if tool in ['phrank', 'phenolyzer', 'GADO', 'phen2gene']:
            columns_fill = columns_fill + [f'{tool}_rank']
            if intersect: columns_fill = columns_fill + [f'{tool}_intersect_rank']
            if len(REVEL_thresh) != 0:
                for thresh in REVEL_thresh:
                    columns_fill = columns_fill + [f'{tool}_rank_REVEL_{thresh}']
            if len(CADD_thresh) != 0:
                for thresh in CADD_thresh:
                    columns_fill = columns_fill + [f'{tool}_rank_CADD_{thresh}']
        else:
            if tool == 'phenoapt':
                if hpo == 'Weight':
                    for k in ['all_hpo_plus_weight', 'only_weight_hpo_weight', 'only_weight_hpo_no_weight']:
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
    return columns_fill

def statistic_to_table(df,case_ids,pwd,filename,hpo,tools,intersect,REVEL_thresh,CADD_thresh):
    columns_fill = output_columns(hpo,tools,intersect,REVEL_thresh,CADD_thresh)
    part_table = pd.DataFrame(columns=columns_fill)
    part_rr_table = pd.DataFrame(columns=columns_fill+['Frequency','Variant_class','Inheritance_ADAR',f'{hpo}_organ_system',f'{hpo}_organ_system_number'])
    case_id_exomiser_tsv_file_dict = get_case_id_integ_file_map(pwd, filename, 'Exomiser', case_ids, hpo)
    case_id_LIRICAL_tsv_file_dict = get_case_id_integ_file_map(pwd, filename, 'LIRICAL', case_ids, hpo)
    part_table_32 = pd.DataFrame(columns=columns_fill)
    part_rr_table_32 = pd.DataFrame(columns=columns_fill+['Frequency','Variant_class','Inheritance_ADAR',f'{hpo}_organ_system',f'{hpo}_organ_system_number'])
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
        print(ensemblid)
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
        part_rr_table.loc[i, 'Frequency'] = df.loc[i,'Frequency']
        part_rr_table.loc[i, 'Variant_class'] = df.loc[i, 'Variant_class']
        part_rr_table.loc[i, 'Inheritance_ADAR'] = df.loc[i, 'Inheritance_ADAR']
        part_rr_table.loc[i, f'{hpo}_organ_system'] =df.loc[i,f'{hpo}_organ_system']
        part_rr_table.loc[i, f'{hpo}_organ_system_number'] = df.loc[i, f'{hpo}_organ_system_number']
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
                        i, f'phenolyzer_rank_REVEL_{thresh}'] = get_rank(symbol, phenolyzer_intersect_gene_name)
            if len(CADD_thresh) != 0:
                for thresh in CADD_thresh:
                    tsv_cadd_filter_dir = f'./{filename}/cadd_{thresh}/{case_id}.txt'
                    tsv_cadd_filter = pd.read_csv(tsv_cadd_filter_dir, sep='\t', header=None)
                    phenolyzer_intersect_gene_name = (
                        phenolyzer_df[phenolyzer_df.Gene.isin(tsv_cadd_filter[0])]).reset_index()
                    phenolyzer_intersect_gene_name = list(phenolyzer_intersect_gene_name['Gene'])
                    part_table.loc[i, f'phenolyzer_rank_CADD_{thresh}'], part_rr_table.loc[
                        i, f'phenolyzer_rank_CADD_{thresh}'] = get_rank(symbol, phenolyzer_intersect_gene_name)

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
            if hpo != 'Weight':
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
                for k in ['all_hpo_plus_weight', 'only_weight_hpo_weight', 'only_weight_hpo_no_weight']:
                    def get_phenoapt_result_dir(k):
                        if k == 'all_hpo_plus_weight':
                            pheno_result_dir = f'{pwd}/{filename}/phenoaptoutput/{case_id}_{hpo}_phenoapt_rank.tsv'
                            return pheno_result_dir
                        if k == 'only_weight_hpo_no_weight':
                            pheno_result_dir = f'{pwd}/{filename}/phenoaptoutput/{case_id}_only_{hpo}_hpo_no_weight_phenoapt_rank.tsv'
                            return pheno_result_dir
                        if k == 'only_weight_hpo_weight':
                            pheno_result_dir = f'{pwd}/{filename}/phenoaptoutput/{case_id}_only_{hpo}_hpo_weight_phenoapt_rank.tsv'
                            return pheno_result_dir

                    pheno_result = pd.read_csv(get_phenoapt_result_dir(k), sep='\t')
                    pheno_rank = list(pheno_result['gene_symbol'])
                    part_table.loc[i, f'{k}_phenoapt_rank'], part_rr_table.loc[i, f'{k}_phenoapt_rank'] = get_rank(
                        symbol, pheno_rank)
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

        # if 'Exomiser13.1.0' in tools:
        #     if case_id in case_id_exomiser1310_tsv_file_dict:
        #         rankdf = pd.read_csv(case_id_exomiser_tsv_file_dict[case_id], sep="\t")
        #         patho_gene_ID_Exo = list(rankdf['ENTREZ_GENE_ID'])
        #         if entrezId in patho_gene_ID_Exo:
        #             rank = patho_gene_ID_Exo.index(entrezId)
        #             rr = 1 / (1 + rank)
        #         else:
        #             rank = 'NA'
        #             rr = 0
        #     else:
        #         rank = 'No_vcf_ranked_by_Exo'
        #         rr = 0
        #     part_table.loc[i, 'Exomiser_rank'] = rank
        #     part_rr_table.loc[i, 'Exomiser_rank'] = rr

        # if 'HANRD' in tools:

        if 'LIRICAL' in tools:
            patho_gene_rank_LIR = {}
            if case_id in case_id_LIRICAL_tsv_file_dict:
                rankdf2 = pd.read_csv(case_id_LIRICAL_tsv_file_dict[case_id], sep="\t")
                patho_gene_ID_LIR = rankdf2['entrezGeneId']
                df_count = pd.DataFrame(enumerate(patho_gene_ID_LIR))
                for s in range(len(df_count)):
                    k = (df_count[df_count[1].isin([df_count.loc[s, 1]])]).reset_index(drop=True)
                    # print(k)
                    patho_gene_rank_LIR[k.loc[0, 1]] = k.loc[0, 0]
                # print(patho_gene_rank_LIR)

            else:
                patho_gene_rank_LIR[entrezGeneId] = 'No_vcf_ranked_by_LIR'
            rank = patho_gene_rank_LIR.get(entrezGeneId, 'NA')
            if rank == 'NA' or rank == 'No_vcf_ranked_by_LIR':
                rr = 0
            else:
                rr = 1 / (rank + 1)
            part_table.loc[i, 'LIRICAL_rank'] = rank
            part_rr_table.loc[i, 'LIRICAL_rank'] = rr

        if 'Phen_gen' in tools:
            part_table.loc[i, 'Phen_gen_rank'], part_rr_table.loc[i, 'Phen_gen_rank'] = phen_gen(case_id, pwd, hpo)

        ## 转化成作图数据

        if part_table.loc[i, 'LIRICAL_rank'] != 'No_vcf_ranked_by_LIR':
            part_table_32.loc[i] = part_table.loc[i]
            part_rr_table_32.loc[i] = part_rr_table.loc[i]
            ## 替换排名
            for column in part_table_32.columns:
                if 'rank' in column:
                    rank_value = part_table_32.loc[i, column]
                    if rank_value != 'NA':
                        if rank_value == 0:
                            part_table_32.loc[i, column] = 'TOP1'
                        if rank_value >= 1 and rank_value <= 4:
                            part_table_32.loc[i, column] = 'TOP5'
                        if rank_value > 4 and rank_value <= 9:
                            part_table_32.loc[i, column] = 'TOP10'
                        if rank_value > 9 and rank_value <= 19:
                            part_table_32.loc[i, column] = 'TOP20'
                        if rank_value > 19 and rank_value <= 49:
                            part_table_32.loc[i, column] = 'TOP50'
                        if rank_value > 49 and rank_value <= 99:
                            part_table_32.loc[i, column] = 'TOP100'
                        if rank_value > 99:
                            part_table_32.loc[i, column] = '>100'
            ## 替换NA
            for tool in tools:
                if tool in ['phrank', 'phenolyzer', 'GADO', 'phen2gene']:
                    if part_table_32.loc[i, f'{tool}_rank'] == 'NA':
                        for column in part_table_32.columns:
                            if tool in str(column):
                                part_table_32.loc[i, column] = 'not ranked'
                    else:
                        if intersect:
                            if part_table_32.loc[i, f'{tool}_intersect_rank'] == 'NA':
                                for column in part_table_32.columns:
                                    if tool in str(column):
                                        if column != f'{tool}_rank':
                                            part_table_32.loc[i, column] = 'not_in_filtered_tsv'
                            else:
                                if len(REVEL_thresh) != 0:
                                    thresh_na = []
                                    for thresh in REVEL_thresh:
                                        # part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}']
                                        if part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}'] == 'NA':
                                            thresh_na = thresh_na + [thresh]
                                    if len(thresh_na) != 0:
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
                                    if len(thresh_na) != 0:
                                        filter_thresh = min(thresh_na)
                                        for thresh in thresh_na:
                                            part_table_32.loc[
                                                i, f'{tool}_rank_CADD_{thresh}'] = f'filtered_by_CADD_{filter_thresh}'

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
                                            if column != f'{tool}_rank':
                                                part_table_32.loc[i, column] = 'not_in_filtered_tsv'
                                else:
                                    if len(REVEL_thresh) != 0:
                                        thresh_na = []
                                        for thresh in REVEL_thresh:
                                            # part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}']
                                            if part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}'] == 'NA':
                                                thresh_na = thresh_na + [thresh]
                                        if len(thresh_na) != 0:
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
                                        if len(thresh_na) != 0:
                                            filter_thresh = min(thresh_na)
                                            for thresh in thresh_na:
                                                part_table_32.loc[
                                                    i, f'{tool}_rank_CADD_{thresh}'] = f'filtered_by_CADD_{filter_thresh}'

                    else:
                        for k in ['all_hpo_plus_weight', 'only_weight_hpo_weight', 'only_weight_hpo_no_weight']:
                            if part_table_32.loc[i, f'{k}_{tool}_rank'] == 'NA':
                                for column in part_table_32.columns:
                                    if tool in str(column):
                                        part_table_32.loc[i, column] = 'not ranked'
                            else:
                                if intersect:
                                    if part_table_32.loc[i, f'{k}_{tool}_intersect_rank'] == 'NA':
                                        for column in part_table_32.columns:
                                            if tool in str(column):
                                                if column != f'{k}_{tool}_rank':
                                                    part_table_32.loc[i, column] = 'not_in_filtered_tsv'
                                    else:
                                        if len(REVEL_thresh) != 0:
                                            thresh_na = []
                                            for thresh in REVEL_thresh:
                                                # part_table_32.loc[i, f'{tool}_rank_REVEL_{thresh}']
                                                if part_table_32.loc[i, f'{k}_{tool}_rank_REVEL_{thresh}'] == 'NA':
                                                    thresh_na = thresh_na + [thresh]
                                            if len(thresh_na) != 0:
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
                                            if len(thresh_na) != 0:
                                                filter_thresh = min(thresh_na)
                                                for thresh in thresh_na:
                                                    part_table_32.loc[
                                                        i, f'{k}_{tool}_rank_CADD_{thresh}'] = f'filtered_by_CADD_{filter_thresh}'
                if tool in ['Exomiser', 'LIRICAL', 'Phen_gen']:
                    if part_table_32.loc[i, f'{tool}_rank'] == 'NA':
                        part_table_32.loc[i, f'{tool}_rank'] = f'not ranked'

    part_table.to_csv(f'{hpo}_{len(tools)}_rank.tsv', sep='\t')
    part_rr_table.to_csv(f'{hpo}_{len(tools)}_rr.tsv', sep='\t')
    part_table_32.to_csv(f'{hpo}_{len(tools)}_rank_with_vcf.tsv', sep='\t')
    part_rr_table_32.to_csv(f'{hpo}_{len(tools)}_rr_with_vcf.tsv', sep='\t')


def file_prepare(cohort_df, pwd, filename,REVEL_thresh=[],CADD_thresh=[],refresh=False):
    if len(REVEL_thresh) != 0:
        genertatevcf_zs_diag_dropbox_tsv(pwd=pwd, filename=filename, refresh=refresh)
        ##检查整个队列中必要的dropbox到vcf以及revel注释的vcf是否已有，没有就生成一个,储存在filename下的文件夹中；可能找不到dropboxtsv，pirint no tsv
    for i in tqdm(range(len(cohort_df))):
        try:
            case_id = cohort_df.loc[i, 'Blood ID']
            if get_case_id_file(case_id, 'sporadic') != 0:
                filedir = get_case_id_file(case_id, 'sporadic')
            else:
                if get_case_id_file(case_id, 'trio') != 0:
                    filedir = get_case_id_file(case_id, 'trio')
                else:
                    print(f'no dropbox tsv file for {case_id}')
                    continue
            variation = pd.read_csv(filedir, sep="\t")
            variation_gene_name = variation['Gene_name']
            variation_gene_list = list(set(variation_gene_name))

            # #检查是否需要产生新的phen2gene分析使用的candidate文件
            candidate_file(df2list=variation_gene_list, pwd=pwd, filename=filename, profile='candidategene', case_id=case_id)
            ##检查candidate gene 的 ensenmble file是否ready
            print('candidategene,ready')
            ensemblidfile(genelist=variation_gene_list, pwd=pwd, filename=filename, case_id=case_id,
                           profile='candidategene_ensemblid', refresh=refresh)
            print('candidategene_ensemblid,ready')

            if len(REVEL_thresh) != 0:
                for thresh in REVEL_thresh:
                    variation_revel = tsv_filtered_by_revel(df_input=variation, case_id=case_id,pwd=pwd, filename=filename,
                                                            thresh=thresh)
                    ##print(variation_revel[:1])
                    ## df_input, case_id, pwd, filename, thresh
                    variation_revel_gene_list = list(set(variation_revel['Gene_name']))
                    candidate_file(variation_revel_gene_list, pwd, filename, case_id, f'revel_{thresh}')
                    print(f'revel_{thresh} ready')
                    ensemblidfile(variation_revel_gene_list, pwd, case_id, filename,
                                   f'revel_{thresh}_ensemblid', refresh=refresh)
                    print(f'revel_{thresh}_ensemblid ready')
                    ##检查是否有按照阈值filter过的genename.txt

            if len(CADD_thresh) != 0:
                for thresh in CADD_thresh:
                    variation_cadd = patho_filter_n(variation, thresh)
                    variation_cadd_gene_list = list(set(variation_cadd['Gene_name']))
                    candidate_file(variation_cadd_gene_list, pwd, filename, case_id, f'cadd_{thresh}')
                    ##df2list, pwd, filename, case_id, profile, refresh
                    print(f'cadd_{thresh} ready')
                    ensemblidfile(variation_cadd_gene_list, pwd, case_id, filename, f'cadd_{thresh}_ensemblid',
                                   refresh=refresh)
                    ##检查是否有按照阈值filter过的genename.txt,refre指的是是否重新查询ensembleid，filterfile是每次都在重新生成
                    ## genelist, pwd, case_id, filename='scoliosis_gVCF_from_zs_updating', profile='candidategene_ensemblid',refresh=False
                    print(f'cadd_{thresh}_ensemblid ready')

        except Exception as e:
            print(e)

def get_rank(pathogene,rank_df):
    if pathogene in rank_df:
        rank = rank_df.index(pathogene)
        rr = 1 / (rank + 1)
    else:
        rank = 'NA'
        rr = 0
    return rank, rr

def load_revel_filter_df(filename, case_id, thresh):
    tsv_revel_filter_dir = f'./{filename}/revel_{thresh}/{case_id}.txt'
    tsv_revel_filter = pd.read_csv(tsv_revel_filter_dir, sep='\t')
    tsv_revel_filer_ensemblid_dir = f'./{filename}/revel_{thresh}_ensemblid/{case_id}.txt'
    tsv_revel_filer_ensemblid = pd.read_csv(tsv_revel_filer_ensemblid_dir, sep='\t')
    return tsv_revel_filter, tsv_revel_filer_ensemblid


def load_cadd_filter_df(filename, case_id, thresh):
    tsv_CADD_filter_dir = f'./{filename}/cadd_{thresh}/{case_id}.txt'
    tsv_CADD_filter = pd.read_csv(tsv_CADD_filter_dir, sep='\t')
    tsv_CADD_filer_ensemblid_dir = f'./{filename}/cadd_{thresh}_ensemblid/{case_id}.txt'
    tsv_CADD_filer_ensemblid = pd.read_csv(tsv_CADD_filer_ensemblid_dir, sep='\t')
    return tsv_CADD_filter, tsv_CADD_filer_ensemblid


def candidate_file(df2list, pwd, filename, case_id, profile):
    if profile not in os.listdir(f'{pwd}/{filename}/'):
        os.system(f'mkdir {pwd}/{filename}/{profile}')
    if 'candidategene' not in os.listdir(f'{pwd}/{filename}/'):
        os.system('mkdir candidategene')
    candidate_gene_list_no_dup = list(set(df2list))
    candidate_gene_list = pd.DataFrame(columns=[0])
    candidate_gene_list[0] = candidate_gene_list_no_dup
    candidate_gene_list.to_csv(f'{pwd}/{filename}/{profile}/{case_id}.txt', index=False, header=False)
    ##每次都存新的gene list


def ensemblidfile(genelist, pwd, case_id, filename='scoliosis_gVCF_from_zs_updating', profile='candidategene_ensemblid',refresh=False):
    if profile not in os.listdir(f'{pwd}/{filename}/'):
        os.system(f'mkdir {pwd}/{filename}/{profile}')
    if 'candidategene_ensemblid' not in os.listdir(f'{pwd}/{filename}/'):
        os.system(f'mkdir {pwd}/{filename}/candidategene_ensemblid')
    if f'{case_id}_index.tsv' not in os.listdir(f'./{filename}/candidategene_ensemblid/'):
        output = pd.DataFrame(columns=[0])
        for gene_name_keys in genelist:
            try:
                id_result = searchensemblid(gene_name_keys)
                output.loc[gene_name_keys] = id_result
            except Exception as e:
                print(e)
        output.to_csv(f'./{filename}/candidategene_ensemblid/{case_id}.txt', index=False, header=False)
        output.to_csv(f'./{filename}/candidategene_ensemblid/{case_id}_index.tsv', sep='\t', index=True, header=False)
    else:
        if refresh:
            output = pd.DataFrame(columns=[0])
            for gene_name_keys in genelist:
                try:
                    id_result = searchensemblid(gene_name_keys)
                    output.loc[gene_name_keys] = id_result
                except Exception as e:
                    print(e)
            output.to_csv(f'./{filename}/candidategene_ensemblid/{case_id}.txt', index=False, header=False)
            output.to_csv(f'./{filename}/candidategene_ensemblid/{case_id}_index.tsv', sep='\t', index=True,
                          header=False)
        if profile != 'candidategene_ensemblid':
            ensemblids_df = pd.read_csv(f'./{filename}/candidategene_ensemblid/{case_id}_index.tsv', sep='\t', header=None)
            ensemblids_df = ensemblids_df.reset_index()
            output_filter = ensemblids_df[ensemblids_df[0].isin(genelist)][1]
            ##print(output_filter[:5])
            output_filter.to_csv(f'./{filename}/{profile}/{case_id}.txt', index=False, header=False)
            ##每次都重新储存了filter后的新ensembleid

def getvcf_from_zs_gz(pwd,filename):
    ## 需要固定的.gz文件名称，需要批量提取的队列名单由xlsx整理找zs获得；单人vcf则直接下载
    os.system(f'{pwd}/{filename}/bash dividevcf.sh')
    print('vcf unziped from .gz')
    ## 从gz到单个vcf的脚本

def read_xlsx(path, sheet):
    wb = load_workbook(path)
    ws = wb[sheet]
    sh = pd.DataFrame(ws.values)
    df = sh.rename(columns=sh.loc[0])[1:].reset_index(drop=True)
    return df


def searchensemblid(genesymbol):
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/symbol/homo_sapiens/{genesymbol}?external_db=HGNC"
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    print(decoded[0]['id'])
    return decoded[0]['id']

import httplib2 as http
def searchHGNCid(genesymbol):
    try:
        from urlparse import urlparse
    except ImportError:
        from urllib.parse import urlparse

    headers = {
        'Accept': 'application/json',
    }

    uri = 'http://rest.genenames.org'
    path = f'/search/symbol:{genesymbol}'

    target = urlparse(uri + path)
    method = 'GET'
    body = ''

    h = http.Http()

    response, content = h.request(
        target.geturl(),
        method,
        body,
        headers)

    if response['status'] == '200':
    # assume that content is a json reply
    # parse content with the json module
        data = json.loads(content)
        #print('Symbol:' + data['response']['docs'][0]['hgnc_id'])
        HGNCid = data['response']['docs'][0]['hgnc_id']
        HGNCid = HGNCid.split(':')[1]
        return HGNCid
    else:
        #print('Error detected: ' + response['status'])
        return 'Error'


# 读取诊断表格
def read_diagnose_xlsx(path):
    wb = xlrd.open_workbook(path)
    sh = wb.sheet_by_name('Sheet1')
    print(sh.nrows)
    print(sh.ncols)
    print(sh.row_values(0))
    data = [sh.row_values(i) for i in range(1, sh.nrows)]
    df = pd.DataFrame(data=data, columns=sh.row_values(0))
    case_ids = df['CaseID']
    print(case_ids)
    return df, case_ids


def read_gvcf_diagnosis_xlsx(path):
    wb = xlrd.open_workbook(path)
    sh = wb.sheet_by_name('Sheet1')
    data = [sh.row_values(i) for i in range(1, sh.nrows)]
    df = pd.DataFrame(data=data, columns=sh.row_values(0))
    return df


# 获取case_id与相应的tsv文件路径的映射关系
def get_case_id_file_map(case_ids):
    case_id_tsv_file_map = {}
    tsv_dir = '/Users/liyaqi/Dropbox/Filtered-SNV-all/Sporadic-WES-All'
    files = os.listdir(tsv_dir)
    for case_id in case_ids:
        for file_name in files:
            if case_id in file_name:
                case_id_tsv_file_map[case_id] = os.path.join(tsv_dir, file_name)
    return case_id_tsv_file_map


def get_case_id_file(case_id, mode):
    if mode == 'sporadic':
        tsv_dir = '/Users/liyaqi/Dropbox/Filtered-SNV-all/Sporadic-WES-All'
        files = os.listdir(tsv_dir)
        for file_name in files:
            if case_id in file_name:
                case_id_tsv_file = os.path.join(tsv_dir, file_name)
                return case_id_tsv_file
        return 0
    else:
        if mode == 'trio':
            tsv_dir = '/Users/liyaqi/Dropbox/Filtered-SNV-all/Trio-WES-all'
            files = os.listdir(tsv_dir)
            for file_name in files:
                if case_id in file_name:
                    case_id_tsv_file = os.path.join(tsv_dir, file_name)
                    return case_id_tsv_file
            return 0


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


def containtype(biao, type):
    data_mut = biao[biao["Mutation_type"].str.contains(type)]
    return data_mut
    # 包含字符串type=“”的rows，Na的单元格自动去掉


def selectsmall(biao, name, thresh):
    data_mut1 = biao[biao[name] != "."]
    data_mut3 = data_mut1[data_mut1[name] != "nan"]  ## df,"Exac"

    def is_float(x):
        try:
            float(x)
        except ValueError:
            return False
        return True

    data_mut12 = data_mut3[data_mut3[name].apply(lambda x: is_float(x))]
    data_mut12[name] = pd.to_numeric(data_mut12[name], errors='raise')
    data_mut2 = (data_mut12[(data_mut12[name]) < thresh].append(biao[biao[name] == "."])).append(
        biao[biao[name] == "nan"])
    return data_mut2


def selectbig(biao, name, thresh):
    data_mut1 = biao[biao[name] != "."]
    data_mut3 = data_mut1[data_mut1[name] != "nan"]  ## df,"Exac"

    def is_float(x):
        try:
            float(x)
        except ValueError:
            return False
        return True

    data_mut12 = data_mut3[data_mut3[name].apply(lambda x: is_float(x))]
    data_mut12[name] = pd.to_numeric(data_mut12[name], errors='raise')
    data_mut2 = (data_mut12[data_mut12[name] > thresh].append(biao[biao[name] == "."])).append(
        biao[biao[name] == "nan"])
    return data_mut2


def patho_filter_n(df, threshold):
    ##data2 = ((containtype(data1, "missense_variant")).append(containtype(data1, "stop_gained"))).append(containtype(data1, "frameshift"))
    data3 = selectsmall(
        selectsmall(selectsmall(selectsmall(df, "ExAC_AF", 0.01), "ExAC_EAS_AF", 0.01), "gnomAD_genome_EAS", 0.01),
        "In_house_freq", 0.05)
    data4 = selectbig(data3, "VAF", 29)
    data5 = selectbig(data4, "CADD", threshold)
    ##data5 = (congrade(data4,"Mutationtatser_prediction","A")).append(data4[data4["Mutationtatser_prediction"].str.contains("D")])
    ##data7 = congrade(data6,"LRT_prediction","D").append(data6[data6["LRT_prediction"].str.contains("U")])
    ##data8 = congrade(data7, "Polyphen2_HDIV_pred", "D").append(data7[data7["Polyphen2_HDIV_pred"].str.contains("P")])
    ##data9 = congrade(data8, "Polyphen2_HVAR_pred", "D").append(data8[data8["Polyphen2_HVAR_pred"].str.contains("P")])
    ##data10 = congrade(data9, "SIFT_prediction", "D")
    return data5


def tsv_filtered_by_revel(df_input, case_id, pwd, filename, thresh):
    ## df_imput 是读出来的dropbox tsv,根据revel注释的vcf文件，生成{case_id}_revel_filter_{REVEL_Thresh}.tsv，也就是可以读取成df的文件，并且输出df_input中符合revel要求的部分
    dir = f'{pwd}/{filename}/dropbox2vcf_revel'
    if f'{case_id}_revel_filter_{thresh}.tsv' not in os.listdir(dir):
        cmd = f'source /Users/liyaqi/opt/anaconda3/etc/profile.d/conda.sh && conda activate vep && filter_vep -i {dir}/{case_id}.vcf --filter "REVEL > {thresh} or not REVEL" -o {dir}/{case_id}_revel_filter_{thresh}.tsv --force_overwrite && bash {dir}/delet2df.sh'
        os.system(cmd)
    vcf_revel_tsv = pd.read_csv(f'{dir}/{case_id}_revel_filter_{thresh}.tsv', sep='\t')
    vcf_revel_tsv['variant_list'] = vcf_revel_tsv['Location'].map(str) + vcf_revel_tsv['Allele']
    for i in range(len(df_input)):
        if df_input.loc[i, 'ALT'] == '-':
            if len(df_input.loc[i, 'REF']) == 1:
                pos = str(df_input.loc[i, 'POS'])
            else:
                pos = str(df_input.loc[i, 'POS']) + '-' + str(int(df_input.loc[i, 'POS']) + len(df_input.loc[i, 'REF']) - 1)
        else:
            if df_input.loc[i, 'REF'] == '-':
                pos = str(int(df_input.loc[i, 'POS']) - 1) + '-' + str(df_input.loc[i, 'POS'])
            else:
                pos = str(df_input.loc[i, 'POS'])
        df_input.loc[i, 'intersect_revel_list'] = str(df_input.loc[i, 'CHR']) + ':' + pos + str(df_input.loc[i, 'ALT'])
    df_result = df_input[df_input['intersect_revel_list'].isin(vcf_revel_tsv['variant_list'])]
    return df_result


def getvcf_from_dropbox_tsv(pwd,case_id, dir, filename):
    df = pd.read_csv(dir, sep="\t")  ##dropbox的tsv
    dfx = df
    for i in range(len(df)):
        ##补充indel的REF和ALT
        if df.at[i, 'ALT'] == '-':  ##deletion
            ##print(i, 'ALT', df.at[i, 'ALT'])
            pos = int(df.at[i, 'POS'])
            chr = df.at[i, 'CHR']
            cmd = f'source /Users/liyaqi/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && samtools faidx /Users/liyaqi/Software/ensembl-vep/Homo_sapiens.GRCh37.dna.primary_assembly.fa {chr}:{pos - 1}-{pos - 1}'
            result = subprocess.getoutput(cmd)
            dfx.at[i, 'ALT'] = result.split('\n')[1]
            dfx.at[i, 'REF'] = result.split('\n')[1] + df.at[i, 'REF']
            dfx.at[i, 'POS'] = pos - 1
        else:
            if df.at[i, 'REF'] == '-':  ## insertion
                ##print(i, 'REF', df.at[i, 'REF'])
                pos = int(df.at[i, 'POS'])
                chr = df.at[i, 'CHR']
                cmd = f'source /Users/liyaqi/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && samtools faidx /Users/liyaqi/Software/ensembl-vep/Homo_sapiens.GRCh37.dna.primary_assembly.fa {chr}:{pos - 1}-{pos - 1}'
                result = subprocess.getoutput(cmd)
                dfx.at[i, 'REF'] = result.split('\n')[1]
                dfx.at[i, 'ALT'] = result.split('\n')[1] + df.at[i, 'ALT']
                dfx.at[i, 'POS'] = pos - 1

    for j in range(len(dfx)):
        if str(dfx.at[j, 'FILTER']) != '-':
            continue
            ## 只把-换成.
        dfx.at[j, 'FILTER'] = '.'

    QUAL = pd.DataFrame({'QUAL': pd.Series([100 for k in range(len(dfx))])})
    INFO = pd.DataFrame({'INFO': pd.Series(['.' for k in range(len(dfx))])})
    ## 假设一些无关紧要的格式信息，完成vcf必须的8列
    # df2 = pd.concat([dfx[['CHR', 'POS', 'Rs_ID', 'REF', 'ALT']], QUAL, dfx['FILTER'], INFO, dfx['FORMAT']], axis=1)
    # df2.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    df2 = pd.concat([dfx[['CHR', 'POS', 'Rs_ID', 'REF', 'ALT']], QUAL, dfx['FILTER'], INFO], axis=1)
    df2.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    df2.to_csv(f"./{filename}/{case_id}.txt", sep="\t", index=False)  ##转化成vcf信息txt,存到对应数据版本file下
    vcfdir = f"{pwd}/{filename}"
    command = f"source /Users/liyaqi/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && awk -f {pwd}/script.awk {vcfdir}/{case_id}.txt | bcftools view -o {vcfdir}/{case_id}.vcf"
    ## awk加一些vcf的表头信息，再bcftools转化完成vcf
    print(command)
    os.system(command)

def genertatevcf_zs_diag_dropbox_tsv(pwd, filename, refresh=False):
    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    ##df = df[2:].reset_index(drop=True)
    print(f"{len(df)}")
    for i in tqdm(range(len(df))):
        case_id = df.loc[i, 'Blood ID']
        if get_case_id_file(case_id, 'sporadic') != 0:
            filedir = get_case_id_file(case_id, 'sporadic')
        else:
            if get_case_id_file(case_id, 'trio') != 0:
                filedir = get_case_id_file(case_id, 'trio')
            else:
                print(f'no dropbox tsv file for {case_id}')
                continue
        if f'{case_id}.vcf' not in os.listdir(f'{pwd}/{filename}/dropbox2vcf'):
            getvcf_from_dropbox_tsv(pwd,case_id, filedir, f'{filename}/dropbox2vcf')
            ##pwd,case_id, dir, filename
        else:
            if refresh:
                getvcf_from_dropbox_tsv(pwd,case_id, filedir, f'{filename}/dropbox2vcf')
        print('vcf done')
        if f'{case_id}.vcf' not in os.listdir(f'{pwd}/{filename}/dropbox2vcf_revel'):
            vcf_revel(pwd,f'{filename}/dropbox2vcf',
                      f'{filename}/dropbox2vcf_revel', case_id=case_id)
        else:
            if refresh:
                vcf_revel(pwd,f'{filename}/dropbox2vcf',
                          f'{filename}/dropbox2vcf_revel', case_id=case_id)
        print('revel done')


def vcf_revel(pwd,inputdir, outputdir, case_id):
    ## 使用前，先设置文件夹读写权限"chmod a+rwx ./output/*" "chmod a+rwx ./input/*"
    cmd_pre = f'chmod a+rwx {pwd}/{inputdir}/* && chmod a+rwx {pwd}/{outputdir}/*'
    os.system(cmd_pre)
    print("right done")
    cmd = f'source /Users/liyaqi/opt/anaconda3/etc/profile.d/conda.sh && conda activate vep && Vep -i {pwd}/{inputdir}/{case_id}.vcf --fork 4 -o {pwd}/{outputdir}/{case_id}.vcf --assembly GRCh37 --cache --dir /Users/liyaqi/Software/ensembl-vep --offline --fasta /Users/liyaqi/Software/ensembl-vep/Homo_sapiens.GRCh37.dna.primary_assembly.fa --plugin REVEL,/Users/liyaqi/Software/ensembl-vep/revel/new_tabbed_revel.tsv.gz --force_overwrite'
    os.system(cmd)



import yaml
def getyml_for_certain_ingeritance(case_id_file, hpo_id_input, dir, outputdir, inheridict, hpo):
    with open(f'./yml/test-analysis-exome.yml', 'r') as f:
        config = yaml.safe_load(f)
        config['analysis']['vcf'] = f'{dir}/vcf/{case_id_file}.vcf'
        config['analysis']['hpoIds'] = hpo_id_input  # add the command as a list for the correct yaml
        config['outputOptions']['outputFormats'] = ['TSV_GENE']
        # config['outputOptions']['outputFormats'] = ['HTML','TSV_GENE','TSV_VARIANT']
        config['outputOptions']['outputPrefix'] = f'{dir}/{outputdir}/{case_id_file}'
        config['analysis']['inheritanceModes'] = inheridict
        print(config['analysis']['hpoIds'])
        print(config['analysis']['inheritanceModes'])
    with open(f'{dir}/yml/{case_id_file}_{hpo}.yml', "a") as f:  # open the file in append mode
        f.truncate(0)
        yaml.dump(config, f)  ##default_flow_style=Faulse

def getyml_for_1310(sequence_id,case_id, hpo_id_input, dir, outputdir, inheridict, hpo, pathogenic = ['REVEL', 'MVP']):
    with open(f'./yml/1310-test-analysis-exome.yml', 'r') as f:
        config = yaml.safe_load(f)
        config['analysis']['vcf'] = f'{dir}/vcf/{sequence_id}.vcf'
        config['analysis']['hpoIds'] = hpo_id_input  # add the command as a list for the correct yaml
        config['analysis']['pathogenicitySources'] = pathogenic
        config['outputOptions']['outputFormats'] = ['TSV_GENE']
        # config['outputOptions']['outputFormats'] = ['HTML','TSV_GENE','TSV_VARIANT']
        config['outputOptions']['outputPrefix'] = f'{dir}/{outputdir}/1310-{case_id}'
        config['inheritanceModes'] = inheridict
        print(config['analysis']['hpoIds'])
        print(config['analysis']['inheritanceModes'])
    with open(f'{dir}/yml/1310-{case_id}_{hpo}.yml', "a") as f:  # open the file in append mode
        f.truncate(0)
        yaml.dump(config, f)  ##default_flow_style=Faulse

def generateyml_zs_diag_zs_gVCF(df,pwd,filename,hpo,refresh=False):
    #df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    ##df = df[36:37].reset_index(drop=True)
    sequence_ids = df['sequence ID']
    print(f"{len(df)}")
    dir_scoliosis_gVCF_from_zs = f"{pwd}/{filename}"
    if 'vcf' not in os.listdir(dir_scoliosis_gVCF_from_zs):
        os.system(f'mkdir {dir_scoliosis_gVCF_from_zs}/vcf')
    sequence_id_file_dict = os.listdir(f'{dir_scoliosis_gVCF_from_zs}/vcf/')
    ##print(sequence_id_file_dict)
    m = 0
    n = 0
    for i in range(len(df)):
        sequence_id = sequence_ids[i]
        sequence_id_file = f'{sequence_id}.vcf'
        case_id = df.loc[i,'Blood ID']
        #Inheritance_ADAR = df.loc[i,'Inheritance_ADAR']
        if sequence_id_file not in sequence_id_file_dict:
            print(f"{sequence_id_file} not in {filename}")
            n = n + 1
            ##去搞vcf吧
            continue
        else:
            print("ok")
            if df.loc[i, hpo][0] == '[':
                hpo_id = df.loc[i, hpo][2:-2]
                hpo_id_input = [str(k) for k in hpo_id.split("', '")]
            else:
                hpo_id = df.loc[i, hpo]
                hpo_id_input = [k for k in hpo_id.split(";")]
            inheri = (str(df.loc[i, 'inheritance'])).split(',')  ##存在输入了多种mode的可能
            inheridict = {}
            for k in inheri:
                j = k.split(':')
                inheridict[f'{j[0]}'] = float(j[1])
            ##inheridict = [k for k in (df.loc[i, 'inheritance']).split(";")]
            print(inheridict)
            if hpo!='hpo_id':
                outputdir = f'{hpo}_Exomiseroutput'
                if outputdir not in os.listdir(dir_scoliosis_gVCF_from_zs):
                    os.system(f'mkdir {dir_scoliosis_gVCF_from_zs}/{outputdir}')
                if (f'{sequence_id}_{hpo}.yml' not in os.listdir(f'{dir_scoliosis_gVCF_from_zs}/yml/')) or refresh:
                    getyml_for_certain_ingeritance(sequence_id, hpo_id_input, dir_scoliosis_gVCF_from_zs, outputdir, inheridict, hpo)
                if f'{hpo}_Exomiseroutput' not in os.listdir(dir_scoliosis_gVCF_from_zs):
                    os.system(f'mkdir {dir_scoliosis_gVCF_from_zs}/{hpo}_Exomiseroutput')
                if (case_id not in " ".join(os.listdir(f'{dir_scoliosis_gVCF_from_zs}/{hpo}_Exomiseroutput'))) or refresh:
                    command = f'cd /Users/liyaqi/Downloads/exomiser-cli-12.1.0 && java -Xms2g -Xmx4g -jar exomiser-cli-12.1.0.jar -analysis {dir_scoliosis_gVCF_from_zs}/yml/{sequence_id}_{hpo}.yml'
                    print(command)
                    os.system(command)
            else:
                if (f'{sequence_id}_{hpo}.yml' not in os.listdir(f'{dir_scoliosis_gVCF_from_zs}/yml/')) or refresh:
                    getyml_for_certain_ingeritance(sequence_id, hpo_id_input, dir_scoliosis_gVCF_from_zs, 'Exomiseroutput',
                                                   inheridict, hpo)
                if (case_id not in " ".join(os.listdir(f'{dir_scoliosis_gVCF_from_zs}/Exomiseroutput'))) or refresh:
                    command = f'cd /Users/liyaqi/Downloads/exomiser-cli-12.1.0 && java -Xms2g -Xmx4g -jar exomiser-cli-12.1.0.jar -analysis {dir_scoliosis_gVCF_from_zs}/yml/{sequence_id}_{hpo}.yml'
                    print(command)
                    os.system(command)
            print(f'{case_id}, NO.{i}, exomiser done')
            m = m + 1
            print(f'{m} case finish exomiser')
            print(f'{n} files were not in {filename}')

## pahtogenicity version合并
def generateyml1310_zs_diag_zs_gVCF(df, pwd, filename, hpo, refresh=False):
    sequence_ids = df['sequence ID']
    print(f"{len(df)}")
    dir_scoliosis_gVCF_from_zs = f"{pwd}/{filename}"
    if 'vcf' not in os.listdir(dir_scoliosis_gVCF_from_zs):
        os.system(f'mkdir {dir_scoliosis_gVCF_from_zs}/vcf')
    sequence_id_file_dict = os.listdir(f'{dir_scoliosis_gVCF_from_zs}/vcf/')
    ##print(sequence_id_file_dict)
    m = 0
    n = 0
    for i in range(len(df)):
        sequence_id = sequence_ids[i]
        sequence_id_file = f'{sequence_id}.vcf'
        case_id = df.loc[i,'Blood ID']
        #Inheritance_ADAR = df.loc[i,'Inheritance_ADAR']
        if sequence_id_file not in sequence_id_file_dict:
            print(f"{sequence_id_file} not in {filename}")
            n = n + 1
            ##去搞vcf吧
            continue
        else:
            print("vcf ok")
            if df.loc[i, hpo][0] == '[':
                hpo_id = df.loc[i, hpo][2:-2]
                hpo_id_input = [str(k) for k in hpo_id.split("', '")]
            else:
                hpo_id = df.loc[i, hpo]
                hpo_id_input = [k for k in hpo_id.split(";")]
            inheri = (str(df.loc[i, 'inheritance'])).split(',')  ##存在输入了多种mode的可能
            inheridict = {}
            for k in inheri:
                j = k.split(':')
                inheridict[f'{j[0]}'] = float(j[1])
            ##inheridict = [k for k in (df.loc[i, 'inheritance']).split(";")]
            print(inheridict)
            if hpo!='hpo_id':
                outputdir = f'{hpo}_Exomiseroutput'
                if outputdir not in os.listdir(dir_scoliosis_gVCF_from_zs):
                    os.system(f'mkdir {dir_scoliosis_gVCF_from_zs}/{outputdir}')
                if (f'1310-{case_id}_{hpo}.yml' not in os.listdir(f'{dir_scoliosis_gVCF_from_zs}/yml/')) or refresh:
                    getyml_for_1310(sequence_id,case_id, hpo_id_input, dir_scoliosis_gVCF_from_zs, outputdir, inheridict, hpo)
                if (f'1310-{case_id}' not in " ".join(os.listdir(f'{dir_scoliosis_gVCF_from_zs}/{hpo}_Exomiseroutput'))) or refresh:
                    command = f'cd /Users/liyaqi/Downloads/exomiser-cli-13.1.0 && java -Xms2g -Xmx4g -jar exomiser-cli-13.1.0.jar -analysis {dir_scoliosis_gVCF_from_zs}/yml/1310-{case_id}_{hpo}.yml'
                    print(command)
                    os.system(command)
            else:
                if (f'1310-{case_id}_{hpo}.yml' not in os.listdir(f'{dir_scoliosis_gVCF_from_zs}/yml/')) or refresh:
                    getyml_for_1310(sequence_id,case_id, hpo_id_input, dir_scoliosis_gVCF_from_zs, 'Exomiseroutput',
                                                   inheridict, hpo)
                if (f'1310-{case_id}' not in " ".join(os.listdir(f'{dir_scoliosis_gVCF_from_zs}/Exomiseroutput'))) or refresh:
                    command = f'cd /Users/liyaqi/Downloads/exomiser-cli-13.1.0 && java -Xms2g -Xmx4g -jar exomiser-cli-13.1.0.jar -analysis {dir_scoliosis_gVCF_from_zs}/yml/1310-{case_id}_{hpo}.yml'
                    print(command)
                    os.system(command)
            print(f'{case_id}, NO.{i}, exomiser 13.1.0 done')
            m = m + 1
            print(f'{m} case finish exomiser 13.1.0')
            print(f'{n} files were not in {filename}')


def get_phenopaket(sequence_id, hpo_id_input, dir, hpo):
    with open('/Users/liyaqi/Downloads/LIRICAL-1.3.4/example.json', 'rb') as f:  # 使用只读模型，并定义名称为f
        params = json.load(f)  # 加载json文件
        pheno_dict = []
        for i in range(len(hpo_id_input)):
            pheno_dict_id = dict()
            pheno_dict_id_2 = dict()
            pheno_dict_id['id'] = hpo_id_input[i]
            pheno_dict_id_2['type'] = pheno_dict_id
            pheno_dict.append(dict(pheno_dict_id_2))
        params["phenotypicFeatures"] = pheno_dict
        params["htsFiles"][0]["uri"] = f'{dir}vcf/{sequence_id}.vcf'
        ##print("params", params)  # 打印
        f.close()  # 关闭json读模式
    with open(f'{dir}phenopacket/{sequence_id}_{hpo}.json', 'w') as r:
        json.dump(params, r)
        # 关闭json写模式
        r.close()


def generate_phenopaket_zs_diag_zs_gVCF(df,pwd,filename,hpo,refresh=False):
    #df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    ##df = df[36:37].reset_index(drop=True)
    sequence_ids = df['sequence ID']
    print(f"{len(df)}")
    dir_scoliosis_gVCF_from_zs = f"{pwd}/{filename}/"
    if 'vcf' not in os.listdir(f'{dir_scoliosis_gVCF_from_zs}'):
        os.system(f'mkdir {dir_scoliosis_gVCF_from_zs}vcf')
    sequence_id_file_dict = os.listdir(f'{dir_scoliosis_gVCF_from_zs}vcf/')
    ##print(sequence_id_file_dict)
    m = 0
    n = 0
    for i in range(len(df)):
        sequence_id = sequence_ids[i]
        sequence_id_file = f'{sequence_id}.vcf'
        case_id = sequence_id.split('-')[0]
        if sequence_id_file not in sequence_id_file_dict:
            print(f"{sequence_id_file} not in {filename}")
            n = n + 1
            ##需要搞到vcf
            continue
        else:
            print("ok")
            if df.loc[i, hpo][0] == '[':
                hpo_id = df.loc[i, hpo][2:-2]
                hpo_id_input = [str(k) for k in hpo_id.split("', '")]
            else:
                hpo_id = df.loc[i, hpo]
                hpo_id_input = [k for k in hpo_id.split(";")]
            if (f'{sequence_id}_{hpo}.json' not in os.listdir()) or refresh:
                get_phenopaket(sequence_id, hpo_id_input, dir_scoliosis_gVCF_from_zs, hpo)
                print(f'{case_id}, NO.{i}, phenopaket done')
            m = m + 1
            if hpo =='hpo_id':
                if (f'{sequence_id}.tsv' not in os.listdir(f'{pwd}/{filename}/LIRICALoutput/')) or refresh:
                    command = f'cd {pwd}/{filename}/LIRICALoutput && java -jar LIRICAL.jar phenopacket -p {dir_scoliosis_gVCF_from_zs}phenopacket/{sequence_id}_{hpo}.json -d /Users/liyaqi/Downloads/LIRICAL-1.3.4/data -e /Users/liyaqi/Downloads/exomiser-cli-12.1.0/data/1909_hg19 -x {sequence_id} --tsv'
                    print(command)
                    os.system(command)
                    ##跑太慢了，所以hpo_id的结果还是用原来的

            else:
                if f'{hpo}_LIRICALoutput' not in os.listdir(f'{pwd}/{filename}'):
                    os.system(f'mkdir {pwd}/{filename}/{hpo}_LIRICALoutput && cp {pwd}/{filename}/LIRICALoutput/LIRICAL.jar {pwd}/{filename}/{hpo}_LIRICALoutput/LIRICAL.jar && cp {pwd}/{filename}/LIRICALoutput/outputtodf.sh {pwd}/{filename}/{hpo}_LIRICALoutput/outputtodf.sh')
                if (f'{sequence_id}.tsv'not in os.listdir(f'{pwd}/{filename}/{hpo}_LIRICALoutput/')) or refresh:
                    command = f'cd {pwd}/{filename}/{hpo}_LIRICALoutput && java -jar LIRICAL.jar phenopacket -p {dir_scoliosis_gVCF_from_zs}phenopacket/{sequence_id}_{hpo}.json -d /Users/liyaqi/Downloads/LIRICAL-1.3.4/data -e /Users/liyaqi/Downloads/exomiser-cli-12.1.0/data/1909_hg19 -x {sequence_id} --tsv'
                    print(command)
                    os.system(command)

        print(f'{m} case finish LIRICAL')
    if hpo == 'hpo_id':
        os.system(f'bash {pwd}/{filename}/LIRICALoutput/outputtodf.sh')
    else:
        os.system(f'bash {pwd}/{filename}/{hpo}_LIRICALoutput/outputtodf.sh')


def generateped():
    df = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx', 'Sheet1')
    # df = df[:1].reset_index(drop=True)
    df_sex = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx', 'gVCF-zs-确诊-WES')
    df_sex = df_sex.rename(index=df_sex['Blood ID'])
    df_sex = dict(df_sex['Sex'])
    for i in range(len(df)):
        case_id = df.loc[i, 'Blood ID']
        ##做ped
        ped_df = (pd.DataFrame(['FAM1', case_id, 'FATHER', 'MOTHER', df_sex.get(case_id), 2]).T)
        ped_df.to_csv(
            f'/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/ped/{case_id}.tsv',
            sep='\t', index=False, header=False)
        cmd = f'mv /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/ped/{case_id}.tsv /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/ped/{case_id}.ped'
        os.system(cmd)
        ##做好了ped for phen-gen

def get_HANRD(df,pwd,filename,hpo):
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
        HANRD_input.loc[i,'HGNC (Casual Gene)'] = df.loc[i, 'HGNC']
        HANRD_input.loc[i,'Orphanet ID (Diagnosed Disorder)'] = '-'
        # HANRD_input.set_index('ID')
    if 'HANRDoutput' not in os.listdir(f'{pwd}/{filename}'):
        os.system('mkdir HANRDoutput')
    if f'{hpo}_gVCF-diagnosed-cohort-2022-HANRD.tsv' in os.listdir(f'{pwd}/{filename}/HANRDoutput/'):
        ##已经有输入file了
        old_input_file = pd.read_csv(f'{pwd}/{filename}/HANRDoutput/{hpo}_gVCF-diagnosed-cohort-2022-HANRD.tsv', sep='\t')
        if old_input_file.equals(HANRD_input):
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
        os.system(f'cd ~/Software/HANRD/gcas && java -Xmx8096M -jar gcas.jar {pwd}/{filename}/HANRDoutput/{hpo}_gVCF-diagnosed-cohort-2022-HANRD.tsv && mv ./data/output/gVCF-diagnosed-cohort-2022-HANRD_genes.out {pwd}/{filename}/HANRDoutput/{hpo}_gVCF-diagnosed-cohort-2022-HANRD_genes.csv')
        print('HANRD done')

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
        if df.loc[i, hpo][0] == '[':
            hpo_id = df.loc[i, hpo][2:-2]
            hpo_id_input = [str(k) for k in hpo_id.split("', '")]
        else:
            hpo_id = df.loc[i, hpo]
            hpo_id_input = [k for k in hpo_id.split(";")]
        #print(hpo_id_input)
        HPO=''
        for k in hpo_id_input:
            HPO = HPO+str(k)+';'
        case_id = df.loc[i,'Blood ID']
        if 'PhenoRankoutput' not in os.listdir(f'{pwd}/{filename}'):
            os.system(f'mkdir {pwd}/{filename}/PhenoRankoutput')
        if f'{case_id}_{hpo}.tsv' not in os.listdir(f'{pwd}/{filename}/PhenoRankoutput') or refresh:
            cmd = f"cd ~/Software/PhenoRank && source ~/opt/anaconda3/etc/profile.d/conda.sh && conda activate phenorank && python run_PhenoRank.py -o {pwd}/{filename}/PhenoRankoutput/{case_id}_{hpo}.tsv -p '{HPO[:-1]}' && exit"
            appscript.app('Terminal').do_script(cmd)
            time.sleep(15)
            print('new phenorank ok')
        else:
            print('phenorank exist')



def GADO(df,pwd='/Users/liyaqi/PycharmProjects/PhenoAptAnalysis',filename='scoliosis_gVCF_from_zs_updating',hpo='hpo_id',refresh=False):
    # df = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx', 'Sheet1')
    # df = df[:1].reset_index(drop=True)
    # if f'totalhpo_{hpo}.txt' in os.listdir(f'{pwd}/{filename}/hpotxt/'):
    #     previous_input_hpo = pd.read_csv(f'./{filename}/hpotxt/totalhpo_{hpo}.txt',sep='\t')
    ##把hpotxt存在txt里
    # else:
    #     previous_input_hpo =pd.DataFrame()
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
    cmd1 = f'cd /Users/liyaqi/Software/GadoCommandline-1.0.1 && java -jar GADO.jar \
     --mode PROCESS \
     --output {pwd}/{filename}/GADOoutput/hpoProcessed.txt \
     --caseHpo {pwd}/{filename}/hpotxt/totalhpo_{hpo}.txt \
     --hpoOntology data/hp.obo \
     --hpoPredictionsInfo data/predictions_auc_bonf.txt'
    cmd2 = f'cd /Users/liyaqi/Software/GadoCommandline-1.0.1 && java -jar GADO.jar \
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
                command_phen2gene = f'cd /Users/liyaqi/Phen2Gene && python3 phen2gene.py -m{HPO_for_Phen2Gene} -out {pwd}/{filename}/phen2geneoutput_nolist && mv {pwd}/{filename}/phen2geneoutput_nolist/output_file.associated_gene_list {pwd}/{filename}/phen2geneoutput_nolist/{case_id}-phen2gene-nolist-output.txt'
                os.system(command_phen2gene)
                print('run phen2gene')
        else:
            if refresh or (f'{case_id}-phen2gene-{profile}-output.txt' not in os.listdir(f'{pwd}/{filename}/phen2geneoutput/')):
                command_phen2gene_intersect = f'cd /Users/liyaqi/Phen2Gene && python3 phen2gene.py -m{HPO_for_Phen2Gene} -out {pwd}/{filename}/phen2geneoutput -l {pwd}/{filename}/{profile}/{case_id}.txt && mv {pwd}/{filename}/phen2geneoutput/output_file.associated_gene_list {pwd}/{filename}/phen2geneoutput/{case_id}-phen2gene-{profile}-output.txt'
                os.system(command_phen2gene_intersect)
                print('run phen2gene')
        ## 有基因列表


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
                cmd_phenolyzer = f'cd /Users/liyaqi/Software/phenolyzer && perl disease_annotation.pl {hpo_txt_dir} -file -prediction -phenotype -logistic -out {pwd}/{filename}/phenolyzeroutput/{case_id}_{hpo} -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 && exit'
                print(cmd_phenolyzer)
                appscript.app('Terminal').do_script(cmd_phenolyzer)
                time.sleep(30)
        else:
            hpo_txt_dir = f'{pwd}/{filename}/hpotxt/{case_id}_{hpo}.txt'
            ##如果有新case，还是用新产生的——hpo_id后缀的hpo.txt
            if (f'{case_id}.' not in ''.join(os.listdir(f'{pwd}/{filename}/phenolyzeroutput/'))) or refresh:
                cmd_phenolyzer = f'cd /Users/liyaqi/Software/phenolyzer && perl disease_annotation.pl {hpo_txt_dir} -file -prediction -phenotype -logistic -out {pwd}/{filename}/phenolyzeroutput/{case_id} -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 && exit'
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
                intrisic_weight_df= pd.read_csv('/Users/liyaqi/Documents/生信/phenoapt/weight.csv')
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
    with open(f'./yml/test-analysis-exome.yml', 'r') as f:
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
        if tool in ['phenoapt','phrank', 'phenolyzer', 'GADO', 'phen2gene']:
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
    ##针对hpo_id做benchmaark的代码
    for tool in tools:
        if tool in ['phenoapt', 'phrank', 'phenolyzer', 'GADO', 'phen2gene']:
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

def Rscript_MRR_matrix(hpo,tools,intersect,REVEL_thresh,CADD_thresh,subgroup):
    ## 作MRR图统计
    columns_fill = output_columns(hpo,tools,intersect,REVEL_thresh,CADD_thresh)
    print(f'copy this class name to R script:\n')
    merge = 'merge('
    # for column in df.columns:
    #     if 'rank' in column:
    #         MRR_script = f"{column}_MRR = aggregate(df${column}, by = list(type=df$Symbol), mean)\nnames({column}_MRR)=c('Symbol','{column}')\n"
    #         print(MRR_script)
    #         merge = merge + f'{column}_MRR,'
    # merge = merge + f"by = '{subgroup}')"
    print(merge)



def collect_variants(pwd,df, filename='scoliosis_gVCF_from_zs_updating'):
    if 'collected_variants.tsv' in os.listdir(f'{pwd}/{filename}'):
        os.system(f'rm -f {pwd}/{filename}/collected_variants.tsv')
    for i in tqdm(range(len(df))):
        symbol = df.loc[i, "Symbol"]
        ensemblid = df.loc[i, 'Ensembl Gene ID']
        print(ensemblid)
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


def readobo(df,hpo,df_organ):
    # Read the taxrank ontology
    url = 'http://purl.obolibrary.org/obo/hp.obo'
    graph = obonet.read_obo(url)
    groupnames = df_organ['Term']
    df[f'{hpo}_organ_system']= ""
    df[f'{hpo}_organ_system_number'] = np.nan
    for j in tqdm(range(len(df))):
        groupname = []
        if df.loc[j, hpo][0] == '[':
            hpo_id = df.loc[j, hpo][2:-2]
            hpo_id_input = [k for k in hpo_id.split("', '")]
        else:
            hpo_id = df.loc[j, hpo]
            hpo_id_input = [k for k in hpo_id.split(";")]
        for hpo_term in hpo_id_input:
            searchgroup = networkx.descendants(graph,hpo_term)
            groupname = groupname+[k for k in groupnames if k in searchgroup]
        groupname = list(set(groupname))
        df.loc[j,f'{hpo}_organ_system'] = "".join(['+'+name for name in groupname])[1:]
        df.loc[j,f'{hpo}_organ_system_number'] = int(len(groupname))
    print(df[f'{hpo}_organ_system'],df[f'{hpo}_organ_system_number'])
    return df
    # # Number of nodes
    # print(len(graph))
    # print(networkx.is_directed_acyclic_graph(graph))
    # print(networkx.descendants(graph,'HP:0000118'))


    # df_organ = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx', ' organ_system')
    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    # for i in df['Symbol']:
    #     print(searchHGNCid(i))