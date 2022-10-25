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


def main_2(tools=[],CADD_thresh=[],REVEL_thresh=[],hpo='hpo_id',intersect=False,statistic = True,refresh=False,newcase_runtool=False,newvcfgz=False,collect_variants_tofile = False,collect_organ = True): ##无精炼HPO版本

    filename = 'scoliosis_gVCF_from_zs_updating'
    pwd = '/Users/liyaqi/PycharmProjects/PhenoAptAnalysis'

    df_organ = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx',' organ_system')
    df_sex = read_xlsx('~/Documents/生信/gvcf-scoliosis-2022.xlsx', 'case')
    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')

    ## 可以按需要修改manual vcf tsv的地址，在get_case_id_file()中输入
    #df = df[7:8].reset_index(drop=True)
    case_ids = df['Blood ID']
    if collect_organ:
        df = readobo(df,hpo,df_organ,save=False)
    print(f"{len(df)}")

    if newcase_runtool:
        ##所有预先需要的文件准备
        if newvcfgz:
            getvcf_from_zs_gz(pwd, filename)
        file_prepare(df, pwd, filename, REVEL_thresh=REVEL_thresh, CADD_thresh=CADD_thresh, refresh=refresh)
        if 'phen2gene' in tools:
            phen2gene_rank(df,pwd=pwd, filename=filename, profile='no_gene_list',refresh=refresh)
            ## phen2gene直接rank版本
            if intersect:
                phen2gene_rank(df,pwd=pwd, filename=filename, profile=f'candidategene',refresh=refresh)
            if len(REVEL_thresh)!=0:
                for thresh in REVEL_thresh:
                    phen2gene_rank(df,pwd=pwd, filename=filename, profile=f'revel_{thresh}',refresh=True)
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
        if 'Exomiser13.1.0' in tools:
            generateyml1310_zs_diag_zs_gVCF(df,pwd,filename,hpo,refresh=refresh)
            print('Exomiser 13.1.0 ready')
            ## 新的hpo，新的case
        if 'HANRD' in tools:
            get_HANRD(df,pwd,filename,hpo,refresh=refresh,refresh_seperate_file=refresh)
            print('HANRD ready')
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
                    phrank_rank(df,pwd,filename, mode='intersect', hpo=hpo, filter=f'revel_{thresh}_ensemblid', refresh=True)
            ## phrank检查是否有漏跑REVEL
            if len(CADD_thresh) != 0:
                for thresh in CADD_thresh:
                    phrank_rank(df,pwd,filename, mode='intersect', hpo=hpo, filter=f'cadd_{thresh}_ensemblid', refresh=refresh)
            ## phrank检查是否有漏跑CADD
            print('phrank ready')

        if 'PhenoRank' in tools:
            get_PhenoRank(df, pwd, filename, hpo,refresh=refresh)
            input('PhenoRank finish with no HPO error? input any character to continue')
            get_PhenoRank(df, pwd, filename, hpo, refresh=False)##只补做上一波HPO有bug的case，不用refresh
            print('PhenoRank ready')


        if 'phenoapt' in tools:
            if hpo != 'Weight':
                phenoapt_rank(df, pwd, filename, hpo, refresh)
            else:
                phenoapt_rank(df, pwd, filename, hpo, refresh=refresh)#只有重点hpo，不加权
                phenoapt_rank(df,pwd,filename,hpo,refresh=refresh, all_hpo_plus_weight=True)#有所有hpo，加权重点hpo
                phenoapt_rank(df, pwd, filename, hpo, refresh=refresh, only_weight_hpo_weight=True)#只有重点hpo，加权重点hpo

        if 'phenomizer' in tools:
            print('manually search on website\nhttps://compbio.charite.de/phenomizer/')
            os.system(f'cd {pwd}/{filename}/phenomizeroutput && bash delet2df.sh')

        if 'phen-gen' in tools:
            generateped(df,pwd,filename,refresh,df_sex)
            print('manually run on PUMCH server and download the rank result')

    if statistic:
        statistic_to_table(df,case_ids,pwd,filename,hpo,tools,intersect,REVEL_thresh,CADD_thresh)

    if collect_variants_tofile:
        collect_variants(pwd, df, filename)

    Rscript_brief_benchmark(tools, intersect, REVEL_thresh, CADD_thresh)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    filename = 'scoliosis_gVCF_from_zs_updating'
    pwd = '/Users/liyaqi/PycharmProjects/PhenoAptAnalysis'
    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    tools = ['phenoapt', 'phenomizer', 'phenolyzer', 'GADO', 'PhenoRank', 'phen2gene', 'phrank', 'HANRD', 'LIRICAL',
                      'Exomiser','Exomiser13.1.0', 'Phen_gen']
    #df = df[17:18].reset_index(drop=True)
    compare_df_for_R(df, tools, pwd, 'wilcoxon', dMRR_heatmap=False)



    ##给df注释器官信息
    #df_organ = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx', ' organ_system')
    # df=pd.read_csv(f'{pwd}/TF_IDF_weight_organ_2.tsv',sep='\t')
    # result=readobo(df, 'hpo_id', df_organ)

    #重复器官的展开
    #TF_IDF_distribution_R()

    ## searchTF_IDF_for_R(df, tools, pwd,filename)
    # subgroup(df,tools,pwd,filename)
    #tools = ['phenomizer', 'phenolyzer', 'GADO', 'PhenoRank', 'phen2gene', 'phrank', 'HANRD', 'LIRICAL',
    #                  'Exomiser','Exomiser13.1.0', 'Phen_gen']
    # compare_df_for_R(df,tools,pwd,'hpo_id_Weight',dMRR_heatmap=True)
    # df0_phenoapt_4_strategy_rr(df,tools,pwd)
    # main_2(tools=['phrank','phenolyzer','GADO','phen2gene','phenoapt','LIRICAL','Exomiser','Phen_gen'],CADD_thresh=[10,15,20],REVEL_thresh=[0.25,0.5,0.75],hpo='hpo_id',intersect=True,refresh=False,newcase=True,statistic=True)
    # Rscript_MRR_matrix(df)
    # collect_variants(pwd,df, filename)
    # generateyml1310_zs_diag_zs_gVCF(df, pwd, filename, 'Weight', refresh=True)
    # get_HANRD(df,pwd,filename,'hpo_id',False,True)
    # pubcasefinder(df, pwd, filename, 'hpo_id', refresh=True)
    # file_prepare(df, pwd, filename,REVEL_thresh=[0.25,0.5,0.75],CADD_thresh=[10,15,20],refresh=False)
    # get_PhenoRank(df, pwd, filename, 'hpo_id')
    # main_2(tools=['phrank','phenolyzer','Exomiser13.1.0','GADO','PhenoRank','phen2gene','phenoapt','LIRICAL','HANRD','Exomiser','Phen_gen','phenomizer'],CADD_thresh=[10,15,20],REVEL_thresh=[0.25,0.5,0.75],hpo='hpo_id',intersect=True,refresh=False,newcase_runtool=False,statistic=True,collect_organ = True)

    ## newcase to run SNV and CNV
    # main_2(tools=['phrank','phenolyzer','Exomiser13.1.0','GADO','PhenoRank','phen2gene','phenoapt','LIRICAL','HANRD','Exomiser','Phen_gen','phenomizer'], hpo='Weight', intersect=True, refresh=False, newcase_runtool=False, statistic=True, collect_organ=True)
    # cnv_FA(pwd,filename,df,True)
    # See PyCharm help at https://www.jetbrains.com/help/pycharm/

