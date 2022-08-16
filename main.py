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


def main_2(tools=[],CADD_thresh=[],REVEL_thresh=[],hpo='hpo_id',intersect=False,statistic = True,refresh=False,newcase_runtool=False,newvcfgz=False,collect_variants_tofile = False): ##无精炼HPO版本

    filename = 'scoliosis_gVCF_from_zs_updating'
    pwd = '/Users/liyaqi/PycharmProjects/PhenoAptAnalysis'

    df_organ = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx',' organ_system')

    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    #df = df[:1].reset_index(drop=True)
    case_ids = df['Blood ID']
    df = readobo(df,hpo,df_organ)
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
            get_HANRD(df)
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

        if 'phenoapt' in tools:
            if hpo != 'Weight':
                phenoapt_rank(df, pwd, filename, hpo, refresh)
            else:
                phenoapt_rank(df, pwd, filename, hpo, refresh=refresh)#只有重点hpo，不加权
                phenoapt_rank(df,pwd,filename,hpo,refresh=refresh, all_hpo_plus_weight=True)#有所有hpo，加权重点hpo
                phenoapt_rank(df, pwd, filename, hpo, refresh=refresh, only_weight_hpo_weight=True)#只有重点hpo，加权重点hpo

    if statistic:
        statistic_to_table(df,case_ids,pwd,filename,hpo,tools,intersect,REVEL_thresh,CADD_thresh)

    if collect_variants_tofile:
        collect_variants(pwd, df, filename)

    Rscript_brief_benchmark(tools, intersect, REVEL_thresh, CADD_thresh)

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


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    filename = 'scoliosis_gVCF_from_zs_updating'
    pwd = '/Users/liyaqi/PycharmProjects/PhenoAptAnalysis'
    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    # df = df[:1].reset_index(drop=True)
    # main_2(tools=['phrank','phenolyzer','GADO','phen2gene','phenoapt','LIRICAL','Exomiser','Phen_gen'],CADD_thresh=[10,15,20],REVEL_thresh=[0.25,0.5,0.75],hpo='hpo_id',intersect=True,refresh=False,newcase=True,statistic=True)
    # Rscript_MRR_matrix(df)
    # collect_variants(pwd,df, filename)
    # generateyml1310_zs_diag_zs_gVCF(df, pwd, filename, 'Weight', refresh=True)
    # get_HANRD(df,pwd,filename,'Weight')
    # pubcasefinder(df, pwd, filename, 'hpo_id', refresh=True)
    get_PhenoRank(df, pwd, filename, 'hpo_id')

    # See PyCharm help at https://www.jetbrains.com/help/pycharm/

