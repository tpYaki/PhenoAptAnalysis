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
from functions import *

def ensemblID(genelist):
    geneinput = ''
    for i in genelist:
        geneinput = geneinput + ','+i
    geneinput = geneinput[1:]
    url = f'https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi?method=db2db&input=genesymbolandsynonyms&inputValues={geneinput}&outputs=ensemblgeneid&taxonId=9606&format=row'
    u = urlre.urlopen(url)
    response = u.read()
    response = json.loads(response)
    output = {}
    for j in range(len(genelist)):
        output[j] = response[j]['Ensembl Gene ID']
    return output

def generatevcf_jeff_diag_dropbox_tsv():
    df, case_ids = read_diagnose_xlsx('/Users/liyaqi/Documents/生信/Inhouse_cohorts_genes_Version_8_MRR_诊断.xlsx')
    print(f"{len(df)}")
    df = df[:1]
    case_id_tsv_file_dict = get_case_id_file_map(case_ids)
    for i in range(len(df)):
        try:
            case_id = case_ids[i]
            if case_id not in case_id_tsv_file_dict:
                continue
            print(f'{case_id}')
            getvcf_from_dropbox_tsv(case_id, case_id_tsv_file_dict[case_id], 'scoliosis_filtered_tsv_dropbox/VCF')
        except Exception as e:
            print(e)

def zs_diag_jeff_gVCF_to_vcf(pwd):  ##gVCF-scolisosis
    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    sequence_ids = df['sequence ID']
    print(f"{len(df)}")
    dir_scoliosis_gVCF_from_jeff = f"/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_jeff/"
    ##df = df[:1]
    ## os.system(f'cd /Users/liyaqi/PycharmProjects/PhenoAptAnalysis && cat scoliosis_2021Sep.refseq.e0.001.i0.001.tsv | awk '{FS="\t";print $4}' > ./scoliosis_gVCF_from_jeff/mutation_id.tsv') ##直接在terminal里输入就可以了
    for i in range(len(df)):
        sequence_id = sequence_ids[i]
        print(f'{i}, {sequence_id}')
        cmd = f'cat {dir_scoliosis_gVCF_from_jeff}mutation_id.tsv |grep -n {sequence_id} | cut -d : -f 1 > {dir_scoliosis_gVCF_from_jeff}realtsv/{sequence_id}_row_number.txt'
        os.system(cmd)
        print('row done')

    ###找到非空文件
    init_list = subprocess.getoutput(f'zsh {dir_scoliosis_gVCF_from_jeff}non_empty.sh')
    init_list = str(init_list).split('\n')
    ## 从scoliosis大表里拉出VCF

    print(len(init_list))
    m = 1
    for sequence_id in init_list:
        case_id = sequence_id.split('-')[0]
        cmd = f'awk -f {dir_scoliosis_gVCF_from_jeff}greprowline.awk {dir_scoliosis_gVCF_from_jeff}realtsv/{sequence_id}_row_number.txt {pwd}scoliosis_2021Sep.refseq.e0.001.i0.001.tsv >  {dir_scoliosis_gVCF_from_jeff}tsvtoVCF/{case_id}.tsv'
        os.system(cmd)
        print(f"{sequence_id}tsv done")
        dfx = pd.read_csv(f'{dir_scoliosis_gVCF_from_jeff}tsvtoVCF/{case_id}.tsv', sep='\t', header=None)
        dfx[['CHR', 'POS', 'REF', 'ALT']] = dfx[0].str.split('_', 0, expand=True)
        dfx_chr = {}
        dfx_genotype = {}
        for i in range(len(dfx)):
            dfx_chr[i] = dfx.at[i, 'CHR'][3:]
            ids = [k for k in dfx.at[i, 3].split(';')]
            genotypes = [k for k in dfx.at[i, 5].split(';')]
            for k in range(len(ids)):
                if ids[k] == sequence_id:
                    dfx_genotype[i] = genotypes[k]
        dfx_chr = pd.Series(dfx_chr)
        dfx_genotype = pd.Series(dfx_genotype)
        QUAL = pd.DataFrame({'QUAL': pd.Series(['.' for k in range(len(dfx))])})
        INFO = pd.DataFrame({'INFO': pd.Series(['.' for k in range(len(dfx))])})
        FILTER = pd.DataFrame({'INFO': pd.Series(['PASS' for k in range(len(dfx))])})
        df2 = pd.concat([dfx_chr, dfx[['POS', 26, 'REF', 'ALT']], QUAL, FILTER, INFO, dfx[4], dfx_genotype], axis=1)
        df2.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
        df2.to_csv(f"./scoliosis_gVCF_from_jeff/realVCF/{case_id}.txt", sep="\t", index=False)
        command = f"source /Users/liyaqi/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && awk -f {dir_scoliosis_gVCF_from_jeff}script.awk {dir_scoliosis_gVCF_from_jeff}realVCF/{case_id}.txt | bcftools view -o {dir_scoliosis_gVCF_from_jeff}realVCF/{case_id}.vcf"
        print(command)
        os.system(command)
        print(f'{case_id}, NO.{m}, vcf done')
        m = m + 1


def generateyml_jeff_diag_dropbox_tsv():
    df, case_ids = read_diagnose_xlsx('/Users/liyaqi/Documents/生信/Inhouse_cohorts_genes_Version_8_MRR_诊断.xlsx')
    print(f"{len(df)}")
    df = df[:1]
    case_id_tsv_file_dict = get_case_id_file_map(case_ids)
    dir = f'./scoliosis_filtered_tsv_dropbox/'
    for i in range(len(df)):
        try:
            case_id = df.loc[i, 'CaseID']
            if case_id not in case_id_tsv_file_dict:
                print(f"{i} not in dropbox")
                continue

            hpo_id = df.loc[i, 'hpo_id']
            hpo_id_input = [k for k in hpo_id.split(";")]

            getyml(case_id, hpo_id_input, dir)  ##用map到的所有文件名case_id_file称生成yml,与main中的case_id一个意思
        except Exception as e:
            print(e)

def generateyml_zs_diag_jeff_gVCF():  ##gVCF-scolisosis
    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    init_list = subprocess.getoutput(
        f'zsh /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_jeff/non_empty.sh')
    init_list = str(init_list).split('\n')
    ##init_list = init_list[:1]
    dfa = pd.DataFrame()
    m = 0
    for sequence_id in init_list:
        for i in range(len(df)):
            if df.at[i, 'sequence ID'] == sequence_id:
                case_id = sequence_id.split('-')[0]
                dfa.at[case_id, 'HPO'] = df.at[i, 'HPO']
                dfa.at[case_id, 'Symbol'] = df.at[i, 'Symbol']
    dfa = dfa.rename_axis('case_id').reset_index(level=0)
    print(len(dfa))
    dfa.to_csv('/Users/liyaqi/Documents/生信/gVCF-2022-确诊-有VCF.csv')  ##dfa是能在scolisis大表中找到VCF条目的ID-HPO-Symbol大表
    dir = f'/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_jeff/'
    for t in range(len(dfa)):
        hpo_id = dfa.loc[t, 'HPO'][2:-2]
        hpo_id_input = [str(k) for k in hpo_id.split("', '")]
        getyml(dfa.at[t, 'case_id'], hpo_id_input, dir)
        case_id = dfa.at[t, 'case_id']
        ##os.system(f'cd /Users/liyaqi/Downloads/exomiser-cli-12.1.0 && java -Xms2g -Xmx4g -jar exomiser-cli-12.1.0.jar -analysis /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_jeff/realyml/{case_id}.yml') ##需要跑exomiser的时候打开
        print(f'{case_id}, NO.{m}, exomiser done')
        m = m + 1
