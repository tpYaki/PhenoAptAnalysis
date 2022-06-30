# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
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


pwd = "/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/"

# def ensemblID(genelist):
#     geneinput = ''
#     for i in genelist:
#         geneinput = geneinput + ','+i
#     geneinput = geneinput[1:]
#     url = f'https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi?method=db2db&input=genesymbolandsynonyms&inputValues={geneinput}&outputs=ensemblgeneid&taxonId=9606&format=row'
#     u = urlre.urlopen(url)
#     response = u.read()
#     response = json.loads(response)
#     output = {}
#     for j in range(len(genelist)):
#         output[j] = response[j]['Ensembl Gene ID']
#     return output

#读取诊断表格
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


#获取case_id与相应的tsv文件路径的映射关系
def get_case_id_file_map(case_ids):
    case_id_tsv_file_map = {}
    tsv_dir = '/Users/liyaqi/Dropbox/Filtered-SNV-all/Sporadic-WES-All'
    files = os.listdir(tsv_dir)
    for case_id in case_ids:
        for file_name in files:
            if case_id in file_name:
                case_id_tsv_file_map[case_id] = os.path.join(tsv_dir, file_name)
    return case_id_tsv_file_map

def get_case_id_file(case_id,mode):
    if mode =='sporadic':
        tsv_dir = '/Users/liyaqi/Dropbox/Filtered-SNV-all/Sporadic-WES-All'
        files = os.listdir(tsv_dir)
        for file_name in files:
            if case_id in file_name:
                case_id_tsv_file = os.path.join(tsv_dir, file_name)
                return case_id_tsv_file
        return 0
    else:
        if mode =='trio':
            tsv_dir = '/Users/liyaqi/Dropbox/Filtered-SNV-all/Trio-WES-all'
            files = os.listdir(tsv_dir)
            for file_name in files:
                if case_id in file_name:
                    case_id_tsv_file = os.path.join(tsv_dir, file_name)
                    return case_id_tsv_file
            return 0


def get_case_id_integ_file_map(integ,case_ids):
    case_id_tsv_file_map = {}
    tsv_dir = f'/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/{integ}output'
    files = os.listdir(tsv_dir)
    for case_id in case_ids:
        for file_name in files:
            if file_name.endswith('.tsv'):
                if case_id in file_name:
                    case_id_tsv_file_map[case_id] = os.path.join(tsv_dir, file_name)
    return case_id_tsv_file_map

def containtype(biao,type):
    data_mut = biao[biao["Mutation_type"].str.contains(type)]
    return data_mut
    # 包含字符串type=“”的rows，Na的单元格自动去掉

def selectsmall(biao,name,thresh):
    data_mut1 = biao[biao[name] != "."]
    data_mut3 = data_mut1[data_mut1[name] != "nan"]## df,"Exac"

    def is_float(x):
        try:
            float(x)
        except ValueError:
            return False
        return True
    data_mut12 = data_mut3[data_mut3[name].apply(lambda x: is_float(x))]
    data_mut12[name] = pd.to_numeric(data_mut12[name], errors='raise')
    data_mut2 = (data_mut12[(data_mut12[name]) < thresh].append(biao[biao[name] == "."])).append(biao[biao[name]=="nan"])
    return data_mut2

def selectbig(biao,name,thresh):
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
    data_mut2 = (data_mut12[data_mut12[name] > thresh].append(biao[biao[name] == "."])).append(biao[biao[name]=="nan"])
    return data_mut2

def patho_filter_n(df,threshold):
    ##data2 = ((containtype(data1, "missense_variant")).append(containtype(data1, "stop_gained"))).append(containtype(data1, "frameshift"))
    data3 = selectsmall(selectsmall(selectsmall(selectsmall(df,"ExAC_AF",0.01),"ExAC_EAS_AF",0.01),"gnomAD_genome_EAS",0.01),"In_house_freq",0.05)
    data4 = selectbig(data3, "VAF", 29)
    data5 = selectbig(data4, "CADD", threshold)
    ##data5 = (congrade(data4,"Mutationtatser_prediction","A")).append(data4[data4["Mutationtatser_prediction"].str.contains("D")])
    ##data7 = congrade(data6,"LRT_prediction","D").append(data6[data6["LRT_prediction"].str.contains("U")])
    ##data8 = congrade(data7, "Polyphen2_HDIV_pred", "D").append(data7[data7["Polyphen2_HDIV_pred"].str.contains("P")])
    ##data9 = congrade(data8, "Polyphen2_HVAR_pred", "D").append(data8[data8["Polyphen2_HVAR_pred"].str.contains("P")])
    ##data10 = congrade(data9, "SIFT_prediction", "D")
    return data5

def getvcf_from_dropbox_tsv(case_id,dir):
    df = pd.read_csv(dir, sep="\t")
    dfx = df
    for i in range(len(df)):
        # if df.at[i, 'ALT'] == '-':
        #     print(i, 'ALT', df.at[i, 'ALT'])
        # else:
        #     if df.at[i, 'REF'] == '-':
        #         print(i, 'REF', df.at[i, 'REF'])
        #     else:
        #         dfx=pd.concat([dfx, df.loc[i]],axis='columns') ##删掉indel

        ##补充indel的REF和ALT
        if df.at[i, 'ALT'] == '-': ##deletion
            print(i, 'ALT', df.at[i, 'ALT'])
            pos=int(df.at[i, 'POS'])
            chr=df.at[i, 'CHR']
            cmd=f'source /Users/liyaqi/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && samtools faidx /Users/liyaqi/Software/ensembl-vep/Homo_sapiens.GRCh37.dna.primary_assembly.fa {chr}:{pos-1}-{pos-1}'
            result = subprocess.getoutput(cmd)
            dfx.at[i, 'ALT']=result.split('\n')[1]
            dfx.at[i, 'REF']=result.split('\n')[1]+df.at[i, 'REF']
            dfx.at[i, 'POS']=pos-1
        else:
            if df.at[i, 'REF'] == '-': ## insertion
                print(i, 'REF', df.at[i, 'REF'])
                pos = int(df.at[i, 'POS'])
                chr = df.at[i, 'CHR']
                cmd = f'source /Users/liyaqi/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && samtools faidx /Users/liyaqi/Software/ensembl-vep/Homo_sapiens.GRCh37.dna.primary_assembly.fa {chr}:{pos - 1}-{pos - 1}'
                result = subprocess.getoutput(cmd)
                dfx.at[i, 'REF'] = result.split('\n')[1]
                dfx.at[i, 'ALT'] = result.split('\n')[1]+df.at[i, 'ALT']
                dfx.at[i, 'POS'] = pos - 1

    for j in range(len(dfx)):
        if str(dfx.at[j,'FILTER']) != '-':
            continue
        dfx.at[j,'FILTER'] = '.'

    QUAL = pd.DataFrame({'QUAL':pd.Series([100 for k in range(len(dfx))])})
    INFO = pd.DataFrame({'INFO': pd.Series(['.' for k in range(len(dfx))])})


    df2 = pd.concat([dfx[['CHR', 'POS', 'Rs_ID', 'REF', 'ALT']],QUAL,dfx['FILTER'],INFO,dfx['FORMAT']],axis=1)
    df2.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT']
    df2.to_csv(f"./scoliosis_filtered_tsv_dropbox/VCF/{case_id}.txt", sep="\t", index=False)
    vcfdir_jeff_diag_dropbox_tsv = f"/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_filtered_tsv_dropbox/"
    command = f"source /Users/liyaqi/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && awk -f {pwd}script.awk {vcfdir_jeff_diag_dropbox_tsv}VCF/{case_id}.txt | bcftools view -o {vcfdir_jeff_diag_dropbox_tsv}VCF/{case_id}.vcf"
    print(command)
    os.system(command)
    print(f" done")
    return True

def generatevcf_jeff_diag_dropbox_tsv():
    df, case_ids= read_diagnose_xlsx('/Users/liyaqi/Documents/生信/Inhouse_cohorts_genes_Version_8_MRR_诊断.xlsx')
    print(f"{len(df)}")
    df = df[:1]
    case_id_tsv_file_dict = get_case_id_file_map(case_ids)
    for i in range(len(df)):
        try:
            case_id = case_ids[i]
            if case_id not in case_id_tsv_file_dict:
                continue
            print(f'{case_id}')
            getvcf_from_dropbox_tsv(case_id, case_id_tsv_file_dict[case_id])
            print('done')
        except Exception as e:
            print(e)

def zs_diag_jeff_gVCF_to_vcf():  ##gVCF-scolisosis
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
        dfx[['CHR', 'POS','REF','ALT']] = dfx[0].str.split('_', 0, expand=True)
        dfx_chr={}
        dfx_genotype = {}
        for i in range(len(dfx)):
            dfx_chr[i]= dfx.at[i,'CHR'][3:]
            ids = [k for k in dfx.at[i, 3].split(';')]
            genotypes = [k for k in dfx.at[i, 5].split(';')]
            for k in range(len(ids)):
                if ids[k] == sequence_id:
                    dfx_genotype[i]=genotypes[k]
        dfx_chr=pd.Series(dfx_chr)
        dfx_genotype = pd.Series(dfx_genotype)
        QUAL = pd.DataFrame({'QUAL': pd.Series(['.' for k in range(len(dfx))])})
        INFO = pd.DataFrame({'INFO': pd.Series(['.' for k in range(len(dfx))])})
        FILTER = pd.DataFrame({'INFO': pd.Series(['PASS' for k in range(len(dfx))])})
        df2 = pd.concat([dfx_chr,dfx[['POS', 26, 'REF', 'ALT']],QUAL,FILTER,INFO,dfx[4],dfx_genotype],axis=1)
        df2.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT','SAMPLE']
        df2.to_csv(f"./scoliosis_gVCF_from_jeff/realVCF/{case_id}.txt", sep="\t", index=False)
        command = f"source /Users/liyaqi/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && awk -f {dir_scoliosis_gVCF_from_jeff}script.awk {dir_scoliosis_gVCF_from_jeff}realVCF/{case_id}.txt | bcftools view -o {dir_scoliosis_gVCF_from_jeff}realVCF/{case_id}.vcf"
        print(command)
        os.system(command)
        print(f'{case_id}, NO.{m}, vcf done')
        m = m + 1

def zs_diag_zs_gVCF_to_vcf():  ##gVCF-scolisosis
    dir_scoliosis_gVCF_from_zs = f"/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/"
    os.system(f'{dir_scoliosis_gVCF_from_zs}bash dividevcf.sh') ##直接在terminal里输入就可以了


import yaml
def getyml(case_id_file, hpo_id_input,dir):
    with open(f'./yml/test-analysis-exome.yml','r') as f:
        config = yaml.safe_load(f)
        config['analysis']['vcf'] = f'{dir}*vcf/{case_id_file}.vcf'
        config['analysis']['hpoIds'] = hpo_id_input # add the command as a list for the correct yaml
        config['outputOptions']['outputFormats']=['TSV_GENE']
        config['outputOptions']['outputPrefix'] = f'{dir}Exomiseroutput/{case_id_file}'
        ##config['inheritanceModes'] = inheri
        print(config['analysis']['hpoIds'])

    with open(f'{dir}*yml/{case_id_file}.yml',"w") as f: # open the file in append mode
        f.truncate(0)
        yaml.dump(config, f) ##default_flow_style=Faulse

def getyml_for_certain_ingeritance(case_id_file, hpo_id_input,dir,inheridict):
    with open(f'./yml/test-analysis-exome.yml','r') as f:
        config = yaml.safe_load(f)
        config['analysis']['vcf'] = f'{dir}vcf/{case_id_file}.vcf'
        config['analysis']['hpoIds'] = hpo_id_input # add the command as a list for the correct yaml
        config['outputOptions']['outputFormats']=['TSV_GENE']
        config['outputOptions']['outputPrefix'] = f'{dir}Exomiseroutput/{case_id_file}'
        config['analysis']['inheritanceModes'] = inheridict
        print(config['analysis']['hpoIds'])
        print(config['analysis']['inheritanceModes'])

    with open(f'{dir}yml/{case_id_file}.yml',"a") as f: # open the file in append mode
        f.truncate(0)
        yaml.dump(config, f) ##default_flow_style=Faulse

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

            getyml(case_id,hpo_id_input,dir) ##用map到的所有文件名case_id_file称生成yml,与main中的case_id一个意思
        except Exception as e:
            print(e)

def generateyml_zs_diag_jeff_gVCF(): ##gVCF-scolisosis
    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    init_list = subprocess.getoutput(f'zsh /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_jeff/non_empty.sh')
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
    dfa.to_csv('/Users/liyaqi/Documents/生信/gVCF-2022-确诊-有VCF.csv') ##dfa是能在scolisis大表中找到VCF条目的ID-HPO-Symbol大表
    dir = f'/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_jeff/'
    for t in range(len(dfa)):
        hpo_id = dfa.loc[t, 'HPO'][2:-2]
        hpo_id_input = [str(k) for k in hpo_id.split("', '")]
        getyml(dfa.at[t,'case_id'], hpo_id_input,dir)
        case_id = dfa.at[t,'case_id']
        ##os.system(f'cd /Users/liyaqi/Downloads/exomiser-cli-12.1.0 && java -Xms2g -Xmx4g -jar exomiser-cli-12.1.0.jar -analysis /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_jeff/realyml/{case_id}.yml') ##需要跑exomiser的时候打开
        print(f'{case_id}, NO.{m}, exomiser done')
        m = m + 1

def generateyml_zs_diag_zs_gVCF():
    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    df = df[36:37].reset_index(drop=True)
    sequence_ids = df['sequence ID']
    print(f"{len(df)}")
    dir_scoliosis_gVCF_from_zs = f"/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/"
    sequence_id_file_dict = os.listdir(f'{dir_scoliosis_gVCF_from_zs}vcf/')
    ##print(sequence_id_file_dict)
    m = 0
    n = 0
    for i in range(len(df)):
        sequence_id = sequence_ids[i]
        sequence_id_file = f'{sequence_id}.vcf'
        case_id = sequence_id.split('-')[0]
        if sequence_id_file not in sequence_id_file_dict:
            print(f"No.{i} not in scoliosis_gVCF_from_zs_updating")
            n = n+1
            continue
        else:
            print("ok")
            if df.loc[i, 'hpo_id'][0] == '[':
                hpo_id = df.loc[i, 'hpo_id'][2:-2]
                hpo_id_input = [str(k) for k in hpo_id.split("', '")]
            else:
                hpo_id = df.loc[i, 'hpo_id']
                hpo_id_input = [k for k in hpo_id.split(";")]
            inheri = (str(df.loc[i, 'inheritance'])).split(',') ##存在输入了多种mode的可能
            inheridict = {}
            for k in inheri:
                j=k.split(':')
                inheridict[f'{j[0]}'] = float(j[1])
            ##inheridict = [k for k in (df.loc[i, 'inheritance']).split(";")]
            print(inheridict)
            getyml_for_certain_ingeritance(sequence_id, hpo_id_input, dir_scoliosis_gVCF_from_zs,inheridict)
            command = f'cd /Users/liyaqi/Downloads/exomiser-cli-12.1.0 && java -Xms2g -Xmx4g -jar exomiser-cli-12.1.0.jar -analysis {dir_scoliosis_gVCF_from_zs}yml/{sequence_id}.yml'
            print(command)
            os.system(command)
            print(f'{case_id}, NO.{i}, yml done')
            m=m+1

    print(f'{m} case finish exomiser')
    print(f'{n} files were not in scoliosis_gVCF_from_zs_updating')


def get_phenopaket(sequence_id, hpo_id_input,dir):
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
        print("params", params)  # 打印
        f.close()  # 关闭json读模式
    with open(f'{dir}phenopacket/{sequence_id}.json', 'w') as r:
        json.dump(params, r)
    # 关闭json写模式
        r.close()

def generate_phenopaket_zs_diag_zs_gVCF():
    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    df = df[36:37].reset_index(drop=True)
    sequence_ids = df['sequence ID']
    print(f"{len(df)}")
    dir_scoliosis_gVCF_from_zs = f"/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/"
    sequence_id_file_dict = os.listdir(f'{dir_scoliosis_gVCF_from_zs}vcf/')
    ##print(sequence_id_file_dict)
    m = 0
    n = 0
    for i in range(len(df)):
        sequence_id = sequence_ids[i]
        sequence_id_file = f'{sequence_id}.vcf'
        case_id = sequence_id.split('-')[0]
        if sequence_id_file not in sequence_id_file_dict:
            print(f"No.{i} not in scoliosis_gVCF_from_zs_updating")
            n = n + 1
            continue
        else:
            print("ok")
            if df.loc[i, 'hpo_id'][0] == '[':
                hpo_id = df.loc[i, 'hpo_id'][2:-2]
                hpo_id_input = [str(k) for k in hpo_id.split("', '")]
            else:
                hpo_id = df.loc[i, 'hpo_id']
                hpo_id_input = [k for k in hpo_id.split(";")]
            get_phenopaket(sequence_id, hpo_id_input, dir_scoliosis_gVCF_from_zs)
            print(f'{case_id}, NO.{i}, phenopaket done')
            m = m + 1
            command = f'cd /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/LIRICALoutput && java -jar LIRICAL.jar phenopacket -p {dir_scoliosis_gVCF_from_zs}phenopacket/{sequence_id}.json -d /Users/liyaqi/Downloads/LIRICAL-1.3.4/data -e /Users/liyaqi/Downloads/exomiser-cli-12.1.0/data/1909_hg19 -x {sequence_id} --tsv && bash outputtodf.sh'
            print(command)
            os.system(command)
        print(f'{m} case finish LIRICAL')
        print(f'{n} files were not in scoliosis_gVCF_from_zs_updating')

def generateped():
    df = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx', 'Sheet1')
    # df = df[:1].reset_index(drop=True)
    df_sex = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx', 'gVCF-zs-确诊-WES')
    df_sex = df_sex.rename(index=df_sex['Blood ID'])
    df_sex = dict(df_sex['Sex'])
    for i in range(len(df)):
        ##做ped
        ped_df = (pd.DataFrame(['FAM1', case_id, 'FATHER', 'MOTHER', df_sex.get(case_id), 2]).T)
        ped_df.to_csv(f'/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/ped/{case_id}.tsv', sep='\t', index=False, header=False)
        cmd = f'mv /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/ped/{case_id}.tsv /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/ped/{case_id}.ped'
        os.system(cmd)
        ##做好了ped

def generatehpotxt_GADO():
    df = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx', 'Sheet1')
    # df = df[:1].reset_index(drop=True)
    hpotxt = {}
    for i in range(len(df)):
        case_id = df['Blood ID'][i]
        if df.loc[i, 'hpo_id'][0] == '[':
            hpo_id = df.loc[i, 'hpo_id'][2:-2]
            hpo_id_input = [k for k in hpo_id.split("', '")]
        else:
            hpo_id = df.loc[i, 'hpo_id']
            hpo_id_input = [k for k in hpo_id.split(";")]
        ##所有case记录在hpotxt{}里
        hpotxt[case_id] = hpo_id_input
    ##把hpotxt存在txt里
    file = open("./scoliosis_gVCF_from_zs_updating/hpotxt/totalhpo.txt", "w")
    for key, value in hpotxt.items():
        hpoline = ''
        for k in value:
            hpoline = hpoline + '\t'+k
        file.write('%s%s\n' % (key,hpoline))
        print(file)
    file.close()


def GADO():
    generatehpotxt_GADO()
    cmd1 = 'cd /Users/liyaqi/Software/GadoCommandline-1.0.1 && java -jar GADO.jar \
 --mode PROCESS \
 --output /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/GADOoutput/hpoProcessed.txt \
 --caseHpo /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/hpotxt/totalhpo.txt \
 --hpoOntology data/hp.obo \
 --hpoPredictionsInfo data/predictions_auc_bonf.txt'
    cmd2 = 'cd /Users/liyaqi/Software/GadoCommandline-1.0.1 && java -jar GADO.jar \
 --mode PRIORITIZE \
 --output /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/GADOoutput/ \
 --caseHpoProcessed /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/GADOoutput/hpoProcessed.txt \
 --genes data/hpo_prediction_genes.txt \
 --hpoPredictions data/genenetwork_bonf_spiked/genenetwork_bonf_spiked'
    os.system(cmd1)
    os.system(cmd2)
    ##按照case_id输出

def hpotxt_phenolyzer():
    df = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx', 'Sheet1')
    # df = df[:1].reset_index(drop=True)
    for i in range(len(df)):
        case_id = df['Blood ID'][i]
        if df.loc[i, 'hpo_id'][0] == '[':
            hpo_id = df.loc[i, 'hpo_id'][2:-2]
            hpo_id_input = [k for k in hpo_id.split("', '")]
        else:
            hpo_id = df.loc[i, 'hpo_id']
            hpo_id_input = [k for k in hpo_id.split(";")]
        ##每个case单独一个txt-phen-gen
        hpo_id_input_df = pd.DataFrame(hpo_id_input)
        hpo_id_input_df.to_csv(
            f'/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/hpotxt/{case_id}.txt', sep = '\n',
            index=False, header=False)
        print(f"{i}txt done")

def phenolyzer():
    ##需要新生成hpotxt for phenolyzer时打开
    hpotxt_phenolyzer()
    ##运行phenolyzer
    df = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx', 'Sheet1')
    # df = df[:1].reset_index(drop=True)
    for i in tqdm(range(len(df))):
        case_id = df['Blood ID'][i]
        hpo_txt_dir = f'/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/hpotxt/{case_id}.txt'
        cmd_phenolyzer = f'cd /Users/liyaqi/Software/phenolyzer && perl disease_annotation.pl {hpo_txt_dir} -file -prediction -phenotype -logistic -out /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/phenolyzeroutput/{case_id} -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 && exit'
        print(cmd_phenolyzer)
        appscript.app('Terminal').do_script(cmd_phenolyzer)
        time.sleep(30)

def phrank_in_pipline(filename,ensemblid,mode,hpo_id_input,case_id):
    DAG = "./phrank/demo/data/hpodag.txt"
    DISEASE_TO_PHENO = "./phrank/demo/data/disease_to_pheno.build127.txt"
    DISEASE_TO_GENE = "./phrank/demo/data/gene_to_disease.build127.txt"
    GENE_TO_PHENO = "./phrank/demo/data/gene_to_pheno.amelie.txt"
    p_hpo = Phrank(DAG, diseaseannotationsfile=DISEASE_TO_PHENO, diseasegenefile=DISEASE_TO_GENE)
    if mode =='intersect':
        patient_genes = pd.read_csv(f'./{filename}/candidategene_ensemblid/{case_id}.txt', header=None)
        phrank_gene_ranking = p_hpo.rank_genes_directly(patient_genes[0], hpo_id_input)
        phrank_gene_ranking = list(pd.DataFrame(phrank_gene_ranking)[1])
        if ensemblid in phrank_gene_ranking:
            Phrank_rank_intersect = phrank_gene_ranking.index(ensemblid)
        else:
            Phrank_rank_intersect = 'NA'
        return Phrank_rank_intersect
    else:
        if mode =='OMIM':
            patient_genes_df = pd.read_csv(f'./{filename}/omim-non-empty.txt', header=None)
            patient_genes = patient_genes_df[0]
            phrank_gene_ranking = p_hpo.rank_genes_directly(patient_genes, hpo_id_input)
            phrank_gene_ranking = list(pd.DataFrame(phrank_gene_ranking)[1])
            if ensemblid in phrank_gene_ranking:
                Phrank_rank = phrank_gene_ranking.index(ensemblid)
            else:
                Phrank_rank = 'NA'
            return Phrank_rank
        else:
            return 'wrong mode input'


def main_1(): ##使用精炼HPO版本
    filename = 'scoliosis_gVCF_from_zs_updating'
    pwd = '/Users/liyaqi/PycharmProjects/PhenoAptAnalysis'
    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    ##df = df[22:23].reset_index(drop=True)
    case_ids = df['Blood ID']
    print(f"{len(df)}")
    final_big_table = pd.DataFrame(
        columns=['CaseID', 'hpo_id', 'Symbol',
                 'phenoapt_rank', 'intersect_rank', 'Patho_rank_CADD_10', 'Patho_rank_CADD_15', 'Patho_rank_CADD_20',
                 'phrank_rank_intersect',
                 'phen2gene','phen2gene_intersect',
                 'Exomiser','LIRICAL'])
    if len(get_case_id_file_map(case_ids))==0:
        case_id_filtered_tsv_file_dict = get_case_id_file_trio(case_ids)
    else:
        case_id_filtered_tsv_file_dict = get_case_id_file(case_ids)
    case_id_exomiser_tsv_file_dict = get_case_id_integ_file_map('Exomiser',case_ids)
    case_id_LIRICAL_tsv_file_dict = get_case_id_integ_file_map('LIRICAL', case_ids)
    for i in range(len(df)):
        # try:
            if df.loc[i, 'hpo_id'][0] == '[':
                hpo_id = df.loc[i, 'hpo_id'][2:-2]
                hpo_id_input = [k for k in hpo_id.split("', '")]
            else:
                hpo_id = df.loc[i, 'hpo_id']
                hpo_id_input = [k for k in hpo_id.split(";")]
            if df.loc[i, 'Weight'][0] == '[':
                weight_HPO_id = df.loc[i, 'Weight'][2:-2]
                weight_HPO_id_input = [k for k in weight_HPO_id.split("', '")]
            else:
                weight_HPO_id = df.loc[i, 'Weight']
                weight_HPO_id_input = [k for k in weight_HPO_id.split(";")]
            symbol = df.loc[i, "Symbol"]
            entrezGeneId = df.loc[i,'entrezGeneId']
            entrezId = df.loc[i,'entrezId']
            ensemblid = df.loc[i,'Ensembl Gene ID']
            print(entrezId)
            case_id = df.loc[i,'Blood ID']

            variation = pd.read_csv(case_id_filtered_tsv_file_dict[case_id], sep="\t")
            variation_gene_name = variation['Gene_name']


            candidate_gene_list = variation_gene_name[~variation_gene_name.duplicated()]
            # #需要产生新的phen2gene分析使用的candidate文件时候打开
            # candidate_gene_list.to_csv(f'./{filename}/candidategene/{case_id}.txt',index=False,header=False)

            ## phrank_intersect
            DAG = "./phrank/demo/data/hpodag.txt"
            DISEASE_TO_PHENO = "./phrank/demo/data/disease_to_pheno.build127.txt"
            DISEASE_TO_GENE = "./phrank/demo/data/gene_to_disease.build127.txt"
            GENE_TO_PHENO = "./phrank/demo/data/gene_to_pheno.amelie.txt"
            p_hpo = Phrank(DAG, diseaseannotationsfile=DISEASE_TO_PHENO, diseasegenefile=DISEASE_TO_GENE)
            output = ensemblID(candidate_gene_list)
            patient_genes = set(output[k] for k in output.keys())
            phrank_gene_ranking = p_hpo.rank_genes_directly(patient_genes, hpo_id_input)
            phrank_gene_ranking = list(pd.DataFrame(phrank_gene_ranking)[1])
            if ensemblid in phrank_gene_ranking:
                Phrank_rank_intersect = phrank_gene_ranking.index(ensemblid)
            else:
                Phrank_rank_intersect = 'NA'

            temp = ''
            for hpo in hpo_id_input:
                temp = temp + ' '+hpo
            HPO_for_Phen2Gene = temp
            # #phen2gene intersect
            command_phen2gene_intersect = f'cd /Users/liyaqi/Phen2Gene && python3 phen2gene.py -m{HPO_for_Phen2Gene} -out {pwd}/{filename}/phen2geneoutput -l {pwd}/{filename}/candidategene/{case_id}.txt && mv {pwd}/{filename}/phen2geneoutput/output_file.associated_gene_list {pwd}/{filename}/phen2geneoutput/{case_id}-phen2gene-output.txt'
            os.system(command_phen2gene_intersect)
            phen2gene_intersect_result = pd.read_csv(f'{pwd}/{filename}/phen2geneoutput/{case_id}-phen2gene-output.txt', sep="\t")
            patho_gene_rank_phen2gene_intersect = {}
            for j, v in enumerate(phen2gene_intersect_result["Gene"]):
                patho_gene_rank_phen2gene_intersect[v] = j

            # #phen2gene
            command_phen2gene = f'cd /Users/liyaqi/Phen2Gene && python3 phen2gene.py -m{HPO_for_Phen2Gene} -out {pwd}/{filename}/phen2geneoutput_nolist && mv {pwd}/{filename}/phen2geneoutput_nolist/output_file.associated_gene_list {pwd}/{filename}/phen2geneoutput_nolist/{case_id}-phen2gene-nolist-output.txt'
            os.system(command_phen2gene)
            phen2gene_result = pd.read_csv(
                f'{pwd}/{filename}/phen2geneoutput_nolist/{case_id}-phen2gene-nolist-output.txt', sep="\t")
            if symbol in list(phen2gene_result["Gene"]):
                patho_gene_rank_phen2gene = list(phen2gene_result["Gene"]).index(symbol)
            else:
                patho_gene_rank_phen2gene = 'NA'


            ## PhenoApt排序
            weight_1 = [1 for k in range(len(hpo_id_input))]
            client = PhenoApt(token='H0pVk00CX07VzkZbdnvHI$24XiU$u9q')
            pheno_result = (client.rank_gene(phenotype=hpo_id_input,weight=weight_1,n=5000)).rank_frame
            print(hpo_id_input)
            print(pheno_result)
            pheno_gene_rank = {}
            for j, v in enumerate(pheno_result["gene_symbol"]):
                pheno_gene_rank[v] = j
            print(pheno_gene_rank)


            # #仅取TSV基因的交集排序
            intersect_gene = pheno_result[pheno_result.gene_symbol.isin(variation_gene_name)]
            intersect_gene_rank = {}
            for j, v in enumerate(intersect_gene["gene_symbol"]):
                intersect_gene_rank[v] = j ##写出rank
            # filtered_result = variation[variation.Gene_name.isin(intersect_gene['gene_symbol'])]
            # filtered_result['pheno_rank'] = [pheno_gene_rank[gene] for gene in filtered_result['Gene_name']]
            # filtered_result['intersect_rank'] = [intersect_gene_rank[gene] for gene in filtered_result['Gene_name']]
            # filtered_result['CaseID'] = [case_id for gene in filtered_result['Gene_name']]
            # filtered_result.to_csv(f"output/{case_id}.csv")

            variation_p_10 = patho_filter_n(variation, 10)
            patho_gene_name_10 = variation_p_10['Gene_name']
            patho_gene_10 = pheno_result[pheno_result.gene_symbol.isin(patho_gene_name_10)]
            patho_gene_rank_10 = {}
            for j, v in enumerate(patho_gene_10["gene_symbol"]):
                patho_gene_rank_10[v] = j

            variation_p_15 = patho_filter_n(variation, 15)
            patho_gene_name_15 = variation_p_15['Gene_name']
            patho_gene_15 = pheno_result[pheno_result.gene_symbol.isin(patho_gene_name_15)]
            patho_gene_rank_15 = {}
            for j, v in enumerate(patho_gene_15["gene_symbol"]):
                patho_gene_rank_15[v] = j

            variation_p_20 = patho_filter_n(variation, 20)
            patho_gene_name_20 = variation_p_20['Gene_name']
            patho_gene_20 = pheno_result[pheno_result.gene_symbol.isin(patho_gene_name_20)]
            patho_gene_rank_20 = {}
            for j, v in enumerate(patho_gene_20["gene_symbol"]):
                patho_gene_rank_20[v] = j

            patho_gene_rank_Exo = {}
            if case_id in list(case_id_exomiser_tsv_file_dict):
                rankdf = pd.read_csv(case_id_exomiser_tsv_file_dict[case_id], sep="\t")
                patho_gene_ID_Exo = rankdf['ENTREZ_GENE_ID']
                for j, v in enumerate(patho_gene_ID_Exo):
                    patho_gene_rank_Exo[v] = j
                print(patho_gene_rank_Exo)
            else:
                patho_gene_rank_Exo[entrezId]= 'No_vcf_ranked_by_Exo'

            patho_gene_rank_LIR = {}
            if case_id in list(case_id_LIRICAL_tsv_file_dict):
                rankdf2 = pd.read_csv(case_id_LIRICAL_tsv_file_dict[case_id], sep="\t")
                patho_gene_ID_LIR = rankdf2['entrezGeneId']
                patho_gene_rank_LIR = {}
                df_count = pd.DataFrame(enumerate(patho_gene_ID_LIR))
                for s in range(len(df_count)):
                    k = (df_count[df_count[1].isin([df_count.loc[s,1]])]).reset_index(drop=True)
                    print(k)
                    patho_gene_rank_LIR[k.loc[0,1]] = k.loc[0,0]
                print(patho_gene_rank_LIR)
            else:
                patho_gene_rank_LIR[entrezGeneId]= 'No_vcf_ranked_by_LIR'

            final_big_table.loc[i] = [case_id, hpo_id, symbol,
                                      pheno_gene_rank.get(symbol, 'NA'), intersect_gene_rank.get(symbol, 'NA'), patho_gene_rank_10.get(symbol, 'NA'), patho_gene_rank_15.get(symbol, 'NA'),patho_gene_rank_20.get(symbol, 'NA'),##写出rank
                                        Phrank_rank_intersect,
                                        patho_gene_rank_phen2gene,patho_gene_rank_phen2gene_intersect.get(symbol, 'NA'),
                                        patho_gene_rank_Exo.get(entrezId, 'NA'), patho_gene_rank_LIR.get(entrezGeneId, 'NA')]
            print(case_id)
            print(f"{i} done")

        # except Exception as e:
        #     print(e)
    print(f"final result length is: {len(final_big_table)}")
    final_big_table.to_csv("phenoapt_integ.csv")

def main_2(): ##无精炼HPO版本
    # ## 新加入case跑LIRICAL
    # generate_phenopaket_zs_diag_zs_gVCF()
    # ## 新加入case跑Exomiser
    # generateyml_zs_diag_zs_gVCF()
    # ## 新加入case跑DAGO
    # GADO()
    # ## 新加入case跑phenolyzer，注意关闭terminal
    # phenolyzer()
    filename = 'scoliosis_gVCF_from_zs_updating'
    pwd = '/Users/liyaqi/PycharmProjects/PhenoAptAnalysis'
    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    ##df = df[44:].reset_index(drop=True)
    case_ids = df['Blood ID']
    print(f"{len(df)}")
    # final_big_table = pd.DataFrame(
    #     columns=['CaseID', 'hpo_id', 'Symbol',
    #              'phenoapt_rank', 'intersect_rank', 'Patho_rank_CADD_10', 'Patho_rank_CADD_15', 'Patho_rank_CADD_20',
    #              'phrank_rank','phrank_rank_intersect',
    #              'phenolyzer_rank','phenolyzer_intersect',
    #              'GADO_rank','GADO_rank_intersect'
    #              'phen2gene','phen2gene_intersect',
    #              'Exomiser','LIRICAL'])
    part_table = pd.DataFrame(
        columns=['CaseID', 'hpo_id', 'Symbol',
                 'phrank_rank','phrank_rank_intersect',
                 'phenolyzer_rank','phenolyzer_intersect',
                 'GADO_rank','GADO_intersect_rank',])
    case_id_exomiser_tsv_file_dict = get_case_id_integ_file_map('Exomiser',case_ids)
    case_id_LIRICAL_tsv_file_dict = get_case_id_integ_file_map('LIRICAL', case_ids)
    for i in range(len(df)):
        # try:
            if df.loc[i, 'hpo_id'][0] == '[':
                hpo_id = df.loc[i, 'hpo_id'][2:-2]
                hpo_id_input = [k for k in hpo_id.split("', '")]
            else:
                hpo_id = df.loc[i, 'hpo_id']
                hpo_id_input = [k for k in hpo_id.split(";")]
            symbol = df.loc[i, "Symbol"]
            entrezGeneId = df.loc[i,'entrezGeneId']
            entrezId = df.loc[i,'entrezId']
            ensemblid = df.loc[i,'Ensembl Gene ID']
            print(entrezId)
            case_id = df.loc[i,'Blood ID']
            if get_case_id_file(case_id,'sporadic')!=0:
                filedir = get_case_id_file(case_id,'sporadic')
            else:
                if get_case_id_file(case_id,'trio')!=0:
                    filedir = get_case_id_file(case_id, 'trio')
                else:
                    print('no tsv file')
                    continue
            variation = pd.read_csv(filedir, sep="\t")
            variation_gene_name = variation['Gene_name']

            candidate_gene_list = variation_gene_name[~variation_gene_name.duplicated()]
            # ##需要产生新的phen2gene分析使用的candidate文件时候打开
            # candidate_gene_list.to_csv(f'./{filename}/candidategene/{case_id}.txt',index=False,header=False)
            # ## 需要产生新的phrank分析的ensemblID版本的candidategene时打开
            # output = {}
            # for gene_name_keys in candidate_gene_list.keys():
            #     try:
            #         output[gene_name_keys] = searchensemblid(candidate_gene_list[gene_name_keys])
            #     except Exception as e:
            #         print(e)
            # pd.DataFrame([output[k] for k in output.keys()]).to_csv(f'./{filename}/candidategene_ensemblid/{case_id}.txt',index=False,header=False)

            ## phrank_intersect
            phrank_rank_intersect = phrank_in_pipline(filename, ensemblid, 'intersect', hpo_id_input, case_id)

            ##phrank
            phrank_rank = phrank_in_pipline(filename, ensemblid, 'OMIM',hpo_id_input,case_id)

            phenolyzer_result = pd.read_csv(
                f'{pwd}/{filename}/phenolyzeroutput/{case_id}.final_gene_list',
                sep='\t')
            ## phenolyzer_rank
            if symbol in list(phenolyzer_result['Gene']):
                phenolyzer_rank = list(phenolyzer_result['Gene']).index(symbol)
            else:
                phenolyzer_rank = 'NA'
            ## phenolyzer_intsersect
            phenolyzer_intersect_gene_name = phenolyzer_result[phenolyzer_result.Gene.isin(variation_gene_name)]
            phenolyzer_rank_intersect = {}
            for j, v in enumerate(phenolyzer_intersect_gene_name["Gene"]):
                phenolyzer_rank_intersect[v] = j  ##写出rank

            ## GADO
            ensemblid_candidate_gene = pd.read_csv(f'./{filename}/candidategene_ensemblid/{case_id}.txt',sep = '\t',header= None)
            GADO_result = pd.read_csv(f'{pwd}/{filename}/GADOoutput/{case_id}.txt',sep = '\t')
            if ensemblid in list(GADO_result['Ensg']):
                GADO_rank = list(GADO_result['Ensg']).index(ensemblid)
            else:
                GADO_rank = 'NA'
            ## GADO intersect
            GADO_intersect_gene_id = GADO_result[GADO_result.Ensg.isin(ensemblid_candidate_gene[0])]
            GADO_rank_intersect = {}
            for j, v in enumerate(GADO_intersect_gene_id["Ensg"]):
                GADO_rank_intersect[v] = j
            print(i)

            # # #phen2gene intersect
            # temp = ''
            # for hpo in hpo_id_input:
            #     temp = temp + ' '+hpo
            # HPO_for_Phen2Gene = temp
            # # command_phen2gene_intersect = f'cd /Users/liyaqi/Phen2Gene && python3 phen2gene.py -m{HPO_for_Phen2Gene} -out {pwd}/{filename}/phen2geneoutput -l {pwd}/{filename}/candidategene/{case_id}.txt && mv {pwd}/{filename}/phen2geneoutput/output_file.associated_gene_list {pwd}/{filename}/phen2geneoutput/{case_id}-phen2gene-output.txt'
            # # os.system(command_phen2gene_intersect)
            # phen2gene_intersect_result = pd.read_csv(f'{pwd}/{filename}/phen2geneoutput/{case_id}-phen2gene-output.txt', sep="\t")
            # patho_gene_rank_phen2gene_intersect = {}
            # for j, v in enumerate(phen2gene_intersect_result["Gene"]):
            #     patho_gene_rank_phen2gene_intersect[v] = j

            # # #phen2gene
            # # command_phen2gene = f'cd /Users/liyaqi/Phen2Gene && python3 phen2gene.py -m{HPO_for_Phen2Gene} -out {pwd}/{filename}/phen2geneoutput_nolist && mv {pwd}/{filename}/phen2geneoutput_nolist/output_file.associated_gene_list {pwd}/{filename}/phen2geneoutput_nolist/{case_id}-phen2gene-nolist-output.txt'
            # # os.system(command_phen2gene)
            # phen2gene_result = pd.read_csv(
            #     f'{pwd}/{filename}/phen2geneoutput_nolist/{case_id}-phen2gene-nolist-output.txt', sep="\t")
            # if symbol in list(phen2gene_result["Gene"]):
            #     patho_gene_rank_phen2gene = list(phen2gene_result["Gene"]).index(symbol)
            # else:
            #     patho_gene_rank_phen2gene = 'NA'


            # ## PhenoApt排序
            # weight_1 = [1 for k in range(len(hpo_id_input))]
            # client = PhenoApt(token='H0pVk00CX07VzkZbdnvHI$24XiU$u9q')
            # pheno_result = (client.rank_gene(phenotype=hpo_id_input,weight=weight_1,n=5000)).rank_frame
            # print(hpo_id_input)
            # print(pheno_result)
            # pheno_gene_rank = {}
            # for j, v in enumerate(pheno_result["gene_symbol"]):
            #     pheno_gene_rank[v] = j
            # print(pheno_gene_rank)
            #
            #
            # # #仅取TSV基因的交集排序
            # intersect_gene = pheno_result[pheno_result.gene_symbol.isin(variation_gene_name)]
            # intersect_gene_rank = {}
            # for j, v in enumerate(intersect_gene["gene_symbol"]):
            #     intersect_gene_rank[v] = j ##写出rank
            # # filtered_result = variation[variation.Gene_name.isin(intersect_gene['gene_symbol'])]
            # # filtered_result['pheno_rank'] = [pheno_gene_rank[gene] for gene in filtered_result['Gene_name']]
            # # filtered_result['intersect_rank'] = [intersect_gene_rank[gene] for gene in filtered_result['Gene_name']]
            # # filtered_result['CaseID'] = [case_id for gene in filtered_result['Gene_name']]
            # # filtered_result.to_csv(f"output/{case_id}.csv")
            #
            # variation_p_10 = patho_filter_n(variation, 10)
            # patho_gene_name_10 = variation_p_10['Gene_name']
            # patho_gene_10 = pheno_result[pheno_result.gene_symbol.isin(patho_gene_name_10)]
            # patho_gene_rank_10 = {}
            # for j, v in enumerate(patho_gene_10["gene_symbol"]):
            #     patho_gene_rank_10[v] = j
            #
            # variation_p_15 = patho_filter_n(variation, 15)
            # patho_gene_name_15 = variation_p_15['Gene_name']
            # patho_gene_15 = pheno_result[pheno_result.gene_symbol.isin(patho_gene_name_15)]
            # patho_gene_rank_15 = {}
            # for j, v in enumerate(patho_gene_15["gene_symbol"]):
            #     patho_gene_rank_15[v] = j
            #
            # variation_p_20 = patho_filter_n(variation, 20)
            # patho_gene_name_20 = variation_p_20['Gene_name']
            # patho_gene_20 = pheno_result[pheno_result.gene_symbol.isin(patho_gene_name_20)]
            # patho_gene_rank_20 = {}
            # for j, v in enumerate(patho_gene_20["gene_symbol"]):
            #     patho_gene_rank_20[v] = j
            #
            # patho_gene_rank_Exo = {}
            # if case_id in list(case_id_exomiser_tsv_file_dict):
            #     rankdf = pd.read_csv(case_id_exomiser_tsv_file_dict[case_id], sep="\t")
            #     patho_gene_ID_Exo = rankdf['ENTREZ_GENE_ID']
            #     for j, v in enumerate(patho_gene_ID_Exo):
            #         patho_gene_rank_Exo[v] = j
            #     print(patho_gene_rank_Exo)
            # else:
            #     patho_gene_rank_Exo[entrezId]= 'No_vcf_ranked_by_Exo'
            #
            # patho_gene_rank_LIR = {}
            # if case_id in list(case_id_LIRICAL_tsv_file_dict):
            #     rankdf2 = pd.read_csv(case_id_LIRICAL_tsv_file_dict[case_id], sep="\t")
            #     patho_gene_ID_LIR = rankdf2['entrezGeneId']
            #     patho_gene_rank_LIR = {}
            #     df_count = pd.DataFrame(enumerate(patho_gene_ID_LIR))
            #     for s in range(len(df_count)):
            #         k = (df_count[df_count[1].isin([df_count.loc[s,1]])]).reset_index(drop=True)
            #         print(k)
            #         patho_gene_rank_LIR[k.loc[0,1]] = k.loc[0,0]
            #     print(patho_gene_rank_LIR)
            # else:
            #     patho_gene_rank_LIR[entrezGeneId]= 'No_vcf_ranked_by_LIR'
            #
            # final_big_table.loc[i] = [case_id, hpo_id, symbol,
            #                           pheno_gene_rank.get(symbol, 'NA'), intersect_gene_rank.get(symbol, 'NA'), patho_gene_rank_10.get(symbol, 'NA'), patho_gene_rank_15.get(symbol, 'NA'),patho_gene_rank_20.get(symbol, 'NA'),##写出rank
            #                             Phrank_rank,Phrank_rank_intersect,
            #                             phenolyzer_rank, phenolyzer_intersect.get(symbol, 'NA'),
            #                             GADO_rank, GADO_rank_intersect.get(ensemblid, 'NA'),
            #                             patho_gene_rank_phen2gene,patho_gene_rank_phen2gene_intersect.get(symbol, 'NA'),
            #                             patho_gene_rank_Exo.get(entrezId, 'NA'), patho_gene_rank_LIR.get(entrezGeneId, 'NA')]
            # print(case_id)
            part_table.loc[i] = [case_id, hpo_id, symbol,
                                 phrank_rank,phrank_rank_intersect,
                                phenolyzer_rank, phenolyzer_rank_intersect.get(symbol, 'NA'),
                                 GADO_rank,GADO_rank_intersect.get(ensemblid, 'NA')]
            print(f"{i} done")

    #     except Exception as e:
    #         print(e)
    # print(f"final result length is: {len(final_big_table)}")
    # final_big_table.to_csv("phenoapt_integ.csv")
    part_table.to_csv("6-30-phrank-phenolyzer-GADO.csv")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main_2()

    # See PyCharm help at https://www.jetbrains.com/help/pycharm/

