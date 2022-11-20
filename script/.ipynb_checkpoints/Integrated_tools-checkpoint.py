import os
import appscript as appscript
import pandas as pd
import xlrd
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
import yaml

def getyml_for_certain_ingeritance(case_id, hpo_id_input,pwd, filename, vcfdir, ymldir,outputdir, inheridict, hpo):
    with open(f'{pwd}/yml/test-analysis-exome.yml', 'r') as f:
        config = yaml.safe_load(f)
        config['analysis']['vcf'] = f'{vcfdir}'
        config['analysis']['hpoIds'] = hpo_id_input  # add the command as a list for the correct yaml
        #config['outputOptions']['outputFormats'] = ['TSV_GENE']
        # config['outputOptions']['outputFormats'] = ['HTML','TSV_GENE','TSV_VARIANT']
        config['outputOptions']['outputPrefix'] = f'{pwd}/{filename}/{outputdir}/{case_id}'
        config['analysis']['inheritanceModes'] = inheridict
        print(config['analysis']['hpoIds'])
        print(config['analysis']['inheritanceModes'])
    with open(f'{pwd}/{filename}/{ymldir}/{case_id}_{hpo}.yml', "a") as r:  # open the file in append mode
        r.truncate(0)
        yaml.dump(config, r)
        r.close()##default_flow_style=Faulse

def getyml_for_1310(case_id, hpo_id_input,pwd,filename, vcfdir,ymldir, outputdir, inheridict, hpo, pathogenic = ['REVEL', 'MVP']):
    with open(f'{pwd}/yml/1310-test-analysis-exome.yml', 'r') as f:
        config = yaml.safe_load(f)
        config['analysis']['vcf'] = f'{vcfdir}'
        config['analysis']['hpoIds'] = hpo_id_input  # add the command as a list for the correct yaml
        config['analysis']['pathogenicitySources'] = pathogenic
        #config['outputOptions']['outputFormats'] = ['TSV_GENE']
        # config['outputOptions']['outputFormats'] = ['HTML','TSV_GENE','TSV_VARIANT']
        config['outputOptions']['outputPrefix'] = f'{pwd}/{filename}/{outputdir}/{case_id}'
        config['inheritanceModes'] = inheridict
        print(config['analysis']['hpoIds'])
        print(config['analysis']['inheritanceModes'])
    with open(f'{pwd}/{filename}/{ymldir}/{case_id}_{hpo}.yml', "a") as r:  # open the file in append mode
        r.truncate(0)
        yaml.dump(config, r)
        r.close()##default_flow_style=Faulse

##不用专门写
def get_intersect_yml_for_1310(case_id, hpo_id_input, dir, outputdir, inheridict, hpo, pathogenic = ['REVEL', 'MVP']):
    with open(f'{dir}/../yml/1310-test-analysis-exome.yml', 'r') as f:
        config = yaml.safe_load(f)
        config['analysis']['vcf'] = f'{dir}/intersect_vcf/{case_id}.vcf'
        config['analysis']['hpoIds'] = hpo_id_input  # add the command as a list for the correct yaml
        config['analysis']['pathogenicitySources'] = pathogenic
        #config['outputOptions']['outputFormats'] = ['TSV_GENE']
        # config['outputOptions']['outputFormats'] = ['HTML','TSV_GENE','TSV_VARIANT']
        config['outputOptions']['outputPrefix'] = f'{dir}/{outputdir}/1310-{case_id}'
        config['inheritanceModes'] = inheridict
        print(config['analysis']['hpoIds'])
        print(config['analysis']['inheritanceModes'])
    with open(f'{dir}/intersect_yml/1310-{case_id}_{hpo}.yml', "a") as r:  # open the file in append mode
        r.truncate(0)
        yaml.dump(config, r)
        r.close()##default_flow_style=Faulse


def generateyml_zs_diag_zs_gVCF(df,vcf_file_dict,pwd,filename,hpo,mode='',refresh=False):
    #df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    ##df = df[36:37].reset_index(drop=True)
    print(f"{len(df)}")
    dir_scoliosis_gVCF_from_zs = f"{pwd}/{filename}"
    ##print(sequence_id_file_dict)
    m = 0
    n = 0
    for i in range(len(df)):
        case_id = df.loc[i,'Blood ID']
        vcfdir = vcf_file_dict[case_id]
        #Inheritance_ADAR = df.loc[i,'Inheritance_ADAR']
        if case_id not in vcf_file_dict.keys():
            print(f"no vcf for {case_id}")
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

            outputdir = f'{hpo}_Exomiseroutput'
            outputdir = outputdir + f'{mode}'
            ymldir = f'Exomiseryml' + f'{mode}'
            if ymldir not in os.listdir(f'{pwd}/{filename}'):
                os.system(f'mkdir {pwd}/{filename}/{ymldir}')
            if outputdir not in os.listdir(dir_scoliosis_gVCF_from_zs):
                os.system(f'mkdir {dir_scoliosis_gVCF_from_zs}/{outputdir}')
            if (case_id not in " ".join(os.listdir(f'{pwd}/{filename}/{outputdir}'))) or refresh:
                if (f'{case_id}_{hpo}.yml' not in os.listdir(f'{pwd}/{filename}/{ymldir}')) or refresh:
                    getyml_for_certain_ingeritance(case_id, hpo_id_input, pwd, filename, vcfdir, ymldir, outputdir,
                                                   inheridict, hpo)
                    print('yml done')
                command = f'cd /Users/liyaqi/Downloads/exomiser-cli-12.1.0 && java -Xms2g -Xmx4g -jar exomiser-cli-12.1.0.jar -analysis {pwd}/{filename}/{ymldir}/{case_id}_{hpo}.yml'
                print(command)
                os.system(command)

            print(f'{case_id}, NO.{i}, exomiser done')
            m = m + 1
            print(f'{m} case finish exomiser')
            print(f'{n} files were not in {filename}')

## pahtogenicity version合并
def generateyml1310_zs_diag_zs_gVCF(df,vcf_file_dict, pwd, filename, hpo, mode='',refresh=False):
    print(f"{len(df)}")
    dir_scoliosis_gVCF_from_zs = f"{pwd}/{filename}"
    ##print(sequence_id_file_dict)
    m = 0
    n = 0
    for i in range(len(df)):
        case_id = df.loc[i, 'Blood ID']
        vcfdir = vcf_file_dict[case_id]
        # Inheritance_ADAR = df.loc[i,'Inheritance_ADAR']
        if case_id not in vcf_file_dict.keys():
            print(f"no vcf for {case_id}")
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

            outputdir = f'{hpo}_Exomiser1310output'+mode
            ymldir = f'Exomiser1310yml'+mode
            if outputdir not in os.listdir(dir_scoliosis_gVCF_from_zs):
                os.system(f'mkdir {dir_scoliosis_gVCF_from_zs}/{outputdir}')
            if ymldir not in os.listdir(dir_scoliosis_gVCF_from_zs):
                os.system(f'mkdir {dir_scoliosis_gVCF_from_zs}/{ymldir}')

            if (f'{case_id}' not in " ".join(os.listdir(f'{dir_scoliosis_gVCF_from_zs}/{outputdir}'))) or refresh:
                if (f'{case_id}_{hpo}.yml' not in os.listdir(f'{dir_scoliosis_gVCF_from_zs}/{ymldir}')) or refresh:
                    getyml_for_1310(case_id,  hpo_id_input,pwd,filename, vcfdir, ymldir, outputdir, inheridict, hpo)
                command = f'cd ~/Downloads/exomiser-cli-13.1.0 && java -Xms2g -Xmx4g -jar exomiser-cli-13.1.0.jar -analysis {pwd}/{filename}/{ymldir}/{case_id}_{hpo}.yml'
                print(command)
                os.system(command)

            print(f'{case_id}, NO.{i}, exomiser 13.1.0 done')
            m = m + 1
            print(f'{m} case finish exomiser 13.1.0')
            print(f'{n} files were not in {filename}')


def get_phenopaket(case_id, hpo_id_input,pwd,filename, phenopaketdir, vcfdir,hpo):
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
        params["htsFiles"][0]["uri"] = f'{vcfdir}'
        ##print("params", params)  # 打印
    with open(f'{pwd}/{filename}/{phenopaketdir}/{case_id}_{hpo}.json', 'a') as r:
        r.truncate(0)
        json.dump(params, r)
        ##设置多种VCF输入的包，还没有安置好储存方式！
        # 关闭json写模式
        r.close()

def generate_phenopaket_zs_diag_zs_gVCF(df,vcf_adress_dict,pwd,filename,hpo,mode='',refresh=False):

    print(f"\n.............\n run LIRICAL\n in {len(df)} cases\n loading..........")
    dir_scoliosis_gVCF_from_zs = f"{pwd}/{filename}/"

    m = 0
    n = 0
    for k in range(len(df)):
        case_id = df.loc[k,'Blood ID']
        if vcf_adress_dict[case_id]==0:
            print(f"no vcf of {case_id}")
            n = n + 1
            ##需要搞到vcf
            continue
        else:
            print(f"{case_id} vcf ok")
            if df.loc[k, hpo][0] == '[':
                hpo_id = df.loc[k, hpo][2:-2]
                hpo_id_input = [str(kt) for kt in hpo_id.split("', '")]
            else:
                hpo_id = df.loc[k, hpo]
                hpo_id_input = [kt for kt in hpo_id.split(";")]
            vcfdir = vcf_adress_dict[case_id]
            phenopaketdir = f'LIRICAL_phenopaket'+mode

            m = m + 1
            outputdir = f'{hpo}_LIRICALoutput'+mode
            print(outputdir)
            if outputdir not in os.listdir(f'{pwd}/{filename}'):
                os.system(f'mkdir {pwd}/{filename}/{outputdir} && cp {pwd}/{filename}/LIRICALoutput/LIRICAL.jar {pwd}/{filename}/{outputdir}/LIRICAL.jar && cp {pwd}/{filename}/LIRICALoutput/outputtodf.sh {pwd}/{filename}/{outputdir}/outputtodf.sh')
            if (f'{case_id}.tsv'not in os.listdir(f'{pwd}/{filename}/{outputdir}/')) or refresh:
                if phenopaketdir not in os.listdir(f'{pwd}/{filename}'):
                    os.system(f'mkdir {pwd}/{filename}/{phenopaketdir}')
                if (f'{case_id}_{hpo}.json' not in os.listdir(f'{pwd}/{filename}/{phenopaketdir}')) or refresh:
                    get_phenopaket(case_id, hpo_id_input, pwd, filename, phenopaketdir, vcfdir, hpo)
                    print(f'{case_id}, NO.{k}, phenopaket done')
                command = f'cd {pwd}/{filename}/{outputdir} && java -jar LIRICAL.jar phenopacket -p {pwd}/{filename}/{phenopaketdir}/{case_id}_{hpo}.json -d ~/Downloads/LIRICAL-1.3.4/data -e ~/Downloads/exomiser-cli-12.1.0/data/1909_hg19 -x {case_id} --tsv'
                print(command)
                os.system(command)

            print(f'{m} case finish LIRICAL')
            os.system(f'bash {pwd}/{filename}/{outputdir}/outputtodf.sh')


def generateped(df,pwd,filename,refresh,df_sex):
    df_sex = df_sex.rename(index=df_sex['Blood ID'])
    df_sex = dict(df_sex['Sex'])
    for i in range(len(df)):
        case_id = df.loc[i, 'Blood ID']
        ##做ped
        if '{case_id}.tsv' not in os.listdir(f'{pwd}/{filename}/ped/') or refresh:
            ped_df = (pd.DataFrame(['FAM1', case_id, 'FATHER', 'MOTHER', df_sex.get(case_id), 2]).T)
            ped_df.to_csv(
                f'{pwd}/{filename}/ped/{case_id}.tsv',
                sep='\t', index=False, header=False)
            cmd = f'mv {pwd}/{filename}/ped/{case_id}.tsv {pwd}/{filename}/ped/{case_id}.ped'
            os.system(cmd)
            ##做好了ped for phen-gen
