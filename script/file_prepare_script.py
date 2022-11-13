import requests

from .process_file_script import get_case_id_file
from tqdm import tqdm
import os
import pandas as pd
import json
import subprocess
import sys



def file_prepare_multiple_source(cohort_df, pwd, filename,intersect,REVEL_thresh=[],CADD_thresh=[],refresh=False):
    ##检查是否有.s.vcf.gz 标准VCF
    combine_vcf_svcf_multiple_source(cohort_df, pwd, filename)
    if intersect:
        genertate_intersect_vcf_and_REVEL_Anno_multiple_source_tsv_vcf(cohort_df, pwd, filename, REVEL_thresh,CADD_thresh,refresh=False)
        ##检查整个队列中必要的dropbox到vcf以及revel\cadd注释的.flt.vcf是否已有，没有就生成一个,储存在filename下的文件夹中；可能找不到dropbox tsv或者manual tsv，pirint no tsv

        for i in tqdm(range(len(cohort_df))):
            try:
                case_id = cohort_df.loc[i, 'Blood ID']
                ## 直接intersect的gene list
                ## ensemblid
                annt_vcf_filedir = f'{pwd}/{filename}/intersect_vcf_reca/{case_id}.vcf'
                idfile_from_annot('Gene', pwd, filename, case_id, annt_vcf_filedir,
                                  profile=f'candidategene_ensemblid', refresh=refresh)
                print(f'candidate ensembl txt file ready')

                ## HGNC_ID
                idfile_from_annot('HGNC_ID', pwd, filename, case_id, annt_vcf_filedir,
                                  profile=f'candidategene_hgncid', refresh=refresh)
                print(f'candidate hgncid txt ready')

                ## SYMBOL
                idfile_from_annot('SYMBOL', pwd, filename, case_id, annt_vcf_filedir,
                                  profile=f'candidategene', refresh=refresh)
                print(f'candidate SYMBOL txt ready')

                if len(REVEL_thresh) != 0:
                    for thresh in REVEL_thresh:
                        thresh_filter_filedir = f'{pwd}/{filename}/intersect_vcf_reca/{case_id}_revel_{thresh}.vcf'
                        vcf_filtered_by_reca(case_id,'REVEL',thresh, annt_vcf_filedir, thresh_filter_filedir,refresh)
                        ## 先把filter的结果和没有过滤的注释后intersect_vcf都存在intersect_vcf_revel;CADD_PHRED

                        ##获得不同形式的txt
                        ## ensemblid
                        idfile_from_annot('Gene', pwd, filename,case_id, thresh_filter_filedir,
                                          profile=f'revel_{thresh}_ensemblid', refresh=refresh)
                        print(f'revel_{thresh} ensembl txt file ready')

                        ## HGNC_ID
                        idfile_from_annot('HGNC_ID', pwd, filename, case_id, thresh_filter_filedir,
                                          profile=f'revel_{thresh}_hgncid', refresh=refresh)
                        print(f'revel_{thresh} hgncid txt ready')

                        ## SYMBOL
                        idfile_from_annot('SYMBOL', pwd, filename, case_id, thresh_filter_filedir,
                                          profile=f'revel_{thresh}', refresh=refresh)
                        print(f'revel_{thresh} SYMBOL txt ready')
                        ##检查是否有按照阈值filter过的genename.txt

                if len(CADD_thresh) != 0:
                    for thresh in CADD_thresh:
                        thresh_filter_filedir = f'{pwd}/{filename}/intersect_vcf_reca/{case_id}_cadd_{thresh}.vcf'
                        vcf_filtered_by_reca(case_id, 'CADD_PHRED', thresh, annt_vcf_filedir, thresh_filter_filedir,refresh)
                        ## 先把filter的结果和没有过滤的注释后intersect_vcf都存在intersect_vcf_revel;CADD_PHRED

                        ##获得不同形式的txt
                        ## ensemblid
                        idfile_from_annot('Gene', pwd, filename, case_id, thresh_filter_filedir,
                                          profile=f'cadd_{thresh}_ensemblid', refresh=refresh)
                        print(f'cadd_{thresh} ensembl txt file ready')

                        ## HGNC_ID
                        idfile_from_annot('HGNC_ID', pwd, filename, case_id, thresh_filter_filedir,
                                          profile=f'cadd_{thresh}_hgncid', refresh=refresh)
                        print(f'cadd_{thresh} hgncid txt ready')

                        ## SYMBOL
                        idfile_from_annot('SYMBOL', pwd, filename, case_id, thresh_filter_filedir,
                                          profile=f'cadd_{thresh}', refresh=refresh)
                        print(f'cadd_{thresh} SYMBOL txt ready')
                        ##检查是否有按照阈值filter过的genename.txt

            except Exception as e:
                print(e)


def getvcf_from_zs_gz(pwd,filename):
    ## 需要固定的.gz文件名称，需要批量提取的队列名单由xlsx整理找zs获得,每次刷新list，但是不会refresh已经分离出来的个人vcf；单人vcf则直接下载
    os.system(f'{pwd}/{filename}/bash dividevcf.sh')
    print('vcf unziped from .gz')
    ## 从gz到单个vcf的脚本

def idfile_from_annot(id_name, pwd, filename, case_id,input_vcfdir,profile='candidategene_hgncid',refresh=False):
    ## id_name = ['HGNC_ID','Gene']
    ## profile = ['cadd_20_hgncid','revel_0.25_ensemblid']
    if profile not in os.listdir(f'{pwd}/{filename}/'):
        os.system(f'mkdir {pwd}/{filename}/{profile}')
    if f'{case_id}.txt' not in os.listdir(f'{pwd}/{filename}/{profile}') or refresh:
        cmd1 = f'echo "$SHELL"'+f"&& bcftools +split-vep {input_vcfdir} -f '%{id_name},'" + f" > {pwd}/{filename}/{profile}/{case_id}.txt"
        print(cmd1)
        os.system(cmd1)
        os.system(f'cd {pwd}/{filename}/{profile} && find . -type f -empty -print -delete')
        print('ok')
    ##每次都重新储存了filter后的新ensembleid


def vcf_filtered_by_reca(case_id,SCORE, thresh, vcfdir, thresh_filter_dir, refresh=False):
    ## profile = ['cadd_20_hgncid','revel_0.25_ensemblid']
    outfile = 'wrong SCORE'
    if SCORE == 'REVEL':
        outfile = f'{case_id}_revel_{thresh}.vcf'
    if SCORE == 'CADD_PHRED':
        outfile = f'{case_id}_cadd_{thresh}.vcf'
    if outfile not in os.listdir(thresh_filter_dir[:-(len(outfile)+1)]) or refresh:
        cmd = f"/bin/zsh -c '"+f'source ~/opt/anaconda3/etc/profile.d/conda.sh && conda activate vep && filter_vep -i {vcfdir} --filter "{SCORE} > {thresh} or not {SCORE}" -o {thresh_filter_dir} --force_overwrite'+"'"
        print(cmd)
        os.system(cmd)
        print('filtered vcf ok')




def combine_vcf_svcf_multiple_source(cohort_df, pwd, filename):
    svcf_search_log = pd.DataFrame(
        columns=['case_id', 'vcf_source_mode', 'vcf_file_address'])
    if 'log' not in os.listdir(f'{pwd}/{filename}'):
        os.system(f'mkdir {pwd}/{filename}/log')
    for k in tqdm(range(len(cohort_df))):
        case_id = cohort_df.loc[k,'Blood ID']
        os.system(f'cd {pwd}/{filename}/manually_downloaded_vcf && find . -type f -empty -print -delete')
        multi_svcf_filedir = get_case_id_file(case_id,'manual',pwd,filename,filetype='.s.vcf')
        if multi_svcf_filedir !=0:
            vcf_source_mode = 'manual/.s.vcf'
        else:
            multi_svcf_filedir = get_case_id_file(case_id,'manual',pwd,filename,filetype='.s.vcf.gz')
            if multi_svcf_filedir !=0:
                vcf_source_mode = 'manual/.s.vcf.gz'
                if '.gz' in multi_svcf_filedir:
                    ##解压
                    strlen = len('.gz')
                    output_unzip = multi_svcf_filedir[:-strlen]
                    cmd = f'cp {multi_svcf_filedir} {multi_svcf_filedir}.backup && bgzip -d {multi_svcf_filedir} && mv {multi_svcf_filedir}.backup {multi_svcf_filedir}'
                    os.system(cmd)
                    multi_svcf_filedir = output_unzip
            else:
                multi_svcf_filedir = get_case_id_file(case_id,'vcfgz',pwd,filename,filetype='.vcf')
                if multi_svcf_filedir!=0:
                    vcf_source_mode = 'vcfgz:/vcf/.flt.vcf'
                else:
                    vcf_source_mode = 'no svcf or vcf anywhere'
        svcf_search_log.loc[k,'case_id']=case_id
        svcf_search_log.loc[k,'vcf_source_mode'] = vcf_source_mode
        svcf_search_log.loc[k,'vcf_file_address'] = multi_svcf_filedir
        os.system(f'cd {pwd}/{filename}/manually_downloaded_vcf && find . -type f -empty -print -delete')
        svcf_search_log.to_csv(f'{pwd}/{filename}/log/standard_vcf_log_file.csv', index=True)
    svcf_search_log_s=svcf_search_log[['case_id','vcf_file_address']].set_index('case_id')
    svcf_search_log_s.to_csv(f'{pwd}/{filename}/log/standard_vcf_adreess_file.csv', index=True)

def getvcf_from_dropbox_tsv(pwd,filename,case_id, dir, outputfilename):
    df_tsv = pd.read_csv(dir, sep="\t")  ##dropbox的tsv
    dfx = df_tsv
    for i in tqdm(range(len(df_tsv))):
        ##补充indel的REF和ALT
        if df_tsv.at[i, 'ALT'] == '-':  ##deletion
            ##print(i, 'ALT', df.at[i, 'ALT'])
            pos = int(df_tsv.at[i, 'POS'])
            chr = df_tsv.at[i, 'CHR']
            cmd = f'/bin/zsh -c "source ~/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && samtools faidx ~/Software/ensembl-vep/Homo_sapiens.GRCh37.dna.primary_assembly.fa {chr}:{pos - 1}-{pos - 1}"'
            print(cmd)
            result = subprocess.getoutput(cmd)
            dfx.at[i, 'ALT'] = result.split('\n')[1]
            dfx.at[i, 'REF'] = result.split('\n')[1] + df_tsv.at[i, 'REF']
            dfx.at[i, 'POS'] = pos - 1
        else:
            if df_tsv.at[i, 'REF'] == '-':  ## insertion
                ##print(i, 'REF', df.at[i, 'REF'])
                pos = int(df_tsv.at[i, 'POS'])
                chr = df_tsv.at[i, 'CHR']
                cmd = f'/bin/zsh -c "source ~/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && samtools faidx ~/Software/ensembl-vep/Homo_sapiens.GRCh37.dna.primary_assembly.fa {chr}:{pos}-{pos}"'
                print(cmd)
                result = subprocess.getoutput(cmd)
                dfx.at[i, 'REF'] = result.split('\n')[1]
                dfx.at[i, 'ALT'] = result.split('\n')[1] + df_tsv.at[i, 'ALT']
                dfx.at[i, 'POS'] = pos

    for j in range(len(dfx)):
        if str(dfx.at[j, 'FILTER']) == '-':
            dfx.at[j, 'FILTER'] = '.'
        else:
            if str(dfx.at[j, 'FILTER']) != 'PASS':
                dfx.at[j, 'FILTER'] = 'FailorLowQ'
            continue
            ## 只把-换成.


    QUAL = pd.DataFrame({'QUAL': pd.Series([100 for k in range(len(dfx))])})
    INFO = pd.DataFrame({'INFO': pd.Series(['.' for k in range(len(dfx))])})
    CHR = 'chr'+dfx['CHR']
    ## 假设一些无关紧要的格式信息，完成vcf必须的8列
    # df2 = pd.concat([dfx[['CHR', 'POS', 'Rs_ID', 'REF', 'ALT']], QUAL, dfx['FILTER'], INFO, dfx['FORMAT']], axis=1)
    # df2.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    df2 = pd.concat([CHR, dfx[['POS', 'Rs_ID', 'REF', 'ALT']], QUAL, dfx['FILTER'], INFO], axis=1)
    df2.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    if f'{outputfilename}' not in os.listdir(f'{pwd}/{filename}'):
        os.system(f'mkdir {pwd}/{filename}/{outputfilename}')
    df2.to_csv(f"{pwd}/{filename}/{outputfilename}/{case_id}.txt", sep="\t",
               index=False)  ##转化成vcf信息txt,存到对应数据版本file下
    vcfdir = f"{pwd}/{filename}/{outputfilename}"
    command = f'/bin/zsh -c "source ~/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && awk -f {pwd}/script.awk {vcfdir}/{case_id}.txt | bcftools view |  bcftools sort -o {vcfdir}/{case_id}.vcf"'
    ## awk加一些vcf的表头信息，再bcftools转化完成vcf
    print(command)
    os.system(command)

def dir_of_Filtered_TSV_MUltiple_source(pwd,filename,case_id):
    ##可以从old_pan来，也可以从manual来，获取地址

    filedir = 'no Filtered TSV anywhere'
    source_mode = 'no Filtered TSV anywhere'
    if get_case_id_file(case_id,'manual',pwd,filename)!=0:
        filedir = get_case_id_file(case_id,'manual',pwd,filename)
        source_mode = 'manual_single'
    else:
        print(f'no manual tsv file for {case_id}')
        if get_case_id_file(case_id, 'sporadic',pwd,filename) != 0:
            filedir = get_case_id_file(case_id, 'sporadic',pwd,filename)
            source_mode = 'sporadic'
        else:
            if get_case_id_file(case_id, 'trio',pwd,filename) != 0:
                filedir = get_case_id_file(case_id, 'trio',pwd,filename)
                source_mode = 'trio'
            else:
                print(f'no old tsv file for {case_id}')
    return filedir,source_mode

def dir_of_Filtered_VCF_MUltiple_source(pwd,filename,case_id,filetype='.vcf'):
    ##可以从old_list_vcf.gz来，也可以从manual来，仅获取原始表（还没有与TSV取交集）地址用
    ## 之后与TSV差不多的过程，分为大表解压和手动下载两种
    filedir = 'no Filtered VCF anywhere'
    source_mode = 'no Filtered VCF anywhere'
    if get_case_id_file(case_id,'manual',pwd,filename,filetype=filetype)!=0:
        filedir = get_case_id_file(case_id,'manual',pwd,filename,filetype=filetype)
        source_mode = 'manual'
    else:
        print(f'no manual {filetype} file for {case_id}')
        if get_case_id_file(case_id, 'vcfgz',pwd,filename,filetype=filetype) != 0:
            filedir = get_case_id_file(case_id, 'vcfgz',pwd,filename,filetype=filetype)
            source_mode = 'big_vcf_gz'
        else:
            print(f'no old {filetype} file for {case_id}')
    return filedir,source_mode

def generatevcf_Filtered_TSV_intersect_Filtered_VCF(df,pwd,filename,refresh=False): ##测试记得refresh！！有5个bug case需要这个
    transform_log = pd.DataFrame(columns=['case_id','tsv_source_mode','tsv_file_address','vcf_source_mode','vcf_file_address','bgzip_tabix','intersect'])
    if 'log' not in os.listdir(f'{pwd}/{filename}'):
        os.system(f'mkdir {pwd}/{filename}/log')
    if 'temp' not in os.listdir(f'{pwd}/{filename}'):
        os.system(f'mkdir {pwd}/{filename}/temp')
    if 'intersect_vcf' not in os.listdir(f'{pwd}/{filename}'):
        os.system(f'mkdir {pwd}/{filename}/intersect_vcf')
    #df=df[:1]
    for k in tqdm(range(len(df))):
        case_id = df.loc[k,'Blood ID']
        if f'{case_id}.vcf' not in os.listdir(f'{pwd}/{filename}/intersect_vcf') or refresh:
            tsv_filedir, tsv_source_mode = dir_of_Filtered_TSV_MUltiple_source(pwd, filename, case_id)
            transform_log.loc[k, 'case_id'] = case_id
            transform_log.loc[k, 'tsv_source_mode'] = tsv_source_mode
            transform_log.loc[k, 'intersect'] = 'ok'
            transform_log.loc[k, 'tsv_file_address'] = tsv_filedir
            print(tsv_filedir, tsv_source_mode)
            if f'{case_id}.vcf' not in os.listdir(f'{pwd}/{filename}/dropbox2vcf') or refresh:
                getvcf_from_dropbox_tsv(pwd,filename,case_id, tsv_filedir, 'dropbox2vcf')
            ##如果是有旧的中间文件了，可以设置整个函数refresh
            ##合成好了TSV中variant的vcf，没有info，存在{pwd}/{filename}/dropbox2vcf，multi就是tsv来源包括manual，并且manual优先,没必要区别不同来源，正式试验都是同一来源了
            ##存在log表格文件中
            ##Filtered VCF 地址 获取
            multi_vcf_filedir,vcf_source_mode = dir_of_Filtered_VCF_MUltiple_source(pwd,filename,case_id)
            print(vcf_source_mode)
            if multi_vcf_filedir=='no Filtered VCF anywhere':
                gzdir=get_case_id_file(case_id, 'manual',pwd,filename,filetype='.flt.vcf.gz')
                if gzdir!=0:
                    os.system(f'gzip -d {gzdir}') ##补漏，可能是下载了没有解压，但文件夹里有其他vcf.gz ;此前处理过的文件，不需要爱再解压一遍
                    multi_vcf_filedir, vcf_source_mode = dir_of_Filtered_VCF_MUltiple_source(pwd, filename, case_id)
            transform_log.loc[k, 'source_mode'] = vcf_source_mode
            transform_log.loc[k, 'vcf_file_address'] = multi_vcf_filedir
            print(multi_vcf_filedir,vcf_source_mode)
            ##原地格式转换
            for file in [f'{pwd}/{filename}/dropbox2vcf/{case_id}.vcf',multi_vcf_filedir]:
                cmd0 = f'/bin/zsh -c "source ~/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && bgzip -c -f {file} > {file[:-4]}.vcf.gz && bcftools sort {file[:-4]}.vcf.gz -o {file[:-4]}.vcf.gz -Ou -T {pwd}/{filename}/temp && tabix -p vcf {file[:-4]}.vcf.gz"'
                print(cmd0)
                try:
                    os.system(cmd0)
                except Exception as e:
                    print(e)
                    transform_log.loc[k, 'bgzip_tabix'] = e
            ##取与Filtered VCF产生的交集
            cmd = f'/bin/zsh -c "source ~/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && bcftools isec -n=2 -w1 {multi_vcf_filedir}.gz {pwd}/{filename}/dropbox2vcf/{case_id}.vcf.gz > {pwd}/{filename}/intersect_vcf/{case_id}.vcf.gz"'
            print(cmd)
            try:
                os.system(cmd)
            except Exception as e:
                print(e)
                transform_log.loc[k, 'intersect'] = e
            transform_log.loc[k, 'intersect_adress'] = f'{pwd}/{filename}/intersect_vcf/{case_id}.vcf'
            transform_log.to_csv(f'{pwd}/{filename}/log/intersect_vcf_log_file.csv', index=True)
            os.system(f'cd {pwd}/{filename}/intersect_vcf && bgzip -d {case_id}.vcf.gz')
    os.system(f'cd {pwd}/{filename}/intersect_vcf && find . -type f -empty -print -delete')


        ## 没有手动下载VCF的患者和手动患者的dropbox2vcf中间文件储存路径一致。

def genertatevcf_zs_diag_dropbox_tsv(df,pwd, filename, refresh=False):
    ##df = df[2:].reset_index(drop=True)
    ##做到流程里的REVEL筛选用的VCF（不完整，需要改成有完整信息的版本_就是所谓原始VCF与F_TSV做interscet的版本）准备和REVEL筛选后的VCF产生
    print(f"{len(df)}")
    for i in tqdm(range(len(df))):
        case_id = df.loc[i, 'Blood ID']
        if f'{case_id}.vcf' not in os.listdir(f'{pwd}/{filename}/dropbox2vcf') or refresh:
            if get_case_id_file(case_id, 'sporadic',pwd,filename)!=0:
                filedir= get_case_id_file(case_id, 'sporadic',pwd,filename)
            else:
                if get_case_id_file(case_id, 'trio',pwd,filename)!=0:
                    filedir = get_case_id_file(case_id, 'trio', pwd, filename)
                else:
                    print(f'no {case_id} tsv in old box')
                    continue
            getvcf_from_dropbox_tsv(pwd,filename,case_id, filedir, 'dropbox2vcf')
            ##pwd,case_id, dir, filename
        print('vcf done')
        if f'{case_id}.vcf' not in os.listdir(f'{pwd}/{filename}/dropbox2vcf_revel'):
            vcf_revel(pwd,filename,'dropbox2vcf',
                      'dropbox2vcf_revel', case_id=case_id)
        else:
            if refresh:
                vcf_revel(pwd,filename,'dropbox2vcf',
                          'dropbox2vcf_revel', case_id=case_id)
        print('revel done')


def vcf_revel(pwd,filename,inputdir, outputdir, case_id):
    ## 使用前，先设置文件夹读写权限"chmod a+rwx ./output/*" "chmod a+rwx ./input/*"
    cmd_pre = f'chmod a+rwx {pwd}/{inputdir}/* && chmod a+rwx {pwd}/{outputdir}/*'
    os.system(cmd_pre)
    print("right done")
    cmd = f'/bin/zsh -c "source ~/opt/anaconda3/etc/profile.d/conda.sh && conda activate vep && Vep -i {pwd}/{filename}/{inputdir}/{case_id}.vcf --fork 4 -o {pwd}/{filename}/{outputdir}/{case_id}.vcf --assembly GRCh37 --cache --dir /Users/liyaqi/Software/ensembl-vep --offline --fasta /Users/liyaqi/Software/ensembl-vep/Homo_sapiens.GRCh37.dna.primary_assembly.fa --plugin REVEL,/Users/liyaqi/Software/ensembl-vep/revel/new_tabbed_revel.tsv.gz --force_overwrite"'
    print(cmd)
    os.system(cmd)

def vcf_reca_VCF_FORMAT(pwd,filename,inputdir, outputdir, case_id):
    ## 使用前，先设置文件夹读写权限"chmod a+rwx ./output/*" "chmod a+rwx ./input/*"
    cmd_pre = f'chmod a+rwx {pwd}/{inputdir}/* && chmod a+rwx {pwd}/{outputdir}/*'
    os.system(cmd_pre)
    print("right done")
    cmd = f'/bin/zsh -c "source ~/opt/anaconda3/etc/profile.d/conda.sh && conda activate vep && Vep -i {pwd}/{filename}/{inputdir}/{case_id}.vcf --fork 4 -o {pwd}/{filename}/{outputdir}/{case_id}.vcf --assembly GRCh37 --cache --dir /Users/liyaqi/Software/ensembl-vep --offline --fasta /Users/liyaqi/Software/ensembl-vep/Homo_sapiens.GRCh37.dna.primary_assembly.fa --plugin REVEL,/Users/liyaqi/Software/ensembl-vep/revel/new_tabbed_revel.tsv.gz --plugin CADD,/Users/liyaqi/Downloads/exomiser-cli-13.1.0/data/cadd/whole_genome_SNVs.tsv.gz,/Users/liyaqi/Downloads/exomiser-cli-13.1.0/data/cadd/InDels.tsv.gz --force_overwrite --vcf"'
    print(cmd)
    os.system(cmd)

def genertate_intersect_vcf_and_REVEL_Anno_multiple_source_tsv_vcf(df,pwd, filename, REVEL_thresh,CADD_thresh, refresh=False):
    ##df = df[2:].reset_index(drop=True)
    ##做到新流程里的REVEL筛选用的VCF（不完整，需要改成有完整信息的版本_就是所谓原始VCF与F_TSV做interscet的版本）准备和REVEL筛选后的VCF产生

    print(f"{len(df)} case running inntersect vcf and revel annotation")
    generatevcf_Filtered_TSV_intersect_Filtered_VCF(df, pwd, filename, refresh=refresh)
    ##pwd,case_id, dir, filename
    print('intersect vcf done')

    ## 使用REVEL注释VCF、CADD，
    if len(REVEL_thresh) != 0 or len(CADD_thresh) != 0:
        if 'intersect_vcf_reca' not in os.listdir(f'{pwd}/{filename}'):
            os.system(f'mkdir {pwd}/{filename}/intersect_vcf_reca')
        for i in tqdm(range(len(df))):
            case_id = df.loc[i, 'Blood ID']
            if f'{case_id}.vcf' not in os.listdir(f'{pwd}/{filename}/intersect_vcf_reca') or refresh:
                vcf_reca_VCF_FORMAT(pwd,filename,'intersect_vcf','intersect_vcf_reca', case_id)
            print('revel cadd VCF FORMAT done')
