import xlrd
import pandas as pd
import os
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
def get_case_id_file_map(case_ids,tsv_dir='/Users/liyaqi/Dropbox/Filtered-SNV-all/Sporadic-WES-All'):
    case_id_tsv_file_map = {}
    files = os.listdir(tsv_dir)
    for case_id in case_ids:
        for file_name in files:
            if case_id in file_name:
                case_id_tsv_file_map[case_id] = os.path.join(tsv_dir, file_name)
    return case_id_tsv_file_map


def get_case_id_file(case_id, mode,pwd,filename,old_dir='/Users/liyaqi/Dropbox/Filtered-SNV-all',manual_dir = '/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/scoliosis_gVCF_from_zs_updating/manually_downloaded_vcf',filetype='.flt.tsv'):
    if mode in ['sporadic','trio']:
        tsv_dir = os.path.join(old_dir,f'{mode}-WES-All')
        files = os.listdir(tsv_dir)
        for file_name in files:
            if case_id in file_name:
                case_id_tsv_file = os.path.join(tsv_dir, file_name)
                return case_id_tsv_file
        return 0
    if mode == 'manual':
        files = os.listdir(manual_dir)
        for file_name in files:
            if case_id in file_name:
                strlen = len(f'{filetype}')
                if f'{filetype}' == file_name[-strlen:]:
                    case_id_tsv_file = os.path.join(manual_dir, file_name)
                    return case_id_tsv_file
        return 0
    if mode == 'vcfgz':
        vcfgz_dir = f'{pwd}/{filename}/vcf'
        files = os.listdir(vcfgz_dir)
        for file_name in files:
            if case_id in file_name:
                strlen = len(f'{filetype}')
                if f'{filetype}' == file_name[-strlen:]:
                    case_id_tsv_file = os.path.join(vcfgz_dir, file_name)
                    return case_id_tsv_file
        return 0
    ##获取一些原始数据manual的vcf和tsv，vcfgz解压来的vcf,old tsv

