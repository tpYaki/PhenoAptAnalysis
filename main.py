# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import os

from ./phenoapt import PhenoApt
client = PhenoApt(token='H0pVk00CX07VzkZbdnvHI$24XiU$u9q')

import xlrd
import pandas as pd
import subprocess
import re
import csv

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.


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

def main():
    df, case_ids = read_diagnose_xlsx('/Users/liyaqi/Documents/生信/Inhouse_cohorts_genes_Version_8_MRR_诊断.xlsx')
    print(f"{len(df)}")
    df = df[:5]
    case_id_tsv_file_dict = get_case_id_file_map(case_ids)
    final_big_table = pd.DataFrame(columns=['CaseID', 'hpo_id', 'phenoapt_rank', 'intersect_rank', 'Symbol'])
    for i in range(len(df)):
        try:
            case_id = df.loc[i, 'CaseID']
            hpo_id = df.loc[i, 'hpo_id']
            symbol = df.loc[i, "Symbol"]
            weight = ','.join([str(1) for i in range(len(hpo_id.split(";")))])
            command = f"phenoapt rank-gene -p '{hpo_id}' -w '{weight}' -n 1000"
            print(command)
            result = subprocess.getoutput(command)
            lines = result.split("\n")
            columns = re.split('\s+', lines[0].strip())
            rows = [re.split("\s+", line.strip()) for line in lines[2:]]
            pheno_result = pd.DataFrame(columns=columns, data=rows)
            pheno_gene_rank = {}
            for j, v in enumerate(pheno_result["gene_symbol"]):
                pheno_gene_rank[v] = j
            print(pheno_result)

            if case_id not in case_id_tsv_file_dict:
                continue
            variation = pd.read_csv(case_id_tsv_file_dict[case_id], sep="\t")
            variation_gene_name = variation['Gene_name']

            #求基因的交集
            intersect_gene = pheno_result[pheno_result.gene_symbol.isin(variation_gene_name)]
            intersect_gene_rank = {}
            for j, v in enumerate(intersect_gene["gene_symbol"]):
                intersect_gene_rank[v] = j
            filtered_result = variation[variation.Gene_name.isin(intersect_gene['gene_symbol'])]
            filtered_result['pheno_rank'] = [pheno_gene_rank[gene] for gene in filtered_result['Gene_name']]
            filtered_result['intersect_rank'] = [intersect_gene_rank[gene] for gene in filtered_result['Gene_name']]
            filtered_result['CaseID'] = [case_id for gene in filtered_result['Gene_name']]
            filtered_result.to_csv(f"output/{case_id}.csv")
            #final_big_table = pd.DataFrame(columns=['CaseID', 'hpo_id', 'phenoapt_rank', 'intersect_rank', 'Symbol'])
            final_big_table.loc[i] = [case_id, hpo_id, pheno_gene_rank.get(symbol, 'NA'),
                                      intersect_gene_rank.get(symbol, 'NA'), symbol]

            print(lines)
            print(case_id)
            print(f"{i} done")
        except Exception as e:
            print(e)
    print(f"final result length is: {len(final_big_table)}")
    final_big_table.to_csv("ground_truth_symbol.csv")



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
