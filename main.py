# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import os

import pandas as pd
import xlrd
from phenoapt import PhenoApt
import numpy as np

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

def getvcf(case_id,dir):
    # command = f"conda activate bcftools \ awk -f script.awk {dir}  |bcftools view -o /Users/liyaqi/PycharmProjects/PhenoAptAnalysis/VCF/{case_id}.vcf \ conda deactivate"
    df = pd.read_csv(dir, sep="\t")
    # df2 = pd.dataframe(columns=["CHROM", "POS", "Rs_ID", "REF", "ALT","QUAL","FILTER","INFO"])
    dfx = pd.DataFrame()
    for i in range(len(df)):
        if df.at[i, 'ALT'] == '-':
            print(i, 'ALT', df.at[i, 'ALT'])
        else:
            if df.at[i, 'REF'] == '-':
                print(i, 'REF', df.at[i, 'REF'])
            else:
                dfx=pd.concat([dfx, df.loc[i]],axis='columns') ##删掉indel
    dfx = (dfx.T).reset_index(drop=True)
    for j in range(len(dfx)):
        if str(dfx.at[j,'FILTER']) != '-':
            continue
        dfx.at[j,'FILTER'] = '.'

    QUAL = pd.DataFrame({'QUAL':pd.Series(['.' for k in range(len(dfx))])})
    INFO = pd.DataFrame({'INFO': pd.Series(['.' for k in range(len(dfx))])})
    pwd = "/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/"
    df2 = pd.concat([dfx[['CHR', 'POS', 'Rs_ID', 'REF', 'ALT']],QUAL,dfx['FILTER'],INFO],axis=1)
    df2.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    df2.to_csv(f"./VCF/{case_id}.txt", sep="\t", index=False)
    command = f"source /Users/liyaqi/opt/anaconda3/etc/profile.d/conda.sh && conda activate bcftools && awk -f {pwd}script.awk {pwd}VCF/{case_id}.txt | bcftools view -o {pwd}VCF/{case_id}.vcf"
    print(command)
    os.system(command)
    print(f" done")
    return True

def generatevcf():
    df, case_ids = read_diagnose_xlsx('/Users/liyaqi/Documents/生信/Inhouse_cohorts_genes_Version_8_MRR_诊断.xlsx')
    print(f"{len(df)}")
    df = df[:1]
    case_id_tsv_file_dict = get_case_id_file_map(case_ids)
    for i in range(len(df)):
        try:
            case_id = df.loc[i, 'CaseID']
            if case_id not in case_id_tsv_file_dict:
                continue
            getvcf(case_id, case_id_tsv_file_dict[case_id])
        except Exception as e:
            print(e)

import yaml
def getyml(case_id_file, hpo_id_input):
    pwd = "/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/"
    with open('./yml/test-analysis-exome.yml','r') as f:
        config = yaml.safe_load(f)
        config['analysis']['vcf'] = f'{pwd}VCF/{case_id_file}.vcf'
        config['analysis']['hpoIds'] = hpo_id_input # add the command as a list for the correct yaml
        config['outputOptions']['outputFormats']=['TSV_GENE']
        config['outputOptions']['outputPrefix'] = f'{pwd}Exomiseroutput/{case_id_file}'
        ## del config['inheritanceModes']   # del the 'hostname' key from config
        print(config['analysis']['hpoIds'])

    with open(f'./yml/{case_id_file}.yml',"w") as f: # open the file in append mode
        f.truncate(0)
        yaml.dump(config, f) ##default_flow_style=Faulse

def generateyml():
    df, case_ids = read_diagnose_xlsx('/Users/liyaqi/Documents/生信/Inhouse_cohorts_genes_Version_8_MRR_诊断.xlsx')
    print(f"{len(df)}")
    df = df[:1]
    case_id_tsv_file_dict = get_case_id_file_map(case_ids)
    for i in range(len(df)):
        try:
            case_id = df.loc[i, 'CaseID']
            if case_id not in case_id_tsv_file_dict:
                print(f"{i} not in dropbox")
                continue

            hpo_id = df.loc[i, 'hpo_id']
            hpo_id_input = [k for k in hpo_id.split(";")]

            getyml(case_id,hpo_id_input) ##用map到的所有文件名case_id_file称生成yml,与main中的case_id一个意思
        except Exception as e:
            print(e)

def main():
    df, case_ids = read_diagnose_xlsx('/Users/liyaqi/Documents/生信/Inhouse_cohorts_genes_Version_8_MRR_诊断.xlsx')
    print(f"{len(df)}")
    ##df = df[:1]
    case_id_tsv_file_dict = get_case_id_file_map(case_ids)
    final_big_table = pd.DataFrame(
        columns=['CaseID', 'hpo_id', 'Symbol', 'phenoapt_rank', 'intersect_rank', 'Patho_rank_CADD_10',
                 'Patho_rank_CADD_15', 'Patho_rank_CADD_20'])

    for i in range(len(df)):
        try:
            case_id = df.loc[i, 'CaseID']
            if case_id not in case_id_tsv_file_dict:
                print(f"{i} not in dropbox")
                continue

            hpo_id = df.loc[i, 'hpo_id']
            hpo_id_input = [k for k in hpo_id.split(";")]
            symbol = df.loc[i, "Symbol"]

            ## PhenoApt排序
            weight_1 = [1 for k in range(len(hpo_id.split(";")))]
            client = PhenoApt(token='H0pVk00CX07VzkZbdnvHI$24XiU$u9q')
            pheno_result = (client.rank_gene(phenotype=hpo_id_input,weight=weight_1,n=5000)).rank_frame
            print(hpo_id_input)
            print(pheno_result)
            pheno_gene_rank = {}
            for j, v in enumerate(pheno_result["gene_symbol"]):
                pheno_gene_rank[v] = j
            print(pheno_gene_rank)

            variation = pd.read_csv(case_id_tsv_file_dict[case_id], sep="\t")
            variation_gene_name = variation['Gene_name']

            #仅取TSV基因的交集排序
            intersect_gene = pheno_result[pheno_result.gene_symbol.isin(variation_gene_name)]
            intersect_gene_rank = {}
            for j, v in enumerate(intersect_gene["gene_symbol"]):
                intersect_gene_rank[v] = j ##写出rank
            filtered_result = variation[variation.Gene_name.isin(intersect_gene['gene_symbol'])]
            filtered_result['pheno_rank'] = [pheno_gene_rank[gene] for gene in filtered_result['Gene_name']]
            filtered_result['intersect_rank'] = [intersect_gene_rank[gene] for gene in filtered_result['Gene_name']]
            filtered_result['CaseID'] = [case_id for gene in filtered_result['Gene_name']]
            filtered_result.to_csv(f"output/{case_id}.csv")

            # #求致病性突变过滤后的基因排序
            # for w in (10,15,20)
            #     variation_p_w = patho_filter_n(variation,w)
            #     patho_gene_name_w = variation_p_w['Gene_name']
            #     patho_gene_w = pheno_result[pheno_result.gene_symbol.isin(patho_gene_name_w)]
            #     patho_gene_rank_w = {}
            #     for j, v in enumerate(patho_gene_w["gene_symbol"]):
            #         patho_gene_rank_w[v] = j
            #     dict(name=f"patho_gene_rank_{w}")=patho_gene_rank_w
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
            final_big_table.loc[i] = [case_id, hpo_id, symbol, pheno_gene_rank.get(symbol, 'NA'),
                                      intersect_gene_rank.get(symbol, 'NA'), patho_gene_rank_10.get(symbol, 'NA'), patho_gene_rank_15.get(symbol, 'NA'),patho_gene_rank_20.get(symbol, 'NA')] ##写出rank


            #final_big_table = pd.DataFrame(columns=['CaseID', 'hpo_id', 'phenoapt_rank', 'intersect_rank', 'Symbol'])

            #print(lines)
            print(case_id)
            print(f"{i} done")
        except Exception as e:
            print(e)
    print(f"final result length is: {len(final_big_table)}")
    final_big_table.to_csv("strategy _compare_3.csv")



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    generateyml()


    # See PyCharm help at https://www.jetbrains.com/help/pycharm/
