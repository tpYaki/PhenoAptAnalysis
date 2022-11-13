# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


from script.process_file_script import *
from script.functions import *
from script.file_prepare_script import *
from script.Integrated_tools import *

def main_2(tools=[], CADD_thresh=[], REVEL_thresh=[], hpo='hpo_id', intersect=False, statistic = True, refresh=False, newcase_runtool=False, newvcfgz=False, collect_variants_tofile = False, collect_organ = True):

    filename = 'scoliosis_gVCF_from_zs_updating'
    pwd = '/Users/liyaqi/PycharmProjects/PhenoAptAnalysis'

    df_organ = read_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx',' organ_system')
    df_sex = read_xlsx('/Users/liyaqi/Documents/生信/gvcf-scoliosis-2022.xlsx', 'case')
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

        file_prepare_multiple_source(df, pwd, filename, intersect,REVEL_thresh=REVEL_thresh, CADD_thresh=CADD_thresh, refresh=refresh)
        ## 准备.s.vcf 或者 .flt.vcf 如果intersect 准备 vcf和转化txt，如果score filter，准备各种预置txt
        svcf_vcf_adress_df = pd.read_csv(f'{pwd}/{filename}/log/standard_vcf_adreess_file.csv').set_index('case_id')
        svcf_vcf_adress_dict = svcf_vcf_adress_df['vcf_file_address'].to_dict()
        if intersect:
            tsv_fil_vcf_adress_df = pd.read_csv(f'{pwd}/{filename}/log/intersect_vcf_log_file.csv').set_index('case_id')
            tsv_fil_vcf_adress_dict = tsv_fil_vcf_adress_df['intersect_adress'].to_dict()

        REVEL_filter_vcf_array = []
        if len(REVEL_thresh) != 0:
            for thresh in REVEL_thresh:
                REVEL_filter_vcf_dict = [f'{pwd}/{filename}/intersect_vcf_reca/{case_id}_revel_{thresh}.vcf' for case_id
                                         in case_ids]
                REVEL_filter_vcf_array.append(REVEL_filter_vcf_dict)
            REVEL_filter_vcf_df = pd.DataFrame(np.array(REVEL_filter_vcf_array).T, index=case_ids, columns=REVEL_thresh)

        CADD_filter_vcf_array = []
        if len(CADD_thresh) != 0:
            for thresh in CADD_thresh:
                CADD_filter_vcf_dict = [f'{pwd}/{filename}/intersect_vcf_reca/{case_id}_cadd_{thresh}.vcf' for case_id
                                         in case_ids]
                CADD_filter_vcf_array.append(CADD_filter_vcf_dict)
            CADD_filter_vcf_df = pd.DataFrame(np.array(CADD_filter_vcf_array).T, index=case_ids, columns=CADD_thresh)

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
        if 'LIRICAL' in tools: ## ready
            generate_phenopaket_zs_diag_zs_gVCF(df,svcf_vcf_adress_dict,pwd,filename,hpo,refresh=refresh)
            if intersect:
                generate_phenopaket_zs_diag_zs_gVCF(df,tsv_fil_vcf_adress_dict,pwd,filename,hpo,'intersect',refresh=refresh)
                if len(CADD_thresh) != 0:
                    for thresh in CADD_thresh:
                        generate_phenopaket_zs_diag_zs_gVCF(df, CADD_filter_vcf_df[thresh].to_dict(), pwd, filename, hpo, f'cadd_{thresh}',
                                           refresh=refresh)
                if len(REVEL_thresh) != 0:
                    for thresh in REVEL_thresh:
                        generate_phenopaket_zs_diag_zs_gVCF(df, REVEL_filter_vcf_df[thresh].to_dict(), pwd, filename, hpo, f'revel_{thresh}',
                                           refresh=refresh)
            ## 新的hpo，新的case
            print('LIRICAL ready')
        if 'Exomiser' in tools: #ready now!
            generateyml_zs_diag_zs_gVCF(df,svcf_vcf_adress_dict,pwd,filename,hpo,refresh=refresh)
            if intersect:
                generateyml_zs_diag_zs_gVCF(df,tsv_fil_vcf_adress_dict,pwd,filename,hpo,'intersect',refresh=refresh)
                if len(CADD_thresh) != 0:
                    for thresh in CADD_thresh:
                        generateyml_zs_diag_zs_gVCF(df, CADD_filter_vcf_df[thresh].to_dict(), pwd, filename, hpo, f'cadd_{thresh}',
                                           refresh=refresh)
                if len(REVEL_thresh) != 0:
                    for thresh in REVEL_thresh:
                        generateyml_zs_diag_zs_gVCF(df, REVEL_filter_vcf_df[thresh].to_dict(), pwd, filename, hpo, f'revel_{thresh}',
                                           refresh=refresh)

            print('Exomiser ready')
            ## 新的hpo，新的case
        if 'Exomiser13.1.0' in tools:
            generateyml1310_zs_diag_zs_gVCF(df,svcf_vcf_adress_dict, pwd, filename, hpo,refresh=refresh)
            if intersect:
                generateyml1310_zs_diag_zs_gVCF(df, tsv_fil_vcf_adress_dict, pwd, filename, hpo, 'intersect', refresh=refresh)
                if len(CADD_thresh) != 0:
                    for thresh in CADD_thresh:
                        generateyml1310_zs_diag_zs_gVCF(df, CADD_filter_vcf_df[thresh].to_dict(), pwd, filename, hpo, f'cadd_{thresh}',
                                           refresh=refresh)
                if len(REVEL_thresh) != 0:
                    for thresh in REVEL_thresh:
                        generateyml1310_zs_diag_zs_gVCF(df, REVEL_filter_vcf_df[thresh].to_dict(), pwd, filename, hpo, f'revel_{thresh}',
                                           refresh=refresh)

            print('Exomiser 13.1.0 ready')
        if 'PhenoGenius' in tools:
            Phenogenius_rank(df, pwd, filename, hpo, refresh)
            print('PhenoGenius ready')
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
        statistic_to_table(df,pwd,filename,hpo,tools,intersect,REVEL_thresh,CADD_thresh,collect_variants_tofile,collect_organ)

    if collect_variants_tofile:
        collect_variants(pwd, df, filename)

    Rscript_brief_benchmark(tools, intersect, REVEL_thresh, CADD_thresh)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    filename = 'scoliosis_gVCF_from_zs_updating'
    pwd = '/Users/liyaqi/PycharmProjects/PhenoAptAnalysis'
    df = read_gvcf_diagnosis_xlsx('/Users/liyaqi/Documents/生信/gVCF-2022-确诊.xlsx')
    tools = ['phenoapt', 'phenomizer', 'phenolyzer', 'GADO', 'PhenoRank', 'phen2gene', 'phrank', 'HANRD', 'LIRICAL',
                      'Exomiser','Exomiser13.1.0','PhenoGenius','phen-gen']

    #compare_df_for_R(df, tools, pwd, 'wilcoxon', dMRR_heatmap=False)
    # file_prepare_multiple_source(df, pwd, filename, intersect=True, REVEL_thresh=[0.25,0.5,0.75],CADD_thresh=[10,15,20], refresh=False)
    ##各种策略 for integrated tools
    main_2(tools=tools,REVEL_thresh=[0.25,0.5,0.75],CADD_thresh=[10,15,20], hpo='hpo_id', intersect=True, statistic = True, refresh=False, newcase_runtool=False, newvcfgz=False, collect_variants_tofile = False, collect_organ = False)

    ##统计
    # statistic_to_table(df, pwd, filename, 'hpo_id', tools=['phenomizer'], intersect=True, REVEL_thresh=[0.25,0.5,0.75], CADD_thresh=[10,15,20],collect_variants_tofile=False,collect_organ=False)
    #'phenoapt', 'GADO', 'phrank', 'phenolyzer'
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

    # pubcasefinder(df, pwd, filename, 'hpo_id', refresh=True)

    ## newcase to run SNV and CNV
    # main_2(tools=['phrank','phenolyzer','Exomiser13.1.0','GADO','PhenoRank','phen2gene','phenoapt','LIRICAL','HANRD','Exomiser','Phen_gen','phenomizer'], hpo='Weight', intersect=True, refresh=False, newcase_runtool=False, statistic=True, collect_organ=True)
    # cnv_FA(pwd,filename,df,True)
    # See PyCharm help at https://www.jetbrains.com/help/pycharm/

