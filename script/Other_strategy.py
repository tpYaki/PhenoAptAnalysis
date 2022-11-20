from script.functions import *

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
    dir = f'../scoliosis_filtered_tsv_dropbox/'
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

def tsv_filtered_by_revel(df_input, case_id, pwd, filename, thresh,vcfdir):
    ## df_imput 是读出来的dropbox tsv,根据revel注释的vcf文件，生成{case_id}_revel_filter_{REVEL_Thresh}.tsv，也就是可以读取成df的文件，并且输出df_input中符合revel要求的部分
    dir = f'{pwd}/{filename}/{vcfdir}'
    if f'{case_id}_revel_filter_{thresh}.tsv' not in os.listdir(dir):
        cmd = f'/bin/zsh -c "source ~/opt/anaconda3/etc/profile.d/conda.sh && conda activate vep && filter_vep -i {dir}/{case_id}.vcf --filter "REVEL > {thresh} or not REVEL" -o {dir}/{case_id}_revel_filter_{thresh}.tsv --force_overwrite && bash {dir}/delet2df.sh"'
        os.system(cmd)
    ##使用原来大范围vcf输出filtered vcf格式的结果，替代tsv vep格式的输出之后，没有解决与tsv intersect的方式问题，vcf没有拼出来variant_list的方式。
    vcf_revel_tsv = pd.read_csv(f'{dir}/{case_id}_revel_filter_{thresh}.tsv', sep='\t')
    vcf_revel_tsv['variant_list'] = vcf_revel_tsv['Location'].map(str) + vcf_revel_tsv['Allele']
    for i in range(len(df_input)): ##是整个Filtered tsv，草
        if df_input.loc[i, 'ALT'] == '-':
            if len(df_input.loc[i, 'REF']) == 1:
                pos = str(df_input.loc[i, 'POS'])
            else:
                pos = str(df_input.loc[i, 'POS']) + '-' + str(int(df_input.loc[i, 'POS']) + len(df_input.loc[i, 'REF']) - 1)
        else:
            if df_input.loc[i, 'REF'] == '-':
                pos = str(int(df_input.loc[i, 'POS']) - 1) + '-' + str(df_input.loc[i, 'POS'])
            else:
                pos = str(df_input.loc[i, 'POS'])
        df_input.loc[i, 'intersect_revel_list'] = str(df_input.loc[i, 'CHR']) + ':' + pos + str(df_input.loc[i, 'ALT'])
    df_result = df_input[df_input['intersect_revel_list'].isin(vcf_revel_tsv['variant_list'])]
    return df_result


def file_prepare(cohort_df, pwd, filename,REVEL_thresh=[],CADD_thresh=[],refresh=False,intersect=False):
    if len(REVEL_thresh) != 0:
        # genertate_intersect_vcf_and_REVEL_Anno_multiple_source_tsv_vcf(cohort_df, pwd, filename, refresh=False)
        genertatevcf_zs_diag_dropbox_tsv(cohort_df,pwd=pwd, filename=filename, refresh=refresh)
        ##加入流程模式开关？
        ##检查整个队列中必要的dropbox到vcf以及revel注释的vcf是否已有，没有就生成一个,储存在filename下的文件夹中；可能找不到dropbox tsv或者manual tsv，pirint no tsv
    for i in tqdm(range(len(cohort_df))):
        try:
            case_id = cohort_df.loc[i, 'Blood ID']
            if get_case_id_file(case_id, 'sporadic') != 0:
                filedir = get_case_id_file(case_id, 'sporadic')
            else:
                if get_case_id_file(case_id, 'trio') != 0:
                    filedir = get_case_id_file(case_id, 'trio')
                else:
                    print(f'no dropbox tsv file for {case_id}')
                    continue
            variation = pd.read_csv(filedir, sep="\t")
            variation_gene_name = variation['Gene_name']
            variation_gene_list = list(set(variation_gene_name))

            if intersect:
                # #检查是否需要产生新的phen2gene分析使用的candidate文件
                candidate_file(df2list=variation_gene_list, pwd=pwd, filename=filename, profile='candidategene', case_id=case_id)
                ##检查candidate gene 的 ensenmble file是否ready
                print('candidategene,ready')
                ensemblidfile(genelist=variation_gene_list, pwd=pwd, filename=filename, case_id=case_id,
                               profile='candidategene_ensemblid', refresh=refresh)
                print('candidategene_ensemblid,ready')
                ##检查candidate gene 的 hgnc file是否ready
                hgncidfile(genelist=variation_gene_list, pwd=pwd, filename=filename, case_id=case_id,
                              profile='candidategene_hgncid', refresh=refresh)
                print('candidategene_hgncid,ready')


            if len(REVEL_thresh) != 0:
                for thresh in REVEL_thresh:
                    variation_revel = tsv_filtered_by_revel(df_input=variation, case_id=case_id,pwd=pwd, filename=filename,
                                                            thresh=thresh,vcfdir='dropbox2vcf_revel')
                    ##print(variation_revel[:1])
                    ## df_input, case_id, pwd, filename, thresh
                    variation_revel_gene_list = list(set(variation_revel['Gene_name']))
                    candidate_file(variation_revel_gene_list, pwd, filename, case_id, f'revel_{thresh}')
                    print(f'revel_{thresh} ready')
                    ensemblidfile(variation_revel_gene_list, pwd, case_id, filename,
                                   f'revel_{thresh}_ensemblid', refresh=refresh)
                    print(f'revel_{thresh}_ensemblid ready')
                    hgncidfile(variation_revel_gene_list, pwd, case_id, filename, f'revel_{thresh}_hgncid',
                               refresh=refresh)
                    print(f'revel_{thresh}_hgncid ready')
                    ##检查是否有按照阈值filter过的genename.txt

            if len(CADD_thresh) != 0:
                for thresh in CADD_thresh:
                    variation_cadd = patho_filter_n(variation, thresh)
                    variation_cadd_gene_list = list(set(variation_cadd['Gene_name']))
                    candidate_file(variation_cadd_gene_list, pwd, filename, case_id, f'cadd_{thresh}')
                    ##df2list, pwd, filename, case_id, profile, refresh
                    print(f'cadd_{thresh} ready')
                    ensemblidfile(variation_cadd_gene_list, pwd, case_id, filename, f'cadd_{thresh}_ensemblid',
                                   refresh=refresh)
                    print(f'cadd_{thresh}_ensemblid ready')
                    hgncidfile(variation_cadd_gene_list, pwd, case_id, filename, f'cadd_{thresh}_hgncid',
                                  refresh=refresh)
                    ##检查是否有按照阈值filter过的genename.txt,refre指的是是否重新查询ensembleid，filterfile是每次都在重新生成
                    ## genelist, pwd, case_id, filename='scoliosis_gVCF_from_zs_updating', profile='candidategene_ensemblid',refresh=False
                    print(f'cadd_{thresh}_hgncid ready')

        except Exception as e:
            print(e)


def candidate_file(df2list, pwd, filename, case_id, profile):
    if profile not in os.listdir(f'{pwd}/{filename}/'):
        os.system(f'mkdir {pwd}/{filename}/{profile}')
    if 'candidategene' not in os.listdir(f'{pwd}/{filename}/'):
        os.system('mkdir candidategene')
    candidate_gene_list_no_dup = list(set(df2list))
    candidate_gene_list = pd.DataFrame(columns=[0])
    candidate_gene_list[0] = candidate_gene_list_no_dup
    candidate_gene_list.to_csv(f'{pwd}/{filename}/{profile}/{case_id}.txt', index=False, header=False)
    ##每次都存新的gene list


def hgncidfile(genelist, pwd, case_id, filename='scoliosis_gVCF_from_zs_updating', profile='candidategene_hgncid',refresh=False):
    if profile not in os.listdir(f'{pwd}/{filename}/'):
        os.system(f'mkdir {pwd}/{filename}/{profile}')
    if 'candidategene_hgncid' not in os.listdir(f'{pwd}/{filename}/'):
        os.system(f'mkdir {pwd}/{filename}/candidategene_hgncid')
    if f'{case_id}_index.tsv' not in os.listdir(f'./{filename}/candidategene_hgncid/'):
        output = pd.DataFrame(columns=[0])
        for gene_name_keys in genelist:
            try:
                id_result = searchHGNCid(gene_name_keys)
                output.loc[gene_name_keys] = id_result
                print(id_result)
            except Exception as e:
                print(e)
        output.to_csv(f'./{filename}/candidategene_hgncid/{case_id}.txt', index=False, header=False)
        output.to_csv(f'./{filename}/candidategene_hgncid/{case_id}_index.tsv', sep='\t', index=True, header=False)
    else:
        if refresh:
            output = pd.DataFrame(columns=[0])
            for gene_name_keys in genelist:
                try:
                    id_result = searchHGNCid(gene_name_keys)
                    output.loc[gene_name_keys] = id_result
                except Exception as e:
                    print(e)
            output.to_csv(f'./{filename}/candidategene_hgncid/{case_id}.txt', index=False, header=False)
            output.to_csv(f'./{filename}/candidategene_hgncid/{case_id}_index.tsv', sep='\t', index=True,
                          header=False)
        if profile != 'candidategene_hgncid':
            hgncids_df = pd.read_csv(f'./{filename}/candidategene_hgncid/{case_id}_index.tsv', sep='\t', header=None)
            hgncids_df = hgncids_df.reset_index()
            output_filter = hgncids_df[hgncids_df[0].isin(genelist)][1]
            ##print(output_filter[:5])
            output_filter.to_csv(f'./{filename}/{profile}/{case_id}.txt', index=False, header=False)
            ##每次都重新储存了filter后的新ensembleid

def ensemblidfile(genelist, pwd, case_id, filename='scoliosis_gVCF_from_zs_updating', profile='candidategene_ensemblid',refresh=False):
    if profile not in os.listdir(f'{pwd}/{filename}/'):
        os.system(f'mkdir {pwd}/{filename}/{profile}')
    if 'candidategene_ensemblid' not in os.listdir(f'{pwd}/{filename}/'):
        os.system(f'mkdir {pwd}/{filename}/candidategene_ensemblid')
    if f'{case_id}_index.tsv' not in os.listdir(f'./{filename}/candidategene_ensemblid/'):
        output = pd.DataFrame(columns=[0])
        for gene_name_keys in genelist:
            try:
                id_result = searchensemblid(gene_name_keys)
                output.loc[gene_name_keys] = id_result
            except Exception as e:
                print(e)
        output.to_csv(f'./{filename}/candidategene_ensemblid/{case_id}.txt', index=False, header=False)
        output.to_csv(f'./{filename}/candidategene_ensemblid/{case_id}_index.tsv', sep='\t', index=True, header=False)
    else:
        if refresh:
            output = pd.DataFrame(columns=[0])
            for gene_name_keys in genelist:
                try:
                    id_result = searchensemblid(gene_name_keys)
                    output.loc[gene_name_keys] = id_result
                except Exception as e:
                    print(e)
            output.to_csv(f'./{filename}/candidategene_ensemblid/{case_id}.txt', index=False, header=False)
            output.to_csv(f'./{filename}/candidategene_ensemblid/{case_id}_index.tsv', sep='\t', index=True,
                          header=False)
        if profile != 'candidategene_ensemblid':
            ensemblids_df = pd.read_csv(f'./{filename}/candidategene_ensemblid/{case_id}_index.tsv', sep='\t', header=None)
            ensemblids_df = ensemblids_df.reset_index()
            output_filter = ensemblids_df[ensemblids_df[0].isin(genelist)][1]
            ##print(output_filter[:5])
            output_filter.to_csv(f'./{filename}/{profile}/{case_id}.txt', index=False, header=False)
            ##每次都重新储存了filter后的新ensembleid


import httplib2 as http
def searchHGNCid(genesymbol):
    try:
        from urlparse import urlparse
    except ImportError:
        from urllib.parse import urlparse

    headers = {
        'Accept': 'application/json',
    }

    uri = 'http://rest.genenames.org'
    path = f'/search/symbol:{genesymbol}'

    target = urlparse(uri + path)
    method = 'GET'
    body = ''

    h = http.Http()

    response, content = h.request(
        target.geturl(),
        method,
        body,
        headers)

    if response['status'] == '200':
    # assume that content is a json reply
    # parse content with the json module
        data = json.loads(content)
        #print('Symbol:' + data['response']['docs'][0]['hgnc_id'])
        HGNCid = data['response']['docs'][0]['hgnc_id']
        HGNCid = HGNCid.split(':')[1]
        return HGNCid
    else:
        #print('Error detected: ' + response['status'])
        return 'Error'


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


###filtered
def load_revel_filter_df(filename, case_id, thresh):
    tsv_revel_filter_dir = f'./{filename}/revel_{thresh}/{case_id}.txt'
    tsv_revel_filter = pd.read_csv(tsv_revel_filter_dir, sep='\t')
    tsv_revel_filer_ensemblid_dir = f'./{filename}/revel_{thresh}_ensemblid/{case_id}.txt'
    tsv_revel_filer_ensemblid = pd.read_csv(tsv_revel_filer_ensemblid_dir, sep='\t')
    return tsv_revel_filter, tsv_revel_filer_ensemblid


def load_cadd_filter_df(filename, case_id, thresh):
    tsv_CADD_filter_dir = f'./{filename}/cadd_{thresh}/{case_id}.txt'
    tsv_CADD_filter = pd.read_csv(tsv_CADD_filter_dir, sep='\t')
    tsv_CADD_filer_ensemblid_dir = f'./{filename}/cadd_{thresh}_ensemblid/{case_id}.txt'
    tsv_CADD_filer_ensemblid = pd.read_csv(tsv_CADD_filer_ensemblid_dir, sep='\t')
    return tsv_CADD_filter, tsv_CADD_filer_ensemblid


def selectsmall(biao, name, thresh):
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
    data_mut2 = (data_mut12[(data_mut12[name]) < thresh].append(biao[biao[name] == "."])).append(
        biao[biao[name] == "nan"])
    return data_mut2

def selectbig(biao, name, thresh):
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
    data_mut2 = (data_mut12[data_mut12[name] > thresh].append(biao[biao[name] == "."])).append(
        biao[biao[name] == "nan"])
    return data_mut2

def patho_filter_n(df, threshold):
    data3 = selectsmall(
        selectsmall(selectsmall(selectsmall(df, "ExAC_AF", 0.01), "ExAC_EAS_AF", 0.01), "gnomAD_genome_EAS", 0.01),
        "In_house_freq", 0.05)
    data4 = selectbig(data3, "VAF", 29)
    data5 = selectbig(data4, "CADD", threshold)
    return data5

