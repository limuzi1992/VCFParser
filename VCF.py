import os
import re
import sys
import glob
import pandas as pd
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import statsmodels.stats.multicomp as ssm



def Get_Options():
    usage = 'Usage: %prog [options] arg'
    version="%prog 1.0"
    parser = OptionParser(usage=usage, version=version)
    parser.add_option('-i', '--input', dest='invcf', metavar='InputVCF', help='Input data from VCF file or folder')
    parser.add_option('-g', '--group', dest='ingroup', metavar='InputGroupFile', help='Input data from group file')
    parser.add_option('-e', '--extract', dest='inextract', metavar='InputExtractionFile', help='Extract a list of mutations')
    parser.add_option('-f', '--filter', dest='filtercond', metavar='FilterCondition', help='Filter the files based on specific conditions (e.g. DP>15)')
    parser.add_option('-a', '--analyze', dest='anatype', metavar='AnalysisType', help='Analyze the files based on analysis types (Plot_Histogram, Plot_Boxplot, Analyze_Genes, Mutation_Mapper)')
    parser.add_option('-w', '--write', dest='dfout', metavar='OutputTextFile', help='Write the dataframe to a text file')
    parser.add_option('-p', '--plot', dest='plotout', metavar='OutputPlotFile', help='Save the plot to a png file')
    (options, args) = parser.parse_args()
    input_vcf = options.invcf
    input_group = options.ingroup
    input_extract = options.inextract
    filter_cond = options.filtercond
    analysis_type = options.anatype
    df_out = options.dfout
    plot_out = options.plotout
    if input_vcf == None:
        print >>sys.stderr, 'Please provide your input folder!'
        sys.exit(1)
    if input_group == None:
        g_flag = False
    else:
        g_flag = True
    if input_extract == None:
        e_flag = False
    else:
        e_flag = True
    if filter_cond == None:
        f_flag = False
    else:
        f_flag = True
    if analysis_type == None:
        a_flag = False
    else:
        a_flag = True
    if df_out == None:
        d_flag = False
    else:
        d_flag = True
    if plot_out == None:
        p_flag = False
    else:
        p_flag = True
    return input_vcf, input_group, g_flag, input_extract, e_flag, filter_cond, f_flag, analysis_type, a_flag, df_out, d_flag, plot_out, p_flag


def Parse_INFO(info, info_field, ann_ids):
    temp = {}
    for i in info: 
        if '=' in i:
            e = i.find('=')
            temp[i[0:e]] = i[e+1:]
        else:
            temp[i] = i
    keys = temp.keys()
    info_dict = {}
    for info_id, info_type, info_num in info_field:
        if info_id in keys:
            if info_num == '.':
                if info_id == 'ANN':
                    anns = temp['ANN'].split(',')[0].split('|')
                    for i in range(0, len(anns)):
                        if anns[i] == '':
                            info_dict[ann_ids[i]] = '.'
                        else:
                            info_dict[ann_ids[i]] = anns[i]
                else:
                    info_value = temp[info_id]
                    info_value = info_value[info_value.find('(')+1: info_value.find(')')]
                    info_dict[info_id] = info_value
            else:
                if info_type == 'Integer':
                    if info_num > 1:
                        values = temp[info_id].split(',')
                        for i in range(0, info_num):
                            info_dict[info_id+'_'+str(i+1)] = int(values[i])
                    else:
                        info_dict[info_id] = int(temp[info_id])
                elif info_type == 'Float':
                    if info_num > 1:
                        values = temp[info_id].split(',')
                        for i in range(0, info_num):
                            info_dict[info_id+'_'+str(i+1)] = float(values[i])
                    else:
                        info_dict[info_id] = float(temp[info_id]) 
                else:
                      info_dict[info_id] = temp[info_id]  
        else:
            if info_num == '.':
                info_dict[info_id] = '.'
            else:
                if info_num > 1:
                    for i in range(0, info_num):
                        info_dict[info_id+'_'+str(i+1)] = '.'
                else:
                    info_dict[info_id] = '.'
    return info_dict


def Parse_VCF(input_vcf):
    df_list = []
    for vcf_file in glob.glob(os.path.join(input_vcf, '*.vcf')):
        print vcf_file
        vcf = open(vcf_file, 'r')
        meta = []
        info_field = []
        ann_oids = []
        info_ids = []
        for line in vcf:
            if line.startswith('##'): 
                meta.append(line) #Get the Meta-information
                if 'INFO' in line: #Get the INFO IDs, Types and Numbers 
                    info_id = line[line.find('ID=')+3:line.find('Number=')-1]
                    info_type = line[line.find('Type=')+5:line.find('Description=')-1]
                    info_num = line[line.find('Number=')+7:line.find('Type=')-1]
                    if info_num == '.': #Get the Functional Annotations 
                        info_field.append((info_id, info_type, info_num)) #Get All INFO fields
                        if info_id == 'ANN': 
                            ann_ids = map(str.strip,line[line.find("'")+1:line.rfind("'")].split('|'))
                        else:
                            ann_oids.append(info_id)
                    else: #Get the INFO IDs Except Annotation Fields
                        info_num = int(info_num)
                        info_field.append((info_id, info_type, info_num)) #Get All INFO fields
                        if info_num > 1:
                            for i in range(0, info_num):
                                info_ids.append(info_id+'_'+str(i+1))
                        else:
                            info_ids.append(info_id) 
            elif line.startswith('#'):
                header = line[1:].strip().split('\t')
                bam = header[-1]
                sample = bam[0: bam.find('.')] #Get the Sample Name
            else:
                data = line.strip().split('\t')
                data_dict = {}
                for i in range(0,7):
                    if i == 1:
                        data_dict[header[i]] = int(data[i])
                    elif i == 5:
                        data_dict[header[i]] = float(data[i])
                    else:
                        data_dict[header[i]] = data[i]
                key = sample+'_'+data[0]+'_'+data[1]+'_'+data[3]+'_'+data[4]
                data_dict['Key'] = key
                if len(data[3]) > len(data[4]):
                    alt_key = sample + '_' + data[0] + '_' + str(int(data[1])+len(data[4]))
                else:
                    alt_key = sample + '_' + data[0] + '_' + data[1]
                data_dict['AlternativeKey'] = alt_key
                format_ids = data[8].split(':')
                format_values = data[9].split(':')
                for i in range(0, len(format_ids)):
                    data_dict[format_ids[i]] = format_values[i]
                data_dict['Sample'] = sample
                info = data[7].split(';')
                info_dict = Parse_INFO(info, info_field, ann_ids)
                data_dict.update(info_dict)
                df_list.append(data_dict) #List of All Mutations
        vcf.close()
    df_cols = ['Sample'] + header[0:6] + format_ids + info_ids + ann_ids + ann_oids + ['Key', 'AlternativeKey']#All Column Names
    df = pd.DataFrame(df_list) #A Dataframe of All Mutations
    df = df[df_cols]
    return df
    

def Add_Key(df, key_flag):
    def Select(row):
        if len(row['REF'])>len(row['ALT']):
            return row['Sample'] + '_' + row['CHROM'] + '_' + str(row['POS']+len(row['ALT']))
        else:
            return row['Sample'] + '_' + row['CHROM'] + '_' + str(row['POS'])
    if key_flag == 'Key':
        df['Key'] = df['Sample'] + '_' + df['CHROM'] + '_' + map(str,df['POS']) + '_' + df['REF'] + '_' + df['ALT']
    if key_flag == 'AlternativeKey':
        df['AlternativeKey'] = df.apply (lambda row: Select (row), axis=1)
    return df
    

def Add_Group(df, input_group):
    group = pd.read_table(input_group, sep='\t')
    group = group.astype(str)
    group_dict = dict(zip(group.iloc[:,0], group.iloc[:,1]))
    df_cols = list(df)
    df_cols = df_cols[0:1] + ['Group'] + df_cols[1:]
    df['Group'] = df['Sample'].apply (lambda sample: group_dict.get(sample, '.'))
    df = df[df_cols]
    return df
    
    
def Extract_List(df, input_extract):
    df_cols = list(df)[0:-2] + ['Allele_Type', 'Score', 'Coverage', 'Key', 'AlternativeKey'] 
    df = df.assign(Key=df['Sample']+'_'+df['CHROM']+'_'+map(str,df['POS'])+'_'+df['REF']+'_'+df['ALT'])
    chroms = []
    for i in range(1,23):
        chroms.append('chr'+str(i))
    chroms = chroms + ['chrX', 'chrY', 'chrM']
    if os.path.isdir(input_extract):
        df_list = []
        for extract_file in glob.glob(os.path.join(input_extract, '*')):
            print extract_file
            sample_dict = {}
            path = extract_file.split('/')[-1]
            sample = path[0:path.find('.')]
            extract = open(extract_file, 'r')
            for line in extract:
                data = line.strip().split('\t')
                key = sample + '_' + data[0] + '_' + data[1] + '_' + data[3] + '_' + data[4] 
                alt_key = sample + '_' + data[0] + '_' + data[1]
                try:
                    select_row = df[df['Key'] == key]
                    row = select_row.iloc[0,]
                    row = row.append(pd.Series([data[5],data[6],int(data[7])] , index=['Allele_Type', 'Score', 'Coverage']))
                    flag = True
                except:
                    try:
                        select_row = df[df['AlternativeKey'] == alt_key]
                        row = select_row.iloc[0,]
                        row = row.append(pd.Series([data[5],data[6],int(data[7])] , index=['Allele_Type', 'Score', 'Coverage']))
                        flag = True
                    except:
                        flag = False
                if flag:
                    if data[0] not in sample_dict:
                        sample_dict[data[0]] = [row]
                    else:
                        sample_dict[data[0]].append(row)
            for chrom in chroms:
                try:
                    sample_list = sample_dict[chrom]
                    sample_list = sorted(sample_list, key=lambda x: x['POS'])
                    df_list.extend(sample_list)
                except:
                    pass
            extract.close()
        df = pd.DataFrame(df_list)
        df = df[df_cols]
        return df       
    
            
#-----------Filter------------------
def create_code(filter_cond):
    split = re.split('&|\|', filter_cond) 
    if len(split)==1:
        c = split[0]
        code = 'df.' + c
    else:
        if '&' in filter_cond:
            c1 = split[0]
            c2 = split[1]
            code = '(df.' + c1 + ')&(df.' + c2 + ')'
        elif '|' in filter_cond:
            c1 = split[0]
            c2 = split[1]
            code = '(df.' + c1 + ')|(df.' + c2 + ')'   
        else:
            print 'Please combine your two conditions using & or |'
    return code
       
def one_filter(filter_cond):
    code = create_code(filter_cond)
    code = 'df = df.loc[' + code + ']' 
    return code

def two_filter(match, filter_conds):
    m1 = match[0]
    m2 = match[1]
    code1 = create_code(m1)
    code2 = create_code(m2)   
    if ')&(' in filter_conds:
        code = 'df = df.loc[(' + code1 + ')&(' + code2 + ')]'
    elif ')|(' in filter_conds:
        code = 'df = df.loc[(' + code1 + ')|(' + code2 + ')]'
    else:
        print 'Please combine your two conditions using & or |'
    return code
         
def attribute_error(c, df):
    match  = re.findall(r'([a-zA-Z_]+)(==|!=|>|>=|<|<=)(.)', c)
    df_cols = list(df)
    old_attri = ''
    err_flag = False
    for i in match:
        attribute = i[0]
        if attribute != old_attri:
            if attribute not in df_cols:
                print '***************ERROR: NO ATTRIBUTE '+"'"+attribute+"'"+'***************'
                err_flag = True
        old_attri = attribute
    return err_flag

def avail_attri(err_flag, df):
    df_cols = list(df)
    if err_flag == True:    
        print 'Available Attributes:'
        df_cols_str = ', '.join(df_cols)
        print df_cols_str
        sys.exit(1)

def Filter_Dataframe(df, filter_conds):
    print filter_conds
    if '.txt' in filter_conds:
        colname = filter_conds.split('/')[-1][:-4]
        filter_df = pd.read_table(filter_conds)
        filter_list = list(filter_df.iloc[:,0])
        filter_output = df.loc[df[colname].isin(filter_list)]
        return filter_output
    else:
        match = re.findall('\((.*?)\)', filter_conds)
        if match==[]:
            code = one_filter(filter_conds)
            try:
                exec code
            except:
                err_flag = attribute_error(filter_conds, df)
                avail_attri(err_flag, df)         
        else:
            code = two_filter(match, filter_conds)
            try:
                exec code
            except:
                err_flag1 = attribute_error(match[0], df)
                err_flag2 = attribute_error(match[1], df)
                err_flag = err_flag1 | err_flag2
                avail_attri(err_flag, df)
    return df


#-----------Output the Result------------------     
def Write_Dataframe(df, df_out):
    print 'Writing...'
    columns = raw_input('Do you want to include all or user-defined columns? (All|Custom)')
    if columns == 'All':
        df.to_csv(df_out, index=None, sep='\t')
    elif columns == 'Custom':
        custom = raw_input('Which columns do you want to include? (e.g. Sample,CHROM,POS,REF,...)')
        customs = map(str.strip,custom.split(','))
        try:
            df = df[customs]
            df.to_csv(df_out, index=None, sep='\t')
        except:
            df_cols = list(df)
            df_cols_str = ', '.join(df_cols)
            print 'Available Attributes:'
            print df_cols_str
            sys.exit(1)        
    else:
        print 'Please choose from All and Custom'


#---------------------------Analyze the VCF Files--------------------------
def Analyze_VCF(df, analysis_type):
    df_cols = list(df)
    if 'Group' in df_cols:
        df_sg = df.groupby(['Sample','Group']).size().reset_index().rename(columns={0:'Frequency'})
        df['Sample_Group'] = df['Sample'] + '_' + df['Group']
        samples = list(df_sg.Sample)
    if analysis_type == 'Plot_Histogram':
        Plot_Histogram(df)
    elif analysis_type == 'Plot_Boxplot':
        Plot_Boxplot(df_sg)
    elif analysis_type == 'Analyze_Genes':
        df_gene = Analyze_Genes(df, df_sg, samples)
        return df_gene
    elif analysis_type == 'Mutation_Mapper':
        df_mapper = Create_Mapper_Input(df)
        return df_mapper
      
def  Plot_Histogram(df):
    if 'Sample_Group' in list(df):
        bar = pd.value_counts(df['Sample_Group'], sort=True)
    else:
        bar = pd.value_counts(df['Sample'], sort=True)
    bar_color = [(x/15.0, 0.5, 0.75) for x in range(len(bar))]
    fig = bar.plot(kind='bar', figsize=(20.0, 20.0), fontsize=15, color=bar_color)
    plt.figtext(.5,.95,'Number of Mutations', fontsize=30, ha='center')
    return fig
    
def Plot_Boxplot(df):
    boxprops = dict(linestyle='-', linewidth=4, color='k')
    medianprops = dict(linestyle='--', linewidth=4, color='k')
    fig = df.boxplot(column='Frequency', by='Group', figsize=(20.0, 20.0), fontsize=20, return_type='dict', boxprops=boxprops, medianprops=medianprops)
    return fig
    
def Get_Gene_Table(df, samples):
    gene_dict = {}
    for row in df.itertuples():
        sample = row[1]
        gene_name = row[42]
        if gene_name not in gene_dict:
            gene_dict[gene_name] = {sample:1}
        else:
            if sample not in gene_dict[gene_name]:
                gene_dict[gene_name][sample]=1
            else:
                gene_dict[gene_name][sample] = gene_dict[gene_name][sample]+1                
    genes = gene_dict.keys()
    df_list = []
    for g in genes:
        freq = gene_dict[g]
        df_dict = {}
        df_dict['Gene'] = g
        for s in samples:
            df_dict[s] = freq.get(s,0)
        df_list.append(df_dict)
    df_gene = pd.DataFrame(df_list)
    df_gene_cols = list(df_gene)
    df_gene_cols.remove('Gene')
    df_gene_cols = ['Gene'] + df_gene_cols
    df_gene = df_gene.loc[:,df_gene_cols]
    return df_gene

def unpaired_t_test(row, g1, g2):
    group1 = row.loc[g1]
    group2 = row.loc[g2]
    t_statistic, p_value = stats.ttest_ind(group1,group2)
    return p_value

def mann_whitney_test(row, g1, g2):
    group1 = row.loc[g1]
    group2 = row.loc[g2]
    u_statistic, p_value = stats.mannwhitneyu(group1,group2,alternative='two-sided')
    return p_value

def fisher_exact_test(row, g1, g2):
    g1m=0
    g1nm=0
    for i in g1:
        if row.loc[i]!=0:
            g1m=g1m+1
        else:
            g1nm=g1nm+1
    g2m=0
    g2nm=0
    for j in g2:
        if row.loc[j]!=0:
            g2m=g2m+1
        else:
            g2nm=g2nm+1
    obs=np.array([[g1m,g2m],[g1nm,g2nm]])
    oddsratio, p_value = stats.fisher_exact(obs)
    return p_value

def Stat_Test(df_gene, df_sg):
    comp_groups = raw_input('Choose the groups you want to compare. (e.g. good/bad or good/no ...)')
    tests = raw_input('Which statistical test you want to do? (e.g. unpaired_t_test or mann_whitney_test or fisher_exact_test or mann_whitney_test&fisher_exact_test ...)')
    group_list = comp_groups.split('/')
    test_list = tests.split('&')
    if len(group_list) == 2:
        group1 = group_list[0]
        group2 = group_list[1]
        list_g1 = list(df_sg.loc[df_sg.Group==group1,'Sample'])
        list_g2 = list(df_sg.loc[df_sg.Group==group2,'Sample'])
        df_gene = df_gene.loc[:, ['Gene'] + list_g1 + list_g2]
        df_gene[group1+'_Sum'] = df_gene.loc[:,list_g1].sum(axis=1)
        df_gene[group2+'_Sum'] = df_gene.loc[:,list_g2].sum(axis=1)
        df_gene = df_gene.loc[(df_gene[group1+'_Sum']!=0) & (df_gene[group2+'_Sum']!=0), :]
        for t in test_list:
            if t == 'unpaired_t_test':
                df_gene['unpaired_t_test'] = df_gene.apply (lambda row: unpaired_t_test(row, list_g1, list_g2), axis=1)
            elif t == 'mann_whitney_test':
                df_gene['mann_whitney_test'] = df_gene.apply (lambda row: mann_whitney_test(row, list_g1, list_g2), axis=1)
            elif t == 'fisher_exact_test':
                df_gene['fisher_exact_test'] = df_gene.apply (lambda row: fisher_exact_test(row, list_g1, list_g2), axis=1)
            else:
                print 'Please choose from unpaired_t_test, mann_whitney_test and fisher_exact_test'
                sys.exit(1)
    p_cor_flag = raw_input('Do you want to perform mutiple testing correction? (Y|N)')
    if p_cor_flag == 'Y':
        p_cor = raw_input('Which correction method do you want to choose? (e.g. bonferroni or fdr_bh or bonferroni&fdr_bh)')
        p_cor = p_cor.split('&')
        p_val = raw_input('Which p value do you want to correct? (e.g. unpaired_t_test or mann_whitney_test or fisher_exact_test')
    for cor in p_cor:
        df_gene[cor] = ssm.MultiComparison(df_gene[p_val], alpha=0.05, method=cor)
    return df_gene

def Analyze_Genes(df, df_sg, samples):
    stat_test = raw_input('Do you want to do statistical test? (Y|N)')
    if stat_test == 'Y':
        df_gene = Get_Gene_Table(df, samples)
        df_gene = Stat_Test(df_gene, df_sg)
    elif stat_test == 'N':
        df_gene = Get_Gene_Table(df, samples)
    else:
        print 'Please choose from Y and N'
        sys.exit(1)
    return df_gene

def Create_Mapper_Input(df):
    genes = raw_input('Please input the genes of interest (MET&IGF1R)')
    if '.txt' in genes:
        genes = pd.read_table(genes, sep='\t')
        genes = list(genes.Gene)
    else:
        genes = genes.split('&')
    df_mapper = df.loc[df['Gene_Name'].isin(genes)]
    group_flag = raw_input('Do you want to select mutations in certain groups? (Y|N)')
    if group_flag == 'Y':
        group = raw_input('Which group do you want to choose? (good|bad|no|good&bad)')
        group = group.split('&')
        df_mapper = df_mapper.loc[df_mapper['Group'].isin(group)]
    elif group_flag == 'N':
        pass
    else:
        print 'Please choose from Y and N'
        sys.exit(1)
    df_mapper = df_mapper[['Gene_Name', 'Sample_Group', 'CHROM', 'POS', 'REF', 'ALT', 'Annotation', 'HGVS.p']]
    def pc(x):
        if x != '.':
            return x[2:]
        else:
            return x
    df_mapper['Protein_Change'] = df_mapper['HGVS.p'].apply(lambda x: pc(x))
    df_mapper['End_Position'] = df_mapper['POS']
    df_mapper = df_mapper[['Gene_Name', 'Sample_Group', 'Protein_Change', 'Annotation', 'CHROM', 'POS', 'End_Position', 'REF', 'ALT']]
    df_mapper.columns = ['Hugo_Symbol', 'Sample_ID', 'Protein_Change', 'Mutation_Type', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Variant_Allele']
    return df_mapper


#------------------------------Main Function------------------------------
def Main():  
    input_vcf, input_group, g_flag, input_extract, e_flag, filter_cond, f_flag, analysis_type, a_flag, df_out, d_flag, plot_out, p_flag = Get_Options()
    if os.path.isdir(input_vcf):
        print 'Importing Data...'
        df = Parse_VCF(input_vcf)
    else:
        print 'Importing Data...'
        print input_vcf
        df = pd.read_table(input_vcf, sep='\t')
        df_cols = list(df)
        if 'Key' not in df_cols:
            df = Add_Key(df, 'Key')
        if 'AlternativeKey' not in df_cols:
            df = Add_Key(df, 'AlternativeKey')
    if g_flag: 
        print 'Adding Group Info...'
        print input_group
        df = Add_Group(df, input_group)
    if e_flag:
        print 'Extracting...'
        df = Extract_List(df, input_extract)
    if f_flag:
        print 'Filtering...'
        df = Filter_Dataframe(df, filter_cond)
    if a_flag:
        print 'Analyzing...' + ' (' + analysis_type + ')'
        if analysis_type == 'Analyze_Genes' or analysis_type == 'Mutation_Mapper':
            df = Analyze_VCF(df, analysis_type)
        elif analysis_type == 'Plot_Histogram' or analysis_type == 'Plot_Boxplot':
            if p_flag:
                print 'Ploting...'
                Analyze_VCF(df, analysis_type)
                plt.savefig(plot_out)
            else:
                print 'Please use -p option to enter the directory where you want to save the figure'
        else:
            print 'Please choose from the four analysis types: Plot_Histogram, Plot_Boxplot, Analyze_Genes, and Mutation_Mapper'
    if d_flag:
        Write_Dataframe(df, df_out)
        
        
Main()
        


             
                

