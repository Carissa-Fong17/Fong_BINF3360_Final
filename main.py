import pandas as pd
from statistics import median

count_matrix = pd.DataFrame()

# get angiosarcoma data ------------------------------------------------------------------------------------------------
asc_dat = pd.DataFrame()

files = ['C:\\Users\\panda\\Downloads\\asc_dat\\asc1.tsv',
         'C:\\Users\\panda\\Downloads\\asc_dat\\asc2.tsv',
         'C:\\Users\\panda\\Downloads\\asc_dat\\asc3.tsv',
         'C:\\Users\\panda\\Downloads\\asc_dat\\asc4.tsv',
         'C:\\Users\\panda\\Downloads\\asc_dat\\asc5.tsv',
         'C:\\Users\\panda\\Downloads\\asc_dat\\asc6.tsv',
         'C:\\Users\\panda\\Downloads\\asc_dat\\asc7.tsv']

for i, file in enumerate(files):
    cur_dat = pd.read_csv(file, sep='\t', header=1, usecols=['gene_id', 'gene_type', 'unstranded'])
    cur_dat.rename(columns={'unstranded': f'asc_{i}'}, inplace=True)

    # keep only rows for protein-coding genes
    cur_dat = cur_dat.loc[cur_dat['gene_type'] == 'protein_coding']
    cur_dat.reset_index(inplace=True, drop=True)
    cur_dat.drop('gene_type', axis=1, inplace=True)

    # trim version numbers from gene IDs
    cur_dat['gene_id'] = [gene_id[:gene_id.index('.')] for gene_id in cur_dat['gene_id']]

    if len(asc_dat) == 0:
        asc_dat['gene_id'] = cur_dat['gene_id']
        asc_dat['asc_0'] = cur_dat['asc_0']
    else:
        asc_dat = pd.merge(asc_dat, cur_dat, on='gene_id', how='outer')

# drop duplicate gene IDs
asc_dat.drop_duplicates(subset='gene_id', inplace=True)

count_matrix['gene_id'] = asc_dat['gene_id']
count_matrix.set_index('gene_id', inplace=True)
count_matrix = pd.merge(count_matrix, asc_dat, on='gene_id')
'''count_matrix['angiosarcoma'] = None
for i in range(len(asc_dat)):
    row = asc_dat.iloc[i, :]
    med_exp = median(row.iloc[1:])
    count_matrix.loc[row['gene_id'], 'angiosarcoma'] = med_exp'''


# get infiltrating duct carcinoma data ---------------------------------------------------------------------------------
idc_dat = pd.DataFrame()

files = ['C:\\Users\\panda\\Downloads\\idc_dat\\idc1.tsv',
         'C:\\Users\\panda\\Downloads\\idc_dat\\idc2.tsv',
         'C:\\Users\\panda\\Downloads\\idc_dat\\idc3.tsv',
         'C:\\Users\\panda\\Downloads\\idc_dat\\idc4.tsv',
         'C:\\Users\\panda\\Downloads\\idc_dat\\idc5.tsv',
         'C:\\Users\\panda\\Downloads\\idc_dat\\idc6.tsv',
         'C:\\Users\\panda\\Downloads\\idc_dat\\idc7.tsv']

for i, file in enumerate(files):
    cur_dat = pd.read_csv(file, sep='\t', header=1, usecols=['gene_id', 'gene_type', 'unstranded'])
    cur_dat.rename(columns={'unstranded': f'idc_{i}'}, inplace=True)

    # keep only rows for protein-coding genes
    cur_dat = cur_dat.loc[cur_dat['gene_type'] == 'protein_coding']
    cur_dat.reset_index(inplace=True, drop=True)
    cur_dat.drop('gene_type', axis=1, inplace=True)

    # trim version numbers from gene IDs
    cur_dat['gene_id'] = [gene_id[:gene_id.index('.')] for gene_id in cur_dat['gene_id']]

    if len(idc_dat) == 0:
        idc_dat['gene_id'] = cur_dat['gene_id']
        idc_dat['idc_0'] = cur_dat['idc_0']
    else:
        idc_dat = pd.merge(idc_dat, cur_dat, on='gene_id', how='outer')

# drop duplicate gene IDs
idc_dat.drop_duplicates(subset='gene_id', inplace=True)

count_matrix = pd.merge(count_matrix, idc_dat, on='gene_id')

'''count_matrix['idc'] = None
for i in range(len(idc_dat)):
    row = idc_dat.iloc[i, :]
    med_exp = median(row.iloc[1:])
    count_matrix.loc[row['gene_id'], 'idc'] = med_exp'''


# get healthy tissue data ----------------------------------------------------------------------------------------------
healthy_dat = pd.read_csv('C:\\Users\\panda\\Downloads\\gene_reads_2017-06-05_v8_breast_mammary_tissue.gct\\healthy.gct',
                          header=2, sep='\t', index_col=0)
healthy_dat.drop('Description', axis=1, inplace=True)
healthy_dat = healthy_dat.iloc[:, :8]
col_names = ['gene_id']
col_names.extend([f'healthy_{i}' for i in range(len(healthy_dat.columns)-1)])
healthy_dat.columns = col_names
#healthy_dat.rename(columns={'Name': 'gene_id'}, inplace=True)

# trim version numbers from gene IDs
healthy_dat['gene_id'] = [gene_id[:gene_id.index('.')] for gene_id in healthy_dat['gene_id']]

# find indices for genes in count matrix
#use_ids = [i for i in range(len(healthy_dat)) if healthy_dat.loc[i, 'gene_id'] in count_matrix['gene_id']]
#healthy_dat = healthy_dat.iloc[use_ids, :].reset_index(drop=True)

# find median expression for every gene
'''healthy_dat['med_exp'] = None
for i in range(len(healthy_dat)):
    row = healthy_dat.iloc[i, :]
    med_exp = median(row.iloc[2:-1])
    healthy_dat.loc[i, 'med_exp'] = med_exp'''
#healthy_dat = healthy_dat[['gene_id', 'med_exp']]
#healthy_dat.rename(columns={'med_exp': 'healthy'}, inplace=True)

# drop duplicate gene IDs
healthy_dat.drop_duplicates(subset='gene_id', inplace=True)

count_matrix = pd.merge(count_matrix, healthy_dat, on='gene_id')
count_matrix.dropna(inplace=True)
count_matrix.set_index('gene_id', inplace=True)
count_matrix.to_csv('C:\\Users\\panda\\Downloads\\deseq_inputs\\count_matrix.csv', sep='\t')


# make coldata dataframe -----------------------------------------------------------------------------------------------
coldata = pd.DataFrame()
condition = []
idx = []
for sample in ['asc', 'idc', 'healthy']:
    for i in range(7):
        condition.append(sample)
        idx.append(f'{sample}_{i}')
coldata['condition'] = condition
coldata['idx'] = idx
coldata.set_index('idx', inplace=True)
coldata.to_csv('C:\\Users\\panda\\Downloads\\deseq_inputs\\coldata.csv', sep='\t')
