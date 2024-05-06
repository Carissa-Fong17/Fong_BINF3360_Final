counts <- read.csv(file = 'deseq_inputs/count_matrix.csv', sep = '\t', 
                  row.names = 'gene_id')
coldata <- read.csv(file = 'deseq_inputs/coldata.csv', sep = '\t', row.names = 1)

# make count and coldata dtaframes for each combination of disease/healthy
# angiosarcoma vs healthy
healthy_asc_counts <- counts[c('asc_0', 'asc_1', 'asc_2', 'asc_3', 'asc_4',
                               'asc_5', 'asc_6', 'healthy_0', 'healthy_1',
                               'healthy_2', 'healthy_3', 'healthy_4', 
                               'healthy_5', 'healthy_6')]
healthy_asc_coldata <- data.frame('condition' = coldata[c(1:7, 15:21), 'condition'],
                                  row.names = rownames(coldata)[c(1:7, 15:21)])


# infiltrating duct carcinoma vs healthy
healthy_idc_counts <- counts[c('idc_0', 'idc_1', 'idc_2', 'idc_3', 'idc_4',
                               'idc_5', 'idc_6', 'healthy_0', 'healthy_1',
                               'healthy_2', 'healthy_3', 'healthy_4', 
                               'healthy_5', 'healthy_6')]
healthy_idc_coldata <- data.frame('condition' = coldata[c(8:21), 'condition'],
                                  row.names = rownames(coldata)[c(8:21)])

# angiosarcoma vs infiltrating duct carcinoma
asc_idc_counts <- counts[c('asc_0', 'asc_1', 'asc_2', 'asc_3', 'asc_4',
                           'asc_5', 'asc_6', 'idc_0', 'idc_1', 'idc_2', 
                           'idc_3', 'idc_4', 'idc_5', 'idc_6')]
asc_idc_coldata <- data.frame('condition' = coldata[c(1:14), 'condition'],
                              row.names = rownames(coldata)[c(1:14)])


# make new dataframe to compare disease (asc and idc) vs healthy
coldata2 <- data.frame(coldata)
coldata2$condition[1:14] <- 'disease'

# make condition columns in coldata factors
healthy_asc_coldata$condition <- factor(healthy_asc_coldata$condition)
healthy_idc_coldata$condition <- factor(healthy_idc_coldata$condition) 
asc_idc_coldata$condition <- factor(asc_idc_coldata$condition)
coldata2$condition = factor(coldata2$condition)


library(DESeq2)

# angiosarcoma vs healthy
dds_healthy_asc <- DESeqDataSetFromMatrix(countData = healthy_asc_counts,
                              colData = healthy_asc_coldata,
                              design = ~ condition)
dds_healthy_asc <- DESeq(dds_healthy_asc)
resultsNames(dds_healthy_asc)
res_asc_healthy <- results(dds_healthy_asc, name='condition_healthy_vs_asc')
res_asc_healthy[order(res_asc_healthy$padj)[1:10],]
sum(res_asc_healthy$padj < 0.01, na.rm = T)
sig_genes <- na.omit(rownames(res_asc_healthy)[res_asc_healthy$padj < 0.01])
write.csv(noquote(sig_genes), 'res_asc_healthy.txt', row.names = F, quote = F)


# infiltrating duct carcinoma vs healthy
dds_healthy_idc <- DESeqDataSetFromMatrix(countData = healthy_idc_counts,
                                          colData = healthy_idc_coldata,
                                          design = ~ condition)
dds_healthy_idc <- DESeq(dds_healthy_idc)
resultsNames(dds_healthy_idc)
res_idc_healthy <- results(dds_healthy_idc, name='condition_idc_vs_healthy')
res_idc_healthy[order(res_idc_healthy$padj, decreasing = T)]
sum(res_idc_healthy$padj < 0.01, na.rm = T)
sig_genes <- na.omit(rownames(res_idc_healthy)[res_idc_healthy$padj < 0.01])
write.csv(noquote(sig_genes), 'res_idc_healthy.txt', row.names = F, quote = F)


# angiosarcoma vs infiltrating duct carcinoma
dds_asc_idc <- DESeqDataSetFromMatrix(countData = asc_idc_counts,
                                      colData = asc_idc_coldata,
                                      design = ~ condition)
dds_asc_idc <- DESeq(dds_asc_idc)
resultsNames(dds_asc_idc)
res_asc_idc <- results(dds_asc_idc, name='condition_idc_vs_asc')
res_asc_idc
sum(res_asc_idc_$padj < 0.01, na.rm = T)
sig_genes <- na.omit(rownames(res_asc_idc)[res_asc_idc$padj < 0.01])
write.csv(noquote(sig_genes), 'res_asc_idc.txt', row.names = F, quote = F)
