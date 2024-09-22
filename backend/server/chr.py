import pyensembl

# 加载 Ensembl 数据，选择人类基因组
data = pyensembl.EnsemblRelease(75)  # Ensembl release for human genome

# 定义染色体和位置
chromosome = "1"
start = 3119542
end = 3120502

# 查找该位置的基因
genes = data.genes_at_locus(contig=chromosome, position=start)

# 提取并输出基因符号
gene_symbols = [gene.gene_name for gene in genes]
print(gene_symbols)
