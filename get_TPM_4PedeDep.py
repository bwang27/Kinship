import pandas as pd 
import numpy as np 
import os

def count_to_tpm(counts_df, gene_lengths):
	"""
	Converts raw gene counts to Transcripts Per Million (TPM) values. 
	Args:
	counts_df (pd.DataFrame): DataFrame with gene counts, where rows are genes and columns are samples.
	gene_lengths (pd.Series): Series containing the length of each gene, index should match counts_df.index.
	Returns:
	pd.DataFrame: DataFrame with TPM values for each gene in each sample.
	"""
	lengths_kb = gene_lengths / 1000

	rpk = counts_df.divide(gene_lengths, axis=0)
	# Calculate TPM (Transcripts Per Million)
	tpm = rpk.divide(rpk.sum(axis=0), axis=1) * 1000000
	return tpm 

genes={}	
counts_data={}
with open ("file.lst","r") as in_fh:
	for line in in_fh:
		fn=line.strip()
		sample=os.path.basename(fn).split(".")[1]
		counts_data[sample]={}
		with open (fn,"r") as in_fh:
			for line in in_fh:
				lst=line.strip().split("\t")
				gene=lst[0]
				if gene.startswith("__"):
					continue
				genes[gene]=1
				count=lst[1]
				counts_data[sample][gene]=count

samples_lst=counts_data.keys()
genes_lst=sorted(genes.keys())

#print("read in done")
count_matrix_fn="PedDep_hg38_RNASEQ_all_count.txt"				
count_matrix_fh=open(count_matrix_fn,"w")
header="gene\t"+("\t").join(samples_lst)+"\n"
count_matrix_fh.write(header)
for gene in genes_lst:
	tmp_count_lst=[]
	for sample in samples_lst:
		if sample not in counts_data:
			continue
		if gene not in counts_data[sample]:
			sample_gene_count="0"
		else:
			sample_gene_count=counts_data[sample][gene]
		tmp_count_lst.append(sample_gene_count)	
	count_matrix_fh.write(gene+"\t"+("\t").join(tmp_count_lst)+"\n")	
count_matrix_fh.close()		

#print("count matrix done")
counts_df = pd.read_csv(count_matrix_fn,sep="\t", index_col=0)
sample_lst=[]
for sample in counts_df.columns:
	if "SJ" in sample:
		sample_lst.append(sample)
counts_df=counts_df[sample_lst]

gene_annotation_data={}
gene_length_data={}
with open ("gencode.v31.primary_assembly.annotation.gtf.gene","r")as in_fh:
	in_fh.readline()
	for line in in_fh:
		lst=line.strip().split("\t")
		gene=lst[1]
		gene_length=int(lst[-1])
		gene_id=lst[0]
		gene_type=lst[2]
		gene_chro=lst[4]
		gene_start=lst[5]
		gene_end=lst[6]
		gene_strand="ND"
		gene_length_data[gene]=[gene_length,gene_id,gene_type,gene_chro,gene_start,gene_end,gene_strand]

drop_genes_lst=[]
selected_genes_lst=[]
selected_gene_length_lst=[]

selected_genes_id_lst=[]
selected_genes_type_lst=[]
selected_genes_chro_lst=[]
selected_genes_start_lst=[]
selected_genes_end_lst=[]
selected_genes_strand_lst=[]

for gene in counts_df.index:
	if gene not in gene_length_data:
		drop_genes_lst.append(gene)
	else:
		gene_length,gene_id,gene_type,gene_chro,gene_start,gene_end,gene_strand=gene_length_data[gene]
		selected_gene_length_lst.append(gene_length)
		selected_genes_lst.append(gene)
		selected_genes_id_lst.append(gene_id)
		selected_genes_type_lst.append(gene_type)
		selected_genes_chro_lst.append(gene_chro)
		selected_genes_start_lst.append(gene_start)
		selected_genes_end_lst.append(gene_end)
		selected_genes_strand_lst.append(gene_strand)
gene_lengths = pd.Series(selected_gene_length_lst, index=selected_genes_lst)  
counts_df=counts_df.loc[selected_genes_lst]
print("counts_df")
print(counts_df)
gene_annotation_df=pd.DataFrame(
	{"gene":selected_genes_lst,
	"geneid":selected_genes_id_lst,
	"gene_type":selected_genes_type_lst,
	"gene_chro":selected_genes_chro_lst,
	"gene_start":selected_genes_start_lst,
	"gene_end":selected_genes_end_lst,
	"gene_strand":selected_genes_strand_lst
	}
)	
gene_annotation_df=gene_annotation_df.set_index("gene")
tpm_result = count_to_tpm(counts_df, gene_lengths)

tpm_result=round(tpm_result,2)
out_df=pd.merge(gene_annotation_df,tpm_result,on="gene",how="inner")
out_df.to_csv('PedDep_hg38_RNASEQ_all_TPM.txt', sep="\t",index=True)
print("done")
