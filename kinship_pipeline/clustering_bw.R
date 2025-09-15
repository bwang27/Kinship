##setwd("R:/geelegrp/projects/PedDep_DataScience/common/WGS/PedDepAcc/Modupe/fromJohn/bam2ai/test_1_v_all_incohort")

#setwd("R:/research/rgs01/home/clusterHome/bwang68/kinship")

#make_mat()

dt<-as.matrix(read.table("thematrix.txt",header=TRUE,sep="\t",row.names=1))

library(pheatmap)

rlt<-pheatmap(dt,cluster_rows = TRUE,cluster_cols = TRUE, clustering_distance_rows = "euclidean",clustering_distance_cols = "euclidean", clustering_method = "complete")

ndt<-dt[rlt$tree_row$order,rlt$tree_col$order]
ndt[ndt>40]<-40
write.table(ndt,"thematrix_clustered.txt",row.names=TRUE,col.names=TRUE,sep="\t",quote=F)


