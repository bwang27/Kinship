
ml R/4.0.0
#CCWD=/research_jude/rgs01_jude/resgen/prod/tartan/index/data/PedDep/PedDepAccelPilot
CCWD=/research/rgs02/collaborations/PedDepAccel/Auto/tartan/index/PedDep/PedDep
QCbam=./obtain_bam_stats.r

removed="removed"


ls ${CCWD}/S*/TRANSCRIPTOME/cicero-post/final_fusions.counts.medals.tab | fgrep -v -f ${removed} | awk -F'/' '{print $NF "\t" $0}'   > cicero_fusion.lst
ls ${CCWD}/S*/TRANSCRIPTOME/cicero-post/final_internal.counts.medals.tab | fgrep -v -f ${removed} | awk -F'/' '{print $NF "\t" $0}'  > cicero_internal.lst


ls ${CCWD}/S*/WHOLE_GENOME/bam/SJ*.bam | fgrep -v -f ${removed}  > bam.WGS.lst
ls ${CCWD}/S*/TRANSCRIPTOME/bam/SJ*.bam | fgrep -v -f ${removed} > bam.RNA.lst
ls ${CCWD}/S*/WHOLE_GENOME/bam/S*.stat.txt | fgrep -v -f ${removed}  > stat.WGS.lst

Rscript ${QCbam} bam.WGS.lst qc.bam.WGS.lst
Rscript ${QCbam} bam.RNA.lst qc.bam.RNA.lst

ls ${CCWD}/S*/WHOLE_GENOME/bam/*.flagstat.txt | fgrep -v -f ${removed} | awk -F'/' '{print $NF "\t" $0}'  | sed 's/.flagstat.txt//' > flagstat.WGS.lst
ls ${CCWD}/S*/TRANSCRIPTOME/bam/*.flagstat.txt | fgrep -v -f ${removed} | awk -F'/' '{print $NF "\t" $0}' | sed 's/.flagstat.txt//' > flagstat.RNA.lst

ls ${CCWD}/S*/WHOLE_GENOME/coverage-post/SJ*.codingexon.cumulative.txt | fgrep -v -f ${removed} | awk -F'/' '{print $NF "\t" $0}' | sed 's/.codingexon.cumulative.txt//' > codingexon.WGS.lst
ls ${CCWD}/S*/TRANSCRIPTOME/coverage-post/SJ*.codingexon.cumulative.txt | fgrep -v -f ${removed}| awk -F'/' '{print $NF "\t" $0}' | sed 's/.codingexon.cumulative.txt//'  > codingexon.RNA.lst

ls  ${CCWD}/S*/WHOLE_GENOME/snv-unpaired/SJ*high*out | fgrep -v -f ${removed} | awk -F'/' '{print $NF "\t" $0}' | sed 's/.high_20.out//' > high20.WGS.lst

ls ${CCWD}/S*/WHOLE_GENOME/coverage/SJ*.bw | fgrep -v -f ${removed} | awk -F'/' '{print $NF "\t" $0}' | sed 's/.bw//'  > bw.WGS.lst
ls ${CCWD}/S*/TRANSCRIPTOME/coverage/SJ*.bw | fgrep -v -f ${removed}| awk -F'/' '{print $NF "\t" $0}' | sed 's/.bw//' > bw.RNA.lst

ls ${CCWD}/S*/WHOLE_GENOME/crest-tumoronly/*predSV | fgrep -v -f ${removed} | awk -F'/' '{print $NF "\t" $0}' | sed 's/.predSV//' > D_predSV.WGS.lst

ls ${CCWD}/S*/WHOLE_GENOME/cnv-unpaired/SJ*txt | fgrep -v -f ${removed} | awk -F'/' '{print $NF "\t" $0}' | sed 's/_Single_CONSERTING_NoCREST_100.txt//' > cnvunpaired_txt.lst
ls ${CCWD}/S*/WHOLE_GENOME/cnv-unpaired/SJ*png | fgrep -v -f ${removed} | awk -F'/' '{print $NF "\t" $0}' | sed 's/_Single_CONSERTING_NoCREST_100.png//' > cnvunpaired_png.lst


cat bam.WGS.lst | awk -F'/' '{print $NF "\t" $0}' | sed 's/.bam//' | sed 's/_/\t/' | awk '{print $1 "_" $2 "\t" $3 "\t" $1 "\t" $2 "\t" $2}' > NameMap_WGS.txt

rm readgroup.txt
touch readgroup.txt
while read line
do
  samtools view -H ${line} | grep ^@RG   | grep SM: >> readgroup.txt
done < bam.WGS.lst


