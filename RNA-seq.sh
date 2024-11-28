cd /data01/wangyf/project2/CyprinusCarpio/15.pop/19.supplementary/3.rna-seq

:<<!
source activate sratoolkit

for d in */; do
    dir="${d%/}"
  # fasterq-dump --threads 30 --split-3 "$dir"/"$dir.sra" -O ./"$dir"
    gzip "$dir"/"$dir"_1.fastq && gzip "$dir"/"$dir"_2.fastq
done
!

genome=/data01/wangyf/project2/CyprinusCarpio/1.data/genome/genome.fa
reseq=/data01/wangyf/project2/CyprinusCarpio/15.pop/19.supplementary/3.rna-seq/1.data
core=8

source activate gatk4_pipeline

for ff in `ls $reseq/*_1.fq.gz | head -17 |tail -1`
do

  f=`basename $ff _1.fq.gz`

 if [[ ! -f ./2.fastqc/"$f"_1_fastqc.zip && ! -f ./2.fastqc/"$f"_2_fastqc.zip  ]];then
  fastqc $reseq/"$f"_1.fq.gz -o ./2.fastqc -t $core
  fastqc $reseq/"$f"_2.fq.gz -o ./2.fastqc -t $core
 fi

 if [[ ! -f ./3.fastp/"$f".json && ! -f ./3.fastp/"$f".html ]];then
  fastp -i $reseq/"$f"_1.fq.gz -I $reseq/"$f"_2.fq.gz -o 3.fastp/"$f"_1.fq.fastp.gz -O 3.fastp/"$f"_2.fq.fastp.gz --json 3.fastp/"$f".json --html  3.fastp/"$f".html -w $core
  fastqc 3.fastp/"$f"_1.fq.fastp.gz -o 3.fastp/ -t $core
  fastqc 3.fastp/"$f"_2.fq.fastp.gz -o 3.fastp/ -t $core
 fi

 if [[ ! -f ./4.mapping/"$f".sam &&  -f ./4.mapping/"$f".sort.bam ]];then
  echo "do nothing"
 elif [[  -f ./4.mapping/"$f".sam  && ! -f ./4.mapping/"$f".sort.bam ]];then
  bwa mem -M -t $core -v 3 -R "@RG\tID:$f\tSM:$f\tLB:$f\tPL:ILLUMINA" $genome 3.fastp/"$f"_1.fq.fastp.gz 3.fastp/"$f"_2.fq.fastp.gz  > ./4.mapping/"$f".sam
 elif [[ ! -f ./4.mapping/"$f".sam &&  ! -f ./4.mapping/"$f".sort.bam ]];then
  bwa mem -M -t $core -v 3 -R "@RG\\tID:$f\\tSM:$f\\tLB:$f\\tPL:ILLUMINA" $genome 3.fastp/"$f"_1.fq.fastp.gz 3.fastp/"$f"_2.fq.fastp.gz  > ./4.mapping/"$f".sam && \
  gatk SortSam -I ./4.mapping/"$f".sam -O ./4.mapping/"$f".sort.bam -SO coordinate --COMPRESSION_LEVEL 4 --TMP_DIR tmp && \
  rm ./4.mapping/"$f".sam
fi

  stringtie -p $core -G /data01/wangyf/project2/CyprinusCarpio/15.pop/7.annovar/new/genome/genome.gtf -o ./5.transcript_quant/"$f".gtf -l "$f"  ./4.mapping/"$f".sort.bam   # Assamble Transcripts for Each Sample
  stringtie --merge -p $core -G /data01/wangyf/project2/CyprinusCarpio/15.pop/7.annovar/new/genome/genome.gtf -o ./5.transcript_quant/stringtie_merged.gtf ./5.transcript_quant/stringtie_mergelist.txt  # 这一步是进入5.transcript_quant目录在本地跑的, Merge Transcripts from All Samples
  stringtie -e -B -p 8 -G /data01/wangyf/project2/CyprinusCarpio/15.pop/7.annovar/new/genome/genome.gtf -o ./5.transcript_quant/2.quant_transcript/"$f"/"$f"_quant.gtf ./4.mapping/"$f".sort.bam   # Quantify
  
done





# 获取raw_counts
# samplelist.txt里面是每个样本的id和文件的完整路径（tab分隔）；id是每一列的表头
# 最后得到gene_count_matrix.csv、transcript_count_matrix.csv
python prepDE.py -i samplelist.txt
# 提取TPM和FPKM可以使用修改后的脚本，最后有用的文件是transcript_count_matrix.csv

