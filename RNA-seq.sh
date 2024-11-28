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
done




############5.Assamble Transcripts for Each Sample############

# 设置要处理的 BAM 文件夹路径
bam_folder="/data01/fanly/03_rnaseq/liver/02.alignment/sorted_bam"

# 循环处理每个 BAM 文件
input_dir="/data01/fanly/03_rnaseq/liver/02.alignment/sorted_bam"
output_dir="/data01/fanly/03_rnaseq/liver/03.transcript_quant/1.assemble_transcript"
cd /data01/fanly/03_rnaseq/liver/02.alignment/sorted_bam

# 循环处理bam文件
for bam_file in ${input_dir}/*.bam; do
        # 提取bam——file前面的名字
        filename=$(basename "$bam_file" .bam)
        # 执行stringtie
        stringtie -p 8 \
		-G /data01/fanly/03_rnaseq/liver/00.ref/qdf_annotation_chr.gtf \
		-o "$output_dir/${filename}.gtf" \
		-l "${filename}" \   #指定转录本的标签，通常与文件名相同
		"$input_dir/${filename}.bam"
done

# Merge Transcripts from All Samples
stringtie \
--merge -p 16 \
-G /data01/fanly/03_rnaseq/liver/00.ref/qdf_annotation_chr.gtf \
-o /data01/fanly/03_rnaseq/liver/03.transcript_quant/1.assemble_transcript/stringtie_merged.gtf \
/data01/fanly/03_rnaseq/liver/03.transcript_quant/1.assemble_transcript/stringtie_mergelist.txt

############6.Quantify############

# 定义输入和输出目录
input_dir="/data01/fanly/03_rnaseq/liver/03.transcript_quant/1.assemble_transcript"
output_dir="/data01/fanly/03_rnaseq/liver/03.transcript_quant/2.quant_transcript"
bam_dir="/data01/fanly/03_rnaseq/liver/02.alignment/sorted_bam/"

# 循环处理每个 GTF 文件
for gtf_file in "$input_dir"/*.gtf; do
    # 提取文件名
    filename=$(basename "${gtf_file}" .gtf)

    # 执行 StringTie
    stringtie -e -B -p 16 \
        -G "${input_dir}/stringtie_merged.gtf" \
        -o "${output_dir}/${filename}_quant.gtf" \
        "${bam_dir}/${filename}.bam"
done

# 获取raw_counts
# samplelist.txt里面是每个样本的id和文件的完整路径（tab分隔）；id是每一列的表头
# 最后得到gene_count_matrix.csv、transcript_count_matrix.csv
python prepDE.py -i samplelist.txt
# 提取TPM和FPKM可以使用修改后的脚本，最后有用的文件是transcript_count_matrix.csv

