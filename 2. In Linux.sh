fasterq-dump -3 'sra_run_id'  
  # parallel fasterq-dump -3 :::: 'srarunidlist.txt'
kallisto index -i 'ref.idx' 'ref.fasta' # ref.fasta = reference genome
mkdir kallisto_quantification
for i in $(cat srarunidlist.txt)
  do kallisto quant -i 'ref.idx' -o kallisto_quantification/${i}/ fastq/${i}_1.fastq fastq/${i}_2.fastq
done


 # creating ref.cdna.csv 
zcat 'ref.cdna.all.fa.gz' | grep '>'| awk '{FS= " "}BEGIN{ print "TXNAME,GENEID"};{print substr($1,2) "," substr($4,6)};' > 'ref.cdna.csv'
