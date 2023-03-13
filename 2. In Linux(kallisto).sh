fasterq-dump -3 'sra_run_id'  
  # parallel fasterq-dump -3 :::: 'srarunidlist.txt'
kallisto index -i 'ref.idx' 'ref.fasta' # ref.fasta = reference genome
mkdir kallisto_quantification

for i in $(cat srarunidlist.txt)  #for making each run id's directory
  do mkdir ${i}
done

#case - paired end

for i in $(cat srarunidlist.txt)
  do kallisto quant -i 'ref.idx' -o kallisto_quantification/${i}/ fastq/${i}_1.fastq fastq/${i}_2.fastq
done

#case - single end

## calculate estimate framgent length, sd
awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("total %d avg=%f stddev=%f\n",n,m,sq/n-m*m);}'  *.fastq

##
for i in $(cat srarunidlist.txt)
  do kallisto quant -i 'ref.idx' -o kallisto_quantification/${i}/ --single -l 200 -s 20 fastq/${i}.fastq
done

 # creating ref.cdna.csv 
zcat 'ref.cdna.all.fa.gz' | grep '>'| awk '{FS= " "}BEGIN{ print "TXNAME,GENEID"};{print substr($1,2) "," substr($4,6)};' > 'ref.cdna.csv'
