fasterq-dump -3 'sra_run_id'  
  # parallel fasterq-dump -3 :::: 'srarunidlist.txt'
kallisto index -i 'ref.idx' 'ref.fasta' # ref.fasta = reference genome
mkdir salmon

cd salmon

for i in $(cat srarunidlist.txt)  #for making each run id's directory
  do mkdir ${i}
done

#making idx
salmon index -t 'ref.fa.gz' -i ref_index

# paired end
for i in $(cat srarunidlist.txt);
do
salmon quant -i ref_index -l A \
         -1 fastq/${i}_1.fastq.gz \
         -2 fastq/${i}_2.fastq.gz \
         -p 8 --validateMappings -o salmon/${i}
done
    # -p option for thread
