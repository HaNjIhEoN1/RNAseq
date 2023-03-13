# RNAseq
using Kallisto, salmon DESeq2

[kallisto manual](https://pachterlab.github.io/kallisto/manual)

[salmon manual](https://combine-lab.github.io/salmon/getting_started/)

## Index
Ⅰ. Preparation process

Ⅱ. In Linux
 
Ⅲ. In R - to DESeq2

Ⅳ. In R - visualize

### Ⅰ. Preparation process
  1. SRA-Toolkit download – https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump 

    1.1. mkdir download -> cd download
    1.2. wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.2/sratoolkit.3.0.2-ubuntu64.tar.gz
    1.3. tar –zxvf sratoolkit.3.0.2-ubuntu64.tar.gz	
    1.4. cd sratoolkit.3.0.2-ubuntu64/bin 
    1.5. pwd 
    1.6. export PATH=$PATH:'$pwd'
      # convert to environment variable
      case1 echo "export PATH=$PATH:$pwd" >> ~/.bashrc >> source ~/.bashrc
      case2 cd  >> vi .bashrc >> add 'export PATH=$PATH:$pwd' in the last line

  2. RNAseq program download
  
    2.1 conda install -c bioconda kallisto
    2.2 conda install -c bioconda salmon

  3. Creating Metadata
  
    3.1 Go to NCBI
    3.2 Search 'Oryza sativa' and filtering (Libarylayout = Paired, LibrarySelection = cDNA, Instrument = Illuminihiseq/ Illumine nova seq 6000/ next seq 500, Organism = Oryza sativa Japonica...etc)
    3.3 click 'send to' > select 'run selector' click 'go' > download(metadata) > open as xlsx file > data select

### Ⅱ. [In Linux](https://github.com/HaNjIhEoN1/RNAseq/blob/main/2.%20In%20Linux.sh)

### Ⅲ. [In R - to DESeq2](https://github.com/HaNjIhEoN1/RNAseq/blob/main/3.%20In%20R%20to%20DESeq2.R)
