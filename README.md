# pythium
Scripts and commands used for the analysis of Pythium spp. genomes



commands used during analysis of the pythium genome. Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/pythium

The following is a summary of the work presented in this Readme:
Data organisation:
  * Preparing data  
Draft Genome assembly
  * Data qc
  * Genome assembly
  * Repeatmasking
  * Gene prediction
  * Functional annotation
Genome analysis
  * Homology between predicted genes & published effectors


#Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.


```bash
  cd /home/groups/harrisonlab/project_files/pythium
  mkdir -p raw_dna/paired/P.violae/HL/F
  mkdir -p raw_dna/paired/P.violae/HL/R
  mkdir -p raw_dna/paired/P.violae/DE/F
  mkdir -p raw_dna/paired/P.violae/DE/R

  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/pythium/150910_M01678_0026_D0HRG
  cp $RawDat/HL-Pviolae_S1_L001_R1_001.fastq.gz raw_dna/paired/P.violae/HL/F/.
  cp $RawDat/HL-Pviolae_S1_L001_R2_001.fastq.gz raw_dna/paired/P.violae/HL/R/.
  cp $RawDat/DE-Pviolae_S2_L001_R1_001.fastq.gz raw_dna/paired/P.violae/DE/F/.
  cp $RawDat/DE-Pviolae_S2_L001_R2_001.fastq.gz raw_dna/paired/P.violae/DE/R/.


  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/pythium/150911_M01678_0027_A5EK9
  cp $RawDat/HL-Pviolae_S1_L001_R1_001.fastq.gz raw_dna/paired/P.violae/HL/F/HL-Pviolae_run2_S1_L001_R1_001.fastq.gz
  cp $RawDat/HL-Pviolae_S1_L001_R2_001.fastq.gz raw_dna/paired/P.violae/HL/R/HL-Pviolae_run2_S1_L001_R2_001.fastq.gz
  cp $RawDat/DE-Pviolae_S2_L001_R1_001.fastq.gz raw_dna/paired/P.violae/DE/F/DE-Pviolae_run2_S2_L001_R1_001.fastq.gz
  cp $RawDat/DE-Pviolae_S2_L001_R2_001.fastq.gz raw_dna/paired/P.violae/DE/R/DE-Pviolae_run2_S2_L001_R2_001.fastq.gz
```


#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash

for RawData in $(ls raw_dna/paired/P.*/*/*/*.fastq.gz); do
echo $RawData;
ProgDir=/home/halesk/git_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
done
```



Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf


```bash
for Strain in $(ls raw_dna/paired/P.*/); do
echo $Strain
IluminaAdapters=/home/halesk/git_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/halesk/git_repos/tools/seq_tools/rna_qc
Read_F=$(ls raw_dna/paired/P.*/$Strain/F/*.fastq.gz | grep -v 'run2')
Read_R=$(ls raw_dna/paired/P.*/$Strain/R/*.fastq.gz | grep -v 'run2')
echo $Read_F
echo $Read_R
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
Read_F=$(ls raw_dna/paired/P.*/$Strain/F/*.fastq.gz | grep 'run2')
Read_R=$(ls raw_dna/paired/P.*/$Strain/R/*.fastq.gz | grep 'run2')
echo $Read_F
echo $Read_R
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
```

Data quality was visualised once again following trimming:

```bash
for RawData in $(ls qc_dna/paired/P.*/*/*/*.fq.gz); do
echo $RawData;
ProgDir=~/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
done
```



kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
for StrainPath in $(ls -d qc_dna/paired/P.*/*); do
echo $StrainPath
Trim_F=$(ls $StrainPath/F/*.fq.gz)
Trim_R=$(ls $StrainPath/R/*.fq.gz)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/dna_qc
qsub $ProgDir/kmc_kmer_counting.sh $Trim_F $Trim_R
done
```

** Estimated Genome Size is: 48490129

** Esimated Coverage is: 42

#Assembly
Assembly was performed using: Velvet / Abyss / Spades

A range of hash lengths were used and the best assembly selected for subsequent analysis


```bash
  for StrainPath in $(ls -d qc_dna/paired/P.*/*); do
    echo $StrainPath
    Trim_F=$(ls $StrainPath/F/*.fq.gz | grep 'run2')
    Trim_R=$(ls $StrainPath/R/*.fq.gz | grep 'run2')
    echo $Trim_F
    echo $Trim_R
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    OutDir=assembly/spades/P.violae/$Strain
    qsub $ProgDir/submit_SPAdes.sh $F_Read $R_Read $OutDir correct 10
  done
```
<!--
Quast

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
	Assembly=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp.fasta
	OutDir=assembly/spades/N.ditissima/R0905_v2/filtered_contigs
	qsub $ProgDir/sub_quast.sh $Assembly $OutDir
```

Assemblies were summarised to allow the best assembly to be determined by eye.

** Assembly stats are:
  * Assembly size:
  * N50:153669
  * N80:
  * N20:
  * Longest contig:687738
  **

As SPADes was run with the option to autodetect a minimum coverage the assembly was assessed to identify the coverage of assembled contigs. This was done using the following command:

```bash
	BestAss=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp.fasta
	cat $BestAss | grep '>' | cut -f6 -d'_' | sort -n | cut -f1 -d '.' | sort -n | uniq -c | less
```

From this it was determined that SPades could not be trusted to set its own minimum threshold for coverage.
In future an option will be be used to set a coverage for spades.
In the meantime contigs with a coverage lower than 10 were filtered out using the following commands:

```bash
	Headers=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_headers.txt
	cat $BestAss | grep '>' | grep -E -v 'cov_.\..*_' > $Headers
	FastaMinCov=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_headers.fasta
	cat $BestAss | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep -A1 -f $Headers | grep -v -E '^\-\-' > $FastaMinCov
```

```bash
	~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants/remove_contaminants.py --inp ../neonectria_ditissima/assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_headers.fasta  --out assembly/spades/N.galligena/R0905_v2/filtered_contigs/contigs_min_500bp_10x_filtered_renamed.fasta  --coord_file editfile.tab

We run Quast again.

	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
	Assembly=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_headers.fasta
	OutDir=assembly/spades/N.ditissima/R0905_v2/contigs_min_500bp_10x_headers
	qsub $ProgDir/sub_quast.sh $Assembly $OutDir -->
