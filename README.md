New ReadMe for going genome assembly and annotation for Pythium isolates

# Pythium
Scripts and commands used for the analysis of Pythium spp. genomes



Commands used during analysis of the Pythium genome. 
Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/pythium

The following is a summary of the work presented in this Readme:
Data organisation:
  * Preparing data  
Draft Genome assembly
  * Data qc
  * Genome assembly
  
  
  #Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.


```bash
  cd /home/groups/harrisonlab/project_files/pythium
  mkdir -p raw_dna/paired/P.violae/HL/F
  mkdir -p raw_dna/paired/P.violae/HL/R
  mkdir -p raw_dna/paired/P.violae/DE/F
  mkdir -p raw_dna/paired/P.violae/DE/R
  mkdir -p raw_dna/paired/P.sulctaum/P67/F
  mkdir -p raw_dna/paired/P.sulctaum/P67/R
  mkdir -p raw_dna/paired/P.intermedium/P107/F
  mkdir -p raw_dna/paired/P.intermedium/P107/R

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

  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160503_M04465_0013_AMLD4
  cp $RawDat/PsulcaumP67_S1_L001_R1_001.fastq.gz raw_dna/paired/P.sulcaum/P67/F/.
  cp $RawDat/PsulcaumP67_S1_L001_R2_001.fastq.gz raw_dna/paired/P.sulcaum/P67/R/.
  cp $RawDat/PintermediumP107_S2_L001_R1_001.fastq.gz raw_dna/paired/P.intermedium/P107/F/.
  cp $RawDat/PintermediumP107_S2_L001_R2_001.fastq.gz raw_dna/paired/P.intermedium/P107/R/.
```

#Data qc
DNA quality checking before trimming. 

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

Echo is just to show/check which data you have called on in the first line (the for loop)
No output directory is set so it will go into the same folder it came out of (A new directory has been created in
/raw_dna/paired/P.violae/HL/R (for example) called  HL-Pviolae_S1_L001_R2_001_fastqc which contains the fastqc files)


```bash

for RawData in $(ls raw_dna/paired/P.*/*/*/*.fastq.gz); do
echo $RawData;
ProgDir=/home/halesk/git_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
done
```
Visualise by navigating to the file in finder (after cluster_mount) and opening the file from finder (/groups, harrisonlab, project files etc.)

(We have visualised (open html file), don't know what it means.)



Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf

Trimming was first performed on all strains that had a single run of data:

```bash
for Strain in $(ls -d raw_dna/paired/*/* | grep -e 'P107' -e 'P67'); do
echo $Strain
IluminaAdapters=/home/halesk/git_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/halesk/git_repos/tools/seq_tools/rna_qc
ReadsF=$(ls $Strain/F/*.fastq*)
ReadsR=$(ls $Strain/R/*.fastq*)
echo $ReadsF
echo $ReadsR
qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IluminaAdapters DNA
done
```


Trimming was then performed on strains that had 2 runs of data:

```bash
for Strain in $(ls -d raw_dna/paired/P.violae/*); do
echo $Strain
IluminaAdapters=/home/halesk/git_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/halesk/git_repos/tools/seq_tools/rna_qc
Read_F=$(ls $Strain/F/*.fastq.gz | grep -v 'run2')
Read_R=$(ls $Strain/R/*.fastq.gz | grep -v 'run2')
echo $Read_F
echo $Read_R
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
Read_F=$(ls $Strain/F/*.fastq.gz | grep 'run2')
Read_R=$(ls $Strain/R/*.fastq.gz | grep 'run2')
echo $Read_F
echo $Read_R
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
```

This goes into a folder in Pythium called qc_dna. The end of the file is trim.fq.gz. 
Have 2 F and 2 R each for HL and DE and 1 F and 1 R each for Int/Sul.


Data quality was visualised once again following trimming:

```bash
for RawData in $(ls qc_dna/paired/P.*/*/*/*.fq.gz); do
echo $RawData;
ProgDir=/home/halesk/git_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
done
```

No output directory is set so it will go into the same folder it came out of (A new directory has been created in
/qc_dna/paired/P.violae/HL/R (for example) called  HL-Pviolae_S1_L001_R2_001_fastqc which contains the fastqc files)

Visualise by navigating to the file in finder (after cluster_mount) and opening the file from finder (/groups, harrisonlab, project files etc.)

(We have visualised (open html file), don't know what it means.)




kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
for StrainPath in $(ls -d qc_dna/paired/P.*/*); do
echo $StrainPath
Trim_F=$(ls $StrainPath/F/*.fq.gz)
Trim_R=$(ls $StrainPath/R/*.fq.gz)
ProgDir=/home/halesk/git_repos/tools/seq_tools/dna_qc
qsub $ProgDir/kmc_kmer_counting.sh $Trim_F $Trim_R
done
```
done the script above, take estimated genome size etc from true txt file


** Estimated Genome Size is: ?

** Esimated Coverage is: ?


#Assembly
Assembly was performed using: Spades

First submitted the isolates with 2 runs for assembly

```bash

for StrainPath in $(ls -d qc_dna/paired/P.violae/*); do
ProgDir=~/git_repos/tools/seq_tools/assemblers/spades/multiple_libraries;
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev);
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev);
F1_Read=$(ls $StrainPath/F/*.fq.gz | grep -v 'run2');
R1_Read=$(ls $StrainPath/R/*.fq.gz | grep -v 'run2');
F2_Read=$(ls $StrainPath/F/*.fq.gz | grep 'run2');
R2_Read=$(ls $StrainPath/R/*.fq.gz | grep 'run2');
OutDir=assembly/spades/$Organism/$Strain;
echo $F1_Read;
echo $R1_Read;
echo $F2_Read;
echo $R2_Read;
qsub $ProgDir/subSpades_2lib.sh $F1_Read $R1_Read $F2_Read $R2_Read $OutDir correct 10
done
```

Then submitted isolates with single runs for assembly:

for StrainPath in $(ls -d qc_dna/paired/*/* | grep -e 'P107' -e 'P67'); do
ProgDir=/home/halesk/git_repos/tools/seq_tools/assemblers/spades;
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev);
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev);
F_Read=$(ls $StrainPath/F/*.fq.gz);
R_Read=$(ls $StrainPath/R/*.fq.gz);
OutDir=assembly/spades/$Organism/$Strain;
echo $F_Read;
echo $R_Read;
qsub $ProgDir/submit_SPAdes.sh $F_Read $R_Read $OutDir correct 10
done

Then I submitted P.intermedium with the submit_dipSPAdes.sh

for StrainPath in $(ls -d qc_dna/paired/*/* | grep -e 'P107'); do
ProgDir=/home/halesk/git_repos/tools/seq_tools/assemblers/spades;
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev);
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev);
F_Read=$(ls $StrainPath/F/*.fq.gz);
R_Read=$(ls $StrainPath/R/*.fq.gz);
OutDir=assembly/spades/$Organism/$Strain;
echo $F_Read;
echo $R_Read;
qsub $ProgDir/submit_SPAdes.sh $F_Read $R_Read $OutDir correct 10
done


Assemblies run

I can't see 2 intermedium assembiles, did one overwrite the other one?

There is a folder called assembly/spades/P.*/*/filtered_contigs and in there is a fasta file.


To generate assembly stats need to run remove contaminants and then quast


#Remove contaminants

```bash
for OutDir in $(ls -d assembly/spades/P.*/*/filtered_contigs); do
echo "$OutDir"
ProgDir=/home/halesk/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
AssFiltered=$OutDir/contigs_min_500bp.fasta
AssRenamed=$OutDir/contigs_min_500bp_renamed.fasta
ls $AssFiltered
printf '.\t.\t.\t.\n' > editfile.tab
$ProgDir/remove_contaminants.py --inp $AssFiltered --out $AssRenamed --coord_file editfile.tab
rm editfile.tab
done
```

Ran the script above, in the same filtered contigs directory as above have another fasta file called renamed.

Then run Quast, this generates assembly stats

#Quast

```bash
ProgDir=/home/halesk/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/*_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```
These report stats can be found in:
less assembly/spades/P.violae/DE/filtered_contigs/report.txt


Didn't have the normal assembly for intermedium, the dipspades most likely wrote over it.

So re-ran with the script below:

for StrainPath in $(ls -d qc_dna/paired/*/* | grep -e 'P107'); do
ProgDir=/home/halesk/git_repos/tools/seq_tools/assemblers/spades;
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev);
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev);
F_Read=$(ls $StrainPath/F/*.fq.gz);
R_Read=$(ls $StrainPath/R/*.fq.gz);
OutDir=assembly/spades/$Organism/$Strain;
echo $F_Read;
echo $R_Read;
qsub $ProgDir/submit_SPAdes.sh $F_Read $R_Read $OutDir correct 10
done

Then re-ran remove contaminants:

```bash
for OutDir in $(ls -d assembly/spades/P.intermedium/*/filtered_contigs); do
echo "$OutDir"
ProgDir=/home/halesk/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
AssFiltered=$OutDir/contigs_min_500bp.fasta
AssRenamed=$OutDir/contigs_min_500bp_renamed.fasta
ls $AssFiltered
printf '.\t.\t.\t.\n' > editfile.tab
$ProgDir/remove_contaminants.py --inp $AssFiltered --out $AssRenamed --coord_file editfile.tab
rm editfile.tab
done
```

Then re-ran quast:

ProgDir=/home/halesk/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/P.intermedium/*/filtered_contigs/*_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done


Then I need to compare the diploid and non diploid genome stats: BUT THE NUMBERS for the assembly I thought was diploid (becuase it was
the last one I ran) was the same as the non diploid I just ran, so either both give exactly the same, or the diploid one didn;t
run first time. I need to re-run diploid one again to check.

For now carry on with rest: so running repeatmasking for all of then in assembly/spades which includes nondip intermedium assembly
Once I've re-run the diploid assembly for intermedium then can run quast etc for that.

Then I submitted P.intermedium with the submit_dipSPAdes.sh I changed the first ls line to take out the grep and send it straight to P.inermedium

for StrainPath in $(ls -d qc_dna/paired/P.intermedium/P107); do
ProgDir=/home/halesk/git_repos/tools/seq_tools/assemblers/spades;
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev);
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev);
F_Read=$(ls $StrainPath/F/*.fq.gz);
R_Read=$(ls $StrainPath/R/*.fq.gz);
OutDir=assembly/spades/$Organism/$Strain;
echo $F_Read;
echo $R_Read;
qsub $ProgDir/submit_dipSPAdes.sh $F_Read $R_Read $OutDir correct 10
done


#Repeat Masking
Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

```bash
ProgDir=/home/halesk/git_repos/tools/seq_tools/repeat_masking
for BestAss in $(ls assembly/spades/*/*/filtered_contigs/*_500bp_renamed.fasta); do
echo $BestAss
qsub $ProgDir/rep_modeling.sh $BestAss
qsub $ProgDir/transposonPSI.sh $BestAss
done
```

This ran, results in ls repeat_masked/P.violae/DE/filtered_contigs_repmask/
Have softmasked and hardmasked.



The number of bases masked by transposonPSI and Repeatmasker were summarised using the following commands:

for RepDir in $(ls -d repeat_masked/P.*/*/filtered_contigs_repmask); do
Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
printf "$Organism\t$Strain\n"
printf "The number of bases masked by RepeatMasker:\t"
sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
printf "The number of bases masked by TransposonPSI:\t"
sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
printf "The total number of masked bases are:\t"
cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
echo
done
Script above gives the following:

P.intermedium	P107
The number of bases masked by RepeatMasker:	9677535
The number of bases masked by TransposonPSI:	2768624
The total number of masked bases are:	10154407

P.sulcatum	P67
The number of bases masked by RepeatMasker:	10693462
The number of bases masked by TransposonPSI:	3545514
The total number of masked bases are:	11801684

P.violae	DE
The number of bases masked by RepeatMasker:	6898545
The number of bases masked by TransposonPSI:	1291867
The total number of masked bases are:	7207639

P.violae	HL
The number of bases masked by RepeatMasker:	7114132
The number of bases masked by TransposonPSI:	1251306
The total number of masked bases are:	7429753


#Pre-gene prediction
Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
This script is for the 2 P.violae's and P.sulcatum as at the time of running the P.intermedium diploid assembly that will go into this directory hadn't run yet

```bash
ProgDir=/home/halesk/git_repos/tools/gene_prediction/cegma
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/*_500bp_renamed.fasta); do
qsub $ProgDir/sub_cegma.sh $Assembly dna
done
```
Results found in Pythium/gene_pred/cegma  -

The script below is for the non diploid P.intermedium assembly. 

```bash
ProgDir=/home/halesk/git_repos/tools/gene_prediction/cegma
for Assembly in $(ls assembly_nondip/P.intermedium/P107/filtered_contigs/*_500bp_renamed.fasta); do
qsub $ProgDir/sub_cegma.sh $Assembly dna
done
```
Outputs were summarised using the commands:

for File in $(ls gene_pred/cegma/*/*/*_dna_cegma.completeness_report); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Species=$(echo $File | rev | cut -f3 -d '/' | rev)
printf "$Species\t$Strain\n"
cat $File | head -n18 | tail -n+4;printf "\n"
done >> gene_pred/cegma/cegma_results_dna_summary.txt
less gene_pred/cegma/cegma_results_dna_summary.txt


Results printed below:
P.intermedium	P107
              #Prots  %Completeness  -  #Total  Average  %Ortho 

  Complete      217       87.50      -   313     1.44     32.26

   Group 1       52       78.79      -    82     1.58     38.46
   Group 2       51       91.07      -    71     1.39     29.41
   Group 3       54       88.52      -    76     1.41     31.48
   Group 4       60       92.31      -    84     1.40     30.00

   Partial      235       94.76      -   364     1.55     37.45

   Group 1       56       84.85      -    93     1.66     42.86
   Group 2       56      100.00      -    82     1.46     33.93
   Group 3       59       96.72      -    92     1.56     37.29
   Group 4       64       98.46      -    97     1.52     35.94

P.sulcatum	P67
              #Prots  %Completeness  -  #Total  Average  %Ortho 

  Complete      189       76.21      -   261     1.38     26.98

   Group 1       50       75.76      -    76     1.52     36.00
   Group 2       43       76.79      -    55     1.28     20.93
   Group 3       45       73.77      -    66     1.47     35.56
   Group 4       51       78.46      -    64     1.25     15.69

   Partial      232       93.55      -   381     1.64     42.67

   Group 1       59       89.39      -    97     1.64     42.37
   Group 2       52       92.86      -    87     1.67     42.31
   Group 3       60       98.36      -    98     1.63     46.67
   Group 4       61       93.85      -    99     1.62     39.34

P.violae	DE
              #Prots  %Completeness  -  #Total  Average  %Ortho 

  Complete      233       93.95      -   343     1.47     32.19

   Group 1       60       90.91      -    79     1.32     25.00
   Group 2       51       91.07      -    79     1.55     37.25
   Group 3       57       93.44      -    84     1.47     33.33
   Group 4       65      100.00      -   101     1.55     33.85

   Partial      239       96.37      -   371     1.55     36.40

   Group 1       61       92.42      -    85     1.39     26.23
   Group 2       53       94.64      -    85     1.60     41.51
   Group 3       60       98.36      -    93     1.55     36.67
   Group 4       65      100.00      -   108     1.66     41.54

P.violae	HL
              #Prots  %Completeness  -  #Total  Average  %Ortho 

  Complete      232       93.55      -   336     1.45     31.03

   Group 1       59       89.39      -    79     1.34     28.81
   Group 2       51       91.07      -    73     1.43     29.41
   Group 3       57       93.44      -    84     1.47     33.33
   Group 4       65      100.00      -   100     1.54     32.31

   Partial      238       95.97      -   362     1.52     34.87

   Group 1       60       90.91      -    84     1.40     28.33
   Group 2       53       94.64      -    79     1.49     33.96
   Group 3       60       98.36      -    92     1.53     36.67
   Group 4       65      100.00      -   107     1.65     40.00
   
   
Still need to run a diploid P.intermedium one when that assembly has run for EVERYTHING from assembly

AT THIS POINT NEED TO DECIDE WHICH ASSEMBLY IS BEST FOR P.INTERMEIDUM AND CARY ON WITH ONLY 1 ASSEMBLY



###Gene Prediction
USING BRAKER

#1 Data organisation

Downloading RNA seq data from NCBI

These commands were used to download SRA files from NCBI
RNA seq data associated with P. ultimatum under X experimental conditions was
downloaded. This data was published in X.

Make a directory first if necessary with species/strain.

```bash
  mkdir -p raw_rna/external/P.ultimum/DAOM_BR144
  fastq-dump -O raw_rna/external/P.ultimum/DAOM_BR144/unpaired SRR058978
  fastq-dump -O raw_rna/external/P.ultimum/DAOM_BR144/unpaired SRR058979
  fastq-dump -O raw_rna/external/P.ultimum/DAOM_BR144/unpaired SRR058980
  fastq-dump -O raw_rna/external/P.ultimum/DAOM_BR144/unpaired SRR059020
  fastq-dump -O raw_rna/external/P.ultimum/DAOM_BR144/unpaired SRR059021
  fastq-dump -O raw_rna/external/P.ultimum/DAOM_BR144/unpaired SRR059022
  fastq-dump -O raw_rna/external/P.ultimum/DAOM_BR144/unpaired SRR059025
  fastq-dump -O raw_rna/external/P.ultimum/DAOM_BR144/unpaired SRR059026
  fastq-dump -O raw_rna/external/P.ultimum/DAOM_BR144/unpaired SRR059027
  fastq-dump -O raw_rna/external/P.ultimum/DAOM_BR144/unpaired SRR236127
  fastq-dump -O raw_rna/external/P.ultimum/DAOM_BR144/unpaired SRR236128
  fastq-dump -O raw_rna/external/P.ultimum/DAOM_BR144/unpaired SRR236130
```
-O is choosing where to put it (Output)
Final number is the accession number from NCBI, under Run


Downloaded data was in an unzipped format. Therefore files were zipped to
reduce the amount they took on the cluster.

```bash
  gzip raw_rna/external/P.ultimum/DAOM_BR144/unpaired/*.fastq
```
Done to here

###2 QC

RNAseq data was trimmed to remove low quality sequence and adapters. As the
rnaseq data represents unpaired libraries, the fastq-mcf wrapper designed for
unpaired reads was used.

 Script below untested!

```bash
    StrainPath=$(raw_rna_external/external/*/*)
    for RnaDat in $(ls $StrainPath/unpaired/*.fastq.gz); do
      IlluminaAdapters=/home/halesk/git_repos/tools/seq_tools/ncbi_adapters.fa
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
      echo $RnaDat
      qsub $ProgDir/dna_qc_fastq-mcf_unpaired.sh $RnaDat $IlluminaAdapters RNA
    done
```
