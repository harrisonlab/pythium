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
qsub $ProgDir/submit_dipSPAdes.sh $F_Read $R_Read $OutDir correct 10
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


Then I need to compare the diploid and non diploid genome stats: the numbers are the same, need to re-run dip one again

For now carry on with rest: so running repeatmasking for all of then in assembly/spades which includes nondip intermedium assembly
Once I've re-run the diploid assembly for intermedium then can run quast etc for that.

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

Take around 24 hours??

