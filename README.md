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
ProgDir=~/halesk/git_repos/tools/seq_tools/dna_qc;
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
ProgDir= /home/halesk/git_repos/emr_repos/tools/seq_tools/dna_qc
qsub $ProgDir/kmc_kmer_counting.sh $Trim_F $Trim_R
done
```

** Estimated Genome Size is: 48490129

** Esimated Coverage is: 42

#Assembly
Assembly was performed using: Spades


This is what I submitted for assembly

```bash

for StrainPath in $(ls -d qc_dna/paired/*/*); do 
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

#Remove contaminants

```bash
for OutDir in $(ls -d assembly/spades/P.violae/HL/filtered_contigs); do
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

#Quast

```bash
ProgDir=/home/halesk/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/*/HL/filtered_contigs/*_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

These report stats can be found in:
less assembly/spades/P.violae/DE/filtered_contigs

#Repeat Masking
Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

```bash
ProgDir=/home/halesk/git_repos/tools/seq_tools/repeat_masking
for BestAss in $(ls assembly/spades/*/HL/filtered_contigs/*_500bp_renamed.fasta); do
echo $BestAss
qsub $ProgDir/rep_modeling.sh $BestAss
qsub $ProgDir/transposonPSI.sh $BestAss
done

If only want to run either HL or DE then change second * to HL or DE
Take around 24 hours


#Pre-gene prediction
Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

```bash
ProgDir=/home/halesk/git_repos/tools/gene_prediction/cegma
for Assembly in $(ls assembly/spades/*/HL/filtered_contigs/*_500bp_renamed.fasta); do
qsub $ProgDir/sub_cegma.sh $Assembly dna
done
```
Results found in Pythium/gene_pred/cegma  - 

###Gene prediction 1  Augustus

Did 2 for AUgustus - 1st with Fusarium - line GeneModel=Fusarium
then second with P.cactorum

1st line for loop copied from above
2/3/4 copied from above, cuts out words from file directory to get strain and organism. 
Then makes out directory
5 where programme is
6 using this data to train it
7 submit programme, the programme, the thing used to train, the thing want to do programme on, 
false = prevent a gene in the opposite orientation being predicted, then where going to


```bash
for Assembly in $(ls assembly/spades/*/HL/filtered_contigs/*_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
OutDir=gene_pred/augustus/$Organism/$Strain/filtered_contigs
ProgDir=/home/halesk/git_repos/tools/gene_prediction/augustus
GeneModel=P.cactorum_10300_masked_braker
echo "ProgDir/submit_augustus.sh $GeneModel $Assembly false $OutDir"
qsub $ProgDir/submit_augustus.sh $GeneModel $Assembly false $OutDir
done
```

###Gene prediction 2 - atg.pl prediction of ORFs

Open reading frame predictions were made using the run_ORF_finder.sh script. This pipeline also 
identifies open reading frames containing Signal peptide sequences and RxLRs. 
This pipeline was run with the following commands:

```bash
ProgDir=/home/halesk/git_repos/tools/gene_prediction/ORF_finder
for Assembly in $(ls assembly/spades/*/HL/filtered_contigs/*_500bp_renamed.fasta); do
qsub $ProgDir/run_ORF_finder.sh $Assembly
done
```


The Gff files from the the ORF finder are not in true Gff3 format. These were corrected using the following commands:

```bash
for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*_ORF.gff); do
Strain=$(echo $ORF_Gff | rev | cut -f2 -d '/' | rev)
Organism=$(echo $ORF_Gff | rev | cut -f3 -d '/' | rev)
ProgDir=~/git_repos/tools/seq_tools/feature_annotation
ORF_Gff_mod=gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF_corrected.gff3
$ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
done
```

run this once got output from ORF above.


###Genomic analysis

#Genes with homology to PHIbase

Predicted gene models were searched against the PHIbase database using tBLASTx.

```bash
ProgDir=/home/halesk/git_repos/tools/pathogen/blast
Query=../../phibase/v3.8/PHI_accessions.fa
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/*_500bp_renamed.fasta); do
qsub $ProgDir/blast_pipe.sh $Query protein $Assembly
done
```

    #Interproscan
    Interproscan was used to give gene models functional annotations
```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/
for Genes in $(ls gene_pred/augustus/*/*/filtered_contigs/*_EMR_singlestrand_aug_out.aa); do
echo $Genes; cat $Genes | grep '>' |wc -l
$ProgDir/sub_interproscan.sh $Genes
done
```
    
    #Signal peptide prediction
    
    
    ### A) From Augustus gene models - Signal peptide

Required programs:
 * SigP
 * biopython


Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
for Proteome in $(ls gene_pred/augustus/*/*/filtered_contigs/*_EMR_singlestrand_aug_out.aa); do
SplitfileDir=/home/halesk/git_repos/tools/seq_tools/feature_annotation/signal_peptides
ProgDir=/home/halesk/git_repos/tools/seq_tools/feature_annotation/signal_peptides
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
SplitDir=gene_pred/Augustus_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_Augustus_preds
$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_Augustus_preds_*); do
Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
while [ $Jobs -ge 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
done
printf "\n"
echo $File
#qsub $ProgDir/pred_sigP.sh $File
qsub $ProgDir/pred_sigP.sh $File signalp-4.1
done
done
```

    
    
Below: what was on Readme before that Andrew did:
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
