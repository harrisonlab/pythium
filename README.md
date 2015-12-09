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

run this once got output from ORF above. - not done?


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

The above didn't work so trying submitting again with commands below


## 5.1 Identifying avirulence genes

Protein sequence of previously characterised oomycete avirulence genes was used in BLAST searches against
assemlies.

```bash
mkdir -p analysis/blast_homology/oomycete_avr_genes/
cp ../idris/analysis/blast_homology/oomycete_avr_genes/appended_oomycete_avr_cds.fasta analysis/blast_homology/oomycete_avr_genes/appended_oomycete_avr_cds.fasta 
ProgDir=/home/halesk/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/oomycete_avr_genes/appended_oomycete_avr_cds.fasta
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/*_500bp_renamed.fasta); do
echo $Assembly
qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
done
```

You will need to run the conversion of these files to gff annotations as was shown in the previous email


Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

2nd line file want to use (one we imported into excel)
3rd line is setting the output directory
4th line - what column called
5th - no of blast hits - do I want every blast hit annotated (some may be very low confidence) so choose what no of blast hits want to make

```bash
ProgDir=/home/halesk/git_repos/tools/pathogen/blast
for BlastHits in $(ls analysis/blast_homology/*/*/*_appended_oomycete_avr_cds.fasta_homologs.csv); do
HitsGff=$(echo $BlastHits | sed 's/csv/gff/g')
Column2=BLAST_homolog
NumHits=5
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
done
```



Following blasting PHIbase to the genome, the hits were filtered by effect on virulence.

The following commands were used to do this:

```bash
for BlastHits in $(ls analysis/blast_homology/*/*/*_PHI_accessions.fa_homologs.csv); do 
OutFile=$(echo $BlastHits | sed 's/.csv/_virulence.csv/g')
paste -d '\t' ../../phibase/v3.8/PHI_headers.csv ../../phibase/v3.8/PHI_virulence.csv $BlastHits | cut -f-3,1185- > $OutFile
cat $OutFile | grep 'contig' | cut -f2 | sort | uniq -c
done
```

This is getting files that are outside of our project directories (in harrison lab) called phibase. 
To turn into a loop  make a loop that captures csv file. 
For output file use sed command - replace .csv with virulence_.csv as earlier or an is written use cut to cut out organisms/strain

results should look like this once phibase run:

  1  
  3 chemistry target
  32 Chemistry target
  8  effector (plant avirulence determinant)
  13 Effector (plant avirulence determinant)
  2 Enhanced antagonism
  8  increased virulence
  5  increased virulence (Hypervirulence)
  2 Increased virulence (hypervirulence)
  21 Increased virulence (Hypervirulence)
  84 Lethal
  13  loss of pathogenicity
  237 Loss of pathogenicity
  9  mixed outcome
  52  mixed outcome
  83 Mixed outcome
  66  reduced virulence
  12 reduced virulence
  696 Reduced virulence
  1 Reduced Virulence
  30  unaffected pathogenicity
  786 Unaffected pathogenicity
  1 Wild-type mutualism
** Blast results of note: **

'Result A'




    ###Interproscan
    Interproscan was used to give gene models functional annotations
```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/
for Genes in $(ls gene_pred/augustus/*/*/filtered_contigs/*_EMR_singlestrand_aug_out.aa); do
echo $Genes; cat $Genes | grep '>' |wc -l
$ProgDir/sub_interproscan.sh $Genes
done
```

Still need to run below once interproscan finished!

 ProgDir=/home/halesk/git_repos/tools/seq_tools/feature_annotation/interproscan
    Genes=gene_pred/augustus/N.ditissima/R0905_v2/R0905_v2_EMR_aug_out.aa
    InterProRaw=gene_pred/interproscan/N.ditissima/R0905_v2/raw
    $ProgDir/append_interpro.sh $Genes $InterProRaw

Already changed first line
Change second line to include this as a for loop
gene_pred/augustus/P.violae/DE/filtered_contigs/DE_EMR_singlestrand_aug_out.aa 

   ### B) SwissProt
Putative ID's were given to genes with homology to SwissProt (the highly curated gene subset of UniProt). IDs were given by through BLASTP searches.

qlogin
ProjDir=/home/groups/harrisonlab/project_files/pythium
cd $ProjDir
for Proteins in $(ls gene_pred/augustus/*/*/filtered_contigs/*_EMR_singlestrand_aug_out.aa); do
Strain=$(echo $ORF_Gff | rev | cut -f3 -d '/' | rev)
Organism=$(echo $ORF_Gff | rev | cut -f4 -d '/' | rev)
OutDir=$ProjDir/gene_pred/uniprot/$Organism/$Strain
mkdir -p $OutDir
blastp \
-db /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot \
-query $Proteins
-out $OutDir/swissprot_v2015_09_hits.tbl  \
-evalue 1e-100 \
-outfmt 6 \
-num_threads 8 \
-num_alignments 10
    
    
    Still need to run 2nd part of interproscan again - on this bit too?
    
    
    
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

    ###Orthology analysis
    
    # For a comparison between 3 Pythium clade 2 isolates (Pult, Paph, Pirr, PvHL, PvDE)


```bash
  ProjDir=/home/groups/harrisonlab/project_files/pythium
  cd $ProjDir
  IsolateAbrv=PvHL_PvDE_Pirr_Pult_Paph
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## Format fasta files

### for P.violae HL
```bash
Taxon_code=PvHL
Fasta_file=gene_pred/augustus/P.violae/HL/filtered_contigs/HL_EMR_singlestrand_aug_out.aa
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.violae DE
```bash
Taxon_code=PvDE
Fasta_file=gene_pred/augustus/P.violae/DE/filtered_contigs/DE_EMR_singlestrand_aug_out.aa
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.irregulare
```bash
Taxon_code=Pirr
Fasta_file=raw_dna/external/P.irregulare/DAOMBR486/pir.maker.proteins.fasta
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.ultimum
```bash
Taxon_code=Pult
Fasta_file=raw_dna/external/P.ultimumvar.ultimum/DAOMBR144/pythium_ultimum_proteins.fasta
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.aphanidermatum
```bash
Taxon_code=Paph
Fasta_file=raw_dna/external/P.aphanidermatum/DAOM_BR444/pag1.maker.proteins.fasta
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```


## Filter proteins into good and poor sets.

```bash
  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```

## Perform an all-vs-all blast of the proteins

```bash
BlastDB=$WorkDir/blastall/$IsolateAbrv.db

makeblastdb -in $Good_proteins_file -dbtype prot -out $BlastDB
BlastOut=$WorkDir/all-vs-all_results.tsv
mkdir -p $WorkDir/splitfiles

SplitDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
$SplitDir/splitfile_500.py --inp_fasta $Good_proteins_file --out_dir $WorkDir/splitfiles --out_base goodProteins

ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology  
for File in $(find $WorkDir/splitfiles); do
Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
done
printf "\n"
echo $File
BlastOut=$(echo $File | sed 's/.fa/.tab/g')
qsub $ProgDir/blast_500.sh $BlastDB $File $BlastOut
done
```

done down to here!
below still calling on variables set above, so may need to set the variable again - if log back into same screen may be ok

## Merge the all-vs-all blast results  
```bash  
  MergeHits="$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
```

## Perform ortholog identification

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts
```

## Plot venn diagrams:

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/venn_diagrams
  $ProgDir/ven_diag_5_way.R --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs

```
  [1] "Pcac (8494)"
  [1] 586
  [1] 126
  [1] "Pcap (7393)"
  [1] 348
  [1] 59
  [1] "Pinf (8079)"
  [1] 601
  [1] 107
  [1] "Ppar (8687)"
  [1] 732
  [1] 95
  [1] "Psoj (7592)"
  [1] 642
  [1] 153
  NULL
```
    
    more on email for above
    
    
    
    ###Downloading RNA seq data from NCBI
    
    The commands used to do this can be found in /pythium/Gene_pred
    
    
    ###The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash

for SplitDir in $(ls -d gene_pred/Augustus_split/*/*/); do
Strain=$(echo $SplitDir | cut -d '/' -f4)
Organism=$(echo $SplitDir | cut -d '/' -f3)
InStringAA=''
InStringNeg=''
InStringTab=''
InStringTxt=''
SigpDir=Augustus_signalp-4.1
for GRP in $(ls -l $SplitDir/*_Augustus_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_Augustus_preds_$GRP""_sp.aa";  
InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_Augustus_preds_$GRP""_sp_neg.aa";  
InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_Augustus_preds_$GRP""_sp.tab";
InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_Augustus_preds_$GRP""_sp.txt";  
done
cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
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
