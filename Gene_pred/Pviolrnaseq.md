###Gene Prediction

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


###2 QC

RNAseq data was trimmed to remove low quality sequence and adapters. As the
rnaseq data represents unpaired libraries, the fastq-mcf wrapper designed for
unpaired reads was used.

Warning: this below is not tested look at Sascha's to check got it right

```bash
  for StrainPath in $(ls -d raw_rna/external/*/*); do
    echo $StrainPath
    for RnaDat in $(ls $StrainPath/unpaired/*.fastq.gz); do
      IlluminaAdapters=/home/halesk/git_repos/tools/seq_tools/ncbi_adapters.fa
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
      echo $RnaDat
      qsub $ProgDir/dna_qc_fastq-mcf_unpaired.sh $RnaDat $IlluminaAdapters RNA
    done
  done
```

Below again is copied, not altered for my stuff

###3 Align reads vs. published Assemblys

Alignments of RNAseq reads were made against each of the Pythium Assembly
assemblies.

3.1) Alignment

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "Aligning RNAseq data against:"
    echo "$Organism  $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    for RnaDat in $(ls $StrainPath/unpaired/*.fq.gz); do
      DatName=$(echo $RnaDat| rev | cut -d '/' -f1 | rev | tr -d 'fq.gz')
      echo "Aligning the following RNAseq file:"
      echo "$DatName"
      OutDir=alignment/$Organism/$Strain/$DatName
      qsub $ProgDir/tophat_alignment_unpaired.sh $Assembly $RnaDat$OutDir
    done
  done
```

The alignment files of RNAseq data against each genome were merged into a single
file that can be passed to a gene prediction program to indicate the location of
aligned RNAseq data against a particular genome.

```bash
	for StrainDir in $(ls -d alignment/*/*); do
		Strain=$(echo $StrainDir | rev | cut -d '/' -f1 | rev)
		ls alignment/*/$Strain/*/accepted_hits.bam > bamlist.txt
		mkdir -p $StrainDir/merged
		bamtools merge -list bamlist.txt -out $StrainDir/merged/accepted_hits_merged.bam
  done
```

4 Run Braker1

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    OutDir=gene_pred/braker/$Organism/$Strain
    AcceptedHits=alignment/$Organism/$Strain/accepted_hits.bam
    GeneModelName="$Organism"_"$Strain"_braker
    qsub $ProgDir/sub_braker.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```

  5 Extract gff and amino acid sequences
```bash
  for File in $(ls gene_pred/braker/*/*/*_braker/augustus.gff); do
    getAnnoFasta.pl $File
    OutDir=$(dirname $File)
    echo "##gff-version 3" > $OutDir/augustus_extracted.gff
    cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
  done
```
