###Gene Prediction

#1 Data organisation

Downloading RNA seq data from NCBI

These commands were used to download SRA files from NCBI


Make a directory first if necessary with species/strain.

```bash

mkdir -p x/x/x

fastq-dump -O raw_rna/external/P.ultimum/DAOM_BR144 SRR059021

-O is choosing where to put it (Output)
Final number is the accession number from NCBI, under Run



###2 QC


Warning: this below is not tested look at Sascha's to check got it right

```
for Strain in $(ls raw_rna/paired/P.*/); do
IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
  qsub /home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA

Below again is copied, not altered for my stuff

###3 Align reads vs. publihsed genomes

Alignments of RNAseq reads were made against the published Genomes using tophat:

3.1) Alignment

  Pinf=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.parsed.fa
  Ppar=assembly/external_group/P.parisitica/310/dna/phytophthora_parasitica_inra_310.i2.scaffolds.genome.parsed.fa
  Pcap=assembly/external_group/P.capsici/LT1534/dna/Phyca11_unmasked_genomic_scaffolds.fasta
  Psoj=assembly/external_group/P.sojae/67593/dna/Phytophthora_sojae.ASM14975v1.26.dna.genome.parsed.fa
  for Genome in $Pinf $Ppar $Pcap $Psoj; do
    Strain=$(echo $Genome| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Genome | rev | cut -d '/' -f4 | rev)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    FileF=qc_rna/raw_rna/genbank/P.cactorum/F/SRR1206032_trim.fq.gz
    FileR=qc_rna/raw_rna/genbank/P.cactorum/R/SRR1206033_trim.fq.gz
    OutDir=alignment/$Organism/$Strain
    qsub $ProgDir/tophat_alignment.sh $Genome $FileF $FileR $OutDir
  done
  
4 Run Braker1
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
  for Assembly in $Pinf $Ppar $Pcap $Psoj; do
    # Jobs=$(qstat | grep 'tophat_ali' | wc -l)
    # while [ $Jobs -gt 0 ]; do
    #   sleep 10
    #   printf "."
    #   Jobs=$(qstat | grep 'tophat_ali' | wc -l)
    # done
    printf "\n"
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    OutDir=gene_pred/braker/$Organism/$Strain
    AcceptedHits=alignment/$Organism/$Strain/accepted_hits.bam
    GeneModelName="$Organism"_"$Strain"_braker
    qsub $ProgDir/sub_braker.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
  
  5 Extract gff and amino acid sequences
  for File in $(ls gene_pred/braker/P.*/*/*_braker/augustus.gff); do
    getAnnoFasta.pl $File
    OutDir=$(dirname $File)
    echo "##gff-version 3" > $OutDir/augustus_extracted.gff
    cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
  done
P.capsici/LT1534 25866
P.infestans/T30-4 86760
P.parisitica/310 20794
P.sojae/67593 34513