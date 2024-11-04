# Genomes

<!-- toc -->

- [`date`](#date)
- [Software](#software)
    * [Install `nwr` and initiate local databases](#install-nwr-and-initiate-local-databases)
    * [Packages managed by Homebrew](#packages-managed-by-homebrew)
    * [Other Packages](#other-packages)
- [Strain info](#strain-info)
    * [NCBI statistics](#ncbi-statistics)
- [Download all valid Bacteria and Archaea genomes](#download-all-valid-bacteria-and-archaea-genomes)
- [Prokaryote groups](#prokaryote-groups)
- [Eukaryote groups](#eukaryote-groups)

<!-- tocstop -->

## `date`

The date of executing `nwr download` is `Thu Oct 10 18:14:27 CST 2024`

## Software

### Install `nwr` and initiate local databases

```shell
brew install wang-q/tap/nwr # 0.5.4 or above
brew install sqlite         # 3.34 or above
#brew link sqlite --force

nwr download
nwr txdb

nwr ardb
nwr ardb --genbank

cd $HOME/.nwr
tar cvfz ncbi.$(date +"%Y%m%d").tar.gz \
    taxdump.tar.gz \
    taxdump.tar.gz.md5 \
    assembly_summary_genbank.txt \
    assembly_summary_refseq.txt

```

### Packages managed by Homebrew

```shell
brew install hmmer
brew install brewsci/bio/easel
brew install mafft
brew install samtools
brew install brewsci/bio/muscle
brew install brewsci/bio/fasttree
brew install brewsci/bio/iqtree2
brew install brewsci/bio/newick-utils
brew install brewsci/bio/trimal
brew install brewsci/bio/mash
brew install mmseqs2

brew install datamash
#cargo install fd-find
brew install wang-q/tap/tsv-utils
brew install wang-q/tap/faops
brew install wang-q/tap/hnsm

brew install librsvg
brew install jq
brew install pup
brew install pigz

# dN/dS
brew install brewsci/bio/clustal-w
brew install brewsci/bio/paml

cpanm Bio::Tools::Run::Alignment::Clustalw
cpanm https://github.com/wang-q/Bio-Tools-Phylo-PAML.git

```

### Other Packages

* Pangenome
    * `PPanGGOLiN` is used in this project. Installation steps can be
      found [here](https://github.com/wang-q/dotfiles/blob/master/others.sh).

## Strain info

* [Bacteria](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2)
* [Archaea](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2157)
* [Fungi](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4751)
* [Oomycota](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4762)

### NCBI statistics

```shell
curl -L "https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=statistics&uncultured=hide&unspecified=hide" |
    pup 'table[bgcolor=#CCCCFF] tr td text{}' |
    grep '\S' |
    paste -d $'\t' - - - - - - |
    head -n 9 |
    mlr --itsv --omd cat

```

| Ranks:        | higher taxa | genus   | species | lower taxa | total   |
|---------------|-------------|---------|---------|------------|---------|
| Archaea       | 730         | 308     | 1,032   | 0          | 2,070   |
| Bacteria      | 7,116       | 5,394   | 26,918  | 968        | 40,396  |
| Eukaryota     | 70,768      | 101,641 | 548,385 | 38,755     | 759,549 |
| Fungi         | 6,755       | 7,801   | 61,101  | 1,612      | 77,269  |
| Metazoa       | 50,879      | 72,702  | 285,716 | 19,303     | 428,600 |
| Viridiplantae | 8,788       | 17,048  | 186,090 | 17,440     | 229,366 |
| Viruses       | 2,302       | 2,852   | 5,815   | 963        | 11,932  |
| All taxa      | 80,947      | 110,196 | 582,136 | 40,686     | 813,965 |

## Download all valid Bacteria and Archaea genomes

* [Trichoderma](./groups/Trichoderma.md) is a good example of familiarizing yourself with the
  processing steps.

* [Bacteria](./Bacteria.md): All genomes of **Bacteria** and **Archaea**, species by species

* [Fungi](./Fungi.md): All genomes of **Fungi**, species by species

## Prokaryote groups

* [Pseudomonas](groups/Pseudomonas.md)
* MTBC
* Tenericutes

## Eukaryote groups

* [Trichoderma](groups/Trichoderma.md)
* [Protistis](groups/Protists.md)
