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

The date of executing `nwr download` is `Fri Nov 15 15:51:01 CST 2024`

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

rm \
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
    rgr md stdin --right 2-6

```

| Ranks:        | higher taxa |   genus | species | lower taxa |   total |
|---------------|------------:|--------:|--------:|-----------:|--------:|
| Archaea       |         742 |     314 |   1,046 |          0 |   2,102 |
| Bacteria      |       7,192 |   5,435 |  27,084 |        972 |  40,683 |
| Eukaryota     |      71,173 | 102,113 | 551,181 |     39,017 | 763,484 |
| Fungi         |       6,810 |   7,842 |  61,751 |      1,619 |  78,022 |
| Metazoa       |      51,179 |  73,101 | 287,542 |     19,570 | 431,392 |
| Viridiplantae |       8,825 |  17,071 | 186,365 |     17,424 | 229,685 |
| Viruses       |       2,328 |   2,853 |   5,815 |        963 |  11,959 |
| All taxa      |      81,466 | 110,716 | 585,112 |     40,952 | 818,246 |

## Download all valid Bacteria and Archaea genomes

* [Trichoderma](./groups/Trichoderma.md) is a good example of familiarizing yourself with the
  processing steps.

* [Bacteria](./Bacteria.md): All genomes of **Bacteria** and **Archaea**, species by species

* [Fungi](./groups/Fungi.md): All genomes of **Fungi**, species by species

## Prokaryote groups

* [Pseudomonas](groups/Pseudomonas.md)
* MTBC
* Tenericutes

## Eukaryote groups

* [Trichoderma](groups/Trichoderma.md)
* [Protistis](groups/Protists.md)
