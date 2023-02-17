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

Thu Feb 9 05:00:03 CST 2023

## Software

### Install `nwr` and initiate local databases

```shell
brew install wang-q/tap/nwr # 0.5.4 or above
brew install sqlite         # 3.34 or above

nwr download
nwr txdb

nwr ardb
nwr ardb --genbank

```

### Packages managed by Homebrew

```shell
brew install hmmer
brew install brewsci/bio/muscle
brew install brewsci/bio/fasttree
brew install brewsci/bio/newick-utils
brew install brewsci/bio/trimal

brew install datamash
brew install miller
brew install wang-q/tap/tsv-utils
brew install wang-q/tap/faops

brew install librsvg
brew install jq
brew install pup

# dN/dS
brew install brewsci/bio/clustal-w
brew install brewsci/bio/paml

cpanm Bio::Tools::Run::Alignment::Clustalw
cpanm https://github.com/wang-q/Bio-Tools-Phylo-PAML.git

```

* Pangenome

    * `PPanGGOLiN` is used in this project. Installation steps can be
      found [here](https://github.com/wang-q/dotfiles/blob/master/others.sh).

### Other Packages

* Pangenome
    * `PPanGGOLiN` is used in this project. Installation steps can be
      found [here](https://github.com/wang-q/dotfiles/blob/master/others.sh).

## Strain info

* [Bacteria](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2)
* [Archaea](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2157)

### NCBI statistics

```shell
curl -L "https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=statistics&uncultured=hide&unspecified=hide" |
    pup 'table[bgcolor=#CCCCFF] tr td text{}' |
    grep '\S' |
    paste -d $'\t' - - - - - - |
    head -n 9 |
    mlr --itsv --omd cat


```

| Ranks:        | higher taxa |   genus | species | lower taxa |   total |
|---------------|------------:|--------:|--------:|-----------:|--------:|
| Archaea       |         593 |     257 |     867 |          0 |   1,717 |
| Bacteria      |       5,783 |   4,973 |  24,503 |        950 |  36,209 |
| Eukaryota     |      66,566 |  98,103 | 512,437 |     36,398 | 713,504 |
| Fungi         |       5,970 |   7,405 |  55,346 |      1,562 |  70,283 |
| Metazoa       |      48,238 |  69,901 | 268,271 |     18,162 | 404,572 |
| Viridiplantae |       8,413 |  16,876 | 174,241 |     16,302 | 215,832 |
| Viruses       |       2,045 |   2,584 |   7,178 |         65 |  11,872 |
| All taxa      |      75,018 | 105,919 | 544,971 |     37,413 | 763,321 |

## Download all valid Bacteria and Archaea genomes

[Bacteria.md](./Bacteria.md): All genomes of *Bacteria* and *Archaea*, species by species

## Prokaryote groups

* [Pseudomonas](groups/Pseudomonas.md)
* MTBC
* Tenericutes

## Eukaryote groups

* [Trichoderma](groups/Trichoderma.md)
* Fungi
* [Protistis](groups/Protists.md)
