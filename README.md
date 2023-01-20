# Genomes

## Detailed steps

* [Bacteria.md](./Bacteria.md): All genomes of *Bacteria* and *Archaea*, species by species

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

brew install librsvg
brew install jq
brew install pup

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
| Archaea       |         592 |     257 |     865 |          0 |   1,714 |
| Bacteria      |       5,752 |   4,960 |  24,441 |        949 |  36,102 |
| Eukaryota     |      66,556 |  97,988 | 511,429 |     36,366 | 712,339 |
| Fungi         |       5,971 |   7,385 |  55,187 |      1,560 |  70,103 |
| Metazoa       |      48,230 |  69,817 | 267,680 |     18,139 | 403,866 |
| Viridiplantae |       8,419 |  16,870 | 174,006 |     16,295 | 215,590 |
| Viruses       |       2,030 |   2,584 |   7,177 |         65 |  11,856 |
| All taxa      |      74,961 | 105,791 | 543,898 |     37,380 | 762,030 |

## Download all valid Bacteria and Archaea genomes

[Bacteria.md](./Bacteria.md)

## Prokaryote groups

* Pseudomonas
* MTBC
* Tenericutes

## Eukaryote groups

* [Trichoderma](groups/Trichoderma.md)
* Fungi
* [Protistis](./groups/Protists.md)
