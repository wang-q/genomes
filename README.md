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
| Archaea       |         540 |     243 |     841 |          0 |   1,624 |
| Bacteria      |       5,464 |   4,816 |  23,955 |        922 |  35,157 |
| Eukaryota     |      65,275 |  96,617 | 500,987 |     35,445 | 698,324 |
| Fungi         |       5,689 |   7,221 |  53,392 |      1,518 |  67,820 |
| Metazoa       |      47,484 |  68,787 | 262,476 |     17,660 | 396,407 |
| Viridiplantae |       8,218 |  16,759 | 170,895 |     15,902 | 211,774 |
| Viruses       |       1,947 |   2,527 |   7,160 |         65 |  11,699 |
| All taxa      |      73,256 | 104,205 | 532,929 |     36,432 | 746,822 |
