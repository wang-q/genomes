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
| Archaea       |         536 |     240 |     830 |          0 |   1,606 |
| Bacteria      |       5,446 |   4,802 |  23,909 |        925 |  35,082 |
| Eukaryota     |      65,206 |  96,518 | 500,097 |     35,417 | 697,238 |
| Fungi         |       5,686 |   7,213 |  53,319 |      1,518 |  67,736 |
| Metazoa       |      47,431 |  68,706 | 262,013 |     17,646 | 395,796 |
| Viridiplantae |       8,205 |  16,756 | 170,584 |     15,888 | 211,433 |
| Viruses       |       1,934 |   2,525 |   7,158 |         65 |  11,682 |
| All taxa      |      73,152 | 104,087 | 531,980 |     36,407 | 745,626 |
