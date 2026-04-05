# Genomes

## `date`

The date `date --utc` of executing `nwr download` is `Sun Apr  5 15:59:45 UTC 2026`

## Software

### Install `nwr` and initiate local databases

```shell
cbp install nwr
cbp install sqlite3 # 3.34 or above

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
    *.dmp \
    taxdump.tar.gz \
    taxdump.tar.gz.md5 \
    assembly_summary_genbank.txt \
    assembly_summary_refseq.txt

```

### Additional Packages

```shell
cbp install hmmer easel
cbp install mafft muscle trimal
cbp install samtools
cbp install fasttree iqtree2
cbp install mash mmseqs

cbp install fd pigz
cbp install jq tectonic
cbp install tva 
cbp install pgr
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
curl -L "https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=statistics&?&unclassified=hide&uncultured=hide" |
    tva from html -q 'table[bgcolor="#CCCCFF"] table[bgcolor="#FFFFFF"] tr td text{}' |
    grep '\S' |
    paste -d $'\t' - - - - - - |
    tva to md --right 2-6

```

| Ranks:        | higher taxa |   genus | species | lower taxa |     total |
| ------------- | ----------: | ------: | ------: | ---------: | --------: |
| Archaea       |           0 |     340 |   1,200 |      2,290 |     2,290 |
| Bacteria      |           0 |   5,782 |  33,615 |     90,218 |    90,218 |
| Eukaryota     |           0 | 104,261 | 631,437 |    804,447 |   804,447 |
| Fungi         |           0 |   8,095 |  74,507 |     88,460 |    88,460 |
| Metazoa       |           0 |  75,546 | 340,416 |    453,240 |   453,240 |
| Viridiplantae |           0 |  16,338 | 198,532 |    237,280 |   237,280 |
| Viruses       |          36 |   3,493 |  14,612 |    200,795 |   201,328 |
| All taxa      |          54 | 113,878 | 700,762 |  1,097,758 | 1,118,224 |

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
