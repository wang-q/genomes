# Oomycota

<!-- TOC -->
* [Oomycota](#oomycota)
  * [Taxon info](#taxon-info)
    * [Symlink](#symlink)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
  * [All assemblies](#all-assemblies)
    * [Extract from `../Protists` and create assembly.tsv](#extract-from-protists-and-create-assemblytsv)
    * [Count `assembly.tsv`](#count-assemblytsv)
  * [MinHash](#minhash)
    * [Condense branches in the minhash tree](#condense-branches-in-the-minhash-tree)
<!-- TOC -->

## Taxon info

* [Oomycota](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4762)

A nice introduction article about Oomycota:

https://doi.org/10.1016/j.cub.2018.05.062

Phylogenomic Reconstruction of the Oomycete Phylogeny Derived from 37 Genomes

https://doi.org/10.1128/msphere.00095-17

### Symlink

```shell
mkdir -p ~/data/Protists/Oomycota
cd ~/data/Protists/Oomycota

rm -fr ASSEMBLY
rm -fr STRAINS

ln -s ../ASSEMBLY ASSEMBLY
ln -s ../STRAINS STRAINS

```

### List all ranks

```shell
nwr member Oomycota |
    grep -v " x " |
    tsv-summarize -H -g rank --count |
    rgr md stdin --num

```

| rank            | count |
|-----------------|------:|
| phylum          |     1 |
| no rank         |    62 |
| family          |    20 |
| genus           |    86 |
| species         |  3349 |
| subspecies      |     3 |
| order           |    11 |
| strain          |    40 |
| isolate         |     3 |
| varietas        |    12 |
| forma           |     4 |
| species group   |     1 |
| forma specialis |     9 |

### Species with assemblies

```shell
cd ~/data/Protists/Oomycota
mkdir -p summary

SPECIES=$(
    nwr member \
        Oomycota\
        -r species |
        sed '1d' |
        cut -f 1 |
        sort |
        uniq
)

for S in $SPECIES; do
    GB=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
            " |
            sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
    )

    CHR=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
                AND assembly_level IN ('Complete Genome', 'Chromosome')
            " |
            sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
    )

    if [[ ${GB} -gt 0 ]]; then
        echo -e "$S\t$GB\t$CHR"
    fi
done |
    nwr append stdin |
    tsv-select -f 1,4,2-3 |
    tsv-sort -k3,3nr -k4,4nr -k2,2 |
    (echo -e 'species_id\tspecies\tGB\tCHR' && cat) \
    > species.count.tsv

cat species.count.tsv |
    tsv-filter -H --ge GB:5 |
    rgr md stdin --num

```

| species_id | species                          | GB | CHR |
|-----------:|----------------------------------|---:|----:|
|      29920 | Phytophthora cactorum            | 28 |   0 |
|     164328 | Phytophthora ramorum             | 28 |   0 |
|      36331 | Globisporangium irregulare       | 22 |   0 |
|       4784 | Phytophthora capsici             | 17 |   1 |
|     114742 | Pythium insidiosum               | 15 |   0 |
|     112090 | Aphanomyces astaci               | 14 |   0 |
|      53985 | Phytophthora fragariae           | 13 |   0 |
|       4792 | Phytophthora nicotianae          | 13 |   0 |
|     325452 | Phytophthora kernoviae           | 12 |   0 |
|      67593 | Phytophthora sojae               | 11 |   9 |
|      65357 | Albugo candida                   | 11 |   0 |
|     529119 | Globisporangium cryptoirregulare | 10 |   0 |
|       4785 | Phytophthora cinnamomi           |  9 |   0 |
|     542832 | Peronospora effusa               |  8 |   1 |
|       4787 | Phytophthora infestans           |  7 |   1 |
|     100861 | Aphanomyces euteiches            |  7 |   0 |
|     186163 | Globisporangium cylindrosporum   |  5 |   0 |
|     129355 | Phytophthora lateralis           |  5 |   0 |
|     129364 | Phytophthora rubi                |  5 |   0 |
|       4781 | Plasmopara halstedii             |  5 |   0 |
|      41045 | Pythium oligandrum               |  5 |   0 |

## All assemblies

### Extract from `../Protists` and create assembly.tsv

```shell
cd ~/data/Protists/Oomycota/

cat ../summary/collect.pass.tsv | # 1702
    nwr restrict Oomycota -f stdin -c 3 | # belongs to Oomycota 475
    tsv-filter -H --le "C:10000" --ge "N50:10000" | # more stringent parameters 389
    keep-header -- sort \
    > summary/collect.pass.tsv

cat ~/Scripts/genomes/assembly/Protists.assembly.tsv |
    tsv-join -H -f summary/collect.pass.tsv -k 1 \
    > summary/assembly.tsv

# biosample.tsv
cp ../summary/attributes.lst summary/

cat ../summary/biosample.tsv |
    grep -Fw -f <(cat summary/collect.pass.tsv | tsv-select -H -f BioSample | sort | uniq) \
    > summary/biosample.tsv

```

### Count `assembly.tsv`

```shell
cd ~/data/Protists/Oomycota/

nwr template summary/assembly.tsv \
    --count \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --fmt

# .lst and .count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    keep-header -- tsv-sort -k1,1 |
    tsv-filter -H --ge 3:2 |
    rgr md stdin --num

```

| item    | count |
|---------|------:|
| strain  |   388 |
| species |   184 |
| genus   |    21 |
| family  |     6 |
| order   |     5 |
| class   |     1 |

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Albugo           |        2 |        9 |
| Aphanomyces      |        5 |       15 |
| Elongisporangium |        4 |        4 |
| Globisporangium  |       48 |       65 |
| Halophytophthora |        2 |        2 |
| Hyaloperonospora |        3 |        6 |
| Lagenidium       |        2 |        3 |
| Peronospora      |        6 |       18 |
| Phytophthora     |       51 |      187 |
| Phytopythium     |       14 |       17 |
| Pilasporangium   |        1 |        3 |
| Plasmopara       |        3 |        7 |
| Pythium          |       34 |       43 |
| Saprolegnia      |        2 |        2 |

## MinHash

```shell
cd ~/data/Protists/Oomycota/

nwr template summary/assembly.tsv \
    --mh \
    --parallel 8 \
    --ani-ab 0.05 \
    --ani-nr 0.005

# Compute assembly sketches
bash MinHash/compute.sh

# Non-redundant strains within species
bash MinHash/nr.sh

find MinHash -name "NR.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/NR.lst
find MinHash -name "redundant.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/redundant.lst
wc -l summary/NR.lst summary/redundant.lst
# 124 summary/NR.lst
# 131 summary/redundant.lst

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#0

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Protists/Oomycota/

nwr template summary/assembly.tsv \
    --mh \
    --parallel 8 \
    --not-in summary/redundant.lst \
    --height 0.4

bash MinHash/dist.sh

```

### Condense branches in the minhash tree

```shell
mkdir -p ~/data/Protists/Oomycota/tree
cd ~/data/Protists/Oomycota/tree

nwr order --nd --an ../MinHash/tree.nwk \
    > minhash.order.newick

nwr pl-condense --map -r order -r family -r genus \
    minhash.order.newick ../MinHash/species.tsv |
    nwr order stdin --nd --an \
    -o minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

nwr topo --bl minhash.condensed.newick | # remove comments
    nwr tex stdin --bl -o Oomycota.minhash.tex

tectonic Oomycota.minhash.tex

```
