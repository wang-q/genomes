# Oomycota

<!-- toc -->

## Taxon info

* [Oomycota](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4762)

A nice introduction article about Oomycota:

https://doi.org/10.1016/j.cub.2018.05.062

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
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

```

| rank       | count |
|------------|------:|
| phylum     |     1 |
| no rank    |    54 |
| family     |    20 |
| genus      |    86 |
| species    |  1219 |
| subspecies |     3 |
| order      |    11 |
| strain     |    39 |
| isolate    |     1 |
| varietas   |    15 |
| forma      |     4 |

### Species with assemblies

```shell
cd ~/data/Protists/Oomycota
mkdir -p summary

SPECIES=$(
    nwr member \
        Oomycota\
        -r species |
        grep -v " sp." |
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
    mlr --itsv --omd cat

```

| species_id | species                 | GB | CHR |
|------------|-------------------------|----|-----|
| 164328     | Phytophthora ramorum    | 28 | 0   |
| 29920      | Phytophthora cactorum   | 21 | 0   |
| 114742     | Pythium insidiosum      | 15 | 0   |
| 112090     | Aphanomyces astaci      | 14 | 0   |
| 4784       | Phytophthora capsici    | 14 | 0   |
| 53985      | Phytophthora fragariae  | 13 | 0   |
| 325452     | Phytophthora kernoviae  | 12 | 0   |
| 65357      | Albugo candida          | 9  | 0   |
| 4792       | Phytophthora parasitica | 9  | 0   |
| 542832     | Peronospora effusa      | 8  | 1   |
| 4785       | Phytophthora cinnamomi  | 8  | 0   |
| 100861     | Aphanomyces euteiches   | 7  | 0   |
| 4787       | Phytophthora infestans  | 6  | 1   |
| 129355     | Phytophthora lateralis  | 5  | 0   |
| 129364     | Phytophthora rubi       | 5  | 0   |
| 4781       | Plasmopara halstedii    | 5  | 0   |
| 41045      | Pythium oligandrum      | 5  | 0   |

## All assemblies

### Extract from `../Protists` and create assembly.tsv

```shell
cd ~/data/Protists/Oomycota/

cat ../summary/collect.pass.tsv |
    tsv-filter -H --str-eq "RefSeq_category:Reference Genome" \
    > summary/collect.pass.tsv

cat ../summary/collect.pass.tsv | # 1243
    nwr restrict Oomycota -f stdin -c 3 | # belongs to Oomycota 401
    tsv-filter -H --le "C:10000" --ge "N50:10000" | # more stringent parameters 331
    sed '1d' |
    sort \
    >> summary/collect.pass.tsv

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
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# .lst and .count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    keep-header -- tsv-sort -k1,1 |
    tsv-filter -H --ge 3:2 |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```

| item    | count |
|---------|------:|
| strain  |   333 |
| species |   172 |
| genus   |    23 |
| family  |     9 |
| order   |     8 |
| class   |     4 |

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Albugo           |        2 |        7 |
| Aphanomyces      |        5 |       15 |
| Elongisporangium |        4 |        4 |
| Globisporangium  |       45 |       49 |
| Halophytophthora |        2 |        2 |
| Hyaloperonospora |        3 |        6 |
| Lagenidium       |        1 |        2 |
| Peronospora      |        6 |       18 |
| Phytophthora     |       43 |      153 |
| Phytopythium     |       14 |       16 |
| Pilasporangium   |        1 |        3 |
| Plasmopara       |        3 |        7 |
| Pythium          |       32 |       40 |
| Saprolegnia      |        2 |        2 |

## MinHash

```shell
cd ~/data/Protists/Oomycota/

nwr template summary/assembly.tsv \
    --mh \
    --parallel 16 \
    --ani-ab 0.05 \
    --ani-nr 0.005 \
    --height 0.4

# Compute assembly sketches
bash MinHash/compute.sh

# Distances within species
bash MinHash/species.sh

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#359

# Non-redundant strains within species
bash MinHash/nr.sh

find MinHash -name "NR.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/NR.lst
wc -l summary/NR.lst
#49

# Distances between all selected sketches, then hierarchical clustering
bash MinHash/dist.sh

```

### Condense branches in the minhash tree

```shell
mkdir -p ~/data/Protists/Oomycota/tree
cd ~/data/Protists/Oomycota/tree

nwr order --nd --an ../MinHash/tree.nwk \
    > minhash.order.newick

nwr pl-condense -r order -r family -r genus --map minhash.order.newick ../Count/species.tsv |
    nwr tex stdin --bl -o minhash.tex

tectonic minhash.tex

```
