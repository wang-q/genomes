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
mkdir -p ~/data/Protists/Oomycota/summary
cd ~/data/Protists/Oomycota/summary

# should have a valid name of genus
nwr member Oomycota -r genus |
    grep -v " sp." |
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#86 genus.list

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS1.tsv

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB1.tsv

wc -l RS*.tsv GB*.tsv
#   8 RS1.tsv
# 197 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     8
#GB1     429

```

## All assemblies

### Extract from `../Protists` and create assembly.tsv

```shell
cd ~/data/Protists/Oomycota/

mkdir -p summary

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

### Count before download

```shell
cd ~/data/Oomycota

nwr template ~/Scripts/genomes/assembly/Oomycota.assembly.tsv \
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
| strain  |   428 |
| species |   198 |
| genus   |    23 |
| family  |     8 |
| order   |     6 |
| class   |     2 |

| genus             | #species | #strains |
|-------------------|---------:|---------:|
| Albugo            |        2 |       10 |
| Aphanomyces       |        5 |       27 |
| Elongisporangium  |        5 |        5 |
| Globisporangium   |       47 |       53 |
| Halophytophthora  |        2 |        2 |
| Hyaloperonospora  |        3 |        7 |
| Lagenidium        |        1 |        3 |
| Paralagenidium    |        1 |        2 |
| Peronospora       |        6 |       18 |
| Phytophthora      |       53 |      191 |
| Phytopythium      |       14 |       18 |
| Pilasporangium    |        1 |        3 |
| Plasmopara        |        4 |       11 |
| Pseudoperonospora |        2 |        2 |
| Pythium           |       43 |       66 |
| Saprolegnia       |        2 |        2 |
| Sclerospora       |        1 |        2 |

### Download and check

```shell
cd ~/data/Oomycota

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Oomycota.assembly.tsv \
    --ass

# Run
bash ASSEMBLY/rsync.sh

# Check md5; create check.lst
# rm ASSEMBLY/check.lst
bash ASSEMBLY/check.sh

# Put the misplaced directory into the right place
#bash ASSEMBLY/reorder.sh
#
# This operation will delete some files in the directory, so please be careful
#cat ASSEMBLY/remove.lst |
#    parallel --no-run-if-empty --linebuffer -k -j 1 '
#        if [[ -e "ASSEMBLY/{}" ]]; then
#            echo Remove {}
#            rm -fr "ASSEMBLY/{}"
#        fi
#    '

find ASSEMBLY/ -name "*_genomic.fna.gz" |
    grep -v "_from_" |
    wc -l
#428

# N50 C S; create n50.tsv and n50.pass.tsv
bash ASSEMBLY/n50.sh 50000 5000 20000000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50,S" --max "C"
#N50_min S_min   C_max
#280942  53131624        4921

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "S:0.1,0.5" --quantile "N50:0.1,0.5"  --quantile "C:0.5,0.9"
#S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
#35693128        50477035.5      8472.2  25917.5 4826    14503.7

# After the above steps are completed, run the following commands.

# Collect; create collect.tsv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```

| #item            | fields | lines |
|------------------|-------:|------:|
| url.tsv          |      3 |   428 |
| check.lst        |      1 |   428 |
| collect.tsv      |     20 |   429 |
| n50.tsv          |      4 |   429 |
| n50.pass.tsv     |      4 |   124 |
| collect.pass.tsv |     23 |   124 |
| pass.lst         |      1 |   123 |
| omit.lst         |      1 |   308 |
| rep.lst          |      1 |    66 |

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Oomycota/ \
    wangq@202.119.37.251:data/Oomycota

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Oomycota/ \
    wangq@58.213.64.36:data/Oomycota

# rsync -avP wangq@202.119.37.251:data/Oomycota/ ~/data/Oomycota

# rsync -avP -e 'ssh -p 8804' wangq@58.213.64.36:data/Oomycota/ ~/data/Oomycota

```

## BioSample

```shell
cd ~/data/Oomycota

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Oomycota.assembly.tsv \
    --bs

# Run this script twice and it will re-download the failed files
bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 10

datamash check < BioSample/biosample.tsv
#424 lines, 36 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## MinHash

```shell
cd ~/data/Oomycota/

nwr template ~/Scripts/genomes/assembly/Oomycota.assembly.tsv \
    --mh \
    --parallel 16 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
    --ani-nr 0.005 \
    --height 0.4

# Compute assembly sketches
bash MinHash/compute.sh

#find MinHash -name "*.msh" -empty | wc -l

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
mkdir -p ~/data/Oomycota/tree
cd ~/data/Oomycota/tree

nwr order --nd --an ../MinHash/tree.nwk \
    > minhash.order.newick

nwr pl-condense -r order -r family -r genus --map minhash.order.newick ../Count/species.tsv |
    nwr tex stdin --bl -o minhash.tex

tectonic minhash.tex

```
