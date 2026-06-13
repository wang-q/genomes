# Archaea

All genomes of *Archaea*

[TOC levels=2-4]: #

  * [Strain info](#strain-info)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
  * [Download all assemblies](#download-all-assemblies)
    * [Create assembly.tsv](#create-assemblytsv)
    * [Count before download](#count-before-download)
    * [Download and check](#download-and-check)
    * [Rsync to hpcc](#rsync-to-hpcc)
  * [BioSample](#biosample)
  * [MinHash](#minhash)
    * [Condense branches in the minhash tree](#condense-branches-in-the-minhash-tree)
  * [Count valid species and strains](#count-valid-species-and-strains)
    * [For *genomic alignments*](#for-genomic-alignments)
    * [For *protein families*](#for-protein-families)
  * [Collect proteins](#collect-proteins)
  * [Phylogenetics with ar53](#phylogenetics-with-ar53)
    * [Find corresponding representative proteins by `hmmsearch`](#find-corresponding-representative-proteins-by-hmmsearch)
    * [Domain related protein sequences](#domain-related-protein-sequences)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
    * [Condense branches in the protein tree](#condense-branches-in-the-protein-tree)


## Strain info

* [Archaea](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2157)

### List all ranks

```shell
nwr member Archaea |
    grep -v " sp." |
    grep -v " x " |
    tva stats -H -g rank --count |
    tva to md --num

```

| rank          | count |
| ------------- | ----: |
| domain        |     1 |
| kingdom       |     4 |
| phylum        |    43 |
| class         |    55 |
| order         |    93 |
| no rank       |   495 |
| species       |  3796 |
| genus         |   340 |
| family        |   132 |
| clade         |    39 |
| strain        |   352 |
| species group |     2 |
| isolate       |     6 |

### Species with assemblies

* 'RefSeq'
    * RS1 - genome_rep: 'Full'
* 'Genbank'
    * GB1 - genome_rep: 'Full'

```shell
mkdir -p ~/data/Archaea/summary
cd ~/data/Archaea/summary

# should have a valid name of genus
nwr member Archaea -r genus |
    grep -v " x " |
    sed '1d' |
    tva sort -n -k 1 \
    > genus.list.tsv

wc -l genus.list.tsv
#340 genus.list

cat genus.list.tsv | tva select -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tva sort -k 2 \
    > RS1.tsv

cat genus.list.tsv | tva select -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tva sort -k 2 \
    > GB1.tsv

wc -l RS*.tsv GB*.tsv
# 1334 RS1.tsv
# 1996 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tva stats --sum 3
        fi
    done
done
# RS1	2871
# GB1	8559

```

## Download all assemblies

### Create assembly.tsv

If a refseq assembly is available, the corresponding genbank one is not downloaded

```shell
cd ~/data/Archaea/summary

echo "
.headers ON
    SELECT
        *
    FROM ar
    WHERE 1=1
        AND species IN ('Saccharomyces cerevisiae')
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tva select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession \
    > raw.tsv

# RS
SPECIES=$(
    cat RS1.tsv |
        tva select -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND genome_rep IN ('Full')
        AND species NOT LIKE '% x %'
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% aff.%'
        AND species NOT LIKE '% cf.%'
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

# Preference for refseq
cat raw.tsv |
    tva select -H -f "assembly_accession" \
    > rs.acc.tsv

# GB
SPECIES=$(
    cat GB1.tsv |
        tva select -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND genome_rep IN ('Full')
        AND species NOT LIKE '% x %'
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% aff.%'
        AND species NOT LIKE '% cf.%'
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tva join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

echo "
    SELECT
        genus || ' sp. ' || infraspecific_name || ' ' || assembly_accession AS name,
        genus || ' sp.', genus, ftp_path, biosample, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND genome_rep IN ('Full')
        AND species NOT LIKE '% x %'
        AND (
            species LIKE '% sp.%'
            OR species LIKE '% aff.%'
            OR species LIKE '% cf.%'
        )
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tva join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

cat raw.tsv |
    tva uniq |
    tva check
#8561 lines, 7 fields

# Create abbr.
cat raw.tsv |
    grep -v '^#' |
    tva uniq |
    tva select -f 1-6 |
    nwr abbr -c "1,2,3" -m 3 --shortsub |
    tva uniq -H -f ftp_path |
    tva uniq -H -f 7 |
    sed '1d' |
    tva select -f 7,4,5,2,6 |
    (echo -e '#name\tftp_path\tbiosample\tspecies\tassembly_level' && cat ) |
    tva filter -H --or --str-in-fld 2:ftp --str-in-fld 2:http |
    tva sort -H -k 4,1 \
    > Archaea.assembly.tsv

tva check < Archaea.assembly.tsv
#8561 lines, 5 fields

# find potential duplicate strains or assemblies
cat Archaea.assembly.tsv |
    tva uniq -f 1 --repeated

cat Archaea.assembly.tsv |
    tva filter --str-not-in-fld 2:ftp

# Edit .assembly.tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Archaea.assembly.tsv
#
# Save the file to another directory to prevent accidentally changing it
# cp Archaea.assembly.tsv ~/Scripts/genomes/assembly

# Cleaning
rm raw*.*sv
```

### Count before download

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class

```shell
cd ~/data/Archaea

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --count \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    tva to md --fmt

# .lst and .count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    tva sort -H -k 1 |
    tva filter -H --ge 3:100 |
    tva to md --fmt
```

| item    | count |
| ------- | ----: |
| strain  | 7,924 |
| species | 1,201 |
| genus   |   326 |
| family  |   125 |
| order   |    84 |
| class   |    52 |

| genus                | #species | #strains |
| -------------------- | -------: | -------: |
| Haloarcula           |       34 |      105 |
| Halorubrum           |       44 |      188 |
| Ignisphaera          |        3 |      117 |
| Metallosphaera       |        8 |      114 |
| Methanobacterium     |       21 |      212 |
| Methanobrevibacter   |       17 |      666 |
| Methanocorpusculum   |        8 |      303 |
| Methanoculleus       |       19 |      191 |
| Methanomethylophilus |        2 |      100 |
| Methanoregula        |        3 |      177 |
| Methanosarcina       |       18 |      240 |
| Methanothrix         |        4 |      321 |
| Nitrososphaera       |        5 |      157 |
| Saccharolobus        |        5 |      117 |
| Thermococcus         |       42 |      219 |

### Download and check

```shell
cd ~/data/Archaea

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --ass

# Run
bash ASSEMBLY/aria2.sh

# Check md5; create check.lst
# rm ASSEMBLY/check.lst
bash ASSEMBLY/check.sh

# Remove failed directories and re-download
bash ASSEMBLY/check.sh 2>&1 |
    grep "checksum failed" |
    sed 's/.*==> //;s/ checksum failed <==//' |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        dir=$(cat ASSEMBLY/url.tsv | tva filter --str-eq "1:{}" | tva select -f 3,1 | tr "\t" "/")
        if [[ -n "$dir" && -e "ASSEMBLY/$dir" ]]; then
            echo Remove ASSEMBLY/$dir
            rm -fr "ASSEMBLY/$dir"
        fi
    '

## Put the misplaced directory into the right place
#bash ASSEMBLY/reorder.sh
#
## This operation will delete some files in the directory, so please be careful
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
#7443

# N50 C S; create n50.tsv and n50.pass.tsv
bash ASSEMBLY/n50.sh 100000 500 500000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tva filter -H --str-in-fld "name:_GCF_" |
    tva stats -H --min "N50" --max "C" --min "S"
# N50_min	C_max	S_min
# 5163	967	449376

cat ASSEMBLY/n50.tsv |
    tva stats -H --quantile "N50:0.1,0.5" --quantile "C:0.5,0.9" --quantile "S:0.1,0.5" |
    tva transpose
# N50_quantile_0.1	5226.4
# N50_quantile_0.5	33717
# C_quantile_0.5	86
# C_quantile_0.9	322
# S_quantile_0.1	949067.8
# S_quantile_0.5	1750667

# After the above steps are completed, run the following commands.

# Collect; create collect.tsv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/
cp ASSEMBLY/omit.lst summary/
cp ASSEMBLY/pass.lst summary/
cp ASSEMBLY/sp.lst summary/
cp ASSEMBLY/rep.lst summary/

cat ASSEMBLY/counts.tsv |
    tva to md --fmt
```

| #item            | fields | lines |
| ---------------- | -----: | ----: |
| url.tsv          |      3 | 8,560 |
| check.lst        |      1 | 8,560 |
| collect.tsv      |     20 | 8,561 |
| n50.pass.tsv     |      4 | 2,426 |
| collect.pass.tsv |     23 | 2,426 |
| pass.lst         |      1 | 2,425 |
| omit.lst         |      1 | 3,364 |
| rep.lst          |      1 |   738 |
| sp.lst           |      1 |   786 |

### Rsync to hpcc

```shell
rsync -avP \
    ~/data/Archaea/ \
    wangq@202.119.37.251:data/Archaea

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Archaea/ \
    wangq@58.213.64.36:data/Archaea

# rsync -avP wangq@202.119.37.251:data/Archaea/ ~/data/Archaea

# rsync -avP -e 'ssh -p 8804' wangq@58.213.64.36:data/Archaea/ ~/data/Archaea

```

## BioSample

```shell
cd ~/data/Archaea

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --bs

# Run this script twice and it will re-download the failed files
bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 50

tva check < BioSample/biosample.tsv
#8548 lines, 112 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## MinHash

```shell
cd ~/data/Archaea

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --mh \
    --parallel 8 \
    --in summary/pass.lst \
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
#  561 summary/NR.lst
#  370 summary/redundant.lst

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#23

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Archaea/

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --mh \
    --parallel 8 \
    --not-in MinHash/abnormal.lst \
    --not-in summary/redundant.lst \
    --height 0.4

bash MinHash/dist.sh
```

### Condense branches in the minhash tree

```shell
mkdir -p ~/data/Archaea/tree
cd ~/data/Archaea/tree

nw_reroot ../MinHash/tree.nwk Saccharom_cere_S288C Saccharom_eub |
    nwr order stdin --nd --an \
    > minhash.reroot.newick

nwr pl-condense --map -r order -r family -r genus \
    minhash.reroot.newick ../MinHash/species.tsv |
    nwr order stdin --nd --an \
    > minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

# svg
nwr topo --bl minhash.condensed.newick | # remove comments
    nw_display -s -b 'visibility:hidden' -w 1200 -v 20 - \
    > Archaea.minhash.svg

```

## Count valid species and strains

### For *genomic alignments*

```shell
cd ~/data/Archaea/

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --rank order --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/order.count.tsv |
    tsv-filter -H --ge "3:20" |
    rgr md stdin --num

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:20" |
    rgr md stdin --num

# Can accept N_COUNT
bash Count/lineage.sh 10

cat Count/lineage.count.tsv |
    rgr md stdin --num

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| item    | count |
|---------|------:|
| strain  |  1437 |
| species |   762 |
| genus   |   237 |
| family  |    79 |
| order   |    52 |
| class   |    28 |

| order                   | #species | #strains |
|-------------------------|---------:|---------:|
| Desulfurococcales       |       21 |       25 |
| Halobacteriales         |      376 |      619 |
| Methanobacteriales      |       34 |      137 |
| Methanococcales         |       15 |       45 |
| Methanomassiliicoccales |        7 |       53 |
| Methanomicrobiales      |       48 |       56 |
| Methanosarcinales       |       46 |       89 |
| Sulfolobales            |       29 |      141 |
| Thermococcales          |       41 |       61 |

| genus                | #species | #strains |
|----------------------|---------:|---------:|
| Haloarcula           |       32 |       57 |
| Halobacterium        |        9 |       30 |
| Haloferax            |       16 |       41 |
| Halorubrum           |       32 |       84 |
| Halorussus           |       16 |       20 |
| Methanobacterium     |       14 |       25 |
| Methanobrevibacter   |       11 |       90 |
| Methanococcus        |        4 |       33 |
| Methanoculleus       |       16 |       21 |
| Methanomethylophilus |        1 |       25 |
| Methanosarcina       |       10 |       37 |
| Natrinema            |       22 |       33 |
| Sulfolobus           |        3 |       81 |
| Thermococcus         |       34 |       41 |

| #family                         | genus                         | species                                       | count |
|---------------------------------|-------------------------------|-----------------------------------------------|------:|
| Candidatus Methanomethylicaceae | Candidatus Methanosuratincola | Candidatus Methanosuratincola petrocarbonis   |    10 |
| Halobacteriaceae                | Halobacterium                 | Halobacterium salinarum                       |    12 |
| Haloferacaceae                  | Haloferax                     | Haloferax volcanii                            |    15 |
|                                 | Halorubrum                    | Halorubrum ezzemoulense                       |    38 |
| Methanobacteriaceae             | Methanobrevibacter            | Methanobrevibacter smithii                    |    64 |
| Methanococcaceae                | Methanococcus                 | Methanococcus maripaludis                     |    27 |
| Methanomassiliicoccaceae        | Methanomassiliicoccus         | Candidatus Methanomassiliicoccus intestinalis |    11 |
| Methanomethylophilaceae         | Methanomethylophilus          | Methanomethylophilus alvi                     |    25 |
| Methanosarcinaceae              | Methanosarcina                | Methanosarcina mazei                          |    13 |
|                                 |                               | Methanosarcina thermophila                    |    12 |
| NA                              | Candidatus Methanarcanum      | Candidatus Methanarcanum hacksteinii          |    12 |
| Sulfolobaceae                   | Metallosphaera                | Metallosphaera sedula                         |    11 |
|                                 | Saccharolobus                 | Saccharolobus solfataricus                    |    14 |
|                                 | Sulfolobus                    | Sulfolobus acidocaldarius                     |    58 |
|                                 |                               | Sulfolobus islandicus                         |    21 |
| Thermococcaceae                 | Pyrococcus                    | Pyrococcus kukulkanii                         |    11 |

### For *protein families*

```shell
cd ~/data/Archaea/

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:10" |
    rgr md stdin --num

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv

```

| item    | count |
|---------|------:|
| strain  |  1356 |
| species |   743 |
| genus   |   225 |
| family  |    75 |
| order   |    48 |
| class   |    26 |

| genus                 | #species | #strains |
|-----------------------|---------:|---------:|
| Haladaptatus          |        7 |       11 |
| Haloarcula            |       32 |       56 |
| Halobacterium         |        9 |       30 |
| Halobaculum           |       10 |       16 |
| Halobellus            |       11 |       17 |
| Haloferax             |       16 |       41 |
| Halogeometricum       |        6 |       11 |
| Halomarina            |        7 |       12 |
| Haloplanus            |        9 |       14 |
| Halorientalis         |        8 |       11 |
| Halorubrum            |       32 |       84 |
| Halorussus            |       16 |       20 |
| Metallosphaera        |        7 |       18 |
| Methanobacterium      |       14 |       24 |
| Methanobrevibacter    |       11 |       63 |
| Methanococcus         |        4 |       31 |
| Methanoculleus        |       16 |       21 |
| Methanohalophilus     |        6 |       14 |
| Methanolobus          |       10 |       11 |
| Methanomassiliicoccus |        2 |       12 |
| Methanomethylophilus  |        1 |       22 |
| Methanosarcina        |       10 |       35 |
| Methanothermobacter   |        5 |       14 |
| Natrinema             |       22 |       33 |
| Natronomonas          |        8 |       10 |
| Natronorubrum         |        9 |       13 |
| Nitrosopumilus        |       10 |       12 |
| Pyrococcus            |        5 |       18 |
| Saccharolobus         |        3 |       19 |
| Sulfolobus            |        3 |       80 |
| Thermococcus          |       34 |       41 |

## Collect proteins

```shell
cd ~/data/Archaea/

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --pro \
    --parallel 8 \
    --in ASSEMBLY/pass.lst \
    --not-in ASSEMBLY/omit.lst

# collect proteins
bash Protein/collect.sh

# clustering
# It may need to be run several times
bash Protein/cluster.sh

rm -fr Protein/tmp/

# info.tsv
bash Protein/info.sh

# counts
bash Protein/count.sh

cat Protein/counts.tsv |
    tsv-summarize -H --count --sum 2-7 |
    sed 's/^count/species/' |
    datamash transpose |
    (echo -e "#item\tcount" && cat) |
    rgr md stdin --fmt

```

| #item      |     count |
|------------|----------:|
| species    |       743 |
| strain_sum |     1,375 |
| total_sum  | 3,926,365 |
| dedup_sum  | 2,793,143 |
| rep_sum    | 2,426,768 |
| fam88_sum  | 2,330,424 |
| fam38_sum  | 1,960,961 |

## Phylogenetics with ar53

```shell
cd ~/data/Archaea/

# The Archaea HMM set
nwr kb ar53 -o HMM
cp HMM/ar53.lst HMM/marker.lst

```

### Find corresponding representative proteins by `hmmsearch`

```shell
cd ~/data/Archaea

cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 \
    > Protein/species-f.tsv

cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ -s Protein/"${SPECIES}"/ar53.tsv ]]; then
        continue
    fi
    if [[ ! -f Protein/"${SPECIES}"/rep_seq.fa.gz ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    cat HMM/marker.lst |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 "
            gzip -dcf Protein/${SPECIES}/rep_seq.fa.gz |
                hmmsearch --cut_nc --noali --notextw HMM/hmm/{}.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), q({}), \$1; '
        " \
        > Protein/${SPECIES}/ar53.tsv
done

cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ ! -s Protein/"${SPECIES}"/ar53.tsv ]]; then
        continue
    fi
    if [[ ! -f Protein/"${SPECIES}"/seq.sqlite ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    nwr seqdb -d Protein/${SPECIES} --rep f3=Protein/${SPECIES}/ar53.tsv

done

```

### Domain related protein sequences

```shell
cd ~/data/Archaea

mkdir -p Domain

# each assembly
cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ ! -f Protein/"${SPECIES}"/seq.sqlite ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    echo "
        SELECT
            seq.name,
            asm.name,
            rep.f3
        FROM asm_seq
        JOIN rep_seq ON asm_seq.seq_id = rep_seq.seq_id
        JOIN seq ON asm_seq.seq_id = seq.id
        JOIN rep ON rep_seq.rep_id = rep.id
        JOIN asm ON asm_seq.asm_id = asm.id
        WHERE 1=1
            AND rep.f3 IS NOT NULL
        ORDER BY
            asm.name,
            rep.f3
        " |
        sqlite3 -tabs Protein/${SPECIES}/seq.sqlite \
        > Protein/${SPECIES}/seq_asm_f3.tsv

    hnsm some Protein/"${SPECIES}"/pro.fa.gz <(
            tsv-select -f 1 Protein/"${SPECIES}"/seq_asm_f3.tsv |
            rgr dedup stdin
        )
done |
    hnsm dedup stdin |
    hnsm gz stdin -o Domain/ar53.fa

fd --full-path "Protein/.+/seq_asm_f3.tsv" -X cat \
    > Domain/seq_asm_f3.tsv

cat Domain/seq_asm_f3.tsv |
    tsv-join -e -d 2 -f summary/redundant.lst -k 1 \
    > Domain/seq_asm_f3.NR.tsv

```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Archaea

# Extract proteins
cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

        mkdir -p Domain/{}

        hnsm some Domain/ar53.fa.gz <(
            cat Domain/seq_asm_f3.tsv |
                tsv-filter --str-eq "3:{}" |
                tsv-select -f 1 |
                rgr dedup stdin
            ) \
            > Domain/{}/{}.pro.fa
    '

# Align each marker
cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo >&2 "==> marker [{}]"
        if [ ! -s Domain/{}/{}.pro.fa ]; then
            exit
        fi
        if [ -s Domain/{}/{}.aln.fa ]; then
            exit
        fi

#        muscle -quiet -in Domain/{}/{}.pro.fa -out Domain/{}/{}.aln.fa
        mafft --auto Domain/{}/{}.pro.fa > Domain/{}/{}.aln.fa
    '

cat HMM/marker.lst |
while read marker; do
    echo >&2 "==> marker [${marker}]"
    if [ ! -s Domain/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sometimes `muscle` can not produce alignments
    if [ ! -s Domain/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # Only NR strains
    # 1 name to many names
    cat Domain/seq_asm_f3.NR.tsv |
        tsv-filter --str-eq "3:${marker}" |
        tsv-select -f 1-2 |
        hnsm replace -s Domain/${marker}/${marker}.aln.fa stdin \
        > Domain/${marker}/${marker}.replace.fa
done

# Concat marker genes
cat HMM/marker.lst |
while read marker; do
    if [ ! -s Domain/${marker}/${marker}.pro.fa ]; then
        continue
    fi
    if [ ! -s Domain/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    cat Domain/${marker}/${marker}.replace.fa

    # empty line for .fas
    echo
done \
    > Domain/ar53.aln.fas

cat Domain/seq_asm_f3.NR.tsv |
    cut -f 2 |
    rgr dedup stdin |
    sort |
    fasops concat Domain/ar53.aln.fas stdin -o Domain/ar53.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Domain/ar53.aln.fa -out Domain/ar53.trim.fa -automated1

hnsm size Domain/ar53.*.fa |
    rgr dedup stdin -f 2 |
    cut -f 2
#45483
#12216

# To make it faster
FastTree -fastest -noml Domain/ar53.trim.fa > Domain/ar53.trim.newick

```

### Condense branches in the protein tree

```shell
cd ~/data/Archaea/tree

nw_reroot  ../Domain/ar53.trim.newick Saccharom_cere_S288C Saccharom_eub |
    nwr order stdin --nd --an \
    > ar53.reroot.newick

nwr pl-condense --map -r order -r family -r genus -r species \
    ar53.reroot.newick ../Count/species.tsv |
    nwr order stdin --nd --an \
    > ar53.condensed.newick

mv condensed.tsv ar53.condense.tsv

# pdf
nwr topo --bl ar53.condensed.newick | # remove comments
    nwr tex stdin --bl -o Archaea.ar53.tex

tectonic Archaea.ar53.tex

```

