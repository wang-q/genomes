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
    tva stats -H --min "N50" --max "C" --min "S" |
    tva transpose
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

fd --full-path "MinHash/.+/NR.lst" -X cat |
    sort |
    uniq \
    > summary/NR.lst
fd --full-path "MinHash/.+/redundant.lst" -X cat |
    sort |
    uniq \
    > summary/redundant.lst
wc -l summary/NR.lst summary/redundant.lst
# 1186 summary/NR.lst
#  629 summary/redundant.lst

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst |
    tva join -e -f summary/sp.lst \
    > MinHash/tmp.lst
mv MinHash/tmp.lst summary/abnormal.lst

wc -l MinHash/abnormal.lst summary/abnormal.lst
# 1429 MinHash/abnormal.lst
#  467 summary/abnormal.lst

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Archaea/

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --mh \
    --parallel 8 \
    --not-in summary/sp.lst \
    --not-in summary/omit.lst \
    --not-in summary/abnormal.lst \
    --not-in summary/redundant.lst \
    --height 0.4

bash MinHash/dist.sh
```

### Condense branches in the minhash tree

```shell
mkdir -p ~/data/Archaea/tree
cd ~/data/Archaea/tree

pgr nwk reroot ../MinHash/tree.nwk -n Saccharom_cere_S288C |
    pgr nwk order stdin --nd --an \
    > minhash.reroot.newick

pgr pl condense --map -t ../Count/strains.taxon.tsv -r 5 -r 4 -r 3 \
    minhash.reroot.newick |
    pgr nwk order stdin --nd --an \
    > minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

# svg
pgr nwk to-svg minhash.condensed.newick \
    > Archaea.minhash.svg
```

## Count valid species and strains

### For *genomic alignments*

```shell
cd ~/data/Archaea/

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --count \
    --in summary/pass.lst \
    --not-in summary/abnormal.lst \
    --rank order --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    tva to md --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/order.count.tsv |
    tva filter -H --ge "3:20" |
    tva to md --num

cat Count/genus.count.tsv |
    tva filter -H --ge "3:20" |
    tva to md --num

# Can accept N_COUNT
bash Count/lineage.sh 15

cat Count/lineage.count.tsv |
    tva to md --num

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv
cp Count/genus.count.tsv summary/genus.genome.tsv
```

| item    | count |
| ------- | ----: |
| strain  |  2093 |
| species |   944 |
| genus   |   259 |
| family  |    90 |
| order   |    60 |
| class   |    36 |

| order                   | #species | #strains |
| ----------------------- | -------: | -------: |
| Desulfurococcales       |       31 |       75 |
| Halobacteriales         |      440 |      766 |
| Methanobacteriales      |       47 |      226 |
| Methanococcales         |       22 |       54 |
| Methanomassiliicoccales |        9 |       61 |
| Methanomicrobiales      |       64 |      111 |
| Methanosarcinales       |       56 |      139 |
| Sulfolobales            |       32 |      213 |
| Thermococcales          |       45 |       75 |
| Thermoproteales         |       20 |       39 |

| genus                | #species | #strains |
| -------------------- | -------: | -------: |
| Haloarcula           |       34 |       68 |
| Halobacterium        |       10 |       36 |
| Halobellus           |       13 |       20 |
| Haloferax            |       17 |       52 |
| Halorubrum           |       33 |       98 |
| Halorussus           |       17 |       24 |
| Ignisphaera          |        3 |       21 |
| Metallosphaera       |        7 |       62 |
| Methanobacterium     |       21 |       64 |
| Methanobrevibacter   |       13 |      108 |
| Methanococcus        |        5 |       34 |
| Methanoculleus       |       19 |       34 |
| Methanomethylophilus |        2 |       28 |
| Methanosarcina       |       16 |       69 |
| Methanothermobacter  |        7 |       38 |
| Natrinema            |       23 |       39 |
| Pyrococcus           |        7 |       23 |
| Saccharolobus        |        5 |       66 |
| Sulfolobus           |        2 |       60 |
| Thermococcus         |       35 |       49 |

| #family                    | genus                   | species                     | count |
| -------------------------- | ----------------------- | --------------------------- | ----: |
| Candidatus Iainarchaeaceae | Candidatus Iainarchaeum | Candidatus Iainarchaeum sp. |    17 |
| Desulfurococcaceae         | Ignisphaera             | Ignisphaera sp.             |    19 |
| Haloferacaceae             | Haloferax               | Haloferax volcanii          |    15 |
|                            | Halorubrum              | Halorubrum ezzemoulense     |    38 |
| Methanobacteriaceae        | Methanobacterium        | Methanobacterium sp.        |    33 |
|                            | Methanobrevibacter      | Methanobrevibacter smithii  |    74 |
|                            | Methanothermobacter     | Methanothermobacter sp.     |    16 |
| Methanococcaceae           | Methanococcus           | Methanococcus maripaludis   |    27 |
| Methanomethylophilaceae    | Methanomethylophilus    | Methanomethylophilus alvi   |    25 |
| Methanosarcinaceae         | Methanosarcina          | Methanosarcina sp.          |    15 |
| Sulfolobaceae              | Metallosphaera          | Metallosphaera sp.          |    44 |
|                            | Saccharolobus           | Saccharolobus islandicus    |    21 |
|                            |                         | Saccharolobus sp.           |    26 |
|                            | Sulfolobus              | Sulfolobus acidocaldarius   |    58 |

### For *protein families*

```shell
cd ~/data/Archaea/

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --count \
    --in summary/pass.lst \
    --not-in summary/abnormal.lst \
    --not-in summary/omit.lst \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    tva to md --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    tva filter -H --ge "3:50" |
    tva to md --num

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv
```

| item    | count |
| ------- | ----: |
| strain  |  1858 |
| species |   908 |
| genus   |   246 |
| family  |    85 |
| order   |    56 |
| class   |    34 |

| genus              | #species | #strains |
| ------------------ | -------: | -------: |
| Haloarcula         |       34 |       67 |
| Haloferax          |       17 |       50 |
| Halorubrum         |       33 |       98 |
| Metallosphaera     |        7 |       62 |
| Methanobrevibacter |       13 |       67 |
| Methanosarcina     |       16 |       57 |
| Saccharolobus      |        5 |       66 |
| Sulfolobus         |        2 |       59 |

## Collect proteins

```shell
cd ~/data/Archaea/

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --pro \
    --parallel 8 \
    --in summary/pass.lst \
    --not-in summary/omit.lst

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
    tva stats -H --count --sum 2-7 |
    sed 's/^count/species/' |
    tva transpose |
    (echo -e "#item\tcount" && cat) |
    tva to md --fmt
```

| #item      |     count |
| ---------- | --------: |
| species    |       915 |
| strain_sum |     1,994 |
| total_sum  | 5,386,858 |
| dedup_sum  | 4,235,285 |
| rep_sum    | 3,361,076 |
| fam88_sum  | 3,069,764 |
| fam38_sum  | 2,459,388 |

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
    tva join -f summary/pass.lst -k 1 |
    tva join -e -f summary/omit.lst -k 1 \
    > Protein/species-f.tsv

cat Protein/species-f.tsv |
    tva select -f 2 |
    tva uniq |
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
    tva select -f 2 |
    tva uniq |
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
    tva select -f 2 |
    tva uniq |
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

    pgr fa some Protein/"${SPECIES}"/pro.fa.gz <(
            tva select -f 1 Protein/"${SPECIES}"/seq_asm_f3.tsv |
            tva uniq
        )
done |
    pgr fa dedup stdin |
    pgr fa gz stdin -o Domain/ar53.fa.gz

fd --full-path "Protein/.+/seq_asm_f3.tsv" -X cat \
    > Domain/seq_asm_f3.tsv

cat Domain/seq_asm_f3.tsv |
    tva join -e -d 2 -f summary/redundant.lst -k 1 \
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

        pgr fa some Domain/ar53.fa.gz <(
            cat Domain/seq_asm_f3.tsv |
                tva filter --str-eq "3:{}" |
                tva select -f 1 |
                tva uniq
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
        tva filter --str-eq "3:${marker}" |
        tva select -f 1-2 |
        pgr fa replace -s Domain/${marker}/${marker}.aln.fa stdin \
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
    tva uniq |
    sort |
    fasops concat Domain/ar53.aln.fas stdin -o Domain/ar53.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Domain/ar53.aln.fa -out Domain/ar53.trim.fa -automated1

pgr fa size Domain/ar53.*.fa |
    tva uniq -f 2 |
    cut -f 2
# 45483
# 11499

# To make it faster
FastTree -fastest -noml Domain/ar53.trim.fa > Domain/ar53.trim.newick
```

### Condense branches in the protein tree

```shell
cd ~/data/Archaea/tree

pgr nwk reroot ../Domain/ar53.trim.newick -n Saccharom_cere_S288C |
    pgr nwk order stdin --nd --an \
    > ar53.reroot.newick

pgr pl condense --map -t ../Count/strains.taxon.tsv -r 5 -r 4 -r 3 -r 2 \
    ar53.reroot.newick |
    pgr nwk order stdin --nd --an \
    > ar53.condensed.newick

mv condensed.tsv ar53.condense.tsv

# svg
pgr nwk to-svg ar53.condensed.newick \
    > Archaea.ar53.svg
```
