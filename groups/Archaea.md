# Archaea

All genomes of *Archaea*

<!-- toc -->

- [Strain info](#strain-info)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
- [Download all assemblies](#download-all-assemblies)
    * [Create assembly.tsv](#create-assemblytsv)
    * [Count before download](#count-before-download)
    * [Download and check](#download-and-check)
    * [Rsync to hpcc](#rsync-to-hpcc)
- [BioSample](#biosample)
- [MinHash](#minhash)
    * [Condense branches in the minhash tree](#condense-branches-in-the-minhash-tree)
- [Count valid species and strains](#count-valid-species-and-strains)
    * [For *genomic alignments*](#for-genomic-alignments)
    * [For *protein families*](#for-protein-families)
- [Collect proteins](#collect-proteins)
- [Phylogenetics with ar53](#phylogenetics-with-ar53)
    * [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)

<!-- tocstop -->

## Strain info

* [Archaea](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2157)

### List all ranks

```shell
nwr member Archaea |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

```

| rank          | count |
|---------------|------:|
| superkingdom  |     1 |
| phylum        |    44 |
| order         |    66 |
| no rank       |   333 |
| species       |  3461 |
| class         |    39 |
| family        |    89 |
| genus         |   265 |
| clade         |    43 |
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
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#265 genus.list

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
#  586 RS1.tsv
#  603 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     1136
#GB1     1380

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
        AND genus IN ('Saccharomyces')
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession \
    > raw.tsv

# RS
SPECIES=$(
    cat RS1.tsv |
        cut -f 1 |
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
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

# Preference for refseq
cat raw.tsv |
    tsv-select -H -f "assembly_accession" \
    > rs.acc.tsv

# GB
SPECIES=$(
    cat GB1.tsv |
        cut -f 1 |
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
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

cat raw.tsv |
    tsv-uniq |
    datamash check
#1382 lines, 7 fields

# Create abbr.
cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
    tsv-select -f 1-6 |
    perl ~/Scripts/genomes/bin/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\tbiosample\tspecies\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[6]}++; # abbr_name
        $seen{$F[6]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\t%s\n}, $F[6], $F[3], $F[4], $F[1], $F[5];
        ' |
    tsv-filter --or --str-in-fld 2:ftp --str-in-fld 2:http |
    keep-header -- tsv-sort -k4,4 -k1,1 \
    > Archaea.assembly.tsv

datamash check < Archaea.assembly.tsv
#1379 lines, 5 fields

# find potential duplicate strains or assemblies
cat Archaea.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Archaea.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

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
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# genus.lst and genus.count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    keep-header -- tsv-sort -k1,1 |
    tsv-filter -H --ge 3:20 |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```

| item    | count |
|---------|------:|
| strain  |  1378 |
| species |   602 |
| genus   |   163 |
| family  |    46 |
| order   |    30 |
| class   |    16 |

| genus              | #species | #strains |
|--------------------|---------:|---------:|
| Haloarcula         |       15 |       29 |
| Halobacterium      |        7 |       24 |
| Haloferax          |       17 |       37 |
| Halorubrum         |       37 |       94 |
| Metallosphaera     |        8 |       20 |
| Methanobacterium   |       13 |       28 |
| Methanobrevibacter |       15 |      137 |
| Methanococcus      |        4 |       30 |
| Methanoculleus     |        8 |       29 |
| Methanosarcina     |       11 |      103 |
| Methanothrix       |        3 |       33 |
| Natrinema          |       22 |       34 |
| Sulfolobus         |        2 |       83 |
| Thermococcus       |       33 |       47 |

### Download and check

```shell
cd ~/data/Archaea

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
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

# N50 C S; create n50.tsv and n50.pass.tsv
bash ASSEMBLY/n50.sh 100000 500 500000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50,S" --max "C"
#N50_min S_min   C_max
#5379    668961  848

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "S:0.1,0.5" --quantile "N50:0.1,0.5"  --quantile "C:0.5,0.9"
#S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
#32255196.8      37316984        98936.2 1332095 167     1276.4

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
| url.tsv          |      3 | 1,378 |
| check.lst        |      1 | 1,378 |
| collect.tsv      |     20 | 1,379 |
| n50.tsv          |      4 | 1,379 |
| n50.pass.tsv     |      4 | 1,000 |
| collect.pass.tsv |     23 | 1,000 |
| pass.lst         |      1 |   999 |
| omit.lst         |      1 |   122 |
| rep.lst          |      1 |   536 |

### Rsync to hpcc

```bash
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

datamash check < BioSample/biosample.tsv
#1373 lines, 41 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## MinHash

```shell
cd ~/data/Archaea

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --mh \
    --parallel 8 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
    --ani-nr 0.005

# Compute assembly sketches
bash MinHash/compute.sh

# Distances within species
bash MinHash/species.sh

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#56

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
#  347 summary/NR.lst
#  236 summary/redundant.lst

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

nw_order -c n ../MinHash/tree.nwk \
    > minhash.sort.newick

# rank::col
ARRAY=(
    'order::5'
    'family::4'
    'genus::3'
#    'species::2'
)

rm minhash.condensed.map
CUR_TREE=minhash.sort.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/genomes/bin/condense_tree.sh ${CUR_TREE} ../Count/strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick minhash.${GROUP_NAME}.newick
    cat condense.map >> minhash.condensed.map

    CUR_TREE=minhash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 minhash.genus.newick |
    rsvg-convert -o Archaea.minhash.png

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
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# .lst and .count.tsv
bash Count/rank.sh

cat Count/order.count.tsv |
    tsv-filter -H --ge "3:20" |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:20" |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

# Can accept N_COUNT
bash Count/lineage.sh 5

cat Count/lineage.count.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| item    | count |
|---------|------:|
| strain  |   943 |
| species |   547 |
| genus   |   158 |
| family  |    45 |
| order   |    29 |
| class   |    15 |

| order              | #species | #strains |
|--------------------|---------:|---------:|
| Halobacteriales    |      122 |      175 |
| Haloferacales      |      102 |      174 |
| Methanobacteriales |       29 |       91 |
| Methanococcales    |       15 |       38 |
| Methanomicrobiales |       34 |       39 |
| Methanosarcinales  |       34 |       65 |
| Natrialbales       |       71 |       89 |
| Sulfolobales       |       25 |      127 |
| Thermococcales     |       39 |       48 |

| genus              | #species | #strains |
|--------------------|---------:|---------:|
| Haloarcula         |       14 |       26 |
| Haloferax          |       16 |       35 |
| Halorubrum         |       26 |       67 |
| Methanobacterium   |       13 |       23 |
| Methanobrevibacter |        9 |       54 |
| Methanococcus      |        4 |       26 |
| Methanosarcina     |       10 |       28 |
| Natrinema          |       20 |       27 |
| Sulfolobus         |        2 |       76 |
| Thermococcus       |       32 |       38 |

| #family             | genus              | species                     | count |
|---------------------|--------------------|-----------------------------|------:|
| Haloarculaceae      | Haloarcula         | Haloarcula hispanica        |     7 |
| Halobacteriaceae    | Halobacterium      | Halobacterium hubeiense     |     8 |
|                     |                    | Halobacterium salinarum     |     7 |
| Haloferacaceae      | Haloferax          | Haloferax volcanii          |     6 |
|                     | Halogeometricum    | Halogeometricum borinquense |     5 |
| Halorubraceae       | Halorubrum         | Halorubrum ezzemoulense     |    38 |
| Methanobacteriaceae | Methanobacterium   | Methanobacterium formicicum |     6 |
|                     | Methanobrevibacter | Methanobrevibacter smithii  |    37 |
|                     |                    | Methanobrevibacter woesei   |     5 |
| Methanococcaceae    | Methanococcus      | Methanococcus maripaludis   |    20 |
| Methanomicrobiaceae | Methanoculleus     | Methanoculleus bourgensis   |     5 |
| Methanosarcinaceae  | Methanosarcina     | Methanosarcina mazei        |    11 |
|                     |                    | Methanosarcina thermophila  |     5 |
| Sulfolobaceae       | Metallosphaera     | Metallosphaera sedula       |    10 |
|                     | Saccharolobus      | Saccharolobus solfataricus  |    13 |
|                     | Sulfolobus         | Sulfolobus acidocaldarius   |    55 |
|                     |                    | Sulfolobus islandicus       |    21 |

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
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:10" |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv

```

| item    | count |
|---------|------:|
| strain  |   928 |
| species |   544 |
| genus   |   157 |
| family  |    45 |
| order   |    29 |
| class   |    15 |

| genus              | #species | #strains |
|--------------------|---------:|---------:|
| Haloarcula         |       14 |       25 |
| Halobacterium      |        5 |       18 |
| Halobellus         |        9 |       11 |
| Haloferax          |       16 |       35 |
| Halomicroarcula    |        8 |       13 |
| Haloplanus         |        8 |       10 |
| Halorubrum         |       26 |       67 |
| Metallosphaera     |        7 |       17 |
| Methanobacterium   |       13 |       23 |
| Methanobrevibacter |        9 |       48 |
| Methanococcus      |        4 |       26 |
| Methanoculleus     |        8 |       13 |
| Methanohalophilus  |        6 |       14 |
| Methanosarcina     |       10 |       28 |
| Natrinema          |       20 |       27 |
| Natronorubrum      |        9 |       13 |
| Saccharolobus      |        3 |       18 |
| Sulfolobus         |        2 |       75 |
| Thermococcus       |       32 |       38 |

## Collect proteins

```shell
cd ~/data/Archaea/

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --pro \
    --in ASSEMBLY/pass.lst \
    --not-in summary/redundant.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst

# collect proteins
bash Protein/collect.sh

cat Protein/counts.tsv |
    mlr --itsv --omd cat

```

| #item                          | count     |
|--------------------------------|-----------|
| Proteins                       | 2,079,140 |
| Unique headers and annotations | 1,854,309 |
| Unique proteins                | 1,843,330 |
| all.replace.fa                 | 2,079,140 |
| all.annotation.tsv             | 2,079,141 |
| all.info.tsv                   | 2,079,141 |

## Phylogenetics with ar53

### Find corresponding proteins by `hmmsearch`

```shell
cd ~/data/Archaea

# The Archaea61 HMM set
nwr kb ar53 -o HMM
cp HMM/ar53.lst HMM/marker.lst

E_VALUE=1e-20

# Find all genes
for marker in $(cat HMM/marker.lst); do
    echo >&2 "==> marker [${marker}]"

    mkdir -p Protein/${marker}

    cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -f summary/redundant.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/hmm/${marker}.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/${marker}/replace.tsv

    echo >&2
done

```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Archaea

cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat Protein/{}/replace.tsv |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
#226     228     230

cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        cat Protein/{}/replace.tsv |
            wc -l
    ' |
    paste - - |
    tsv-filter --invert --ge 2:150 --le 2:300 |
    cut -f 1 \
    > Protein/marker.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process
cat HMM/marker.lst |
    grep -v -Fx -f Protein/marker.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

        cat Protein/{}/replace.tsv \
            > Protein/{}/{}.replace.tsv

        faops some Protein/all.uniq.fa.gz <(
            cat Protein/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > Protein/{}/{}.pro.fa
    '

# Align each markers with muscle
cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo >&2 "==> marker [{}]"
        if [ ! -s Protein/{}/{}.pro.fa ]; then
            exit
        fi
        if [ -s Protein/{}/{}.aln.fa ]; then
            exit
        fi

        muscle -quiet -in Protein/{}/{}.pro.fa -out Protein/{}/{}.aln.fa
    '

for marker in $(cat HMM/marker.lst); do
    echo >&2 "==> marker [${marker}]"
    if [ ! -s Protein/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sometimes `muscle` can not produce alignments
    if [ ! -s Protein/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # 1 name to many names
    cat Protein/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s Protein/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > Protein/${marker}/${marker}.replace.fa
done

# Concat marker genes
for marker in $(cat HMM/marker.lst); do
    if [ ! -s Protein/${marker}/${marker}.pro.fa ]; then
        continue
    fi
    if [ ! -s Protein/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # sequences in one line
    faops filter -l 0 Protein/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > Protein/ar53.aln.fas

cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -f summary/redundant.lst -k 1 |
    tsv-join -e -f MinHash/abnormal.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
    cut -f 1 |
    fasops concat Protein/ar53.aln.fas stdin -o Protein/ar53.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Protein/ar53.aln.fa -out Protein/ar53.trim.fa -automated1

faops size Protein/ar53.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#19766
#11266

# To make it faster
FastTree -fastest -noml Protein/ar53.trim.fa > Protein/ar53.trim.newick

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 Protein/ar53.trim.newick |
    rsvg-convert -o tree/Archaea.ar53.png

```
