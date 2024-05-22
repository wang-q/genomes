# *Bacillus*

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

<!-- tocstop -->

## Strain info

* [Bacillus](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1386)
* [Paenibacillus](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=44249)

A recent phylogenomics [study](https://doi.org/10.1016/j.syapm.2018.10.007) has altered our
understanding of the order Bacillales.

![The order Bacillales](../images/1-s2.0-S0723202018303291-gr3.png)

We include the following families:

* *Bacillaceae*
* *Paenibacillaceae*
* *Sporolactobacillaceae*
* *Thermoactinomycetaceae*
* *Alicyclobacillaceae*

### List all ranks

```shell
mkdir -p ~/data/Bacillus
cd ~/data/Bacillus

nwr member Bacillaceae Paenibacillaceae Sporolactobacillaceae Thermoactinomycetaceae Alicyclobacillaceae |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

nwr member Bacillaceae Paenibacillaceae Sporolactobacillaceae Thermoactinomycetaceae Alicyclobacillaceae -r "species group" -r "species subgroup" |
    tsv-select -f 1-3 |
    keep-header -- tsv-sort -k3,3 -k2,2 |
    mlr --itsv --omd cat

```

| rank             | count |
|------------------|------:|
| family           |     5 |
| genus            |   186 |
| species          | 44915 |
| subspecies       |    42 |
| no rank          |   254 |
| strain           |   751 |
| species group    |     5 |
| species subgroup |     2 |
| biotype          |     1 |

| #tax_id | sci_name                          | rank             |
|---------|-----------------------------------|------------------|
| 1792192 | Bacillus altitudinis complex      | species group    |
| 86661   | Bacillus cereus group             | species group    |
| 653685  | Bacillus subtilis group           | species group    |
| 1505648 | Geobacillus thermoleovorans group | species group    |
| 2044880 | Paenibacillus sonchi group        | species group    |
| 1938374 | Bacillus amyloliquefaciens group  | species subgroup |
| 653388  | Bacillus mojavensis subgroup      | species subgroup |

### Species with assemblies

* 'RefSeq'
* 'Genbank'

```shell
mkdir -p ~/data/Bacillus/summary
cd ~/data/Bacillus/summary

# should have a valid name of genus
nwr member \
    Bacillaceae Paenibacillaceae Sporolactobacillaceae Thermoactinomycetaceae Alicyclobacillaceae \
    -r genus |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#186 genus.list.tsv

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
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB1.tsv

wc -l RS*.tsv GB*.tsv
#  3328 RS1.tsv
#  3690 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     13439
#GB1     17569

```

## Download all assemblies

### Create assembly.tsv

If a refseq assembly is available, the corresponding genbank one is not downloaded

```shell
cd ~/data/Bacillus/summary

echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus IN ('Bacillus', 'Staphylococcus')
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
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% x %'
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

echo "
    SELECT
        genus || ' sp. ' || infraspecific_name || ' ' || assembly_accession AS name,
        genus || ' sp. ', genus, ftp_path, biosample, assembly_level,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species LIKE '% sp.%'
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
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% x %'
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

cat raw.tsv |
    tsv-uniq |
    datamash check
#17520 lines, 7 fields

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
    > Bacillus.assembly.tsv

datamash check < Bacillus.assembly.tsv
#17511 lines, 5 fields

# find potential duplicate strains or assemblies
cat Bacillus.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Bacillus.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

# Edit .assembly.tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Bacillus.assembly.tsv
#
# Save the file to another directory to prevent accidentally changing it
# cp Bacillus.assembly.tsv ~/Scripts/genomes/assembly

# Cleaning
rm raw*.*sv

```

### Count before download

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class

```shell
cd ~/data/Bacillus

nwr template ~/Scripts/genomes/assembly/Bacillus.assembly.tsv \
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
    tsv-filter -H --ge 3:50 |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```

| item    | count |
|---------|------:|
| strain  | 17409 |
| species |  3594 |
| genus   |   169 |
| family  |     6 |
| order   |     1 |
| class   |     1 |

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Alicyclobacillus |       30 |       65 |
| Anoxybacillus    |       35 |      105 |
| Bacillus         |     1413 |    12180 |
| Brevibacillus    |       66 |      248 |
| Cytobacillus     |       39 |      131 |
| Geobacillus      |       61 |      177 |
| Halobacillus     |       35 |       55 |
| Heyndrickxia     |       23 |      176 |
| Lysinibacillus   |      107 |      350 |
| Metabacillus     |       27 |       55 |
| Neobacillus      |       45 |      112 |
| Niallia          |       24 |       69 |
| Oceanobacillus   |       54 |       87 |
| Paenibacillus    |      806 |     1453 |
| Peribacillus     |       40 |      171 |
| Priestia         |       28 |      571 |
| Rossellomorea    |       15 |       63 |
| Shouchella       |       11 |       73 |
| Virgibacillus    |       47 |      109 |

### Download and check

```shell
cd ~/data/Bacillus

nwr template ~/Scripts/genomes/assembly/Bacillus.assembly.tsv \
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
bash ASSEMBLY/n50.sh 50000 500 500000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50,S" --max "C"
#N50_min S_min   C_max
#0       0       1976

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "S:0.1,0.5" --quantile "N50:0.1,0.5"  --quantile "C:0.5,0.9"
#S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
#3748962.2       5271792 61791.8 316786  57      262

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

| #item            | fields |  lines |
|------------------|-------:|-------:|
| url.tsv          |      3 | 12,787 |
| check.lst        |      1 | 12,787 |
| collect.tsv      |     20 | 12,788 |
| n50.tsv          |      4 | 12,788 |
| n50.pass.tsv     |      4 | 11,764 |
| collect.pass.tsv |     23 | 11,764 |
| pass.lst         |      1 | 11,763 |
| omit.lst         |      1 |    431 |
| rep.lst          |      1 |  1,067 |

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Bacillus/ \
    wangq@202.119.37.251:data/Bacillus

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Bacillus/ \
    wangq@58.213.64.36:data/Bacillus

# rsync -avP wangq@202.119.37.251:data/Bacillus/ ~/data/Bacillus

# rsync -avP -e 'ssh -p 8804' wangq@58.213.64.36:data/Bacillus/ ~/data/Bacillus

```

## BioSample

```shell
cd ~/data/Bacillus

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Bacillus.assembly.tsv \
    --bs

# Run this script twice and it will re-download the failed files
bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 50

datamash check < BioSample/biosample.tsv
#12766 lines, 84 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## MinHash

```shell
cd ~/data/Bacillus

nwr template ~/Scripts/genomes/assembly/Bacillus.assembly.tsv \
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
#540

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
#  4072 summary/NR.lst
#  5475 summary/redundant.lst

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Bacillus/

nwr template ~/Scripts/genomes/assembly/Bacillus.assembly.tsv \
    --mh \
    --parallel 8 \
    --not-in MinHash/abnormal.lst \
    --not-in summary/redundant.lst \
    --height 0.4

bash MinHash/dist.sh

```

### Condense branches in the minhash tree

```shell
mkdir -p ~/data/Bacillus/tree
cd ~/data/Bacillus/tree

nwr order --an --nd ../MinHash/tree.nwk -o minhash.sort.newick

nwr pl-condense --map -r genus \
    minhash.sort.newick ../Count/species.tsv \
    -o minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

# png
nw_display -s -b 'visibility:hidden' -w 1000 -v 10 minhash.condensed.newick |
    rsvg-convert -o Bacillus.minhash.png

```

## Count valid species and strains

### For *genomic alignments*

```shell
cd ~/data/Bacillus/

nwr template ~/Scripts/genomes/assembly/Bacillus.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --rank family --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# .lst and .count.tsv
bash Count/rank.sh

cat Count/family.count.tsv |
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
| strain  | 11156 |
| species |  2536 |
| genus   |   155 |
| family  |     5 |
| order   |     1 |
| class   |     1 |

| family                 | #species | #strains |
|------------------------|---------:|---------:|
| Alicyclobacillaceae    |       42 |       55 |
| Bacillaceae            |     1750 |     9862 |
| Paenibacillaceae       |      666 |     1136 |
| Sporolactobacillaceae  |       27 |       35 |
| Thermoactinomycetaceae |       51 |       68 |

| genus              | #species | #strains |
|--------------------|---------:|---------:|
| Alicyclobacillus   |       25 |       32 |
| Alkalihalobacillus |       27 |       68 |
| Anoxybacillus      |       27 |       45 |
| Bacillus           |      973 |     8128 |
| Brevibacillus      |       58 |      153 |
| Caldibacillus      |        5 |       20 |
| Cohnella           |       36 |       39 |
| Cytobacillus       |       18 |       64 |
| Fictibacillus      |       19 |       23 |
| Geobacillus        |       47 |       88 |
| Gracilibacillus    |       22 |       25 |
| Halobacillus       |       32 |       43 |
| Lysinibacillus     |       75 |      183 |
| Metabacillus       |       22 |       34 |
| Neobacillus        |       28 |       36 |
| Niallia            |       10 |       27 |
| Oceanobacillus     |       37 |       55 |
| Paenibacillus      |      534 |      892 |
| Parageobacillus    |        9 |       31 |
| Peribacillus       |       24 |       90 |
| Priestia           |       18 |      332 |
| Psychrobacillus    |       18 |       22 |
| Rossellomorea      |        9 |       25 |
| Sporolactobacillus |       17 |       24 |
| Thermoactinomyces  |       17 |       26 |
| Virgibacillus      |       37 |       74 |
| Weizmannia         |        7 |       54 |

| #family                | genus              | species                             | count |
|------------------------|--------------------|-------------------------------------|------:|
| Alicyclobacillaceae    | Kyrpidia           | Kyrpidia spormannii                 |     5 |
| Bacillaceae            | Alkalihalobacillus | Alkalihalobacillus clausii          |    38 |
|                        | Anoxybacillus      | Anoxybacillus flavithermus          |     6 |
|                        | Bacillus           | Bacillus aerophilus                 |     7 |
|                        |                    | Bacillus albus                      |    15 |
|                        |                    | Bacillus altitudinis                |   148 |
|                        |                    | Bacillus amyloliquefaciens          |   171 |
|                        |                    | Bacillus anthracis                  |   677 |
|                        |                    | Bacillus atrophaeus                 |    88 |
|                        |                    | Bacillus badius                     |     9 |
|                        |                    | Bacillus bombysepticus              |     5 |
|                        |                    | Bacillus cereus                     |  1695 |
|                        |                    | Bacillus cytotoxicus                |    38 |
|                        |                    | Bacillus glycinifermentans          |     9 |
|                        |                    | Bacillus halotolerans               |    50 |
|                        |                    | Bacillus haynesii                   |    93 |
|                        |                    | Bacillus inaquosorum                |   100 |
|                        |                    | Bacillus infantis                   |     6 |
|                        |                    | Bacillus intestinalis               |     7 |
|                        |                    | Bacillus licheniformis              |   263 |
|                        |                    | Bacillus methanolicus               |     5 |
|                        |                    | Bacillus mobilis                    |    25 |
|                        |                    | Bacillus mojavensis                 |    11 |
|                        |                    | Bacillus mycoides                   |   135 |
|                        |                    | Bacillus nitratireducens            |     9 |
|                        |                    | Bacillus pacificus                  |    46 |
|                        |                    | Bacillus paralicheniformis          |   130 |
|                        |                    | Bacillus paranthracis               |   130 |
|                        |                    | Bacillus pseudomycoides             |    99 |
|                        |                    | Bacillus pumilus                    |   202 |
|                        |                    | Bacillus safensis                   |   148 |
|                        |                    | Bacillus siamensis                  |    13 |
|                        |                    | Bacillus sonorensis                 |    23 |
|                        |                    | Bacillus sp. WP8                    |     6 |
|                        |                    | Bacillus spizizenii                 |   126 |
|                        |                    | Bacillus stercoris                  |     6 |
|                        |                    | Bacillus stratosphericus            |    11 |
|                        |                    | Bacillus subtilis                   |   721 |
|                        |                    | Bacillus swezeyi                    |     5 |
|                        |                    | Bacillus thuringiensis              |   785 |
|                        |                    | Bacillus toyonensis                 |   279 |
|                        |                    | Bacillus tropicus                   |    34 |
|                        |                    | Bacillus vallismortis               |    11 |
|                        |                    | Bacillus velezensis                 |   613 |
|                        |                    | Bacillus wiedmannii                 |   194 |
|                        | Caldibacillus      | Caldibacillus debilis               |     5 |
|                        |                    | Caldibacillus thermoamylovorans     |    10 |
|                        | Cytobacillus       | Cytobacillus firmus                 |    26 |
|                        |                    | Cytobacillus kochii                 |     5 |
|                        |                    | Cytobacillus oceanisediminis        |    10 |
|                        | Geobacillus        | Geobacillus kaustophilus            |     6 |
|                        |                    | Geobacillus stearothermophilus      |     8 |
|                        |                    | Geobacillus thermodenitrificans     |     8 |
|                        |                    | Geobacillus thermoleovorans         |    11 |
|                        | Heyndrickxia       | Heyndrickxia oleronia               |    12 |
|                        | Lysinibacillus     | Lysinibacillus boronitolerans       |     6 |
|                        |                    | Lysinibacillus capsici              |    14 |
|                        |                    | Lysinibacillus fusiformis           |    42 |
|                        |                    | Lysinibacillus sphaericus           |    34 |
|                        | Niallia            | Niallia circulans                   |    14 |
|                        |                    | Niallia taxi                        |     5 |
|                        | Parageobacillus    | Parageobacillus caldoxylosilyticus  |     8 |
|                        |                    | Parageobacillus thermoglucosidasius |    13 |
|                        | Peribacillus       | Peribacillus butanolivorans         |     5 |
|                        |                    | Peribacillus frigoritolerans        |    40 |
|                        |                    | Peribacillus simplex                |    19 |
|                        | Priestia           | Priestia aryabhattai                |    58 |
|                        |                    | Priestia endophytica                |    11 |
|                        |                    | Priestia filamentosa                |     6 |
|                        |                    | Priestia flexa                      |    23 |
|                        |                    | Priestia megaterium                 |   218 |
|                        | Rossellomorea      | Rossellomorea marisflavi            |     9 |
|                        |                    | Rossellomorea vietnamensis          |     5 |
|                        | Salipaludibacillus | Salipaludibacillus agaradhaerens    |     5 |
|                        | Siminovitchia      | Siminovitchia fortis                |     9 |
|                        | Terribacillus      | Terribacillus saccharophilus        |     7 |
|                        | Virgibacillus      | Virgibacillus halodenitrificans     |    12 |
|                        |                    | Virgibacillus pantothenticus        |    16 |
|                        | Weizmannia         | Weizmannia coagulans                |    43 |
| Paenibacillaceae       | Aneurinibacillus   | Aneurinibacillus migulanus          |     8 |
|                        | Brevibacillus      | Brevibacillus agri                  |     8 |
|                        |                    | Brevibacillus borstelensis          |    15 |
|                        |                    | Brevibacillus brevis                |    14 |
|                        |                    | Brevibacillus formosus              |     6 |
|                        |                    | Brevibacillus laterosporus          |    38 |
|                        |                    | Brevibacillus parabrevis            |     7 |
|                        | Paenibacillus      | Paenibacillus alvei                 |    16 |
|                        |                    | Paenibacillus amylolyticus          |     9 |
|                        |                    | Paenibacillus apiarius              |     9 |
|                        |                    | Paenibacillus dendritiformis        |     8 |
|                        |                    | Paenibacillus elgii                 |     8 |
|                        |                    | Paenibacillus glucanolyticus        |     6 |
|                        |                    | Paenibacillus ihbetae               |     8 |
|                        |                    | Paenibacillus lactis                |     5 |
|                        |                    | Paenibacillus larvae                |    32 |
|                        |                    | Paenibacillus lautus                |     9 |
|                        |                    | Paenibacillus macerans              |     9 |
|                        |                    | Paenibacillus mucilaginosus         |     5 |
|                        |                    | Paenibacillus odorifer              |    33 |
|                        |                    | Paenibacillus peoriae               |     8 |
|                        |                    | Paenibacillus polymyxa              |    89 |
|                        |                    | Paenibacillus silvae                |     5 |
|                        |                    | Paenibacillus thiaminolyticus       |    23 |
|                        |                    | Paenibacillus xylanexedens          |    10 |
| Sporolactobacillaceae  | Sporolactobacillus | Sporolactobacillus terrae           |     7 |
| Thermoactinomycetaceae | Thermoactinomyces  | Thermoactinomyces vulgaris          |     7 |

### For *protein families*

```shell
cd ~/data/Bacillus/

nwr template ~/Scripts/genomes/assembly/Bacillus.assembly.tsv \
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
| strain  | 10871 |
| species |  2470 |
| genus   |   155 |
| family  |     5 |
| order   |     1 |
| class   |     1 |

| genus              | #species | #strains |
|--------------------|---------:|---------:|
| Alicyclobacillus   |       25 |       32 |
| Alkalihalobacillus |       26 |       66 |
| Aneurinibacillus   |        8 |       18 |
| Anoxybacillus      |       26 |       44 |
| Bacillus           |      938 |     7917 |
| Brevibacillus      |       58 |      152 |
| Caldibacillus      |        5 |       20 |
| Cohnella           |       36 |       39 |
| Cytobacillus       |       18 |       64 |
| Domibacillus       |       10 |       13 |
| Fictibacillus      |       19 |       23 |
| Geobacillus        |       45 |       86 |
| Gracilibacillus    |       19 |       22 |
| Halobacillus       |       32 |       43 |
| Heyndrickxia       |        3 |       16 |
| Lederbergia        |        9 |       15 |
| Lentibacillus      |       14 |       15 |
| Lysinibacillus     |       67 |      159 |
| Mesobacillus       |       12 |       19 |
| Metabacillus       |       22 |       34 |
| Neobacillus        |       28 |       36 |
| Niallia            |       10 |       26 |
| Oceanobacillus     |       36 |       54 |
| Ornithinibacillus  |       11 |       11 |
| Paenibacillus      |      526 |      876 |
| Parageobacillus    |        9 |       31 |
| Paraliobacillus    |        7 |       10 |
| Peribacillus       |       24 |       88 |
| Pontibacillus      |        8 |       10 |
| Priestia           |       18 |      323 |
| Psychrobacillus    |       18 |       22 |
| Rossellomorea      |        9 |       25 |
| Saccharibacillus   |       11 |       14 |
| Salipaludibacillus |        7 |       12 |
| Siminovitchia      |        5 |       16 |
| Sporolactobacillus |       17 |       24 |
| Terribacillus      |        6 |       13 |
| Thermoactinomyces  |       16 |       25 |
| Virgibacillus      |       37 |       74 |
| Weizmannia         |        7 |       54 |
