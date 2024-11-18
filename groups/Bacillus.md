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
* *Planococcaceae*
* *Pasteuriaceae*

* *Desulfuribacillaceae* as closet outgroups

According to a recent [publication](https://doi.org/10.1007/s10482-023-01857-6), the families below
are found to be more closely related to *Staphylococcaceae*:

- *Abyssicoccaceae*
- *Gemellaceae*
- *Listeriaceae*
- *Salinicoccaceae*

### List all ranks

```shell
mkdir -p ~/data/Bacillus
cd ~/data/Bacillus

nwr member Bacillales -r family |
    keep-header -- tsv-sort -k2,2 |
    rgr md stdin --num

nwr member \
    Bacillaceae Paenibacillaceae \
    Sporolactobacillaceae Thermoactinomycetaceae Alicyclobacillaceae \
    Planococcaceae Pasteuriaceae \
    Desulfuribacillaceae |
    tsv-summarize -H -g 3 --count |
    rgr md stdin --fmt

nwr member \
    Bacillaceae Paenibacillaceae \
    Sporolactobacillaceae Thermoactinomycetaceae Alicyclobacillaceae \
    Planococcaceae Pasteuriaceae \
    Desulfuribacillaceae \
    -r "species group" -r "species subgroup" |
    tsv-select -f 1-3 |
    keep-header -- tsv-sort -k3,3 -k2,2 |
    rgr md stdin --num

```

| #tax_id | sci_name                            | rank   | division |
|--------:|-------------------------------------|--------|----------|
| 3076164 | Abyssicoccaceae                     | family | Bacteria |
|  186823 | Alicyclobacillaceae                 | family | Bacteria |
| 3120669 | Anoxybacillaceae                    | family | Bacteria |
|  186817 | Bacillaceae                         | family | Bacteria |
|  539003 | Bacillales Family X. Incertae Sedis | family | Bacteria |
|  539738 | Gemellaceae                         | family | Bacteria |
|  186820 | Listeriaceae                        | family | Bacteria |
|  186822 | Paenibacillaceae                    | family | Bacteria |
|  538998 | Pasteuriaceae                       | family | Bacteria |
|  186818 | Planococcaceae                      | family | Bacteria |
| 3076165 | Salinicoccaceae                     | family | Bacteria |
|  186821 | Sporolactobacillaceae               | family | Bacteria |
|   90964 | Staphylococcaceae                   | family | Bacteria |
|  186824 | Thermoactinomycetaceae              | family | Bacteria |

| rank             |  count |
|------------------|-------:|
| family           |      8 |
| genus            |    212 |
| species          | 47,021 |
| subspecies       |     42 |
| no rank          |    279 |
| strain           |    769 |
| species group    |      5 |
| species subgroup |      2 |
| biotype          |      1 |

| #tax_id | sci_name                          | rank             |
|--------:|-----------------------------------|------------------|
| 1792192 | Bacillus altitudinis complex      | species group    |
|   86661 | Bacillus cereus group             | species group    |
|  653685 | Bacillus subtilis group           | species group    |
| 1505648 | Geobacillus thermoleovorans group | species group    |
| 2044880 | Paenibacillus sonchi group        | species group    |
| 1938374 | Bacillus amyloliquefaciens group  | species subgroup |
|  653388 | Bacillus mojavensis subgroup      | species subgroup |

### Species with assemblies

* 'RefSeq'
* 'Genbank'

```shell
mkdir -p ~/data/Bacillus/summary
cd ~/data/Bacillus/summary

# should have a valid name of genus
nwr member \
    Bacillaceae Paenibacillaceae \
    Sporolactobacillaceae Thermoactinomycetaceae Alicyclobacillaceae \
    Planococcaceae Pasteuriaceae \
    Desulfuribacillaceae \
    -r genus |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#212 genus.list.tsv

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
#  3940 RS1.tsv
#  4332 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     15147
#GB1     20472

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
        AND genus IN ('Bacillus', 'Staphylococcus', 'Listeria')
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-join -H -d assembly_accession -f ~/Scripts/genomes/assembly/Bacteria.reference.tsv -k assembly_accession |
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
        genus || ' sp.', genus, ftp_path, biosample, assembly_level,
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

echo "
    SELECT
        genus || ' sp. ' || infraspecific_name || ' ' || assembly_accession AS name,
        genus || ' sp.', genus, ftp_path, biosample, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species LIKE '% sp.%'
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

cat raw.tsv |
    rgr dedup stdin |
    datamash check
#20478 lines, 7 fields

# Create abbr.
cat raw.tsv |
    grep -v '^#' |
    rgr dedup stdin |
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
#20466 lines, 5 fields

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
    rgr md stdin --fmt

# .lst and .count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    keep-header -- tsv-sort -k1,1 |
    tsv-filter -H --ge 3:50 |
    rgr md stdin --num

```

| item    |  count |
|---------|-------:|
| strain  | 18,912 |
| species |  1,541 |
| genus   |    199 |
| family  |     10 |
| order   |      2 |
| class   |      2 |

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Alicyclobacillus |       28 |       72 |
| Anoxybacillus    |       18 |      113 |
| Bacillus         |      136 |    12019 |
| Brevibacillus    |       30 |      270 |
| Cohnella         |       37 |       66 |
| Cytobacillus     |       19 |      154 |
| Fictibacillus    |       16 |       56 |
| Geobacillus      |       17 |      187 |
| Halobacillus     |       25 |       62 |
| Heyndrickxia     |       13 |      196 |
| Kurthia          |       10 |       52 |
| Lysinibacillus   |       28 |      438 |
| Metabacillus     |       22 |       62 |
| Neobacillus      |       27 |      127 |
| Niallia          |        8 |       77 |
| Oceanobacillus   |       35 |      105 |
| Paenibacillus    |      316 |     1975 |
| Peribacillus     |       19 |      259 |
| Planococcus      |       29 |       81 |
| Priestia         |       10 |      654 |
| Psychrobacillus  |       10 |       56 |
| Rossellomorea    |        7 |       78 |
| Shouchella       |       12 |       77 |
| Solibacillus     |       10 |       52 |
| Sporosarcina     |       27 |      133 |
| Virgibacillus    |       34 |      124 |

### Download and check

```shell
cd ~/data/Bacillus

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Bacillus.assembly.tsv \
    --ass

# Run
bash ASSEMBLY/rsync.sh

# Check md5; create check.lst
# rm ASSEMBLY/check.lst
bash ASSEMBLY/check.sh

## Put the misplaced directories into the right ones
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

# N50 C S; create n50.tsv and n50.pass.tsv
bash ASSEMBLY/n50.sh 50000 500 500000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50" --max "C" --min "S"
#N50_min S_min   C_max
#5966    2222169 1976

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "N50:0.1,0.5" --quantile "C:0.5,0.9" --quantile "S:0.1,0.5" |
    datamash transpose
#S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
#3752678.5       5249709 58246.7 315670.5        55      274

# After the above steps are completed, run the following commands.

# Collect; create collect.tsv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    rgr md stdin --fmt

```

| #item            | fields |  lines |
|------------------|-------:|-------:|
| url.tsv          |      3 | 17,513 |
| check.lst        |      1 | 17,513 |
| collect.tsv      |     20 | 17,514 |
| n50.tsv          |      4 | 17,514 |
| n50.pass.tsv     |      4 | 15,967 |
| collect.pass.tsv |     23 | 15,967 |
| pass.lst         |      1 | 15,966 |
| omit.lst         |      1 |    661 |
| rep.lst          |      1 |  1,135 |
| sp.lst           |      1 |  2,400 |

### Rsync to hpcc

```shell
rsync -avP \
    ~/data/Bacillus/ \
    wangq@202.119.37.251:data/Bacillus

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Bacillus/ \
    wangq@58.213.64.36:data/Bacillus

# rsync -avP wangq@202.119.37.251:data/Bacillus/ ~/data/Bacillus

```

For such huge collections, we can rsync files inside ASSEMBLY/ in parallel.

```shell
# Copy Directory Structure
rsync -avP \
    -f"+ */" -f"- *" \
    ~/data/Bacillus/ASSEMBLY/ \
    wangq@202.119.37.251:data/Bacillus/ASSEMBLY

# Transfer species directories in parallel
cat ~/data/Bacillus/ASSEMBLY/url.tsv |
    tsv-select -f 3 |
    rgr dedup stdin |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo -e "\n==> {}"
        rsync -avP \
            ~/data/Bacillus/ASSEMBLY/{}/
            wangq@202.119.37.251:data/Bacillus/ASSEMBLY/{} \
    '

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
#17960 lines, 91 fields

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
#  7016 summary/NR.lst
#  8202 summary/redundant.lst

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#1266

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Bacillus/

nwr template ~/Scripts/genomes/assembly/Bacillus.assembly.tsv \
    --mh \
    --parallel 8 \
    --in ASSEMBLY/rep.lst \
    --not-in ASSEMBLY/sp.lst \
    --not-in MinHash/abnormal.lst \
    --not-in summary/redundant.lst \
    --height 0.4

bash MinHash/dist.sh

```

### Condense branches in the minhash tree

```shell
mkdir -p ~/data/Bacillus/tree
cd ~/data/Bacillus/tree

nw_reroot ../MinHash/tree.nwk Desu_alkalia_AHT28_GCF_001730225_1 Desu_sti_MLFW_2_GCF_001742305_1 |
    nwr order stdin --nd --an \
    > minhash.reroot.newick

nwr pl-condense --map -r family -r genus \
    minhash.reroot.newick ../MinHash/species.tsv |
    nwr order stdin --nd --an \
    > minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

# svg
nwr topo --bl minhash.condensed.newick | # remove comments
    nw_display -s -b 'visibility:hidden' -w 1200 -v 20 - \
    > Bacillus.minhash.svg

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
    rgr md stdin --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/family.count.tsv |
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
| strain  | 13637 |
| species |  1209 |
| genus   |   167 |
| family  |     8 |
| order   |     2 |
| class   |     2 |

| family                 | #species | #strains |
|------------------------|---------:|---------:|
| Alicyclobacillaceae    |       39 |       60 |
| Bacillaceae            |      732 |    12068 |
| Paenibacillaceae       |      370 |     1404 |
| Sporolactobacillaceae  |       23 |       33 |
| Thermoactinomycetaceae |       41 |       68 |

| genus              | #species | #strains |
|--------------------|---------:|---------:|
| Aeribacillus       |        3 |       23 |
| Alicyclobacillus   |       22 |       37 |
| Aneurinibacillus   |        7 |       32 |
| Anoxybacillus      |       13 |       58 |
| Bacillus           |      130 |     9791 |
| Brevibacillus      |       30 |      208 |
| Cohnella           |       31 |       35 |
| Cytobacillus       |       19 |       90 |
| Fictibacillus      |       12 |       24 |
| Geobacillus        |       17 |      106 |
| Gracilibacillus    |       20 |       25 |
| Halobacillus       |       24 |       38 |
| Heyndrickxia       |       13 |      110 |
| Lysinibacillus     |       27 |      259 |
| Mesobacillus       |       10 |       20 |
| Metabacillus       |       20 |       39 |
| Neobacillus        |       24 |       47 |
| Niallia            |        6 |       52 |
| Oceanobacillus     |       30 |       65 |
| Paenibacillus      |      278 |     1100 |
| Parageobacillus    |        7 |       36 |
| Peribacillus       |       18 |      136 |
| Priestia           |       10 |      496 |
| Psychrobacillus    |       10 |       27 |
| Rossellomorea      |        7 |       29 |
| Schinkia           |        1 |       24 |
| Shouchella         |       11 |       68 |
| Sporolactobacillus |       13 |       22 |
| Thermoactinomyces  |        5 |       23 |
| Virgibacillus      |       25 |       75 |

| #family                | genus                 | species                             | count |
|------------------------|-----------------------|-------------------------------------|------:|
| Bacillaceae            | Aeribacillus          | Aeribacillus sp.                    |    11 |
|                        | Anoxybacillus         | Anoxybacillus sp.                   |    20 |
|                        | Bacillus              | Bacillus albus                      |    28 |
|                        |                       | Bacillus altitudinis                |   204 |
|                        |                       | Bacillus amyloliquefaciens          |   232 |
|                        |                       | Bacillus anthracis                  |   755 |
|                        |                       | Bacillus atrophaeus                 |   166 |
|                        |                       | Bacillus badius                     |    12 |
|                        |                       | Bacillus bombysepticus              |    11 |
|                        |                       | Bacillus cereus                     |  2167 |
|                        |                       | Bacillus cytotoxicus                |    39 |
|                        |                       | Bacillus glycinifermentans          |    15 |
|                        |                       | Bacillus halotolerans               |    69 |
|                        |                       | Bacillus haynesii                   |   125 |
|                        |                       | Bacillus inaquosorum                |   133 |
|                        |                       | Bacillus infantis                   |    10 |
|                        |                       | Bacillus licheniformis              |   347 |
|                        |                       | Bacillus mobilis                    |    41 |
|                        |                       | Bacillus mojavensis                 |    41 |
|                        |                       | Bacillus mycoides                   |   168 |
|                        |                       | Bacillus nitratireducens            |    14 |
|                        |                       | Bacillus pacificus                  |    83 |
|                        |                       | Bacillus paralicheniformis          |   212 |
|                        |                       | Bacillus paramycoides               |    16 |
|                        |                       | Bacillus paranthracis               |   265 |
|                        |                       | Bacillus pseudomycoides             |   111 |
|                        |                       | Bacillus pumilus                    |   252 |
|                        |                       | Bacillus safensis                   |   219 |
|                        |                       | Bacillus siamensis                  |    21 |
|                        |                       | Bacillus sonorensis                 |    36 |
|                        |                       | Bacillus spizizenii                 |   181 |
|                        |                       | Bacillus stercoris                  |    16 |
|                        |                       | Bacillus stratosphericus            |    11 |
|                        |                       | Bacillus subtilis                   |  1011 |
|                        |                       | Bacillus swezeyi                    |    14 |
|                        |                       | Bacillus thuringiensis              |  1146 |
|                        |                       | Bacillus toyonensis                 |   313 |
|                        |                       | Bacillus tropicus                   |    52 |
|                        |                       | Bacillus vallismortis               |    14 |
|                        |                       | Bacillus velezensis                 |   851 |
|                        |                       | Bacillus wiedmannii                 |   216 |
|                        | Caldibacillus         | Caldibacillus thermoamylovorans     |    11 |
|                        | Caldifermentibacillus | Caldifermentibacillus hisashii      |    16 |
|                        | Cytobacillus          | Cytobacillus firmus                 |    35 |
|                        |                       | Cytobacillus oceanisediminis        |    10 |
|                        | Geobacillus           | Geobacillus sp.                     |    37 |
|                        |                       | Geobacillus stearothermophilus      |    16 |
|                        |                       | Geobacillus thermodenitrificans     |    11 |
|                        |                       | Geobacillus thermoleovorans         |    11 |
|                        | Halalkalibacterium    | Halalkalibacterium halodurans       |    17 |
|                        | Heyndrickxia          | Heyndrickxia coagulans              |    56 |
|                        |                       | Heyndrickxia oleronia               |    13 |
|                        |                       | Heyndrickxia sp.                    |    14 |
|                        | Lysinibacillus        | Lysinibacillus capsici              |    29 |
|                        |                       | Lysinibacillus fusiformis           |    65 |
|                        |                       | Lysinibacillus sp.                  |    64 |
|                        |                       | Lysinibacillus sphaericus           |    44 |
|                        | Niallia               | Niallia circulans                   |    21 |
|                        |                       | Niallia sp.                         |    13 |
|                        |                       | Niallia taxi                        |    14 |
|                        | Oceanobacillus        | Oceanobacillus sp.                  |    11 |
|                        | Parageobacillus       | Parageobacillus thermoglucosidasius |    17 |
|                        | Peribacillus          | Peribacillus frigoritolerans        |    61 |
|                        |                       | Peribacillus simplex                |    29 |
|                        |                       | Peribacillus sp.                    |    14 |
|                        | Priestia              | Priestia aryabhattai                |    89 |
|                        |                       | Priestia endophytica                |    12 |
|                        |                       | Priestia filamentosa                |    10 |
|                        |                       | Priestia flexa                      |    23 |
|                        |                       | Priestia megaterium                 |   337 |
|                        |                       | Priestia sp.                        |    16 |
|                        | Psychrobacillus       | Psychrobacillus sp.                 |    14 |
|                        | Rossellomorea         | Rossellomorea marisflavi            |    14 |
|                        | Schinkia              | Schinkia azotoformans               |    24 |
|                        | Shouchella            | Shouchella clausii                  |    53 |
|                        | Terribacillus         | Terribacillus saccharophilus        |    10 |
|                        | Virgibacillus         | Virgibacillus halodenitrificans     |    13 |
|                        |                       | Virgibacillus pantothenticus        |    22 |
| Paenibacillaceae       | Aneurinibacillus      | Aneurinibacillus migulanus          |    11 |
|                        | Brevibacillus         | Brevibacillus agri                  |    24 |
|                        |                       | Brevibacillus borstelensis          |    24 |
|                        |                       | Brevibacillus brevis                |    15 |
|                        |                       | Brevibacillus formosus              |    10 |
|                        |                       | Brevibacillus laterosporus          |    49 |
|                        |                       | Brevibacillus parabrevis            |    11 |
|                        |                       | Brevibacillus sp.                   |    19 |
|                        | Paenibacillus         | Paenibacillus alvei                 |    17 |
|                        |                       | Paenibacillus amylolyticus          |    11 |
|                        |                       | Paenibacillus apiarius              |    11 |
|                        |                       | Paenibacillus larvae                |    58 |
|                        |                       | Paenibacillus lautus                |    13 |
|                        |                       | Paenibacillus odorifer              |    36 |
|                        |                       | Paenibacillus peoriae               |    12 |
|                        |                       | Paenibacillus polymyxa              |   110 |
|                        |                       | Paenibacillus sp.                   |   367 |
|                        |                       | Paenibacillus thiaminolyticus       |    26 |
|                        |                       | Paenibacillus xylanexedens          |    10 |
| Thermoactinomycetaceae | Thermoactinomyces     | Thermoactinomyces sp.               |    10 |

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
| strain  | 13297 |
| species |  1200 |
| genus   |   167 |
| family  |     8 |
| order   |     2 |
| class   |     2 |

| genus                 | #species | #strains |
|-----------------------|---------:|---------:|
| Aeribacillus          |        3 |       23 |
| Alicyclobacillus      |       22 |       37 |
| Alkalihalobacillus    |        8 |       12 |
| Aneurinibacillus      |        7 |       32 |
| Anoxybacillus         |       13 |       57 |
| Bacillus              |      128 |     9556 |
| Brevibacillus         |       30 |      205 |
| Caldibacillus         |        2 |       17 |
| Caldifermentibacillus |        1 |       16 |
| Cohnella              |       31 |       35 |
| Cytobacillus          |       19 |       89 |
| Domibacillus          |        9 |       15 |
| Fictibacillus         |       12 |       24 |
| Geobacillus           |       17 |      104 |
| Gracilibacillus       |       19 |       23 |
| Halalkalibacter       |       10 |       14 |
| Halobacillus          |       24 |       38 |
| Heyndrickxia          |       13 |      110 |
| Lederbergia           |        9 |       17 |
| Lentibacillus         |       12 |       13 |
| Lysinibacillus        |       27 |      211 |
| Mesobacillus          |       10 |       20 |
| Metabacillus          |       20 |       39 |
| Neobacillus           |       24 |       47 |
| Niallia               |        6 |       52 |
| Oceanobacillus        |       30 |       63 |
| Ornithinibacillus     |       11 |       12 |
| Paenibacillus         |      275 |     1077 |
| Parageobacillus       |        7 |       36 |
| Peribacillus          |       18 |      134 |
| Priestia              |       10 |      484 |
| Psychrobacillus       |       10 |       27 |
| Rossellomorea         |        7 |       29 |
| Salimicrobium         |        7 |       10 |
| Salipaludibacillus    |        5 |       10 |
| Schinkia              |        1 |       24 |
| Shouchella            |       11 |       67 |
| Siminovitchia         |        6 |       18 |
| Sporolactobacillus    |       13 |       22 |
| Sutcliffiella         |        6 |       10 |
| Terribacillus         |        4 |       15 |
| Thermoactinomyces     |        5 |       23 |
| Virgibacillus         |       25 |       74 |

## Collect proteins

```shell
cd ~/data/Bacillus/

nwr template ~/Scripts/genomes/assembly/Bacillus.assembly.tsv \
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

| #item      | count      |
|------------|------------|
| species    | 1,394      |
| strain_sum | 15,882     |
| total_sum  | 76,715,515 |
| dedup_sum  | 36,090,356 |
| rep_sum    | 12,688,227 |

## Phylogenetics with bac120

### Find corresponding proteins by `hmmsearch`

```shell
cd ~/data/Bacillus

# The HMM set
nwr kb bac120 -o HMM
cp HMM/bac120.lst HMM/marker.lst

E_VALUE=1e-20

# Find all genes
for marker in $(cat HMM/marker.lst); do
    echo >&2 "==> marker [${marker}]"

    mkdir -p Domain/${marker}

    cat Protein/species-f.tsv |
        tsv-join -f ASSEMBLY/rep.lst -k 1 |
        tsv-join -e -f ASSEMBLY/sp.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/hmm/${marker}.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\t%s\n), \$1, {1}, {2}; '
        " \
        > Domain/${marker}/replace.tsv

    echo >&2
done

```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Bacillus

cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat Domain/{}/replace.tsv |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
#1231    1242    1985.25

cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        cat Domain/{}/replace.tsv |
            wc -l
    ' |
    paste - - |
    tsv-filter --invert --ge 2:1000 --le 2:1400 |
    cut -f 1 \
    > Domain/marker.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process
cat Protein/species-f.tsv |
    tsv-join -f ASSEMBLY/rep.lst -k 1 |
    tsv-join -e -f ASSEMBLY/sp.lst -k 1 |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ ! -f Protein/"${SPECIES}"/pro.fa.gz ]]; then
        continue
    fi

    cat Protein/"${SPECIES}"/pro.fa.gz
done \
    > Domain/all.uniq.fa.gz

cat HMM/marker.lst |
    grep -v -Fx -f Domain/marker.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

        cat Domain/{}/replace.tsv \
            > Domain/{}/{}.replace.tsv

        faops some Domain/all.uniq.fa.gz <(
            cat Domain/{}/{}.replace.tsv |
                cut -f 1 |
                rgr dedup stdin
            ) stdout \
            > Domain/{}/{}.pro.fa
    '

# Align each markers with muscle
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

for marker in $(cat HMM/marker.lst); do
    echo >&2 "==> marker [${marker}]"
    if [ ! -s Domain/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sometimes `muscle` can not produce alignments
    if [ ! -s Domain/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # 1 name to many names
    cat Domain/${marker}/${marker}.replace.tsv |
        tsv-select -f 1-2 |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s Domain/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > Domain/${marker}/${marker}.replace.fa
done

# Concat marker genes
for marker in $(cat HMM/marker.lst); do
    if [ ! -s Domain/${marker}/${marker}.pro.fa ]; then
        continue
    fi
    if [ ! -s Domain/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # sequences in one line
    faops filter -l 0 Domain/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > Domain/bac120.aln.fas

cat Protein/species-f.tsv |
    tsv-join -f ASSEMBLY/rep.lst -k 1 |
    tsv-join -e -f ASSEMBLY/sp.lst -k 1 |
    cut -f 1 |
    fasops concat Domain/bac120.aln.fas stdin -o Domain/bac120.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Domain/bac120.aln.fa -out Domain/bac120.trim.fa -automated1

faops size Domain/bac120.*.fa |
    rgr dedup stdin -f 2 |
    cut -f 2
#34569
#19138

# To make it faster
FastTree -fastest -noml Domain/bac120.trim.fa > Domain/bac120.trim.newick

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 Domain/bac120.trim.newick |
    rsvg-convert -o tree/Bacillus.bac120.png

```

### Condense branches in the protein tree

```shell
cd ~/data/Bacillus/tree

nw_reroot  ../Domain/bac120.trim.newick Desu_alkalia_AHT28_GCF_001730225_1 Desu_sti_MLFW_2_GCF_001742305_1 |
    nwr order stdin --nd --an \
    > bac120.reroot.newick

nwr pl-condense --map -r order -r family -r genus \
    bac120.reroot.newick ../Count/species.tsv |
    nwr order stdin --nd --an \
    -o bac120.condensed.newick

mv condensed.tsv bac120.condense.tsv

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 bac120.condensed.newick |
    rsvg-convert -o Bacillus.bac120.png

```
