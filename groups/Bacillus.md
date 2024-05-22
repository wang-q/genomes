# *Bacillus*

<!-- toc -->

- [Strain info](#strain-info)
    * [Symlink](#symlink)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
- [Download all assemblies](#download-all-assemblies)
    * [Create assembly.tsv](#create-assemblytsv)
    * [Count before download](#count-before-download)
    * [Download and check](#download-and-check)

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

### Symlink

```shell
mkdir -p ~/data/Bacillus
cd ~/data/Bacillus

rm -fr ASSEMBLY
rm -fr STRAINS

ln -s ../Bacteria/ASSEMBLY ASSEMBLY
ln -s ../Bacteria/STRAINS STRAINS

```

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
| genus            |   182 |
| species          | 43978 |
| no rank          |   245 |
| strain           |   751 |
| subspecies       |    39 |
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
#182 genus.list.tsv

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
#  2482 RS1.tsv
#  2796 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     10128
#GB1     12794

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
        AND genus IN ('Bacillus')
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
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

cat raw.tsv |
    tsv-uniq |
    datamash check
#12796 lines, 7 fields

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
#12788 lines, 5 fields

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
| strain  | 12704 |
| species |  2721 |
| genus   |   159 |
| family  |     5 |
| order   |     1 |
| class   |     1 |

| genus              | #species | #strains |
|--------------------|---------:|---------:|
| Alicyclobacillus   |       30 |       56 |
| Alkalihalobacillus |       28 |       78 |
| Anoxybacillus      |       34 |       90 |
| Bacillus           |     1045 |     9016 |
| Brevibacillus      |       60 |      170 |
| Cytobacillus       |       18 |       87 |
| Geobacillus        |       55 |      117 |
| Halobacillus       |       33 |       52 |
| Heyndrickxia       |        3 |       53 |
| Lysinibacillus     |       79 |      222 |
| Neobacillus        |       33 |       58 |
| Oceanobacillus     |       40 |       62 |
| Paenibacillus      |      561 |     1029 |
| Peribacillus       |       25 |      108 |
| Priestia           |       18 |      352 |
| Virgibacillus      |       42 |       87 |
| Weizmannia         |        7 |       69 |

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
#1373 lines, 41 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```
