# Plants

### List all ranks

```shell
for phylum in Viridiplantae Rhodophyta Glaucophyta; do
    nwr member ${phylum} |
        grep -v " sp." |
        grep -v " x " |
        tsv-summarize -H -g rank --count |
        rgr md stdin --num
    echo
done

```

| rank          |  count |
|---------------|-------:|
| kingdom       |      1 |
| phylum        |      3 |
| class         |     38 |
| order         |    189 |
| family        |    924 |
| genus         |  17030 |
| no rank       |   5218 |
| species       | 186344 |
| strain        |     40 |
| subphylum     |      1 |
| subspecies    |   8152 |
| varietas      |   8590 |
| subclass      |     28 |
| clade         |    152 |
| forma         |    398 |
| subfamily     |    276 |
| suborder      |     24 |
| subgenus      |    163 |
| section       |    494 |
| tribe         |    856 |
| subtribe      |    469 |
| subsection    |     41 |
| morph         |      5 |
| series        |      4 |
| isolate       |      3 |
| genotype      |     14 |
| species group |      5 |
| superorder    |      2 |

| rank       | count |
|------------|------:|
| phylum     |     1 |
| genus      |   809 |
| no rank    |   451 |
| class      |     5 |
| order      |    43 |
| family     |   119 |
| species    |  3719 |
| subclass   |     5 |
| subspecies |     7 |
| varietas   |    30 |
| forma      |    13 |
| clade      |     2 |
| subfamily  |    14 |
| tribe      |     8 |
| strain     |     3 |

| rank    | count |
|---------|------:|
| class   |     1 |
| no rank |     2 |
| genus   |     4 |
| species |    15 |
| order   |     2 |
| family  |     3 |

### Species with assemblies

```shell
mkdir -p ~/data/Plants/summary
cd ~/data/Plants/summary

# should have a valid name of genus
nwr member Viridiplantae Rhodophyta Glaucophyta -r genus |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#17873 genus.list

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
#  182 RS1.tsv
# 2050 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     182
#GB1     4753

```

## Download all assemblies

### Create assembly.tsv

If a refseq assembly is available, the corresponding genbank one is not downloaded

```shell
cd ~/data/Plants/summary

# RS
SPECIES=$(
    cat RS1.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
.header ON
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% x %'
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
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% x %'
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

cat raw.tsv |
    tsv-uniq |
    datamash check
#4764 lines, 7 fields

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
    > Plants.assembly.tsv

datamash check < Plants.assembly.tsv
#4466 lines, 5 fields

# find potential duplicate strains or assemblies
cat Plants.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Plants.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

# Edit .assembly.tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Plants.assembly.tsv
#
# Save the file to another directory to prevent accidentally changing it
# cp Plants.assembly.tsv ~/Scripts/genomes/assembly

# Cleaning
rm raw*.*sv

```

### Count before download

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class

```shell
cd ~/data/Plants

nwr template ~/Scripts/genomes/assembly/Plants.assembly.tsv \
    --count \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --num

# genus.lst and genus.count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    keep-header -- tsv-sort -k1,1 |
    tsv-filter -H --ge 3:20 |
    rgr md stdin --num

```

| item    | count |
|---------|------:|
| strain  |  4464 |
| species |  2027 |
| genus   |   988 |
| family  |   292 |
| order   |   113 |
| class   |    32 |

| genus        | #species | #strains |
|--------------|---------:|---------:|
| Arabidopsis  |        5 |      158 |
| Arachis      |        6 |       27 |
| Bathycoccus  |        2 |       23 |
| Beta         |        5 |      127 |
| Brassica     |        9 |       71 |
| Capsicum     |        5 |       31 |
| Carex        |       14 |       22 |
| Chlorella    |        5 |       22 |
| Citrullus    |        7 |       35 |
| Citrus       |       15 |       48 |
| Cucumis      |        5 |       53 |
| Erythroxylum |       62 |       63 |
| Eucalyptus   |       34 |       41 |
| Fraxinus     |       23 |       31 |
| Glycine      |        4 |       29 |
| Gossypium    |       29 |       61 |
| Hordeum      |        4 |      189 |
| Juglans      |        9 |       22 |
| Linum        |        5 |       27 |
| Malus        |       11 |       61 |
| Micromonas   |        2 |       31 |
| Nicotiana    |       13 |       20 |
| Oryza        |       14 |      160 |
| Panicum      |        3 |       40 |
| Populus      |       11 |       22 |
| Prunus       |       17 |       47 |
| Quercus      |       21 |       34 |
| Raphanus     |        2 |       20 |
| Salix        |       11 |       27 |
| Solanum      |       54 |      198 |
| Sorghum      |        1 |       21 |
| Trifolium    |        9 |       25 |
| Triticum     |        6 |       45 |
| Vigna        |        9 |       31 |
| Vitis        |        7 |       39 |
| Zea          |        2 |      116 |

### Download and check

```shell
cd ~/data/Plants

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Plants.assembly.tsv \
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
#5302    449376  848

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "S:0.1,0.5" --quantile "N50:0.1,0.5"  --quantile "C:0.5,0.9"
#S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
#1220388 2242415 11520   239729.5        27      217.5

# After the above steps are completed, run the following commands.

# Collect; create collect.tsv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    rgr md stdin --right 2-3

```
