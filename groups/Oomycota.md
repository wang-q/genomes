# Oomycota

<!-- toc -->

## Taxon info

* [Oomycota](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4762)

A nice introduction article about Oomycota:

https://doi.org/10.1016/j.cub.2018.05.062

### List all ranks

```shell
nwr member Oomycota |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat

```

| rank       | count |
|------------|-------|
| phylum     | 1     |
| no rank    | 55    |
| family     | 20    |
| genus      | 85    |
| species    | 1216  |
| subspecies | 3     |
| order      | 11    |
| strain     | 39    |
| isolate    | 1     |
| varietas   | 15    |
| forma      | 4     |

### Species with assemblies

```shell
mkdir -p ~/data/Oomycota/summary
cd ~/data/Oomycota/summary

# should have a valid name of genus
nwr member Oomycota -r genus |
    grep -v -i "Candidatus " |
    grep -v -i "candidate " |
    grep -v " sp." |
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#85 genus.list

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
# 194 GB1.tsv

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
#GB1     415

```

## Download all assemblies

### Create assembly.tsv

If a refseq assembly is available, the corresponding genbank one is not downloaded

```shell
cd ~/data/Oomycota/summary

cat ~/Scripts/genomes/assembly/Fungi.reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,assembly_level,assembly_accession \
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
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level,
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
    cut -f 6 \
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
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 6 -e \
    >> raw.tsv

cat raw.tsv |
    tsv-uniq |
    datamash check
#417 lines, 6 fields

# Create abbr.
cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
    tsv-select -f 1-5 |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[5]}++; # abbr_name
        $seen{$F[5]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $F[5], $F[3], $F[1], $F[4];
        ' |
    tsv-filter --or --str-in-fld 2:ftp --str-in-fld 2:http |
    keep-header -- sort -k3,3 -k1,1 \
    > Oomycota.assembly.tsv

datamash check < Oomycota.assembly.tsv
#417 lines, 4 fields

# find potential duplicate strains or assemblies
cat Oomycota.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Oomycota.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Oomycota.assembly.tsv
# cp Oomycota.assembly.tsv ~/Scripts/genomes/assembly

# Comment out unneeded strains

# Cleaning
rm raw*.*sv

```

### Rsync and check

```shell
cd ~/data/Oomycota

nwr assembly ~/Scripts/genomes/assembly/Oomycota.assembly.tsv \
    -o ASSEMBLY

# Remove dirs not in the list
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    tr "/" "\t" |
    cut -f 2 |
    tsv-join --exclude -k 1 -f ASSEMBLY/url.tsv -d 1 |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        echo Remove {}
        rm -fr ASSEMBLY/{}
    '

# Run
bash ASSEMBLY/rsync.sh

# Check md5; create check.lst
bash ASSEMBLY/check.sh

# N50 C S; create n50.tsv and n50.pass.csv
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
#35704092        49127533.5      8042    24984.5 5018    14664

# Temporary files, possibly caused by an interrupted rsync process
find ASSEMBLY/ -type f -name ".*" > ASSEMBLY/temp.list

cat ASSEMBLY/temp.list |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        if [[ -f {} ]]; then
            echo Remove {}
            rm {}
        fi
    '

# Collect; create collect.csv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.csv summary/

cat ASSEMBLY/counts.tsv |
    mlr --itsv --omd cat

```

| #item            | count |
|------------------|-------|
| url.tsv          | 416   |
| check.lst        | 416   |
| collect.csv      | 417   |
| n50.tsv          | 416   |
| n50.pass.csv     | 94    |
| omit.lst         | 305   |
| collect.pass.csv | 94    |

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
