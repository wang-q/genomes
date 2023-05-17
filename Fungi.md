# Fungi

Download all genomes and analyze representative strains.

<!-- toc -->

## Taxon info

* [Fungi](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4751)

### List all ranks

```shell
nwr member Fungi |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat

```

| rank          | count |
|---------------|------:|
| kingdom       |     1 |
| no rank       |  4577 |
| species       | 64236 |
| subkingdom    |     1 |
| class         |    64 |
| order         |   240 |
| family        |   901 |
| genus         |  7410 |
| phylum        |    10 |
| subphylum     |    14 |
| strain        |  2237 |
| varietas      |  1053 |
| subspecies    |   161 |
| forma         |   207 |
| isolate       |    29 |
| subclass      |    19 |
| clade         |    26 |
| suborder      |    24 |
| subfamily     |    17 |
| subgenus      |    10 |
| section       |    37 |
| species group |     1 |
| tribe         |     3 |
| superfamily   |     1 |
| morph         |     2 |

### Species with assemblies

In the vast majority of fungal species, only one genome was selected for refseq.

* 'RefSeq'
    * RS1 - assembly_level: 'Complete Genome', 'Chromosome'
    * RS2 - genome_rep: 'Full'
* 'Genbank'
    * '>= 20 genomes'
        * GB1 - With strain ID; assembly_level: 'Complete Genome', 'Chromosome'
        * GB2 - assembly_level: 'Complete Genome', 'Chromosome'
        * GB3 - NOT 'contig'; genome_rep: 'Full'
    * '>= 2 genomes'
        * GB4 - assembly_level: 'Complete Genome', 'Chromosome'

```shell
mkdir -p ~/data/Fungi/summary
cd ~/data/Fungi/summary

# should have a valid name of genus
nwr member Fungi -r genus |
    grep -v -i "Candidatus " |
    grep -v -i "candidate " |
    grep -v " sp." |
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#7410 genus.list

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
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
            AND assembly_level IN ('Complete Genome', 'Chromosome')
            AND genome_rep IN ('Full') -- fully representative
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS1.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
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
            AND genome_rep IN ('Full') -- fully representative
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS2.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
    echo "
        SELECT
            species_id,
            species,
            COUNT(DISTINCT tax_id) AS count -- with strain ID
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND assembly_level IN ('Complete Genome', 'Chromosome')
            AND genome_rep IN ('Full') -- fully representative
        GROUP BY species_id
        HAVING count >= 20
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB1.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
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
            AND assembly_level IN ('Complete Genome', 'Chromosome')
            AND genome_rep IN ('Full') -- fully representative
        GROUP BY species_id
        HAVING count >= 20
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB2.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
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
            AND assembly_level IN ('Complete Genome', 'Chromosome', 'Scaffold')
            AND genome_rep IN ('Full') -- fully representative
        GROUP BY species_id
        HAVING count >= 20
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB3.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
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
            AND assembly_level IN ('Complete Genome', 'Chromosome')
        GROUP BY species_id
        HAVING count >= 2
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB4.tsv

wc -l RS*.tsv GB*.tsv
#   81 RS1.tsv
#  487 RS2.tsv
#    1 GB1.tsv
#    4 GB2.tsv
#   51 GB3.tsv
#   82 GB4.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e ${C}${N}.tsv ]; then
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#83
#492
#109
#750
#4473
#1082

```

* The names of some genera are abnormal

```shell
cd ~/data/Fungi/summary

cat RS*.tsv GB*.tsv |
    cut -f 1,2 |
    tsv-uniq |
    grep "\[" |
    nwr append stdin -r genus
#498019	[Candida] auris	Candida/Metschnikowiaceae
#1231522	[Candida] duobushaemulonis	Candida/Metschnikowiaceae
#45357	[Candida] haemuloni	Candida/Metschnikowiaceae
#418784	[Candida] pseudohaemulonii	Candida/Metschnikowiaceae
#561895	[Candida] subhashii	Spathaspora
#5477	[Candida] boidinii	Ogataea
#45354	[Candida] intermedia	Candida/Metschnikowiaceae

```

### Model organisms

There is only one model genome in Fungi, Saccharomyces cerevisiae S288C

```shell
cd ~/data/Fungi/summary

GENUS=$(
    cat genus.list.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > reference.tsv

```

## Download all assemblies

### Create assembly.tsv

Five levels:

* 'RefSeq'
    * RS2 - genome_rep: 'Full'
* 'Genbank'
    * '>= 20 genomes'
        * GB1 - With strain ID; assembly_level: 'Complete Genome', 'Chromosome'
        * GB2 - assembly_level: 'Complete Genome', 'Chromosome'
        * GB3 - NOT 'contig'; genome_rep: 'Full'
    * '>= 2 genomes'
        * GB4 - assembly_level: 'Complete Genome', 'Chromosome'

If a refseq assembly is available, the corresponding genbank one is not downloaded

```shell
cd ~/data/Fungi/summary

cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,assembly_level,assembly_accession \
    > raw.tsv

# RS2
SPECIES=$(
    cat RS2.tsv |
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
        AND genome_rep IN ('Full') -- fully representative
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

# Avoid refseq in genbank
cat raw.tsv |
    cut -f 6 \
    > rs.acc.tsv

# GB1
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
        AND tax_id NOT IN ($SPECIES) -- With strain ID
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 6 -e \
    >> raw.tsv

# GB2
SPECIES=$(
    cat GB2.tsv |
        tsv-join -f GB1.tsv -k 1 -e |
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
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 6 -e \
    >> raw.tsv

# GB3
SPECIES=$(
    cat GB3.tsv |
        tsv-join -f GB1.tsv -k 1 -e |
        tsv-join -f GB2.tsv -k 1 -e |
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
        AND assembly_level IN ('Complete Genome', 'Chromosome', 'Scaffold')
        AND genome_rep IN ('Full') -- fully representative
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 6 -e \
    >> raw.tsv

# GB4
SPECIES=$(
    cat GB4.tsv |
        tsv-join -f GB1.tsv -k 1 -e |
        tsv-join -f GB2.tsv -k 1 -e |
        tsv-join -f GB3.tsv -k 1 -e |
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
        AND assembly_level IN ('Complete Genome', 'Chromosome', 'Scaffold')
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 6 -e \
    >> raw.tsv

datamash check < raw.tsv
#3602 lines, 6 fields

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
    keep-header -- sort -k3,3 -k1,1 \
    > Fungi.assembly.tsv

datamash check < Fungi.assembly.tsv
#3601 lines, 4 fields

# find potential duplicate strains or assemblies
cat Fungi.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Fungi.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

cat Fungi.assembly.tsv |
    tsv-filter --or --str-in-fld 1:genomosp --str-in-fld 1:genomovar

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Fungi.assembly.tsv
# cp Fungi.assembly.tsv ~/Scripts/genomes/assembly

# Comment out unneeded strains

# Cleaning
rm raw*.*sv

```

### rsync and check

```shell
cd ~/data/Fungi

cat ~/Scripts/genomes/assembly/Fungi.assembly.tsv |
    tsv-filter -v --str-in-fld 2:http

nwr assembly ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
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

# Check md5
# rm ASSEMBLY/check.list
bash ASSEMBLY/check.sh

# Collect
bash ASSEMBLY/collect.sh

# Temporary files, possibly caused by an interrupted rsync process
find ASSEMBLY/ -type f -name ".*" > ASSEMBLY/temp.list

cat ASSEMBLY/temp.list |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        if [[ -f {} ]]; then
            echo Remove {}
            rm {}
        fi
    '

```

### Check N50 of assemblies

```shell
cd ~/data/Fungi

cat ASSEMBLY/collect.csv | # head |
    sed '1d' |
    tsv-select -d "," -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        1>&2 echo "==> {}"

        find ASSEMBLY/{} -type f -name "*_genomic.fna.gz" |
            grep -v "_from_" | # exclude CDS and rna
            xargs cat |
            faops n50 -C -S stdin |
            (echo -e "name\t{}" && cat) |
            datamash transpose
    ' |
    tsv-uniq | # keep the first header
    tee ASSEMBLY/n50.tsv

cat ASSEMBLY/n50.tsv |
    tsv-filter \
        -H --or \
        --le 4:100 \
        --ge 2:100000 |
    tsv-filter -H --ge 3:1000000 |
    tr "\t" "," \
    > ASSEMBLY/n50.pass.csv

wc -l ASSEMBLY/n50* ASSEMBLY/collect.csv
#   35159 ASSEMBLY/n50.pass.csv
#   38660 ASSEMBLY/n50.tsv
#   38660 ASSEMBLY/collect.csv

# Omit strains without protein annotations
for STRAIN in $(cat ASSEMBLY/n50.pass.csv | cut -d, -f 1); do
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_protein.faa.gz" > /dev/null; then
        echo ${STRAIN}
    fi
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_cds_from_genomic.fna.gz" > /dev/null; then
        echo ${STRAIN}
    fi
done |
    tsv-uniq
# All OK

tsv-join \
    ASSEMBLY/collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv \
    > summary/collect.pass.csv

wc -l summary/collect.pass.csv
#35159 summary/collect.pass.csv

```

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Fungi/ \
    wangq@202.119.37.251:data/Fungi

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Fungi/ \
    wangq@58.213.64.36:data/Fungi

# rsync -avP wangq@202.119.37.251:data/Fungi/ ~/data/Fungi

# rsync -avP -e "ssh -T -c chacha20-poly1305@openssh.com -o Compression=no -x" \
#   wangq@202.119.37.251:data/Fungi/ASSEMBLY/ ~/data/Fungi/ASSEMBLY

```
