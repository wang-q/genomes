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
    1. assembly_level: 'Complete Genome', 'Chromosome'
    2. genome_rep: 'Full'
* 'Genbank'
    * '>= 20 genomes'
        1. With strain ID; assembly_level: 'Complete Genome', 'Chromosome'
        2. With strain ID; genome_rep: 'Full'
        3. assembly_level: 'Complete Genome', 'Chromosome'; genome_rep: 'Full'
        4. genome_rep: 'Full'
    * '>= 2 genomes'
        5. assembly_level: 'Complete Genome', 'Chromosome'

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
            COUNT(DISTINCT tax_id) AS count -- with strain ID
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
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
            AND assembly_level IN ('Complete Genome', 'Chromosome')
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
            AND genome_rep IN ('Full') -- fully representative
        GROUP BY species_id
        HAVING count >= 20
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB4.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND assembly_level IN ('Complete Genome', 'Chromosome')
            AND genome_rep IN ('Full') -- fully representative
        GROUP BY species_id
        HAVING count >= 2
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB5.tsv

wc -l RS*.tsv GB*.tsv
#      80 RS1.tsv
#     482 RS2.tsv
#       1 GB1.tsv
#       5 GB2.tsv
#       4 GB3.tsv
#      82 GB4.tsv
#      83 GB5.tsv

for C in RS GB; do
    for N in $(seq 1 1 5); do
        if [ -e ${C}${N}.tsv ]; then
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#82
#487
#109
#294
#750
#6731
#1084

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
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite \
    > reference.tsv

```
