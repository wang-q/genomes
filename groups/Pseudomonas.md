# *Pseudomonas* HGT

## Strain info

* [Pseudomonas](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=286)
* [Acinetobacter](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=469)

According to a recent [paper](https://doi.org/10.1128/mSystems.00543-20), there are some order-level
changes in Gammaproteobacteria. We include both old and new orders.

* Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
* New ones: Moraxellales, Kangiellales, and Pseudomonadales

### List all ranks

```shell
mkdir -p ~/data/Pseudomonas
cd ~/data/Pseudomonas

nwr member Pseudomonas |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

nwr member Acinetobacter |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

nwr member Pseudomonas Acinetobacter -r "species group" -r "species subgroup" |
    tsv-select -f 1-3 |
    keep-header -- tsv-sort -k3,3 -k2,2 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    mlr --itsv --omd cat

```

| rank             | count |
|------------------|-------|
| genus            | 1     |
| species          | 425   |
| strain           | 724   |
| subspecies       | 13    |
| no rank          | 121   |
| species group    | 5     |
| species subgroup | 4     |
| isolate          | 1     |

| rank             | count |
|------------------|-------|
| genus            | 1     |
| species          | 117   |
| species group    | 2     |
| species subgroup | 3     |
| strain           | 1110  |
| no rank          | 2     |
| subspecies       | 1     |
| isolate          | 2     |

| #tax_id | sci_name                                 | rank             |
|---------|------------------------------------------|------------------|
| 2839056 | A. Taxon 24                              | species group    |
| 909768  | A. calcoaceticus/baumannii complex       | species group    |
| 136841  | P. aeruginosa group                      | species group    |
| 136842  | P. chlororaphis group                    | species group    |
| 136843  | P. fluorescens group                     | species group    |
| 136845  | P. putida group                          | species group    |
| 136849  | P. syringae group                        | species group    |
| 2839060 | A. Taxon 24C                             | species subgroup |
| 2839057 | A. Taxon 24D                             | species subgroup |
| 2839061 | A. Taxon 24E                             | species subgroup |
| 627141  | P. nitroreducens/multiresinivorans group | species subgroup |
| 1232139 | P. oleovorans/pseudoalcaligenes group    | species subgroup |
| 251695  | P. syringae group genomosp. 1            | species subgroup |
| 251698  | P. syringae group genomosp. 2            | species subgroup |


### Species with assemblies

Check also the order Pseudomonadales.

* Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
* New ones: Moraxellales, Kangiellales, and Pseudomonadales

```shell
mkdir -p ~/data/Pseudomonas/summary
cd ~/data/Pseudomonas/summary

SPECIES=$(
    nwr member \
        Cellvibrionales Oceanospirillales Alteromonadales \
        Moraxellales Kangiellales Pseudomonadales \
        -r species |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        grep -v " sp." |
        grep -v -E "\bbacterium\b" |
        grep -v -E "\bsymbiont\b" |
        sed '1d' |
        cut -f 1 |
        sort |
        uniq
)

for S in $SPECIES; do
    RS=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    CHR=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
                AND assembly_level IN ('Complete Genome', 'Chromosome')
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    if [[ ${RS} -gt 0 ]]; then
        echo -e "$S\t$RS\t$CHR"
    fi
done |
    nwr append stdin |
    tsv-select -f 1,4,2-3 |
    tsv-sort -k3,3nr -k4,4nr -k2,2 |
    (echo -e 'species_id\tspecies\tRS\tCHR' && cat) \
    > species.count.tsv

cat species.count.tsv |
    tsv-filter -H --ge CHR:5 |
    tsv-filter -H --invert --str-in-fld species:Pseudomonas --lt RS:30 |
    tsv-filter -H --invert --str-in-fld species:Acinetobacter --lt RS:30 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    mlr --itsv --omd cat

```

| species_id | species                         | RS   | CHR |
|------------|---------------------------------|------|-----|
| 287        | P. aeruginosa                   | 7458 | 640 |
| 470        | A. baumannii                    | 7191 | 528 |
| 33069      | P. viridiflava                  | 1541 | 9   |
| 317        | P. syringae                     | 612  | 50  |
| 48296      | A. pittii                       | 385  | 35  |
| 294        | P. fluorescens                  | 262  | 40  |
| 303        | P. putida                       | 223  | 67  |
| 480        | Moraxella catarrhalis           | 209  | 16  |
| 106654     | A. nosocomialis                 | 173  | 11  |
| 316        | Stutzerimonas stutzeri          | 141  | 28  |
| 29438      | P. savastanoi                   | 120  | 5   |
| 38313      | Shewanella algae                | 112  | 24  |
| 587753     | P. chlororaphis                 | 107  | 61  |
| 756892     | A. indicus                      | 102  | 21  |
| 47877      | P. amygdali                     | 90   | 10  |
| 380021     | P. protegens                    | 83   | 27  |
| 47880      | P. fulva                        | 75   | 5   |
| 40215      | A. junii                        | 71   | 11  |
| 108980     | A. ursingii                     | 68   | 6   |
| 106648     | A. bereziniae                   | 63   | 5   |
| 1530123    | A. seifertii                    | 60   | 25  |
| 40216      | A. radioresistens               | 60   | 5   |
| 29430      | A. haemolyticus                 | 57   | 14  |
| 62322      | Shewanella baltica              | 48   | 12  |
| 40214      | A. johnsonii                    | 47   | 19  |
| 28090      | A. lwoffii                      | 47   | 11  |
| 76759      | P. monteilii                    | 47   | 9   |
| 296        | P. fragi                        | 40   | 7   |
| 476        | Moraxella bovis                 | 39   | 37  |
| 202956     | A. towneri                      | 32   | 8   |
| 930166     | P. brassicacearum               | 31   | 9   |
| 78327      | P. mosselii                     | 30   | 6   |
| 332186     | Shewanella xiamenensis          | 30   | 5   |
| 28108      | Alteromonas macleodii           | 26   | 16  |
| 43657      | Pseudoalteromonas luteoviolacea | 25   | 5   |
| 34062      | Moraxella osloensis             | 23   | 10  |
| 43662      | Pseudoalteromonas piscicida     | 19   | 6   |
| 314275     | Alteromonas mediterranea        | 17   | 16  |
| 24         | Shewanella putrefaciens         | 15   | 9   |
| 2968969    | Stutzerimonas frequens          | 15   | 5   |
| 44935      | Halomonas venusta               | 10   | 7   |
| 386891     | Moraxella bovoculi              | 9    | 7   |
| 1697053    | Thiopseudomonas alkaliphila     | 7    | 7   |

## Download all assemblies

### Create assembly.tsv

```shell
cd ~/data/Pseudomonas/summary

cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    head -n 1 \
    > Pseudomonas.assembly.tsv

# Pseudomonas aeruginosa PAO1 is in the reference list
cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    grep -F -f <(
        cat ~/Scripts/genomes/assembly/Bacteria.reference.tsv |
            tsv-select -H -f assembly_accession
        ) \
    >> Pseudomonas.assembly.tsv

# Pseudomonadales and close relatives
# Species with 2 or more complete genomes were retained.
cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    sed '1d' |
    tsv-join -d 3 \
        -k 2 -f <(
            cat species.count.tsv |
                tsv-filter -H --ge CHR:2 |
                sed '1d'
        )\
    >> Pseudomonas.assembly.tsv

# Also includes representative strains of Gammaproteobacteria.
# families not in our orders
FAMILY=$(
    nwr member Gammaproteobacteria -r family |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        sed '1d' |
        cut -f 1 |
        nwr append stdin -r order |
        tsv-filter --str-ne 2:"Cellvibrionales" --str-ne 2:"Oceanospirillales" --str-ne 2:"Alteromonadales" |
        tsv-filter --str-ne 2:"Moraxellales" --str-ne 2:"Kangiellales" --str-ne 2:"Pseudomonadales" |
        tsv-filter --str-ne 2:"NA" |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        assembly_accession
    FROM ar
    WHERE 1=1
        AND family_id IN ($FAMILY)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND refseq_category IN ('representative genome')
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    grep -v -i "symbiont " |
    tsv-filter --str-not-in-fld 1:"[" \
    > tmp.tsv

cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    grep -F -f tmp.tsv \
    >> Pseudomonas.assembly.tsv

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Pseudomonas.assembly.tsv
# cp Pseudomonas.assembly.tsv ~/Scripts/genomes/assembly

# Comment out unneeded strains

# Cleaning
rm tmp*.*sv

```

### rsync and check

```shell
cd ~/data/Pseudomonas

cat ~/Scripts/genomes/assembly/Pseudomonas.assembly.tsv |
    tsv-filter -v --str-in-fld 2:http

nwr assembly ~/Scripts/genomes/assembly/Pseudomonas.assembly.tsv \
    -o ASSEMBLY

# Copy downloaded files
cat ASSEMBLY/url.tsv |
    tsv-join -f ASSEMBLY/check.list -k 1 -e |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -d ASSEMBLY/{1} ]; then
            if [ -d ~/data/Bacteria/ASSEMBLY/{1} ]; then
                echo >&2 "==> {1}"
                cp -r ~/data/Bacteria/ASSEMBLY/{1} ASSEMBLY/{1}
            fi
        fi
    '

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
proxychains4 bash ASSEMBLY/rsync.sh

# md5
# rm ASSEMBLY/check.list
bash ASSEMBLY/check.sh

# collect
bash ASSEMBLY/collect.sh

```
