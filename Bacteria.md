# Bacteria

All genomes of *Bacteria*, species by species.

Download all genomes and analyze representative strains.

<!-- toc -->

- [Taxon info](#taxon-info)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
    * [Model organisms](#model-organisms)
- [Download all assemblies](#download-all-assemblies)
    * [Create assembly.tsv](#create-assemblytsv)
    * [Count before download](#count-before-download)
    * [Download and check](#download-and-check)
    * [Remove unnecessary files](#remove-unnecessary-files)
    * [Rsync to hpcc](#rsync-to-hpcc)
- [BioSample](#biosample)
- [Early divergence of Bacteria](#early-divergence-of-bacteria)
    * [ReRoot](#reroot)
- [MinHash](#minhash)
    * [Condense branches in the minhash tree](#condense-branches-in-the-minhash-tree)
- [Count valid species and strains](#count-valid-species-and-strains)
    * [For *genomic alignments* and *protein
      families*](#for-genomic-alignments-and-protein-families)
- [Collect proteins](#collect-proteins)
- [Phylogenetics with bac120](#phylogenetics-with-bac120)
    * [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
    * [Condense branches in the protein tree](#condense-branches-in-the-protein-tree)
- [InterProScan on all proteins of representative and typical strains](#interproscan-on-all-proteins-of-representative-and-typical-strains)

<!-- tocstop -->

## Taxon info

* [Bacteria](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2)

### List all ranks

```shell
nwr member Bacteria |
    grep -v " sp." |
    grep -v " x " |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

```

| rank             |  count |
|------------------|-------:|
| superkingdom     |      1 |
| phylum           |    181 |
| family           |    961 |
| no rank          |   7824 |
| genus            |   5409 |
| class            |    189 |
| order            |    407 |
| species          | 114752 |
| strain           |  42007 |
| subclass         |      6 |
| subspecies       |    714 |
| suborder         |      7 |
| clade            |    135 |
| isolate          |    456 |
| varietas         |     24 |
| forma            |      4 |
| species group    |     98 |
| species subgroup |     31 |
| biotype          |      7 |
| serotype         |    255 |
| serogroup        |    152 |
| subgenus         |      1 |
| tribe            |      2 |
| pathogroup       |      5 |
| subfamily        |      1 |

### Species with assemblies

* RefSeq
    * '>= 200 genomes'
        * RS1 - assembly_level: 'Complete Genome', 'Chromosome'
        * RS2 - genome_rep: 'Full'
    * '>= 2 genomes'
        * RS3 - assembly_level: 'Complete Genome', 'Chromosome'
        * RS4 - genome_rep: 'Full'
* Genbank
    * '>= 200 genomes'
        * GB1 - assembly_level: 'Complete Genome', 'Chromosome'
        * GB2 - genome_rep: 'Full'
    * '>= 2 genomes'
        * GB3 - assembly_level: 'Complete Genome', 'Chromosome'
        * GB4 - genome_rep: 'Full'

```shell
mkdir -p ~/data/Bacteria/summary
cd ~/data/Bacteria/summary

# should have a valid name of genus
# NCBI division Bacteria includes Bacteria and Archaea
nwr member Bacteria -r genus |
    grep -v " sp." |
    grep -v " x " |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#5409 genus.list

QUERY="
    SELECT
        species_id,
        species,
        COUNT(*) AS count
    FROM ar
    WHERE 1=1
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% x %'
        REPLACE_1
        REPLACE_2
    GROUP BY species_id
    REPLACE_3
"

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo -e "${QUERY}" |
        sed "s/REPLACE_1/AND genus_id = ${RANK_ID}/" |
        sed "s/REPLACE_2/AND assembly_level IN ('Complete Genome', 'Chromosome')/" |
        sed "s/REPLACE_3/HAVING count >= 200/" |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS1.tsv

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo -e "${QUERY}" |
        sed "s/REPLACE_1/AND genus_id = ${RANK_ID}/" |
        sed "s/REPLACE_2/AND genome_rep IN ('Full')/" |
        sed "s/REPLACE_3/HAVING count >= 200/" |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS2.tsv

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo -e "${QUERY}" |
        sed "s/REPLACE_1/AND genus_id = ${RANK_ID}/" |
        sed "s/REPLACE_2/AND assembly_level IN ('Complete Genome', 'Chromosome')/" |
        sed "s/REPLACE_3/HAVING count >= 2/" |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS3.tsv

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo -e "${QUERY}" |
        sed "s/REPLACE_1/AND genus_id = ${RANK_ID}/" |
        sed "s/REPLACE_2/AND genome_rep IN ('Full')/" |
        sed "s/REPLACE_3/HAVING count >= 2/" |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS4.tsv

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo -e "${QUERY}" |
        sed "s/REPLACE_1/AND genus_id = ${RANK_ID}/" |
        sed "s/REPLACE_2/AND assembly_level IN ('Complete Genome', 'Chromosome')/" |
        sed "s/REPLACE_3/HAVING count >= 200/" |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB1.tsv

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo -e "${QUERY}" |
        sed "s/REPLACE_1/AND genus_id = ${RANK_ID}/" |
        sed "s/REPLACE_2/AND genome_rep IN ('Full')/" |
        sed "s/REPLACE_3/HAVING count >= 200/" |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB2.tsv

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo -e "${QUERY}" |
        sed "s/REPLACE_1/AND genus_id = ${RANK_ID}/" |
        sed "s/REPLACE_2/AND assembly_level IN ('Complete Genome', 'Chromosome')/" |
        sed "s/REPLACE_3/HAVING count >= 2/" |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB3.tsv

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo -e "${QUERY}" |
        sed "s/REPLACE_1/AND genus_id = ${RANK_ID}/" |
        sed "s/REPLACE_2/AND genome_rep IN ('Full')/" |
        sed "s/REPLACE_3/HAVING count >= 2/" |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB4.tsv

wc -l RS*.tsv GB*.tsv
#    21 RS1.tsv
#   171 RS2.tsv
#  2135 RS3.tsv
#  8628 RS4.tsv
#    28 GB1.tsv
#   250 GB2.tsv
#  2272 GB3.tsv
# 10880 GB4.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     17429
#RS2     236784
#RS3     39114
#RS4     317947
#GB1     26689
#GB2     1686322
#GB3     50411
#GB4     1789020

cat RS2.tsv |
    tsv-join -f RS1.tsv -k 1 -e |
    tsv-summarize --sum 3
#76451

cat RS3.tsv |
    tsv-join -f RS1.tsv -k 1 -e |
    tsv-join -f RS2.tsv -k 1 -e |
    tsv-summarize --sum 3
#12682

cat RS4.tsv |
    tsv-join -f RS1.tsv -k 1 -e |
    tsv-join -f RS2.tsv -k 1 -e |
    tsv-join -f RS3.tsv -k 1 -e |
    tsv-summarize --sum 3
#30760

cat RS4.tsv |
    tsv-join -f RS1.tsv -k 1 -e |
    tsv-join -f RS2.tsv -k 1 -e |
    tsv-summarize --sum 3
#81163

```

* Some species are divided into several separate species

```shell
cd ~/data/Bacteria/summary

nwr member Pseudomonas -r species |
    grep -v " sp." |
    grep -E "syringae|genomosp"
#251701  Pseudomonas syringae group genomosp. 3  species Bacteria
#251699  Pseudomonas syringae group genomosp. 7  species Bacteria
#317659  Pseudomonas syringae pv. coryli species Bacteria
#317     Pseudomonas syringae    species Bacteria

echo "
    SELECT
        *
    FROM ar
    WHERE 1=1
        AND species_id IN (251701,251699,317659)
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    GROUP BY species_id
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite

```

### Model organisms

```shell
cd ~/data/Bacteria/summary

GENUS=$(
    cat genus.list.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
.headers ON

    SELECT
        refseq_category,
        COUNT(*) AS count
    FROM ar
    WHERE 1=1
    GROUP BY refseq_category
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
#refseq_category count
#na      352870
#reference genome        22415

echo "
.headers ON

    SELECT
        ar.*
    FROM (
        SELECT
            species_id,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND species_id != 0
            AND genus_id != 0
        GROUP BY species_id
        HAVING count >= 1000
        ) AS subquery
    JOIN ar ON ar.species_id = subquery.species_id
    WHERE 1=1
        AND ar.refseq_category IN ('reference genome')
        AND ar.tax_id != ar.species_id
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > reference.tsv

cat reference.tsv |
    sed '1s/^/#/' |
    nwr append stdin -r phylum -r class |
    tsv-select -H -f 1-3,phylum,class |
    keep-header -- sort -t$'\t' -k4,4 -k5,5 -k2,2 |
    sed 's/strain=//' |
    parallel --col-sep "\t" -j 1 '
        if [[ "{4}" == "Proteobacteria" || "{4}" == "Pseudomonadota" ]]; then
            printf "%s\t%s\t%s\t%s\n" {1} {2} {3} {5}
        else
            printf "%s\t%s\t%s\t%s\n" {1} {2} {3} {4}
        fi
    ' |
    mlr --itsv --omd cat

cp reference.tsv ~/Scripts/genomes/assembly/Bacteria.reference.tsv

```

| #tax_id | organism_name                                                    | infraspecific_name                | phylum              |
|---------|------------------------------------------------------------------|-----------------------------------|---------------------|
| 83332   | Mycobacterium tuberculosis H37Rv                                 | H37Rv                             | Actinomycetota      |
| 561007  | Mycobacteroides abscessus ATCC 19977                             | ATCC 19977                        | Actinomycetota      |
| 1169293 | Enterococcus faecalis EnGen0336                                  | T5                                | Bacillota           |
| 169963  | Listeria monocytogenes EGD-e                                     | EGD-e                             | Bacillota           |
| 93061   | Staphylococcus aureus subsp. aureus NCTC 8325                    | NCTC 8325                         | Bacillota           |
| 568814  | Streptococcus suis BM407                                         | BM407                             | Bacillota           |
| 192222  | Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819      | NCTC 11168                        | Campylobacterota    |
| 386585  | Escherichia coli O157:H7 str. Sakai                              | Sakai substr. RIMD 0509952        | Gammaproteobacteria |
| 511145  | Escherichia coli str. K-12 substr. MG1655                        | K-12 substr. MG1655               | Gammaproteobacteria |
| 1125630 | Klebsiella pneumoniae subsp. pneumoniae HS11286                  | HS11286                           | Gammaproteobacteria |
| 2590157 | Klebsiella variicola subsp. variicola                            | F2R9T                             | Gammaproteobacteria |
| 529507  | Proteus mirabilis HI4320                                         | HI4320                            | Gammaproteobacteria |
| 208964  | Pseudomonas aeruginosa PAO1                                      | PAO1                              | Gammaproteobacteria |
| 99287   | Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 | LT2                               | Gammaproteobacteria |
| 223926  | Vibrio parahaemolyticus RIMD 2210633                             | RIMD 2210633 substr. RIMD 2210633 | Gammaproteobacteria |

## Download all assemblies

### Create assembly.tsv

* RefSeq
    * '>= 2 genomes'
        * RS4 - genome_rep: 'Full'

```shell
cd ~/data/Bacteria/summary

cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession \
    > raw.tsv

# RS4
SPECIES=$(
    cat RS4.tsv |
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

datamash check < raw.tsv
#317981 lines, 7 fields

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
    > Bacteria.assembly.tsv

datamash check < Bacteria.assembly.tsv
#317862 lines, 5 fields

# find potential duplicate strains or assemblies
cat Bacteria.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Bacteria.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

# Edit .assembly.tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Bacteria.assembly.tsv
#
# Save the file to another directory to prevent accidentally changing it
# cp Bacteria.assembly.tsv ~/Scripts/genomes/assembly

# Cleaning
rm raw*.*sv

```

### Count before download

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class

```shell
cd ~/data/Bacteria

nwr template ~/Scripts/genomes/assembly/Bacteria.assembly.tsv \
    --count \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# .lst and .count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    keep-header -- tsv-sort -k1,1 |
    tsv-filter -H --ge 3:1000 |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```

| item    |  count |
|---------|-------:|
| strain  | 317859 |
| species |   8626 |
| genus   |   2022 |
| family  |    483 |
| order   |    203 |
| class   |     88 |

| genus               | #species | #strains |
|---------------------|---------:|---------:|
| Acinetobacter       |       71 |    12036 |
| Aeromonas           |       25 |     1584 |
| Bacillus            |       79 |     7990 |
| Bacteroides         |       47 |     3468 |
| Bifidobacterium     |       56 |     2669 |
| Bordetella          |       18 |     1353 |
| Brucella            |       28 |     1290 |
| Burkholderia        |       38 |     5033 |
| Campylobacter       |       41 |     5678 |
| Citrobacter         |       18 |     2014 |
| Clostridioides      |        1 |     2978 |
| Clostridium         |       66 |     2492 |
| Corynebacterium     |      126 |     1820 |
| Enterobacter        |       27 |     6117 |
| Enterococcus        |       49 |     9454 |
| Escherichia         |        6 |    40072 |
| Francisella         |       11 |     1076 |
| Haemophilus         |       11 |     1095 |
| Helicobacter        |       35 |     3154 |
| Klebsiella          |       14 |    25445 |
| Lactiplantibacillus |        8 |     1217 |
| Lactobacillus       |       39 |     1972 |
| Lactococcus         |       16 |     1000 |
| Legionella          |       50 |     1222 |
| Listeria            |       19 |     5921 |
| Mycobacterium       |       77 |     8042 |
| Mycobacteroides     |        6 |     2115 |
| Neisseria           |       35 |     3455 |
| Phocaeicola         |       12 |     1007 |
| Proteus             |       10 |     1461 |
| Pseudomonas         |      260 |    15476 |
| Salmonella          |        2 |    14276 |
| Serratia            |       22 |     2060 |
| Shigella            |        4 |     2609 |
| Staphylococcus      |       56 |    22998 |
| Stenotrophomonas    |       20 |     1003 |
| Streptococcus       |       84 |    19487 |
| Streptomyces        |      377 |     2578 |
| Vibrio              |      110 |     7349 |
| Xanthomonas         |       38 |     2920 |
| Yersinia            |       26 |     1629 |

### Download and check

* A standardized archaeal taxonomy for the Genome Taxonomy Database.
  Nat Microbiol 6, 946–959 (2021).
    * Only genomes comprising ≤200 contigs with an N50 of ≥20 kb and with CheckM completeness and
      contamination estimates of ≥95% and ≤5%, respectively, were considered.

```shell
cd ~/data/Bacteria

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Bacteria.assembly.tsv \
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

find ASSEMBLY/ -name "*_genomic.fna.gz" |
    grep -v "_from_" |
    wc -l
#163583
find ASSEMBLY -type f -name "*_protein.faa.gz" -size -4096c


# N50 C S; create n50.tsv and n50.pass.tsv
bash ASSEMBLY/n50.sh 20000 500 100000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50,S" --max "C"
#N50_min S_min   C_max
#0       0       1973

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "S:0.1,0.5" --quantile "N50:0.1,0.5"  --quantile "C:0.5,0.9"
#S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
#1986616 4098739 40286.6 250927  47      210

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

| #item            | fields |   lines |
|------------------|-------:|--------:|
| url.tsv          |      3 | 163,585 |
| check.lst        |      1 | 163,585 |
| collect.tsv      |     20 | 163,586 |
| n50.tsv          |      4 | 163,586 |
| n50.pass.tsv     |      4 | 157,846 |
| collect.pass.tsv |     23 | 157,846 |
| pass.lst         |      1 | 157,845 |
| omit.lst         |      1 |       2 |
| rep.lst          |      1 |   7,943 |
| sp.lst           |      1 |      44 |

### Remove unnecessary files

```shell
cd ~/data/Bacteria

find ASSEMBLY -type d -name "*assembly_structure" | xargs rm -fr
find ASSEMBLY -type f -name "annotation_hashes.txt" | xargs rm -f
find ASSEMBLY -type f -name "*_feature_table.txt.gz" | xargs rm -f
find ASSEMBLY -type f -name "*_genomic_gaps.txt.gz" | xargs rm -f
find ASSEMBLY -type f -name "*_genomic.gtf.gz" | xargs rm -f
find ASSEMBLY -type f -name "*_protein.gpff.gz" | xargs rm -f
find ASSEMBLY -type f -name "*_translated_cds.faa.gz" | xargs rm -f
find ASSEMBLY -type f -name "*_wgsmaster.gbff.gz" | xargs rm -f

```

### Rsync to hpcc

```shell
rsync -avP \
    ~/data/Bacteria/ \
    wangq@202.119.37.251:data/Bacteria

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Bacteria/ \
    wangq@58.213.64.36:data/Bacteria

# back
rsync -avP \
    wangq@202.119.37.251:data/Bacteria/ \
    ~/data/Bacteria

rsync -avP \
    -e 'ssh -p 8804' \
    wangq@58.213.64.36:data/Bacteria/ \
    ~/data/Bacteria

```

For such huge collections, we can rsync files inside ASSEMBLY/ in parallel.

```shell
# Copy Directory Structure
rsync -avP \
    -f"+ */" -f"- *" \
    ~/data/Bacteria/ASSEMBLY/ \
    wangq@202.119.37.251:data/Bacteria/ASSEMBLY

# Transfer species directories in parallel
cat ~/data/Bacteria/ASSEMBLY/url.tsv |
    tsv-select -f 3 |
    tsv-uniq |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo -e "\n==> {}"
        rsync -avP \
            ~/data/Bacteria/ASSEMBLY/{}/
            wangq@202.119.37.251:data/Bacteria/ASSEMBLY/{} \
    '

# Copy Directory Structure
rsync -avP \
    -e 'ssh -p 8804' \
    -f"+ */" -f"- *" \
    ~/data/Bacteria/ASSEMBLY/ \
    wangq@58.213.64.36:data/Bacteria/ASSEMBLY

# Transfer species directories in parallel
cat ~/data/Bacteria/ASSEMBLY/url.tsv |
    tsv-select -f 3 |
    tsv-uniq |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo -e "\n==> {}"
        rsync -avP \
            -e "ssh -p 8804" \
            ~/data/Bacteria/ASSEMBLY/{}/
            wangq@58.213.64.36:data/Bacteria/ASSEMBLY/{} \
    '

```

## BioSample

```shell
cd ~/data/Bacteria

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Bacteria.assembly.tsv \
    --bs

# Run this script twice and it will re-download the failed files
bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 500

datamash check < BioSample/biosample.tsv
# 71588 lines, 64 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

# Save spaces
cd ~/data/Bacteria
cat BioSample/sample.tsv |
    tsv-select -f 3 |
    tsv-uniq |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        cd BioSample
        if [[ ! -d {} ]]; then
            exit
        fi
        if [[ ! -s {}.tar ]]; then
            echo -e "==> {}"
            tar -cvf {}.tar --remove-files {}
        fi
    '

```

## Early divergence of Bacteria

![molbiolevolmsn247f02_ht.jpeg](images%2Fmolbiolevolmsn247f02_ht.jpeg)

### ReRoot

Thermotogota and Aquificota are basal taxa of bacteria.

* For the latter steps, use the following two as the outgroups
    * Thermot_petr_RKU_1_GCF_000016785_1
    * Hyd_thermophilus_TK_6_GCF_000164905_1

```shell
cd ~/data/Bacteria

cat summary/collect.pass.tsv |
    tsv-select -f 1,3 |
    sed '1d' |
    grep -v -Fw -f ASSEMBLY/omit.lst |
    nwr append stdin -c 2 -r species -r phylum |
    tsv-filter --or \
        --str-eq "4:Thermotogota" \
        --str-eq "4:Aquificota" |
    tsv-select -f 1,3,4 |
    tsv-sort -k3,3 -k1,1 |
    tsv-summarize -g 3,2 --count
#Aquificota      Hydrogenobacter thermophilus    2
#Thermotogota    Fervidobacterium pennivorans    3
#Thermotogota    Pseudothermotoga hypogea        2
#Thermotogota    Thermosipho melanesiensis       2
#Thermotogota    Thermotoga maritima     6
#Thermotogota    Thermotoga petrophila   2

cat summary/collect.pass.tsv |
    tsv-filter -H --not-blank RefSeq_category |
    tsv-filter -H --or \
        --str-in-fld "2:Hydrogenobacter" \
        --str-in-fld "2:Fervidobacterium" \
        --str-in-fld "2:Pseudothermotoga" \
        --str-in-fld "2:Thermosipho" \
        --str-in-fld "2:Thermotoga" |
    grep -v -Fw -f ASSEMBLY/omit.lst |
    tsv-select -H -f "name,Assembly_level,Assembly_method,Genome_coverage,Sequencing_technology"

```

## MinHash

```shell
cd ~/data/Bacteria/

nwr template ~/Scripts/genomes/assembly/Bacteria.assembly.tsv \
    --mh \
    --parallel 16 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
    --ani-nr 0.005

# Compute assembly sketches
bash MinHash/compute.sh

#find MinHash -name "*.msh" -empty | wc -l

# Distances within species
bash MinHash/species.sh

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#1477

# Non-redundant strains within species
bash MinHash/nr.sh

find MinHash -name "mash.dist.tsv" -size +0 | wc -l
#1743

find MinHash -name "redundant.lst" -size +0 | wc -l
#990

find MinHash -name "redundant.lst" -empty | wc -l
#753

find MinHash -name "NR.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/NR.lst
wc -l summary/NR.lst
#22805

# All representative should be in NR
cat ASSEMBLY/rep.lst |
    grep -v -F -f summary/NR.lst

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Bacteria/

nwr template ~/Scripts/genomes/assembly/Bacteria.assembly.tsv \
    --mh \
    --parallel 16 \
    --in ASSEMBLY/rep.lst \
    --not-in ASSEMBLY/omit.lst \
    --not-in MinHash/abnormal.lst \
    --height 0.4

bash MinHash/dist.sh

```

### Condense branches in the minhash tree

```shell
mkdir -p ~/data/Bacteria/tree
cd ~/data/Bacteria/tree

nwr reroot ../MinHash/tree.nwk -n Thermot_petr_RKU_1_GCF_000016785_1 -n Hyd_thermophilus_TK_6_GCF_000164905_1 |
    nwr order stdin --nd --an \
    > minhash.reroot.newick

nwr pl-condense -r class -r order -r family -r genus -r species \
    minhash.reroot.newick ../Count/species.tsv --map \
    -o minhash.condensed.newick

mv condensed.tsv minhash.condense.tsv

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 minhash.condensed.newick |
    rsvg-convert -o Bacteria.minhash.png

```

## Count valid species and strains

### For *genomic alignments* and *protein families*

```shell
cd ~/data/Bacteria/

nwr template ~/Scripts/genomes/assembly/Bacteria.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst \
    --rank order --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv
bash Count/strains.sh

# .lst and .count.tsv
bash Count/rank.sh

cat Count/order.count.tsv |
    tsv-filter -H --ge "3:500" |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:500" |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

# Can accept N_COUNT
bash Count/lineage.sh 300

cat Count/lineage.count.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| order             | #species | #strains |
|-------------------|---------:|---------:|
| Aeromonadales     |       11 |      688 |
| Bacillales        |      138 |     7693 |
| Bacteroidales     |       44 |     2296 |
| Bifidobacteriales |       14 |     1656 |
| Burkholderiales   |       95 |     3527 |
| Campylobacterales |       40 |     2614 |
| Enterobacterales  |      148 |    13596 |
| Eubacteriales     |       60 |     1861 |
| Flavobacteriales  |       46 |      847 |
| Hyphomicrobiales  |       74 |     1450 |
| Lactobacillales   |      148 |     8118 |
| Moraxellales      |       28 |     1517 |
| Mycobacteriales   |       94 |     3604 |
| Pasteurellales    |       25 |     1644 |
| Pseudomonadales   |       91 |     3410 |
| Spirochaetales    |       23 |      516 |
| Thiotrichales     |       12 |      955 |
| Vibrionales       |       46 |     1366 |
| Xanthomonadales   |       32 |     1948 |

| genus              | #species | #strains |
|--------------------|---------:|---------:|
| Acinetobacter      |       22 |     1265 |
| Aeromonas          |       11 |      688 |
| Bacillus           |       32 |     2846 |
| Bacteroides        |       15 |     1391 |
| Bifidobacterium    |       12 |     1577 |
| Bordetella         |        9 |      807 |
| Burkholderia       |       26 |     2198 |
| Campylobacter      |       26 |     2032 |
| Citrobacter        |       11 |      502 |
| Clostridium        |       19 |     1315 |
| Corynebacterium    |       34 |      692 |
| Enterobacter       |       10 |     1545 |
| Enterococcus       |       12 |     1291 |
| Escherichia        |        4 |     3316 |
| Francisella        |        9 |      871 |
| Haemophilus        |        6 |      876 |
| Klebsiella         |       10 |     3567 |
| Lacticaseibacillus |        5 |      566 |
| Lactobacillus      |       15 |      596 |
| Listeria           |        6 |      740 |
| Mycobacterium      |       16 |      822 |
| Mycobacteroides    |        3 |     1900 |
| Pseudomonas        |       81 |     3262 |
| Salmonella         |        2 |     1557 |
| Staphylococcus     |       35 |     3446 |
| Stenotrophomonas   |        5 |      597 |
| Streptococcus      |       37 |     3193 |
| Vibrio             |       40 |     1332 |
| Xanthomonas        |       18 |     1130 |

| #family              | genus            | species                         | count |
|----------------------|------------------|---------------------------------|------:|
| Alcaligenaceae       | Bordetella       | Bordetella pertussis            |   593 |
| Bacillaceae          | Bacillus         | Bacillus subtilis               |   302 |
|                      |                  | Bacillus velezensis             |   306 |
| Bacteroidaceae       | Bacteroides      | Bacteroides fragilis            |   355 |
|                      |                  | Bacteroides uniformis           |   323 |
| Bifidobacteriaceae   | Bifidobacterium  | Bifidobacterium longum          |   559 |
| Burkholderiaceae     | Burkholderia     | Burkholderia cenocepacia        |   418 |
|                      |                  | Burkholderia multivorans        |   454 |
| Campylobacteraceae   | Campylobacter    | Campylobacter coli              |  1204 |
| Clostridiaceae       | Clostridium      | Clostridium botulinum           |   399 |
|                      |                  | Clostridium perfringens         |   523 |
| Corynebacteriaceae   | Corynebacterium  | Corynebacterium diphtheriae     |   381 |
| Enterobacteriaceae   | Cronobacter      | Cronobacter sakazakii           |   415 |
|                      | Escherichia      | Escherichia coli                |  2967 |
|                      | Klebsiella       | Klebsiella aerogenes            |   355 |
|                      |                  | Klebsiella michiganensis        |   356 |
|                      |                  | Klebsiella pneumoniae           |  1718 |
|                      |                  | Klebsiella variicola            |   619 |
|                      | Salmonella       | Salmonella enterica             |  1545 |
| Enterococcaceae      | Enterococcus     | Enterococcus faecium            |   319 |
| Francisellaceae      | Francisella      | Francisella tularensis          |   829 |
| Helicobacteraceae    | Helicobacter     | Helicobacter pylori             |   406 |
| Listeriaceae         | Listeria         | Listeria monocytogenes          |   360 |
| Moraxellaceae        | Acinetobacter    | Acinetobacter baumannii         |   601 |
|                      |                  | Acinetobacter pittii            |   347 |
| Mycobacteriaceae     | Mycobacterium    | Mycobacterium tuberculosis      |   571 |
|                      | Mycobacteroides  | Mycobacteroides abscessus       |  1887 |
| Pasteurellaceae      | Haemophilus      | Haemophilus influenzae          |   823 |
| Propionibacteriaceae | Cutibacterium    | Cutibacterium acnes             |   417 |
| Pseudomonadaceae     | Pseudomonas      | Pseudomonas aeruginosa          |   706 |
|                      |                  | Pseudomonas syringae            |   351 |
|                      |                  | Pseudomonas viridiflava         |  1367 |
| Rhizobiaceae         | Rhizobium        | Rhizobium leguminosarum         |   440 |
| Staphylococcaceae    | Staphylococcus   | Staphylococcus aureus           |  1159 |
|                      |                  | Staphylococcus haemolyticus     |   450 |
|                      |                  | Staphylococcus pseudintermedius |   401 |
| Streptococcaceae     | Lactococcus      | Lactococcus lactis              |   307 |
|                      | Streptococcus    | Streptococcus equi              |   495 |
| Xanthomonadaceae     | Stenotrophomonas | Stenotrophomonas maltophilia    |   587 |

## InterProScan on all proteins of representative and typical strains

To be filled by other projects

```shell
mkdir -p ~/data/Bacteria/STRAINS

```
