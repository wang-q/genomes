# Bacteria

All genomes of *Bacteria*, species by species.

Download all genomes and analyze representative strains.

<!-- toc -->

- [Strain info](#strain-info)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
    * [Model organisms](#model-organisms)
- [Download all assemblies](#download-all-assemblies)
    * [Create assembly.tsv](#create-assemblytsv)
    * [rsync and check](#rsync-and-check)
    * [Check N50 of assemblies](#check-n50-of-assemblies)
    * [Rsync to hpcc](#rsync-to-hpcc)
- [BioSample](#biosample)
- [Count species and strains](#count-species-and-strains)
    * [Order](#order)
    * [Genus](#genus)
- [MinHash](#minhash)
    * [Compute MinHash](#compute-minhash)
    * [Raw phylo-tree of representative assemblies](#raw-phylo-tree-of-representative-assemblies)
    * [Tweak the mash tree](#tweak-the-mash-tree)
- [Non-redundant strains within species](#non-redundant-strains-within-species)
    * [Genus](#genus-1)
- [Collect proteins](#collect-proteins)
    * [`all.pro.fa`](#allprofa)
    * [`all.replace.fa`](#allreplacefa)
    * [`all.info.tsv`](#allinfotsv)
- [Phylogenetics with bac120](#phylogenetics-with-bac120)
    * [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
    * [Tweak the concat tree](#tweak-the-concat-tree)
- [InterProScan on all proteins of representative and typical strains](#interproscan-on-all-proteins-of-representative-and-typical-strains)
- [Early divergence of Bacteria](#early-divergence-of-bacteria)
    * [Terrabacteria group](#terrabacteria-group)
    * [Proteobacteria](#proteobacteria)

<!-- tocstop -->

## Strain info

* [Bacteria](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2)

### List all ranks

```shell
nwr member Bacteria |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

```

| rank             |  count |
|------------------|-------:|
| superkingdom     |      1 |
| phylum           |    179 |
| no rank          |   6704 |
| species          | 111042 |
| class            |    150 |
| order            |    321 |
| family           |    797 |
| genus            |   5053 |
| strain           |  41351 |
| subclass         |      6 |
| subspecies       |    718 |
| suborder         |      7 |
| clade            |    134 |
| isolate          |    455 |
| varietas         |     24 |
| forma            |      4 |
| species group    |     99 |
| species subgroup |     31 |
| biotype          |      7 |
| serotype         |    258 |
| serogroup        |    145 |
| subgenus         |      1 |
| tribe            |      2 |
| pathogroup       |      5 |
| subfamily        |      1 |

### Species with assemblies

* RefSeq
    * '>= 100 genomes'
        * RS1 - assembly_level: 'Complete Genome', 'Chromosome'
        * RS2 - genome_rep: 'Full'
    * '>= 2 genomes'
        * RS3 - assembly_level: 'Complete Genome', 'Chromosome'
        * RS4 - genome_rep: 'Full'
* Genbank
    * '>= 100 genomes'
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
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#4470 genus.list

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
        sed "s/REPLACE_3/HAVING count >= 100/" |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS1.tsv

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo -e "${QUERY}" |
        sed "s/REPLACE_1/AND genus_id = ${RANK_ID}/" |
        sed "s/REPLACE_2/AND genome_rep IN ('Full')/" |
        sed "s/REPLACE_3/HAVING count >= 100/" |
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
        sed "s/REPLACE_3/HAVING count >= 100/" |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB1.tsv

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo -e "${QUERY}" |
        sed "s/REPLACE_1/AND genus_id = ${RANK_ID}/" |
        sed "s/REPLACE_2/AND genome_rep IN ('Full')/" |
        sed "s/REPLACE_3/HAVING count >= 100/" |
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
#    48 RS1.tsv
#   230 RS2.tsv
#  1741 RS3.tsv
#  6799 RS4.tsv
#    53 GB1.tsv
#   289 GB2.tsv
#  1822 GB3.tsv
#  7611 GB4.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     16893
#RS2     204296
#RS3     30510
#RS4     252538
#GB1     22591
#GB2     1244156
#GB3     37300
#GB4     1299256

cat RS2.tsv |
    tsv-join -f RS1.tsv -k 1 -e |
    tsv-summarize --sum 3
#46739

cat RS3.tsv |
    tsv-join -f RS1.tsv -k 1 -e |
    tsv-join -f RS2.tsv -k 1 -e |
    tsv-summarize --sum 3
#8084

cat RS4.tsv |
    tsv-join -f RS1.tsv -k 1 -e |
    tsv-join -f RS2.tsv -k 1 -e |
    tsv-join -f RS3.tsv -k 1 -e |
    tsv-summarize --sum 3
#20130

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
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > reference.tsv

cat reference.tsv |
    sed '1s/^/#/' |
    nwr append stdin -r phylum -r class |
    tsv-select -H -f 1,2,phylum,class |
    parallel --col-sep "\t" -j 1 '
        if [[ "{3}" == "Proteobacteria" || "{3}" == "Pseudomonadota" ]]; then
            printf "%s\t%s\t%s\n" {1} {2} {4}
        else
            printf "%s\t%s\t%s\n" {1} {2} {3}
        fi
    ' |
    mlr --itsv --omd cat

cp reference.tsv ~/Scripts/genomes/assembly/Bacteria.reference.tsv

```

| #tax_id | organism_name                                                    | phylum              |
|---------|------------------------------------------------------------------|---------------------|
| 565050  | Caulobacter vibrioides NA1000                                    | Alphaproteobacteria |
| 192222  | Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819      | Campylobacterota    |
| 208964  | Pseudomonas aeruginosa PAO1                                      | Gammaproteobacteria |
| 871585  | Acinetobacter pittii PHEA-2                                      | Gammaproteobacteria |
| 511145  | Escherichia coli str. K-12 substr. MG1655                        | Gammaproteobacteria |
| 386585  | Escherichia coli O157:H7 str. Sakai                              | Gammaproteobacteria |
| 1125630 | Klebsiella pneumoniae subsp. pneumoniae HS11286                  | Gammaproteobacteria |
| 99287   | Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 | Gammaproteobacteria |
| 198214  | Shigella flexneri 2a str. 301                                    | Gammaproteobacteria |
| 227377  | Coxiella burnetii RSA 493                                        | Gammaproteobacteria |
| 272561  | Chlamydia trachomatis D/UW-3/CX                                  | Chlamydiota         |
| 93061   | Staphylococcus aureus subsp. aureus NCTC 8325                    | Bacillota           |
| 224308  | Bacillus subtilis subsp. subtilis str. 168                       | Bacillota           |
| 169963  | Listeria monocytogenes EGD-e                                     | Bacillota           |
| 83332   | Mycobacterium tuberculosis H37Rv                                 | Actinomycetota      |

## Download all assemblies

### Create assembly.tsv

Three levels:

* RefSeq
    * '>= 100 genomes'
        * RS1 - assembly_level: 'Complete Genome', 'Chromosome'
        * RS2 - genome_rep: 'Full'
    * '>= 2 genomes'
        * RS3 - assembly_level: 'Complete Genome', 'Chromosome'

```shell
cd ~/data/Bacteria/summary

cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession \
    > raw.tsv

# RS1
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
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

# RS2
SPECIES=$(
    cat RS2.tsv |
        tsv-join -f RS1.tsv -k 1 -e |
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

# RS3
SPECIES=$(
    cat RS3.tsv |
        tsv-join -f RS1.tsv -k 1 -e |
        tsv-join -f RS2.tsv -k 1 -e |
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
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

datamash check < raw.tsv
#71732 lines, 7 fields

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
#71686 lines, 5 fields

# find potential duplicate strains or assemblies
cat Bacteria.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Bacteria.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

cat Bacteria.assembly.tsv |
    tsv-filter --or --str-in-fld 1:genomosp --str-in-fld 1:genomovar

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
    tsv-filter -H --ge 3:500 |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```

| item    | count |
|---------|------:|
| strain  | 71685 |
| species |  1745 |
| genus   |   525 |
| family  |   204 |
| order   |    95 |
| class   |    44 |

| genus               | #species | #strains |
|---------------------|---------:|---------:|
| Acinetobacter       |       22 |     1340 |
| Aeromonas           |       11 |      828 |
| Bacillus            |       32 |     3129 |
| Bacteroides         |       15 |     1583 |
| Bifidobacterium     |       12 |     1599 |
| Bordetella          |        9 |      807 |
| Burkholderia        |       26 |     2462 |
| Campylobacter       |       27 |     2111 |
| Citrobacter         |       11 |      532 |
| Clostridium         |       20 |     1463 |
| Corynebacterium     |       34 |      697 |
| Enterobacter        |       10 |     1692 |
| Enterococcus        |       12 |     1327 |
| Escherichia         |        4 |     3358 |
| Francisella         |        9 |      959 |
| Haemophilus         |        6 |      898 |
| Klebsiella          |       10 |     3808 |
| Lacticaseibacillus  |        5 |      623 |
| Lactobacillus       |       15 |      951 |
| Lactococcus         |        8 |      508 |
| Limosilactobacillus |        5 |      523 |
| Listeria            |        6 |      748 |
| Mycobacterium       |       16 |      938 |
| Mycobacteroides     |        3 |     1920 |
| Pseudomonas         |       89 |     4079 |
| Rhizobium           |        9 |      538 |
| Salmonella          |        2 |     1564 |
| Shigella            |        4 |     1800 |
| Staphylococcus      |       35 |     3593 |
| Stenotrophomonas    |        5 |      758 |
| Streptococcus       |       38 |     3372 |
| Vibrio              |       41 |     1548 |
| Xanthomonas         |       18 |     1346 |

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
#

# N50 C S; create n50.tsv and n50.pass.tsv
bash ASSEMBLY/n50.sh 100000 200 1000000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50,S" --max "C"
#N50_min S_min   C_max
#32179   2187595 2016

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "S:0.1,0.5" --quantile "N50:0.1,0.5"  --quantile "C:0.5,0.9"
#S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
#11671954.4      32505046        26704.8 368791.5        349.5   4717.6

# After the above steps are completed, run the following commands.

# Collect; create collect.tsv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

wc -l ASSEMBLY/n50* ASSEMBLY/collect.csv
#   35159 ASSEMBLY/n50.pass.csv
#   38660 ASSEMBLY/n50.tsv
#   38660 ASSEMBLY/collect.csv

```

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Bacteria/ \
    wangq@202.119.37.251:data/Bacteria

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Bacteria/ \
    wangq@58.213.64.36:data/Bacteria

# back
rsync -avP \
    -e 'ssh -p 8804' \
    wangq@58.213.64.36:data/Bacteria/ \
    ~/data/Bacteria

```

For such huge collections, we can rsync files inside ASSEMBLY/ in parallel.

```shell
# Copy Directory Structure
rsync -avP \
    -e 'ssh -p 8804' \
    -f"+ */" -f"- *" \
    wangq@58.213.64.36:data/Bacteria/ASSEMBLY/ \
    ~/data/Bacteria/ASSEMBLY

# Transfer species directories in parallel
cat ~/data/Bacteria/ASSEMBLY/url.tsv |
    tsv-select -f 3 |
    tsv-uniq |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo -e "\n==> {}"
        rsync -avP \
            -e "ssh -p 8804" \
            wangq@58.213.64.36:data/Bacteria/ASSEMBLY/{}/ \
            ~/data/Bacteria/ASSEMBLY/{}
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
#

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## Count species and strains

```shell
cd ~/data/Bacteria

echo "
    SELECT
        COUNT(DISTINCT species_id)
    FROM ar
    WHERE 1=1
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
#19169

echo "
    SELECT
        COUNT(DISTINCT species_id)
    FROM ar
    WHERE 1=1
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
#6580

echo "
    SELECT
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '%symbiont %'
        AND refseq_category IN ('reference genome', 'representative genome')
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-filter --str-not-in-fld 1:"[" \
    > summary/assembly_accession.lst

wc -l summary/assembly_accession.lst
#5275 assembly_accession.lst

cat ASSEMBLY/collect.csv |
    tsv-select -H -d, -f Taxid |
    sed '1d' |
    tsv-uniq |
    nwr append stdin -r species |
    tsv-select -f 2 |
    tsv-uniq |
    wc -l
#1689

cat summary/collect.pass.csv |
    tsv-select -H -d, -f Taxid |
    sed '1d' |
    tsv-uniq |
    nwr append stdin -r species |
    tsv-select -f 2 |
    tsv-uniq |
    wc -l
#1656

cat summary/collect.pass.csv |
    tsv-filter -H -d, --or \
        --istr-eq "RefSeq_category:reference genome" --istr-eq "RefSeq_category:representative genome" |
    grep -v "symbiont " |
    tsv-select -H -d, -f name |
    sed '1d' \
    > summary/representative.lst

wc -l summary/representative.lst
#1396 summary/representative.lst

cat summary/collect.pass.csv |
    sed -e '1d' |
    tr "," "\t" |
    tsv-select -f 1,3 |
    nwr append stdin -c 2 -r species -r genus -r family -r order \
    > summary/strains.taxon.tsv

cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    tsv-summarize -H -g 3 --count |
    tsv-filter -H --gt 2:200 |
    mlr --itsv --omd cat

```

| organism                      | count |
|-------------------------------|------:|
| Acinetobacter baumannii       |   531 |
| Bacillus cereus               |   413 |
| Bacillus subtilis             |   261 |
| Bacillus toyonensis           |   238 |
| Bacillus velezensis           |   278 |
| Bacteroides ovatus            |   201 |
| Bifidobacterium longum        |   224 |
| Bordetella pertussis          |   593 |
| Campylobacter coli            |   240 |
| Campylobacter jejuni          |   284 |
| Cronobacter sakazakii         |   208 |
| Enterobacter hormaechei       |   270 |
| Enterococcus faecium          |   313 |
| Erwinia amylovora             |   211 |
| Escherichia coli              |  2739 |
| Helicobacter pylori           |   395 |
| Klebsiella michiganensis      |   211 |
| Klebsiella pneumoniae         |  1611 |
| Klebsiella variicola          |   252 |
| Lactiplantibacillus plantarum |   211 |
| Limosilactobacillus reuteri   |   249 |
| Listeria monocytogenes        |   354 |
| Mycobacterium tuberculosis    |   546 |
| Mycobacteroides abscessus     |   946 |
| Pseudomonas aeruginosa        |   648 |
| Pseudomonas syringae          |   257 |
| Rhizobium leguminosarum       |   372 |
| Salmonella enterica           |  1377 |
| Shigella sonnei               |  1113 |
| Staphylococcus aureus         |  1124 |
| Staphylococcus haemolyticus   |   221 |
| Stenotrophomonas maltophilia  |   373 |
| Streptococcus pneumoniae      |   208 |
| Streptococcus pyogenes        |   263 |

## MinHash

### Compute MinHash

```shell
mkdir -p ~/data/Bacteria/mash
cd ~/data/Bacteria/mash

# Remove .msh files not in the list
find . -maxdepth 1 -mindepth 1 -type f -name "*.msh" |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        basename {} .msh
    ' |
    tsv-join --exclude -k 1 -f ../ASSEMBLY/url.tsv -d 1 |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        echo Remove {}
        rm {}.msh
    '

# Compute MinHash
cat ../ASSEMBLY/url.tsv |
    cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        if [[ -e {}.msh ]]; then
            exit
        fi

        2>&1 echo "==> {}"

        find ../ASSEMBLY/{} -name "*_genomic.fna.gz" |
            grep -v "_from_" |
            xargs cat |
            mash sketch -k 21 -s 100000 -p 2 - -I "{}" -o {}
    '

```

### Raw phylo-tree of representative assemblies

```shell
mkdir -p ~/data/Bacteria/tree
cd ~/data/Bacteria/tree

mash triangle -E -p 8 -l <(
    cat ../summary/representative.lst |
        parallel --no-run-if-empty --linebuffer -k -j 1 '
            if [[ -e ../mash/{}.msh ]]; then
                echo "../mash/{}.msh"
            fi
        '
    ) \
    > mash.dist.tsv

# Fill matrix with lower triangle
tsv-select -f 1-3 mash.dist.tsv |
    (tsv-select -f 2,1,3 mash.dist.tsv && cat) |
    (
        cut -f 1 mash.dist.tsv |
            tsv-uniq |
            parallel -j 1 --keep-order 'echo -e "{}\t{}\t0"' &&
        cat
    ) \
    > mash.dist_full.tsv

cat mash.dist_full.tsv |
    Rscript -e '
        library(readr);
        library(tidyr);
        library(ape);
        pair_dist <- read_tsv(file("stdin"), col_names=F);
        tmp <- pair_dist %>%
            pivot_wider( names_from = X2, values_from = X3, values_fill = list(X3 = 1.0) )
        tmp <- as.matrix(tmp)
        mat <- tmp[,-1]
        rownames(mat) <- tmp[,1]

        dist_mat <- as.dist(mat)
        clusters <- hclust(dist_mat, method = "ward.D2")
        tree <- as.phylo(clusters)
        write.tree(phy=tree, file="tree.nwk")

        group <- cutree(clusters, h=0.4) # k=5
        groups <- as.data.frame(group)
        groups$ids <- rownames(groups)
        rownames(groups) <- NULL
        groups <- groups[order(groups$group), ]
        write_tsv(groups, "groups.tsv")
    '

```

### Tweak the mash tree

```shell
cd ~/data/Bacteria/tree

# rank::col
ARRAY=(
    'order::6'
    'family::5'
#    'genus::4'
    'species::3'
)

rm mash.condensed.map
CUR_TREE=tree.nwk

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../summary/strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick mash.${GROUP_NAME}.newick
    cat condense.map >> mash.condensed.map

    CUR_TREE=mash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 mash.species.newick |
    rsvg-convert -o Bacteria.mash.png

```

## Non-redundant strains within species

This [paper](https://doi.org/10.1038/s41467-018-07641-9) showed that >95% intra-species and and <83%
inter-species ANI values.

```shell
cd ~/data/Bacteria

mkdir NR

cat summary/strains.taxon.tsv |
    tsv-summarize -H -g 3 --count |
    tsv-filter -H --ge 2:2 \
    > NR/species.tsv

tsv-summarize NR/species.tsv -H --count --sum count
#count   count_sum
#1651    35135

# each species
cat NR/species.tsv | sed '1d' | tsv-select -f 1 | #head -n 10 |
while read SPECIES; do
    1>&2 echo "==> ${SPECIES}"

    SPECIES_=$(
        echo "${SPECIES}" |
            tr " " "_"
    )

    mkdir -p NR/${SPECIES_}

    cat summary/strains.taxon.tsv |
        tsv-filter --str-eq "3:${SPECIES}" |
        tsv-select -f 1 \
        > "NR/${SPECIES_}/assembly.lst"

    1>&2 echo "    mash distances"

    if [[ ! -f "NR/${SPECIES_}/mash.dist.tsv" ]]; then
        mash triangle -E -p 8 -l <(
            cat "NR/${SPECIES_}/assembly.lst" |
                parallel --no-run-if-empty --linebuffer -k -j 1 '
                    if [[ -e mash/{}.msh ]]; then
                        echo "mash/{}.msh"
                    fi
                '
            ) \
            > "NR/${SPECIES_}/mash.dist.tsv"
    fi

    1>&2 echo "    List NR"

    cat "NR/${SPECIES_}/mash.dist.tsv" |
        tsv-filter --ff-str-ne 1:2 --le 3:0.01 \
        > "NR/${SPECIES_}/redundant.dist.tsv"

    cat "NR/${SPECIES_}/redundant.dist.tsv" |
        perl -nla -F"\t" -MGraph::Undirected -e '
            BEGIN {
                our $g = Graph::Undirected->new;
            }

            $g->add_edge($F[0], $F[1]);

            END {
                for my $cc ( $g->connected_components ) {
                    print join qq{\t}, sort @{$cc};
                }
            }
        ' \
        > "NR/${SPECIES_}/connected_components.tsv"

    cat "NR/${SPECIES_}/connected_components.tsv" |
        perl -nla -MPath::Tiny -F"\t" -e '
            BEGIN {
                our %rep = map { ($_, 1) } path(q{summary/representative.lst})->lines({chomp => 1});
            }

            # Representative strains are preferred
            if ( grep { $rep{$_} } @F ) {
                @F = grep { ! $rep{$_} } @F
            }
            else {
                shift @F;
            }
            printf qq{%s\n}, $_ for @F;
            ' \
        > "NR/${SPECIES_}/redundant.lst"

    cat "NR/${SPECIES_}/assembly.lst" |
        tsv-join --exclude -f "NR/${SPECIES_}/redundant.lst" \
        > "NR/${SPECIES_}/NR.lst"

done

find NR -name "redundant.lst" -size +0 | wc -l
#1093

find NR -name "redundant.lst" -empty | wc -l
#564

find NR -name "NR.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/NR.lst

wc -l summary/NR.lst
#8475 summary/NR.lst

```

### Genus

```shell
cd ~/data/Bacteria

cat summary/genus.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(
            cat summary/collect.pass.csv |
                sed "1d" |
                grep -F -w -f <( cat summary/NR.lst summary/representative.lst | sort | uniq ) |
                tsv-select -d, -f 3 |
                nwr append stdin -r genus -r species |
                grep {} |
                tsv-select -f 1,3 |
                tsv-uniq |
                wc -l
        )

        n_strains=$(
            cat summary/collect.pass.csv |
                sed "1d" |
                grep -F -w -f <( cat summary/NR.lst summary/representative.lst | sort | uniq ) |
                tsv-select -d, -f 3 |
                nwr append stdin -r genus |
                grep {} |
                wc -l
        )

        printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' |
    nwr append stdin --id |
    tsv-select -f 5,4,2,3 |
    tsv-sort -k2,2 |
    tsv-filter --ge 4:50 |
    (echo -e '#tax_id\tgenus\t#species\t#strains' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | genus               | #species | #strains |
|---------|---------------------|----------|----------|
| 469     | Acinetobacter       | 37       | 244      |
| 642     | Aeromonas           | 14       | 166      |
| 357     | Agrobacterium       | 14       | 79       |
| 1386    | Bacillus            | 97       | 359      |
| 816     | Bacteroides         | 30       | 207      |
| 1678    | Bifidobacterium     | 46       | 275      |
| 32008   | Burkholderia        | 33       | 123      |
| 194     | Campylobacter       | 54       | 106      |
| 544     | Citrobacter         | 15       | 91       |
| 1485    | Clostridium         | 48       | 123      |
| 102106  | Collinsella         | 3        | 54       |
| 1716    | Corynebacterium     | 51       | 96       |
| 547     | Enterobacter        | 18       | 247      |
| 1350    | Enterococcus        | 18       | 82       |
| 561     | Escherichia         | 20       | 181      |
| 724     | Haemophilus         | 12       | 80       |
| 209     | Helicobacter        | 55       | 289      |
| 570     | Klebsiella          | 18       | 121      |
| 1578    | Lactobacillus       | 35       | 100      |
| 1357    | Lactococcus         | 12       | 52       |
| 2742598 | Limosilactobacillus | 7        | 51       |
| 286     | Pseudomonas         | 161      | 511      |
| 379     | Rhizobium           | 31       | 159      |
| 1827    | Rhodococcus         | 18       | 54       |
| 590     | Salmonella          | 37       | 73       |
| 613     | Serratia            | 14       | 69       |
| 1279    | Staphylococcus      | 41       | 100      |
| 40323   | Stenotrophomonas    | 9        | 132      |
| 1301    | Streptococcus       | 70       | 370      |
| 1883    | Streptomyces        | 56       | 87       |
| 662     | Vibrio              | 57       | 364      |
| 338     | Xanthomonas         | 55       | 83       |

## Collect proteins

### `all.pro.fa`

```shell
cd ~/data/Bacteria

mkdir -p PROTEINS

for STRAIN in $(cat summary/representative.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz
done |
    pigz -p4 \
    > PROTEINS/all.pro.fa.gz

gzip -dcf PROTEINS/all.pro.fa.gz |
    perl -nl -e '
        BEGIN { our %seen; our $h; }

        if (/^>/) {
            $h = (split(" ", $_))[0];
            $seen{$h}++;
            $_ = $h;
        }
        print if $seen{$h} == 1;
    ' |
    pigz -p4 \
    > PROTEINS/all.uniq.fa.gz

# counting proteins
gzip -dcf PROTEINS/all.pro.fa.gz |
    grep "^>" |
    wc -l |
    numfmt --to=si
#5.2M

gzip -dcf PROTEINS/all.pro.fa.gz |
    grep "^>" |
    tsv-uniq |
    wc -l |
    numfmt --to=si
#5.1M

# annotations may be different
gzip -dcf PROTEINS/all.uniq.fa.gz |
    grep "^>" |
    wc -l |
    numfmt --to=si
#5.1M

# ribonuclease
gzip -dcf PROTEINS/all.pro.fa.gz |
    grep "ribonuclease" |
    grep -v "deoxyribonuclease" |
    perl -nl -e 's/^>\w+\.\d+\s+//g; print' |
    perl -nl -e 's/\s+\[.+?\]$//g; print' |
    perl -nl -e 's/MULTISPECIES: //g; print' |
    sort |
    uniq -c |
    sort -nr

```

### `all.replace.fa`

```shell
cd ~/data/Bacteria

rm PROTEINS/all.strain.tsv PROTEINS/all.replace.fa.gz
for STRAIN in $(cat summary/representative.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
        grep "^>" |
        cut -d" " -f 1 |
        sed "s/^>//" |
        STRAIN=${STRAIN} perl -nl -e '
            $n = $_;
            $s = $n;
            $s =~ s/\.\d+//;
            printf qq{%s\t%s_%s\t%s\n}, $n, $ENV{STRAIN}, $s, $ENV{STRAIN};
        ' \
    > PROTEINS/${STRAIN}.replace.tsv

    cut -f 2,3 PROTEINS/${STRAIN}.replace.tsv >> PROTEINS/all.strain.tsv

    faops replace -s \
        ASSEMBLY/${STRAIN}/*_protein.faa.gz \
        <(cut -f 1,2 PROTEINS/${STRAIN}.replace.tsv) \
        stdout |
        pigz -p4 \
        >> PROTEINS/all.replace.fa.gz

    rm PROTEINS/${STRAIN}.replace.tsv
done

gzip -dcf PROTEINS/all.replace.fa.gz |
    grep "^>" |
    wc -l |
    numfmt --to=si
#5.2M

(echo -e "#name\tstrain" && cat PROTEINS/all.strain.tsv)  \
    > temp &&
    mv temp PROTEINS/all.strain.tsv

faops size PROTEINS/all.replace.fa.gz > PROTEINS/all.replace.sizes

(echo -e "#name\tsize" && cat PROTEINS/all.replace.sizes) > PROTEINS/all.size.tsv

rm PROTEINS/all.replace.sizes

```

### `all.info.tsv`

```shell
cd ~/data/Bacteria

for STRAIN in $(cat summary/representative.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
        grep "^>" |
        sed "s/^>//" |
        perl -nl -e '/\[.+\[/ and s/\[/\(/; print' |
        perl -nl -e '/\].+\]/ and s/\]/\)/; print' |
        perl -nl -e 's/\s+\[.+?\]$//g; print' |
        perl -nl -e 's/MULTISPECIES: //g; print' |
        STRAIN=${STRAIN} perl -nl -e '
            /^(\w+)\.\d+\s+(.+)$/ or next;
            printf qq{%s_%s\t%s\n}, $ENV{STRAIN}, $1, $2;
        '
done \
    > PROTEINS/all.annotation.tsv

cat PROTEINS/all.annotation.tsv |
    wc -l |
    numfmt --to=si
#5.2M

(echo -e "#name\tannotation" && cat PROTEINS/all.annotation.tsv) \
    > temp &&
    mv temp PROTEINS/all.annotation.tsv

# check differences
cat PROTEINS/all.size.tsv |
    grep -F -f <(cut -f 1 PROTEINS/all.annotation.tsv) -v

tsv-join \
    PROTEINS/all.strain.tsv \
    --data-fields 1 \
    -f PROTEINS/all.size.tsv \
    --key-fields 1 \
    --append-fields 2 \
    > PROTEINS/all.strain_size.tsv

tsv-join \
    PROTEINS/all.strain_size.tsv \
    --data-fields 1 \
    -f PROTEINS/all.annotation.tsv \
    --key-fields 1 \
    --append-fields 2 \
    > PROTEINS/all.info.tsv

cat PROTEINS/all.info.tsv |
    wc -l |
    numfmt --to=si
#5.2M

```

## Phylogenetics with bac120

### Find corresponding proteins by `hmmsearch`

* Download HMM models as described in [`HMM.md`](HMM.md)

* The `E_VALUE` was manually adjusted to 1e-20 to reach a balance between sensitivity and
  speciality.

```shell
E_VALUE=1e-20

cd ~/data/Bacteria

# Find all genes
for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    >&2 echo "==> marker [${marker}]"

    mkdir -p PROTEINS/${marker}

    cat summary/representative.lst |
        parallel --no-run-if-empty --linebuffer -k -j 8 "
            gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/HMM/bac120/HMM/${marker}.HMM - |
                grep '>>' |
                perl -nl -e ' m{>>\s+(\S+)} and printf qq{%s\t%s\n}, \$1, {}; '
        " \
        > PROTEINS/${marker}/replace.tsv

    >&2 echo

done

```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Bacteria

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat PROTEINS/{}/replace.tsv |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
# 1385    1429    2168.25

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        cat PROTEINS/{}/replace.tsv |
            wc -l
    ' |
    paste - - |
    tsv-filter --invert --ge 2:1000 --le 2:2000 |
    cut -f 1 \
    > PROTEINS/bac120.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process
cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    grep -v -Fx -f PROTEINS/bac120.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        cat PROTEINS/{}/replace.tsv \
            > PROTEINS/{}/{}.replace.tsv

        faops some PROTEINS/all.uniq.fa.gz <(
            cat PROTEINS/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > PROTEINS/{}/{}.pro.fa
    '

# Align each markers with muscle
cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        >&2 echo "==> marker [{}]"
        if [ ! -s PROTEINS/{}/{}.pro.fa ]; then
            exit
        fi
        if [ -s PROTEINS/{}/{}.aln.fa ]; then
            exit
        fi

        muscle -quiet -in PROTEINS/{}/{}.pro.fa -out PROTEINS/{}/{}.aln.fa
    '

for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    >&2 echo "==> marker [${marker}]"
    if [ ! -s PROTEINS/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sometimes `muscle` can not produce alignments
    if [ ! -s PROTEINS/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # 1 name to many names
    cat PROTEINS/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s PROTEINS/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > PROTEINS/${marker}/${marker}.replace.fa
done

# Concat marker genes
for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    if [ ! -s PROTEINS/${marker}/${marker}.pro.fa ]; then
        continue
    fi
    if [ ! -s PROTEINS/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # sequences in one line
    faops filter -l 0 PROTEINS/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > PROTEINS/bac120.aln.fas

fasops concat PROTEINS/bac120.aln.fas summary/representative.lst -o PROTEINS/bac120.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in PROTEINS/bac120.aln.fa -out PROTEINS/bac120.trim.fa -automated1

faops size PROTEINS/bac120.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#94432
#15688

# To make it faster
FastTree -fastest -noml PROTEINS/bac120.trim.fa > PROTEINS/bac120.trim.newick

```

### Tweak the concat tree

```shell
cd ~/data/Bacteria/tree

cp ../PROTEINS/bac120.trim.newick .

# rank::col
ARRAY=(
    'order::6'
    'family::5'
#    'genus::4'
    'species::3'
)

rm bac120.condensed.map
CUR_TREE=bac120.trim.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../summary/strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick bac120.${GROUP_NAME}.newick
    cat condense.map >> bac120.condensed.map

    CUR_TREE=bac120.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 bac120.species.newick |
    rsvg-convert -o Bacteria.bac120.png

```

## InterProScan on all proteins of representative and typical strains

To be filled by other projects

```shell
mkdir -p ~/data/Bacteria/STRAINS

```

## Early divergence of Bacteria

![molbiolevolmsn247f02_ht.jpeg](images%2Fmolbiolevolmsn247f02_ht.jpeg)

### Terrabacteria group

Bacillota == Firmicutes

Actinomycetota == Actinobacteria

```shell
mkdir -p ~/data/Bacteria/Terrabacteria
cd ~/data/Bacteria/Terrabacteria

FAMILY=$(
    nwr member "Terrabacteria group" -r family |
        sed '1d' |
        cut -f 1 |
        nwr append stdin -r order |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND family_id IN ($FAMILY)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    grep -v -i "symbiont " |
    tsv-filter --str-not-in-fld 1:"[" \
    > tmp.tsv

cat ../ASSEMBLY/collect.csv |
    grep -F -f <(cut -f 2 tmp.tsv) |
    grep -F -f <(cat ../summary/NR.lst ../summary/representative.lst | sort | uniq) |
    tsv-select -d, -f 1,3 |
    tr "," "\t" |
    nwr append stdin -c 2 -r species -r genus -r family -r order |
    sed 's/Pseudomonas paraeruginosa/Pseudomonas aeruginosa/g' \
    > strains.taxon.tsv

wc -l strains.taxon.tsv
#3106 strains.taxon.tsv

```

### Proteobacteria

Pseudomonadota == Proteobacteria

```shell
mkdir -p ~/data/Bacteria/Proteobacteria
cd ~/data/Bacteria/Proteobacteria

FAMILY=$(
    nwr member Proteobacteria -r family |
        sed '1d' |
        cut -f 1 |
        nwr append stdin -r order |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND family_id IN ($FAMILY)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    grep -v -i "symbiont " |
    tsv-filter --str-not-in-fld 1:"[" \
    > tmp.tsv

cat ../ASSEMBLY/collect.csv |
    grep -F -f <(cut -f 2 tmp.tsv) |
    grep -F -f <(cat ../summary/NR.lst ../summary/representative.lst | sort | uniq) |
    tsv-select -d, -f 1,3 |
    tr "," "\t" |
    nwr append stdin -c 2 -r species -r genus -r family -r order \
    > strains.taxon.tsv

wc -l strains.taxon.tsv
#4905 strains.taxon.tsv

```
