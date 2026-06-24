# *Escherichia*

[TOC levels=2-4]: #
  * [Strain info](#strain-info)
    * [Symlink](#symlink)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
  * [All assemblies](#all-assemblies)
    * [Extract from `../Bacteria` and create assembly.tsv](#extract-from-bacteria-and-create-assemblytsv)
    * [Count `assembly.tsv`](#count-assemblytsv)
  * [MinHash](#minhash)
  * [Count valid species and strains](#count-valid-species-and-strains)
    * [For *genomic alignments* and *protein families*](#for-genomic-alignments-and-protein-families)
    * [Count strains - Genus](#count-strains---genus)
  * [Collect proteins](#collect-proteins)
  * [Phylogenetics with bac120](#phylogenetics-with-bac120)
    * [Find corresponding representative proteins by `hmmsearch`](#find-corresponding-representative-proteins-by-hmmsearch)
    * [Domain related protein sequences](#domain-related-protein-sequences)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
    * [Condense branches in the protein tree](#condense-branches-in-the-protein-tree)

## Strain info

* [Escherichia](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=561)
* [Shigella](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=620)
* [Salmonella](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=590)
* [Klebsiella](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=570)

### List all ranks

```shell
mkdir -p ~/data/Escherichia
cd ~/data/Escherichia

nwr member Enterobacterales -r family |
    grep -v -i "Candidatus" |
    tva sort -H -k 2 |
    tva to md --num

nwr member Enterobacterales |
    tva stats -H -g 3 --count |
    tva to md --num

nwr member Enterobacterales -r "species group" -r "species subgroup" |
    tva select -f 1-3 |
    tva sort -H -k 3,2 |
    tva to md --num

```

| #tax_id | sci_name            | rank   | division |
| ------: | ------------------- | ------ | -------- |
| 2812006 | Bruguierivoracaceae | family | Bacteria |
| 1903416 | Budviciaceae        | family | Bacteria |
|     543 | Enterobacteriaceae  | family | Bacteria |
| 1903409 | Erwiniaceae         | family | Bacteria |
| 3120698 | Gallaecimonadaceae  | family | Bacteria |
| 1903412 | Hafniaceae          | family | Bacteria |
| 1903414 | Morganellaceae      | family | Bacteria |
| 1903410 | Pectobacteriaceae   | family | Bacteria |
| 1771359 | Thorselliaceae      | family | Bacteria |
| 1903411 | Yersiniaceae        | family | Bacteria |

| rank             | count |
| ---------------- | ----: |
| clade            |     5 |
| family           |    11 |
| forma specialis  |   563 |
| genus            |   138 |
| isolate          |    35 |
| no rank          |  2297 |
| order            |     1 |
| serogroup        |   117 |
| serotype         |    51 |
| species          | 30823 |
| species group    |     5 |
| species subgroup |     6 |
| strain           |  6770 |
| subspecies       |    76 |
| varietas         |     1 |

| #tax_id | sci_name                             | rank             |
| ------: | ------------------------------------ | ---------------- |
| 1344959 | Citrobacter freundii complex         | species group    |
|  354276 | Enterobacter cloacae complex         | species group    |
| 3390273 | Klebsiella pneumoniae complex        | species group    |
| 1654067 | Pantoea agglomerans group            | species group    |
| 1649845 | Yersinia pseudotuberculosis complex  | species group    |
| 2302426 | Enterobacter cloacae complex clade K | species subgroup |
| 2302427 | Enterobacter cloacae complex clade L | species subgroup |
| 2302428 | Enterobacter cloacae complex clade N | species subgroup |
| 2302429 | Enterobacter cloacae complex clade O | species subgroup |
| 2302430 | Enterobacter cloacae complex clade P | species subgroup |
| 2302431 | Enterobacter cloacae complex clade S | species subgroup |

### Count assemblies - Genus

```bash
mkdir -p ~/data/Escherichia/summary
cd ~/data/Escherichia/summary

GENUS=$(
    nwr member Enterobacterales -r genus |
        sed '1d' |
        tva select -f 1 |
        sort |
        uniq
)

for S in $GENUS; do
    RS=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND genus_id = $S
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    CHR=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND genus_id = $S
                AND assembly_level IN ('Complete Genome', 'Chromosome')
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    if [[ ${RS} -gt 0 ]]; then
        echo -e "$S\t$RS\t$CHR"
    fi
done |
    nwr append stdin |
    tva select -f 1,4,2-3 |
    tva sort -k 2 |
    tva sort -n -r -k 3,4 |
    (echo -e 'genus_id\tgenus\tRS\tCHR' && cat) \
    > genus.RS.tsv

cat genus.RS.tsv |
    tva filter -H --ge CHR:20 |
    tva to md --num

```

| genus_id | genus                     |    RS |  CHR |
| -------: | ------------------------- | ----: | ---: |
|      561 | Escherichia               | 53725 | 5621 |
|      570 | Klebsiella                | 37224 | 4427 |
|      590 | Salmonella                | 20563 | 2407 |
|      547 | Enterobacter              |  9600 | 1082 |
|      544 | Citrobacter               |  3485 |  498 |
|      613 | Serratia                  |  3119 |  410 |
|      620 | Shigella                  |  3088 |  234 |
|      629 | Yersinia                  |  2258 |  446 |
|      583 | Proteus                   |  1898 |  242 |
|    53335 | Pantoea                   |  1388 |  178 |
|      586 | Providencia               |  1341 |  175 |
|   413496 | Cronobacter               |   804 |   66 |
|   122277 | Pectobacterium            |   688 |  134 |
|      551 | Erwinia                   |   591 |  186 |
|      581 | Morganella                |   548 |   85 |
|   204037 | Dickeya                   |   270 |  118 |
|      635 | Edwardsiella              |   262 |  214 |
|      626 | Xenorhabdus               |   219 |   20 |
|    83654 | Leclercia                 |   176 |   41 |
|    32199 | Buchnera                  |   149 |  142 |
|      568 | Hafnia                    |   139 |   28 |
|  1330547 | Kosakonia                 |   135 |   40 |
|    84565 | Sodalis                   |    41 |   28 |
|   203804 | Candidatus Blochmanniella |    26 |   26 |

### Species with assemblies

* 'RefSeq'
* The number of strains in GenBank is simply too large to be used

```shell
mkdir -p ~/data/Escherichia/summary
cd ~/data/Escherichia/summary

# should have a valid name of genus
nwr member Enterobacterales Pasteurellales -r genus |
    sed '1d' |
    tva sort -n -k 1 \
    > genus.list.tsv

wc -l genus.list.tsv
#176 genus.list.tsv

cat genus.list.tsv | tva select -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tva sort -k 2 \
    > RS1.tsv

cat genus.list.tsv | tva select -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tva sort -k 2 \
    > GB1.tsv

wc -l RS*.tsv GB*.tsv
# 6890 RS1.tsv
# 7324 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tva stats --sum 3
        fi
    done
done
# RS1	147611
# GB1	1274712
```

## Download all assemblies

### Create assembly.tsv

```bash
cd ~/data/Escherichia/summary

echo "
.headers ON
    SELECT
        *
    FROM ar
    WHERE 1=1
        AND species IN ('Escherichia coli', 'Klebsiella pneumoniae', 'Salmonella enterica', 'Shigella sonnei')
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tva select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession \
    > raw.tsv

# RS
SPECIES=$(
    cat RS1.tsv |
        tva select -f 1 |
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

cat raw.tsv |
    tva uniq |
    tva check
#142113 lines, 7 fields

# Create abbr.
cat raw.tsv |
    grep -v '^#' |
    tva uniq |
    tva select -f 1-6 |
    nwr abbr -c "1,2,3" -m 3 --shortsub |
    tva uniq -H -f ftp_path |
    tva uniq -H -f 7 |
    sed '1d' |
    tva select -f 7,4,5,2,6 |
    (echo -e '#name\tftp_path\tbiosample\tspecies\tassembly_level' && cat ) |
    tva filter -H --or --str-in-fld 2:ftp --str-in-fld 2:http |
    tva sort -H -k 4,1 \
    > Escherichia.assembly.tsv

tva check < Escherichia.assembly.tsv
#142107 lines, 5 fields

# find potential duplicate strains or assemblies
cat Escherichia.assembly.tsv |
    tva uniq -f 1 --repeated

cat Escherichia.assembly.tsv |
    tva filter --str-not-in-fld 2:ftp

# Edit .assembly.tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Escherichia.assembly.tsv
#
# Save the file to another directory to prevent accidentally changing it
# cp Escherichia.assembly.tsv ~/Scripts/genomes/assembly

# Cleaning
rm raw*.*sv
```

### Count before download

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class

```shell
cd ~/data/Escherichia

nwr template ~/Scripts/genomes/assembly/Escherichia.assembly.tsv \
    --count \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    tva to md --fmt

# .lst and .count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    tva sort -H -k 1 |
    tva filter -H --ge 3:100 |
    tva to md --fmt
```

| item    |   count |
| ------- | ------: |
| strain  | 141,658 |
| species |     786 |
| genus   |     144 |
| family  |      12 |
| order   |       2 |
| class   |       1 |

| genus           | #species | #strains |
| --------------- | -------: | -------: |
| Actinobacillus  |       17 |      271 |
| Aggregatibacter |        5 |      248 |
| Atlantibacter   |        3 |      130 |
| Avibacterium    |        6 |      218 |
| Buchnera        |        1 |      136 |
| Citrobacter     |       22 |    3,447 |
| Cronobacter     |        8 |      795 |
| Dickeya         |       14 |      257 |
| Edwardsiella    |        5 |      250 |
| Enterobacter    |       31 |    9,338 |
| Erwinia         |       28 |      591 |
| Escherichia     |        7 |   50,793 |
| Gallibacterium  |        9 |      124 |
| Glaesserella    |        3 |      346 |
| Haemophilus     |       13 |    1,276 |
| Hafnia          |        4 |      129 |
| Klebsiella      |       20 |   36,253 |
| Kluyvera        |       11 |      114 |
| Kosakonia       |       13 |      135 |
| Leclercia       |        5 |      172 |
| Lonsdalea       |        4 |      119 |
| Mannheimia      |       11 |      413 |
| Morganella      |        2 |      529 |
| Pantoea         |       49 |    1,370 |
| Pasteurella     |        8 |      801 |
| Pectobacterium  |       25 |      682 |
| Photorhabdus    |       26 |      176 |
| Proteus         |       13 |    1,787 |
| Providencia     |       19 |    1,306 |
| Rahnella        |       16 |      120 |
| Salmonella      |        3 |   19,762 |
| Serratia        |       25 |    2,767 |
| Shigella        |        5 |    3,076 |
| Xenorhabdus     |       37 |      213 |
| Yersinia        |       27 |    2,197 |

### Download and check

```shell
cd ~/data/Escherichia

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Escherichia.assembly.tsv \
    --ass

# Run
bash ASSEMBLY/aria2.sh

# Check md5; create check.lst
# rm ASSEMBLY/check.lst
bash ASSEMBLY/check.sh

# Remove failed directories and re-download
bash ASSEMBLY/check.sh 2>&1 |
    grep "checksum failed" |
    sed 's/.*==> //;s/ checksum failed <==//' |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        dir=$(cat ASSEMBLY/url.tsv | tva filter --str-eq "1:{}" | tva select -f 3,1 | tr "\t" "/")
        if [[ -n "$dir" && -e "ASSEMBLY/$dir" ]]; then
            echo Remove ASSEMBLY/$dir
            rm -fr "ASSEMBLY/$dir"
        fi
    '

# # Put the misplaced directory into the right place
# bash ASSEMBLY/reorder.sh

# # This operation will delete some files in the directory, so please be careful
# cat ASSEMBLY/remove.lst |
#    parallel --no-run-if-empty --linebuffer -k -j 1 '
#        if [[ -e "ASSEMBLY/{}" ]]; then
#            echo Remove {}
#            rm -fr "ASSEMBLY/{}"
#        fi
#    '

find ASSEMBLY/ -name "*_genomic.fna.gz" |
    grep -v "_from_" |
    wc -l
#26830

# N50 C S; create n50.tsv and n50.pass.tsv
bash ASSEMBLY/n50.sh 50000 500 500000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tva filter -H --str-in-fld "name:_GCF_" |
    tva stats -H --min "N50" --max "C" --min "S" |
    tva transpose
# N50_min	5041
# C_max	1987
# S_min	231767

cat ASSEMBLY/n50.tsv |
    tva stats -H --quantile "N50:0.1,0.5" --quantile "C:0.5,0.9" --quantile "S:0.1,0.5" |
    tva transpose
# N50_quantile_0.1	78792
# N50_quantile_0.5	215003.5
# C_quantile_0.5	91
# C_quantile_0.9	271
# S_quantile_0.1	4595768.6
# S_quantile_0.5	5096136.5

# After the above steps are completed, run the following commands.

# Collect; create collect.tsv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/
cp ASSEMBLY/omit.lst summary/
cp ASSEMBLY/pass.lst summary/
cp ASSEMBLY/sp.lst summary/
cp ASSEMBLY/rep.lst summary/

cat ASSEMBLY/counts.tsv |
    tva to md --fmt
```

| #item            | fields |   lines |
| ---------------- | -----: | ------: |
| url.tsv          |      3 | 142,106 |
| check.lst        |      1 | 142,106 |
| collect.tsv      |     20 | 142,107 |
| n50.pass.tsv     |      4 | 132,573 |
| collect.pass.tsv |     23 | 132,573 |
| pass.lst         |      1 | 132,572 |
| omit.lst         |      1 |      82 |
| rep.lst          |      1 |     685 |
| sp.lst           |      1 |   6,165 |

### Rsync to hpcc

```shell
rsync -avP \
    ~/data/Escherichia/ \
    wangq@202.119.37.251:data/Escherichia

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Escherichia/ \
    wangq@58.213.64.36:data/Escherichia

rsync -avP \
    ~/data/Escherichia/ \
    wangq@192.168.31.209:/share/data/Escherichia

# rsync -avP wangq@202.119.37.251:data/Escherichia/ ~/data/Escherichia
```

## BioSample

```shell
cd ~/data/Escherichia

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Escherichia.assembly.tsv \
    --bs

# Run this script twice and it will re-download the failed files
bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 50

tva check < BioSample/biosample.tsv
#141734 lines, 240 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/
```

## MinHash

```bash
cd ~/data/Escherichia/

# relaxed thresholds
nwr template ~/Scripts/genomes/assembly/Escherichia.assembly.tsv \
    --mh \
    --parallel 8 \
    --in summary/pass.lst \
    --ani-ab 0.12 \
    --ani-nr 0.005

# Compute assembly sketches
bash MinHash/compute.sh

# Non-redundant strains within species
bash MinHash/nr.sh

fd --full-path "MinHash/.+/NR.lst" -X cat |
    sort |
    uniq \
    > summary/NR.lst
fd --full-path "MinHash/.+/redundant.lst" -X cat |
    sort |
    uniq \
    > summary/redundant.lst
wc -l summary/NR.lst summary/redundant.lst
#    15574 summary/NR.lst
#   116785 summary/redundant.lst

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst |
    tva join -e -f summary/sp.lst \
    > MinHash/tmp.lst
mv MinHash/tmp.lst summary/abnormal.lst

wc -l MinHash/abnormal.lst summary/abnormal.lst
#  281 MinHash/abnormal.lst
#  143 summary/abnormal.lst

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Escherichia/

nwr template ~/Scripts/genomes/assembly/Escherichia.assembly.tsv \
    --mh \
    --parallel 8 \
    --in summary/rep.lst \
    --not-in summary/sp.lst \
    --not-in summary/abnormal.lst \
    --not-in summary/redundant.lst \
    --height 0.4

bash MinHash/dist.sh
```

### Condense branches in the minhash tree

```shell
mkdir -p ~/data/Escherichia/tree
cd ~/data/Escherichia/tree

pgr nwk reroot ../MinHash/tree.nwk -n Wig_glossinidia_GCF_000247565_1 -n Pas_multo_FDAARGOS_218_GCF_002073255_2 |
    pgr nwk order stdin --nd --an \
    > minhash.reroot.newick

pgr pl condense --map -t ../Count/strains.taxon.tsv -r 4 -r 3 \
    minhash.reroot.newick |
    pgr nwk order stdin --nd --an \
    > minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

# svg
pgr nwk to-svg minhash.condensed.newick \
    > Escherichia.minhash.svg
```

## Count valid species and strains

### For *genomic alignments* and *protein families*

```bash
cd ~/data/Escherichia/

nwr template ~/Scripts/genomes/assembly/Escherichia.assembly.tsv \
    --count \
    --in summary/pass.lst \
    --not-in summary/abnormal.lst \
    --rank family --rank genus \
    --lineage genus

# strains.taxon.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    tva to md --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/family.count.tsv |
    tva filter -H --ge "3:50" |
    tva to md --num

cat Count/genus.count.tsv |
    tva filter -H --ge "3:50" |
    tva to md --num

# Can accept N_COUNT
bash Count/lineage.sh 100

cat Count/lineage.count.tsv |
    tva to md --num

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv
cp Count/genus.count.tsv summary/genus.genome.tsv
```

| item    |  count |
| ------- | -----: |
| strain  | 132046 |
| species |    764 |
| genus   |    140 |
| family  |     12 |
| order   |      2 |
| class   |      1 |

| family             | #species | #strains |
| ------------------ | -------: | -------: |
| Enterobacteriaceae |      241 |   116714 |
| Erwiniaceae        |      108 |     2047 |
| Hafniaceae         |       11 |      377 |
| Morganellaceae     |      108 |     3703 |
| Pasteurellaceae    |      118 |     3542 |
| Pectobacteriaceae  |       63 |     1091 |
| Yersiniaceae       |       89 |     4520 |

| genus           | #species | #strains |
| --------------- | -------: | -------: |
| Actinobacillus  |       17 |      248 |
| Aggregatibacter |        5 |      212 |
| Atlantibacter   |        3 |      120 |
| Avibacterium    |        6 |      195 |
| Brenneria       |       13 |       54 |
| Citrobacter     |       22 |     3342 |
| Cronobacter     |        8 |      760 |
| Dickeya         |       14 |      239 |
| Edwardsiella    |        5 |      247 |
| Enterobacter    |       31 |     9013 |
| Erwinia         |       27 |      576 |
| Escherichia     |        7 |    47911 |
| Gallibacterium  |        8 |       98 |
| Glaesserella    |        3 |      157 |
| Haemophilus     |       13 |     1197 |
| Hafnia          |        4 |      120 |
| Histophilus     |        1 |       60 |
| Klebsiella      |       20 |    35082 |
| Kluyvera        |       11 |      102 |
| Kosakonia       |       13 |      134 |
| Leclercia       |        5 |      159 |
| Lelliottia      |        6 |       94 |
| Lonsdalea       |        4 |      109 |
| Mannheimia      |       11 |      403 |
| Morganella      |        2 |      512 |
| Pantoea         |       49 |     1315 |
| Pasteurella     |        8 |      783 |
| Pectobacterium  |       25 |      675 |
| Photorhabdus    |       25 |      153 |
| Plesiomonas     |        2 |       61 |
| Proteus         |       12 |     1625 |
| Providencia     |       19 |     1200 |
| Rahnella        |       16 |      115 |
| Salmonella      |        3 |    19105 |
| Serratia        |       25 |     2614 |
| Shigella        |        5 |      439 |
| Xenorhabdus     |       35 |      152 |
| Yersinia        |       27 |     1703 |

| #genus          | species                               | count |
| --------------- | ------------------------------------- | ----: |
| Actinobacillus  | Actinobacillus pleuropneumoniae       |   179 |
| Aggregatibacter | Aggregatibacter actinomycetemcomitans |   159 |
| Avibacterium    | Avibacterium paragallinarum           |   115 |
| Citrobacter     | Citrobacter amalonaticus              |   108 |
|                 | Citrobacter braakii                   |   254 |
|                 | Citrobacter freundii                  |  1628 |
|                 | Citrobacter koseri                    |   225 |
|                 | Citrobacter portucalensis             |   382 |
|                 | Citrobacter sp.                       |   423 |
| Cronobacter     | Cronobacter sakazakii                 |   541 |
| Enterobacter    | Enterobacter asburiae                 |   605 |
|                 | Enterobacter bugandensis              |   256 |
|                 | Enterobacter cloacae                  |   599 |
|                 | Enterobacter hormaechei               |  4828 |
|                 | Enterobacter kobei                    |   512 |
|                 | Enterobacter ludwigii                 |   339 |
|                 | Enterobacter roggenkampii             |   683 |
|                 | Enterobacter sp.                      |   768 |
| Erwinia         | Erwinia amylovora                     |   325 |
|                 | Erwinia aphidicola                    |   103 |
| Escherichia     | Escherichia albertii                  |   297 |
|                 | Escherichia coli                      | 46947 |
|                 | Escherichia fergusonii                |   196 |
|                 | Escherichia marmotae                  |   134 |
|                 | Escherichia sp.                       |   306 |
| Glaesserella    | Glaesserella parasuis                 |   150 |
| Haemophilus     | Haemophilus influenzae                |   942 |
| Klebsiella      | Klebsiella aerogenes                  |   766 |
|                 | Klebsiella grimontii                  |   243 |
|                 | Klebsiella michiganensis              |   778 |
|                 | Klebsiella ornithinolytica            |   269 |
|                 | Klebsiella oxytoca                    |   572 |
|                 | Klebsiella planticola                 |   113 |
|                 | Klebsiella pneumoniae                 | 28402 |
|                 | Klebsiella quasipneumoniae            |  1767 |
|                 | Klebsiella sp.                        |   760 |
|                 | Klebsiella variicola                  |  1203 |
| Mannheimia      | Mannheimia haemolytica                |   366 |
| Morganella      | Morganella morganii                   |   504 |
| Pantoea         | Pantoea agglomerans                   |   352 |
|                 | Pantoea ananatis                      |   303 |
|                 | Pantoea sp.                           |   266 |
| Pasteurella     | Pasteurella multocida                 |   689 |
| Pectobacterium  | Pectobacterium brasiliense            |   164 |
|                 | Pectobacterium versatile              |   110 |
| Proteus         | Proteus mirabilis                     |  1389 |
| Providencia     | Providencia huaxiensis                |   102 |
|                 | Providencia rettgeri                  |   200 |
|                 | Providencia sp.                       |   422 |
|                 | Providencia stuartii                  |   201 |
| Salmonella      | Salmonella enterica                   | 16759 |
|                 | Salmonella sp.                        |  2314 |
| Serratia        | Serratia bockelmannii                 |   113 |
|                 | Serratia fonticola                    |   155 |
|                 | Serratia marcescens                   |  1477 |
|                 | Serratia nevei                        |   188 |
|                 | Serratia ureilytica                   |   207 |
| Shigella        | Shigella flexneri                     |   201 |
|                 | Shigella sonnei                       |   134 |
| Yersinia        | Yersinia enterocolitica               |   692 |
|                 | Yersinia pestis                       |   416 |
|                 | Yersinia pseudotuberculosis           |   125 |
|                 | Yersinia ruckeri                      |   169 |

### For *protein families*

```shell
cd ~/data/Escherichia/

nwr template ~/Scripts/genomes/assembly/Escherichia.assembly.tsv \
    --count \
    --in summary/pass.lst \
    --not-in summary/abnormal.lst \
    --not-in summary/omit.lst \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    tva to md --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    tva filter -H --ge "3:50" |
    tva to md --num

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv
```

| item    |  count |
| ------- | -----: |
| strain  | 132046 |
| species |    764 |
| genus   |    140 |
| family  |     12 |
| order   |      2 |
| class   |      1 |

| genus           | #species | #strains |
| --------------- | -------: | -------: |
| Actinobacillus  |       17 |      248 |
| Aggregatibacter |        5 |      212 |
| Atlantibacter   |        3 |      120 |
| Avibacterium    |        6 |      195 |
| Brenneria       |       13 |       54 |
| Citrobacter     |       22 |     3342 |
| Cronobacter     |        8 |      760 |
| Dickeya         |       14 |      239 |
| Edwardsiella    |        5 |      247 |
| Enterobacter    |       31 |     9013 |
| Erwinia         |       27 |      576 |
| Escherichia     |        7 |    47911 |
| Gallibacterium  |        8 |       98 |
| Glaesserella    |        3 |      157 |
| Haemophilus     |       13 |     1197 |
| Hafnia          |        4 |      120 |
| Histophilus     |        1 |       60 |
| Klebsiella      |       20 |    35082 |
| Kluyvera        |       11 |      102 |
| Kosakonia       |       13 |      134 |
| Leclercia       |        5 |      159 |
| Lelliottia      |        6 |       94 |
| Lonsdalea       |        4 |      109 |
| Mannheimia      |       11 |      403 |
| Morganella      |        2 |      512 |
| Pantoea         |       49 |     1315 |
| Pasteurella     |        8 |      783 |
| Pectobacterium  |       25 |      675 |
| Photorhabdus    |       25 |      153 |
| Plesiomonas     |        2 |       61 |
| Proteus         |       12 |     1625 |
| Providencia     |       19 |     1200 |
| Rahnella        |       16 |      115 |
| Salmonella      |        3 |    19105 |
| Serratia        |       25 |     2614 |
| Shigella        |        5 |      439 |
| Xenorhabdus     |       35 |      152 |
| Yersinia        |       27 |     1703 |

## Collect proteins

```bash
cd ~/data/Escherichia/

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Escherichia.assembly.tsv \
    --pro \
    --parallel 8 \
    --in summary/pass.lst \
    --not-in summary/omit.lst

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
    tva stats -H --count --sum 2-7 |
    sed 's/^count/species/' |
    tva transpose |
    (echo -e "#item\tcount" && cat) |
    tva to md --fmt
```

| #item      |       count |
|------------|------------:|
| species    |         170 |
| strain_sum |      91,356 |
| total_sum  | 440,702,715 |
| dedup_sum  |  15,979,554 |
| rep_sum    |   4,086,492 |
| fam88_sum  |   2,453,485 |
| fam38_sum  |   1,752,300 |

## Phylogenetics with bac120

```bash
cd ~/data/Escherichia/

# The Bacteria HMM set
nwr kb bac120 -o HMM
cp HMM/bac120.lst HMM/marker.lst

```

### Find corresponding representative proteins by `hmmsearch`

```bash
cd ~/data/Escherichia

cat Protein/species.tsv |
    tva join -f summary/pass.lst -k 1 |
    tva join -e -f summary/omit.lst -k 1 \
    > Protein/species-f.tsv

cat Protein/species-f.tsv |
    tva select -f 2 |
    tva uniq |
while read SPECIES; do
    if [[ -s Protein/"${SPECIES}"/bac120.tsv ]]; then
        continue
    fi
    if [[ ! -f Protein/"${SPECIES}"/rep_seq.fa.gz ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    cat HMM/marker.lst |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 "
            gzip -dcf Protein/${SPECIES}/rep_seq.fa.gz |
                hmmsearch --cut_nc --noali --notextw HMM/hmm/{}.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), q({}), \$1; '
        " \
        > Protein/${SPECIES}/bac120.tsv
done

cat Protein/species-f.tsv |
    tva select -f 2 |
    tva uniq |
while read SPECIES; do
    if [[ ! -s Protein/"${SPECIES}"/bac120.tsv ]]; then
        continue
    fi
    if [[ ! -f Protein/"${SPECIES}"/seq.sqlite ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    nwr seqdb -d Protein/${SPECIES} --rep f3=Protein/${SPECIES}/bac120.tsv

done
```

### Domain related protein sequences

```bash
cd ~/data/Escherichia

mkdir -p Domain

# each assembly
cat Protein/species-f.tsv |
    tva select -f 2 |
    tva uniq |
while read SPECIES; do
    if [[ ! -f Protein/"${SPECIES}"/seq.sqlite ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    echo "
        SELECT
            seq.name,
            asm.name,
            rep.f3
        FROM asm_seq
        JOIN rep_seq ON asm_seq.seq_id = rep_seq.seq_id
        JOIN seq ON asm_seq.seq_id = seq.id
        JOIN rep ON rep_seq.rep_id = rep.id
        JOIN asm ON asm_seq.asm_id = asm.id
        WHERE 1=1
            AND rep.f3 IS NOT NULL
        ORDER BY
            asm.name,
            rep.f3
        " |
        sqlite3 -tabs Protein/${SPECIES}/seq.sqlite \
        > Protein/${SPECIES}/seq_asm_f3.tsv

    hnsm some Protein/"${SPECIES}"/pro.fa.gz <(
            tsv-select -f 1 Protein/"${SPECIES}"/seq_asm_f3.tsv |
            rgr dedup stdin
        )
done |
    pgr fa dedup stdin |
    pgr fa gz stdin -o Domain/bac120.fa

fd --full-path "Protein/.+/seq_asm_f3.tsv" -X cat \
    > Domain/seq_asm_f3.tsv

cat Domain/seq_asm_f3.tsv |
    tva join -e -d 2 -f summary/redundant.lst -k 1 |
    tva join -e -d 2 -f summary/sp.lst -k 1 \
    > Domain/seq_asm_f3.NR.tsv
```

### Align and concat marker genes to create species tree

```bash
cd ~/data/Escherichia

# Extract proteins
cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

        mkdir -p Domain/{}

        pgr fa some Domain/bac120.fa.gz <(
            cat Domain/seq_asm_f3.tsv |
                tva filter --str-eq "3:{}" |
                tva select -f 1 |
                tva uniq
            ) \
            > Domain/{}/{}.pro.fa
    '

# Align each marker
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

cat HMM/marker.lst |
while read marker; do
    echo >&2 "==> marker [${marker}]"
    if [ ! -s Domain/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sometimes `muscle` can not produce alignments
    if [ ! -s Domain/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # Only NR strains
    # 1 name to many names
    cat Domain/seq_asm_f3.NR.tsv |
        tva filter --str-eq "3:${marker}" |
        tva select -f 1-2 |
        pgr fa replace -s Domain/${marker}/${marker}.aln.fa stdin \
        > Domain/${marker}/${marker}.replace.fa
done

# Concat marker genes
cat HMM/marker.lst |
while read marker; do
    if [ ! -s Domain/${marker}/${marker}.pro.fa ]; then
        continue
    fi
    if [ ! -s Domain/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    cat Domain/${marker}/${marker}.replace.fa

    # empty line for .fas
    echo
done \
    > Domain/bac120.aln.fas

cat Domain/seq_asm_f3.NR.tsv |
    cut -f 2 |
    tva uniq |
    sort |
    fasops concat Domain/bac120.aln.fas stdin -o Domain/bac120.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Domain/bac120.aln.fa -out Domain/bac120.trim.fa -automated1

pgr fa size Domain/bac120.*.fa |
    tva uniq -f 2 |
    cut -f 2
#71124
#46096

# To make it faster
FastTree -fastest -noml Domain/bac120.trim.fa > Domain/bac120.trim.newick
```

### Condense branches in the protein tree

```shell
mkdir -p ~/data/Escherichia/tree
cd ~/data/Escherichia/tree

cat ../Domain/bac120.trim.newick |
    nwr order stdin --nd --an \
    > bac120.order.newick

nwr pl-condense --map -r species \
    bac120.order.newick ../Count/species.tsv |
    nwr order stdin --nd --an \
    > bac120.condensed.newick

mv condensed.tsv bac120.condense.tsv

# svg
nwr topo --bl bac120.condensed.newick | # remove comments
    nw_display -s -b 'visibility:hidden' -w 1200 -v 20 - \
    > Escherichia.bac120.svg

```
