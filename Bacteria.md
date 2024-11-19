# Bacteria

All genomes of *Bacteria*, species by species.

Download all genomes and analyze representative strains.

<!-- TOC -->

* [Bacteria](#bacteria)
    * [Taxon info](#taxon-info)
        * [List all ranks](#list-all-ranks)
        * [Species with assemblies](#species-with-assemblies)
        * [Model organisms](#model-organisms)
    * [Download all assemblies](#download-all-assemblies)
        * [Create assembly.tsv](#create-assemblytsv)
        * [Count before download](#count-before-download)
        * [Download and check](#download-and-check)
        * [Remove unnecessary files](#remove-unnecessary-files)
        * [Rsync to hpcc](#rsync-to-hpcc)
    * [BioSample](#biosample)
    * [Early divergence of Bacteria](#early-divergence-of-bacteria)
        * [ReRoot](#reroot)
    * [MinHash](#minhash)
        * [Condense branches in the minhash tree](#condense-branches-in-the-minhash-tree)
    * [Count valid species and strains](#count-valid-species-and-strains)
        * [For *genomic alignments* and *protein
          families*](#for-genomic-alignments-and-protein-families)
    * [InterProScan on all proteins of representative and typical strains](#interproscan-on-all-proteins-of-representative-and-typical-strains)

<!-- TOC -->

## Taxon info

* [Bacteria](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2)

### List all ranks

```shell
nwr member Bacteria |
    grep -v " sp." |
    grep -v " x " |
    tsv-summarize -H -g rank --count |
    rgr md stdin --fmt

```

| rank             |   count |
|------------------|--------:|
| superkingdom     |       1 |
| kingdom          |       2 |
| phylum           |     181 |
| class            |     191 |
| no rank          |   7,899 |
| species          | 114,945 |
| genus            |   5,438 |
| order            |     415 |
| family           |     977 |
| strain           |  42,083 |
| clade            |     134 |
| isolate          |     456 |
| subspecies       |     712 |
| serogroup        |     152 |
| species group    |      98 |
| subgenus         |       1 |
| tribe            |       2 |
| species subgroup |      30 |
| varietas         |      24 |
| serotype         |     255 |
| biotype          |       7 |
| pathogroup       |       5 |
| subfamily        |       1 |
| suborder         |       7 |
| subclass         |       6 |
| forma            |       4 |

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
#5438 genus.list

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
#   175 RS2.tsv
#  2176 RS3.tsv
#  9291 RS4.tsv
#    28 GB1.tsv
#   252 GB2.tsv
#  2328 GB3.tsv
# 11005 GB4.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     17750
#RS2     241357
#RS3     40091
#RS4     327053
#GB1     27564
#GB2     1777182
#GB3     51991
#GB4     1883056

cat RS2.tsv |
    tsv-join -f RS1.tsv -k 1 -e |
    tsv-summarize --sum 3
#78557

cat RS3.tsv |
    tsv-join -f RS1.tsv -k 1 -e |
    tsv-join -f RS2.tsv -k 1 -e |
    tsv-summarize --sum 3
#13016

cat RS4.tsv |
    tsv-join -f RS1.tsv -k 1 -e |
    tsv-join -f RS2.tsv -k 1 -e |
    tsv-join -f RS3.tsv -k 1 -e |
    tsv-summarize --sum 3
#32941

cat RS4.tsv |
    tsv-join -f RS1.tsv -k 1 -e |
    tsv-join -f RS2.tsv -k 1 -e |
    tsv-summarize --sum 3
#85696

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
#na      365512
#reference genome        22429

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
        HAVING count >= 500
        ) AS subquery
    JOIN ar ON ar.species_id = subquery.species_id
    WHERE 1=1
        AND ar.refseq_category IN ('reference genome')
        AND ar.tax_id != ar.species_id
    ORDER BY ar.organism_name
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > reference.tsv

cat reference.tsv |
    sed '1s/^/#/' |
    nwr append stdin -r phylum -r class |
    tsv-select -H -f 1-2,phylum,class |
    keep-header -- sort -t$'\t' -k3,3 -k4,4 -k2,2 |
    sed 's/strain=//' |
    parallel --col-sep "\t" -j 1 '
        if [[ "{3}" == "Proteobacteria" || "{3}" == "Pseudomonadota" ]]; then
            printf "%s\t%s\t%s\n" {1} {2} {4}
        else
            printf "%s\t%s\t%s\n" {1} {2} {3}
        fi
    ' |
    rgr md stdin

cp reference.tsv ~/Scripts/genomes/assembly/Bacteria.reference.tsv

```

| #tax_id | organism_name                                                    | phylum              |
|---------|------------------------------------------------------------------|---------------------|
| 565042  | Bifidobacterium longum subsp. longum JCM 1217                    | Actinomycetota      |
| 83332   | Mycobacterium tuberculosis H37Rv                                 | Actinomycetota      |
| 561007  | Mycobacteroides abscessus ATCC 19977                             | Actinomycetota      |
| 224308  | Bacillus subtilis subsp. subtilis str. 168                       | Bacillota           |
| 527031  | Bacillus thuringiensis serovar berliner ATCC 10792               | Bacillota           |
| 326423  | Bacillus velezensis FZB42                                        | Bacillota           |
| 1169293 | Enterococcus faecalis EnGen0336                                  | Bacillota           |
| 169963  | Listeria monocytogenes EGD-e                                     | Bacillota           |
| 93061   | Staphylococcus aureus subsp. aureus NCTC 8325                    | Bacillota           |
| 40041   | Streptococcus equi subsp. zooepidemicus                          | Bacillota           |
| 568814  | Streptococcus suis BM407                                         | Bacillota           |
| 272559  | Bacteroides fragilis NCTC 9343                                   | Bacteroidota        |
| 435590  | Phocaeicola vulgatus ATCC 8482                                   | Bacteroidota        |
| 192222  | Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819      | Campylobacterota    |
| 224914  | Brucella melitensis bv. 1 str. 16M                               | Alphaproteobacteria |
| 871585  | Acinetobacter pittii PHEA-2                                      | Gammaproteobacteria |
| 386585  | Escherichia coli O157:H7 str. Sakai                              | Gammaproteobacteria |
| 511145  | Escherichia coli str. K-12 substr. MG1655                        | Gammaproteobacteria |
| 1450527 | Francisella tularensis subsp. novicida D9876                     | Gammaproteobacteria |
| 1125630 | Klebsiella pneumoniae subsp. pneumoniae HS11286                  | Gammaproteobacteria |
| 2590157 | Klebsiella variicola subsp. variicola                            | Gammaproteobacteria |
| 529507  | Proteus mirabilis HI4320                                         | Gammaproteobacteria |
| 208964  | Pseudomonas aeruginosa PAO1                                      | Gammaproteobacteria |
| 99287   | Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 | Gammaproteobacteria |
| 198214  | Shigella flexneri 2a str. 301                                    | Gammaproteobacteria |
| 223926  | Vibrio parahaemolyticus RIMD 2210633                             | Gammaproteobacteria |
| 359385  | Xanthomonas campestris pv. raphani                               | Gammaproteobacteria |
| 1035377 | Yersinia pestis A1122                                            | Gammaproteobacteria |

## Download all assemblies

### Create assembly.tsv

* RefSeq
    * '>= 2 genomes'
        * RS4 - genome_rep: 'Full'

```shell
cd ~/data/Bacteria/summary

cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession,infraspecific_name |
    sed 's/strain=//' |
    keep-header -- perl -nla -F'\t' -e '$F[0] = qq($F[1] $F[7]); print join qq(\t), @F;' |
    tsv-select -f 1-7 \
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

cat raw.tsv |
    rgr dedup stdin |
    datamash check
#327082 lines, 7 fields

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
    > Bacteria.assembly.tsv

datamash check < Bacteria.assembly.tsv
#326966 lines, 5 fields

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
    rgr md stdin --fmt

# .lst and .count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    keep-header -- tsv-sort -k1,1 |
    tsv-filter -H --ge 3:1000 |
    rgr md stdin --num

```

| item    |   count |
|---------|--------:|
| strain  | 326,963 |
| species |   9,289 |
| genus   |   2,111 |
| family  |     488 |
| order   |     204 |
| class   |      88 |

| genus               | #species | #strains |
|---------------------|---------:|---------:|
| Acinetobacter       |       73 |    12299 |
| Aeromonas           |       26 |     1640 |
| Bacillus            |       81 |     8175 |
| Bacteroides         |       47 |     3487 |
| Bifidobacterium     |       59 |     2690 |
| Bordetella          |       18 |     1378 |
| Brucella            |       28 |     1300 |
| Burkholderia        |       38 |     5042 |
| Campylobacter       |       42 |     6021 |
| Citrobacter         |       18 |     2073 |
| Clostridioides      |        1 |     3010 |
| Clostridium         |       67 |     2579 |
| Corynebacterium     |      130 |     1855 |
| Enterobacter        |       27 |     6289 |
| Enterococcus        |       52 |     9606 |
| Escherichia         |        6 |    40739 |
| Francisella         |       11 |     1077 |
| Haemophilus         |       11 |     1097 |
| Helicobacter        |       37 |     3213 |
| Klebsiella          |       14 |    26229 |
| Lactiplantibacillus |       14 |     1239 |
| Lactobacillus       |       39 |     2001 |
| Lactococcus         |       17 |     1029 |
| Legionella          |       51 |     1234 |
| Listeria            |       19 |     5990 |
| Mycobacterium       |       77 |     8053 |
| Mycobacteroides     |        6 |     2114 |
| Neisseria           |       36 |     3464 |
| Phocaeicola         |       12 |     1012 |
| Proteus             |       10 |     1494 |
| Pseudomonas         |      247 |    15607 |
| Salmonella          |        2 |    14234 |
| Serratia            |       22 |     2092 |
| Shigella            |        4 |     2641 |
| Staphylococcus      |       56 |    23311 |
| Stenotrophomonas    |       21 |     1011 |
| Streptococcus       |       88 |    19802 |
| Streptomyces        |      425 |     4414 |
| Vibrio              |      114 |     7560 |
| Xanthomonas         |       38 |     2952 |
| Yersinia            |       26 |     1637 |

### Download and check

* A standardized archaeal taxonomy for the Genome Taxonomy Database. Nat Microbiol 6, 946–959
  (2021).
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

#cat ASSEMBLY/url.tsv |
#    tsv-join -f ASSEMBLY/check.lst -k 1 -e |
#    tsv-select -f 3 |
#    sed 's/_/\t/' |
#    tsv-select -f 1 |
#    tsv-summarize --group-by 1 --count |
#    sort |
#    tsv-filter --ge 2:100
#
#scp ~/data/Bacteria/ASSEMBLY/check.lst wangq@192.168.31.163:data/Bacteria/ASSEMBLY/

# from another machine
#bash ASSEMBLY/rsync.sh Klebsiella
#rsync -avP \
#    -e "ssh -T -c aes128-gcm@openssh.com -o Compression=no -x" \
#    wangq@192.168.31.163:data/Bacteria/ASSEMBLY/ ~/data/Bacteria/ASSEMBLY \
#    --exclude="check.lst" \
#    --exclude="url.tsv"

## Put the misplaced directories into the right ones
#bash ASSEMBLY/reorder.sh
##
## This operation will delete some files in the directory, so please be careful
#cat ASSEMBLY/remove.lst |
#    parallel --no-run-if-empty --linebuffer -k -j 1 '
#        if [[ -e "ASSEMBLY/{}" ]]; then
#            echo Remove {}
#            rm -fr "ASSEMBLY/{}"
#        fi
#    '

find ASSEMBLY/ -name "*_genomic.gff.gz" | wc -l
#326965
find ASSEMBLY -type f -name "*_protein.faa.gz" -size -4096c

#fd _genomic.gff.gz D:\data\Bacteria\ASSEMBLY -tf | Measure-Object -Line
#fd _genomic.gff.gz ~/data/Bacteria/ASSEMBLY -tf | wc -l

#fd _protein.faa.gz D:\data\Bacteria\ASSEMBLY -tf --size +10k | Measure-Object -Line

# N50 C S; create n50.tsv and n50.pass.tsv
bash ASSEMBLY/n50.sh 20000 500 100000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50" --max "C" --min "S"
#N50_min C_max   S_min
#5007    1996    107786

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "N50:0.1,0.5" --quantile "C:0.5,0.9" --quantile "S:0.1,0.5" |
    datamash transpose
#N50_pct10       49582.4
#N50_pct50       199874
#C_pct50 65
#C_pct90 229
#S_pct10 2070817.4
#S_pct50 4476404

# After the above steps are completed, run the following commands.

# Collect; create collect.tsv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    rgr md stdin --fmt

```

| #item            | fields |   lines |
|------------------|-------:|--------:|
| url.tsv          |      3 | 326,965 |
| check.lst        |      1 | 326,965 |
| collect.tsv      |     20 | 326,966 |
| n50.tsv          |      4 | 326,966 |
| n50.pass.tsv     |      4 | 316,617 |
| collect.pass.tsv |     23 | 316,617 |
| pass.lst         |      1 | 316,616 |
| omit.lst         |      0 |       0 |
| rep.lst          |      1 |   9,040 |
| sp.lst           |      1 |      47 |

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

# Windows
# fd assembly_stats.txt D:\data\Bacteria\ASSEMBLY -tf | rm -force

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
    rgr dedup stdin |
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
    rgr dedup stdin |
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

# fd 'SAM.*\.txt' D:\data\Bacteria\BioSample -tf | Measure-Object -Line
# fd 'SAM.*\.txt' ~/data/Bacteria/BioSample -tf | wc -l

# Ignore rare attributes
bash BioSample/collect.sh 500

datamash check < BioSample/biosample.tsv
# 317193 lines, 132 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

# Save spaces
cd ~/data/Bacteria
cat BioSample/sample.tsv |
    tsv-select -f 3 |
    rgr dedup stdin |
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

# Restore dir
cd ~/data/Bacteria
cat BioSample/sample.tsv |
    tsv-select -f 3 |
    rgr dedup stdin |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        cd BioSample
        if [[ -d {} ]]; then
            exit
        fi
        if [[ ! -f {}.tar ]]; then
            exit
        fi

        echo -e "==> {}"
        tar -xvf {}.tar
        rm {}.tar
    '

```

## Early divergence of Bacteria

![molbiolevolmsn247f02_ht.jpeg](images%2Fmolbiolevolmsn247f02_ht.jpeg)

### ReRoot

Thermotogota and Aquificota are basal taxa of bacteria.

* For the latter steps, use the following two as the outgroups
    * Thermot_petrop_RKU_1_GCF_000016785_1
    * Hydrogenob_thermophilus_TK_6_GCF_000164905_1

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
#Aquificota      Desulfurobacterium thermolithotrophum   2
#Aquificota      Hydrogenobacter thermophilus    5
#Aquificota      Sulfurihydrogenibium azorense   2
#Thermotogota    Fervidobacterium changbaicum    2
#Thermotogota    Fervidobacterium pennivorans    3
#Thermotogota    Geotoga petraea 3
#Thermotogota    Kosmotoga pacifica      3
#Thermotogota    Mesotoga prima  16
#Thermotogota    Oceanotoga teriensis    2
#Thermotogota    Petrotoga olearia       2
#Thermotogota    Petrotoga sibirica      2
#Thermotogota    Pseudothermotoga elfii  2
#Thermotogota    Pseudothermotoga hypogea        3
#Thermotogota    Pseudothermotoga lettingae      2
#Thermotogota    Thermosipho africanus   3
#Thermotogota    Thermosipho melanesiensis       7
#Thermotogota    Thermotoga maritima     6
#Thermotogota    Thermotoga neapolitana  2
#Thermotogota    Thermotoga petrophila   3

cat summary/collect.pass.tsv |
    tsv-filter -H --not-blank RefSeq_category |
    tsv-filter -H --or \
        --str-in-fld "2:Desulfurobacterium" \
        --str-in-fld "2:Hydrogenobacter" \
        --str-in-fld "2:Sulfurihydrogenibium" \
        --str-in-fld "2:Fervidobacterium" \
        --str-in-fld "2:Pseudothermotoga" \
        --str-in-fld "2:Thermosipho" \
        --str-in-fld "2:Thermotoga" |
    grep -v -Fw -f ASSEMBLY/omit.lst |
    tsv-select -H -f "#name,Assembly_level,Assembly_method,Genome_coverage,Sequencing_technology"

```

### Groups

* Bacillota
    * Bacillota == Firmicutes
    * Mycoplasmatota == Tenericutes
    * Clostridia belongs to Bacillota

* Terrabacteria group
    * Actinomycetota == Actinobacteria
    * Deinococcota == Deinococcus-Thermus
    * Cyanobacteriota == Cyanobacteria
    * Chloroflexi belongs to Chloroflexota

* Pseudomonadota
    * Pseudomonadota == Proteobacteria

* FCB group

* The rest

```shell
nwr member Bacillota Mycoplasmatota -r family

nwr member "Terrabacteria group" -r family |
    sed '1d' |
    nwr restrict Bacillota --exclude -f stdin -c 1 |
    nwr restrict Mycoplasmatota --exclude -f stdin -c 1

nwr member Pseudomonadota -r family

nwr member "FCB group" -r family

nwr member Bacteria -r family |
    nwr restrict "Terrabacteria group" --exclude -f stdin -c 1 |
    nwr restrict Pseudomonadota --exclude -f stdin -c 1 |
    nwr restrict "FCB group" --exclude -f stdin -c 1

```

## MinHash

```shell
cd ~/data/Bacteria/

nwr template ~/Scripts/genomes/assembly/Bacteria.assembly.tsv \
    --mh \
    --parallel 8 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.10 \
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
#   72756 summary/NR.lst
#  243658 summary/redundant.lst

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#1669

# All representative should be in NR
cat ASSEMBLY/rep.lst |
    grep -v -F -f summary/NR.lst

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Bacteria/

nwr template ~/Scripts/genomes/assembly/Bacteria.assembly.tsv \
    --mh \
    --parallel 8 \
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

nw_reroot ../MinHash/tree.nwk Thermot_petrophila_RKU_1_GCF_000016785_1 Hydrogenob_thermophilus_TK_6_GCF_000164905_1 |
    nwr order stdin --nd --an \
    > minhash.reroot.newick

nwr pl-condense --map -r order -r family -r genus \
    minhash.reroot.newick ../MinHash/species.tsv |
    nwr order stdin --nd --an \
    > minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

# svg
nwr topo --bl minhash.condensed.newick | # remove comments
    nw_display -s -b 'visibility:hidden' -w 1200 -v 20 - \
    > Bacteria.minhash.svg

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

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/order.count.tsv |
    tsv-filter -H --ge "3:1000" |
    rgr md stdin --num

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:1000" |
    rgr md stdin --num

# Can accept N_COUNT
bash Count/lineage.sh 1000

cat Count/lineage.count.tsv |
    rgr md stdin --num

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| item    |  count |
|---------|-------:|
| strain  | 314945 |
| species |   9260 |
| genus   |   2105 |
| family  |    485 |
| order   |    202 |
| class   |     87 |

| order                | #species | #strains |
|----------------------|---------:|---------:|
| Aeromonadales        |       36 |     1636 |
| Alteromonadales      |      166 |     1140 |
| Bacillales           |      713 |    40770 |
| Bacteroidales        |      199 |     7500 |
| Bifidobacteriales    |       72 |     2820 |
| Burkholderiales      |      412 |     9051 |
| Campylobacterales    |      115 |     9450 |
| Enterobacterales     |      407 |   100235 |
| Erysipelotrichales   |       49 |     1090 |
| Eubacteriales        |      204 |     3508 |
| Flavobacteriales     |      413 |     2832 |
| Hyphomicrobiales     |      493 |     4738 |
| Kitasatosporales     |      456 |     4175 |
| Lachnospirales       |      153 |     2268 |
| Lactobacillales      |      518 |    38445 |
| Legionellales        |       57 |     1402 |
| Lysobacterales       |      131 |     4185 |
| Micrococcales        |      522 |     2681 |
| Moraxellales         |      118 |    12376 |
| Mycobacteriales      |      467 |    13529 |
| Mycoplasmoidales     |       56 |     1288 |
| Neisseriales         |       94 |     3712 |
| Pasteurellales       |       72 |     3046 |
| Peptostreptococcales |       43 |     3169 |
| Pseudomonadales      |      330 |    15521 |
| Rhodobacterales      |      233 |     1316 |
| Thiotrichales        |       36 |     1123 |
| Vibrionales          |      159 |     7681 |

| genus               | #species | #strains |
|---------------------|---------:|---------:|
| Acinetobacter       |       73 |    11922 |
| Aeromonas           |       26 |     1603 |
| Bacillus            |       81 |     7937 |
| Bacteroides         |       47 |     3292 |
| Bifidobacterium     |       59 |     2672 |
| Bordetella          |       18 |     1131 |
| Brucella            |       28 |     1287 |
| Burkholderia        |       38 |     4976 |
| Campylobacter       |       42 |     5863 |
| Citrobacter         |       18 |     2032 |
| Clostridioides      |        1 |     2790 |
| Clostridium         |       67 |     2503 |
| Corynebacterium     |      130 |     1796 |
| Enterobacter        |       27 |     6168 |
| Enterococcus        |       52 |     9410 |
| Escherichia         |        6 |    39441 |
| Haemophilus         |       11 |     1071 |
| Helicobacter        |       37 |     3100 |
| Klebsiella          |       14 |    25695 |
| Lactiplantibacillus |       14 |     1237 |
| Lactobacillus       |       38 |     1812 |
| Lactococcus         |       17 |     1016 |
| Legionella          |       51 |     1230 |
| Listeria            |       19 |     5964 |
| Mycobacterium       |       77 |     7892 |
| Mycobacteroides     |        6 |     2108 |
| Neisseria           |       36 |     3216 |
| Proteus             |       10 |     1415 |
| Pseudomonas         |      247 |    14776 |
| Salmonella          |        2 |    14009 |
| Serratia            |       22 |     2026 |
| Shigella            |        4 |     2291 |
| Staphylococcus      |       56 |    22892 |
| Streptococcus       |       87 |    19424 |
| Streptomyces        |      423 |     4034 |
| Vibrio              |      114 |     7110 |
| Xanthomonas         |       38 |     2764 |
| Yersinia            |       26 |     1457 |

| #family               | genus               | species                       | count |
|-----------------------|---------------------|-------------------------------|------:|
| Bacillaceae           | Bacillus            | Bacillus cereus               |  1136 |
| Burkholderiaceae      | Burkholderia        | Burkholderia pseudomallei     |  1853 |
| Campylobacteraceae    | Campylobacter       | Campylobacter coli            |  1603 |
|                       |                     | Campylobacter jejuni          |  3156 |
| Enterobacteriaceae    | Citrobacter         | Citrobacter freundii          |  1083 |
|                       | Enterobacter        | Enterobacter hormaechei       |  3697 |
|                       | Escherichia         | Escherichia coli              | 38841 |
|                       | Klebsiella          | Klebsiella pneumoniae         | 21442 |
|                       |                     | Klebsiella quasipneumoniae    |  1262 |
|                       | Salmonella          | Salmonella enterica           | 13978 |
|                       | Shigella            | Shigella sonnei               |  1449 |
| Enterococcaceae       | Enterococcus        | Enterococcus faecalis         |  3728 |
|                       |                     | Enterococcus faecium          |  3929 |
| Helicobacteraceae     | Helicobacter        | Helicobacter pylori           |  2448 |
| Lactobacillaceae      | Lactiplantibacillus | Lactiplantibacillus plantarum |  1074 |
| Legionellaceae        | Legionella          | Legionella pneumophila        |  1034 |
| Listeriaceae          | Listeria            | Listeria monocytogenes        |  5209 |
| Moraxellaceae         | Acinetobacter       | Acinetobacter baumannii       |  9470 |
| Morganellaceae        | Proteus             | Proteus mirabilis             |  1189 |
| Mycobacteriaceae      | Mycobacterium       | Mycobacterium tuberculosis    |  7104 |
|                       | Mycobacteroides     | Mycobacteroides abscessus     |  2001 |
| Neisseriaceae         | Neisseria           | Neisseria gonorrhoeae         |  1231 |
|                       |                     | Neisseria meningitidis        |  1603 |
| Peptostreptococcaceae | Clostridioides      | Clostridioides difficile      |  2790 |
| Pseudomonadaceae      | Pseudomonas         | Pseudomonas aeruginosa        |  9706 |
|                       |                     | Pseudomonas viridiflava       |  1417 |
| Staphylococcaceae     | Staphylococcus      | Staphylococcus aureus         | 17054 |
|                       |                     | Staphylococcus epidermidis    |  1623 |
| Streptococcaceae      | Streptococcus       | Streptococcus agalactiae      |  1891 |
|                       |                     | Streptococcus pneumoniae      |  8851 |
|                       |                     | Streptococcus pyogenes        |  2426 |
|                       |                     | Streptococcus suis            |  2346 |
| Vibrionaceae          | Vibrio              | Vibrio cholerae               |  1880 |
|                       |                     | Vibrio parahaemolyticus       |  2312 |
| Yersiniaceae          | Serratia            | Serratia marcescens           |  1227 |

## Collect proteins

```shell
cd ~/data/Bacteria/

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Bacteria.assembly.tsv \
    --pro \
    --parallel 8 \
    --in ASSEMBLY/pass.lst \
    --not-in ASSEMBLY/omit.lst

# collect proteins
bash Protein/collect.sh

find ~/data/Bacteria/Protein -name "pro.fa.gz" -exec stat --format="%s" {} + |
    perl -nl -MNumber::Format -e '
        BEGIN { $total = 0 }
        $total += $_;
        END { print Number::Format::format_bytes($total, precision => 1) }
        '

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

| #item      |         count |
|------------|--------------:|
| species    |         9,261 |
| strain_sum |       316,616 |
| total_sum  | 1,225,154,555 |
| dedup_sum  |   164,288,277 |
| rep_sum    |    72,138,906 |
| fam88_sum  |    57,316,703 |
| fam38_sum  |    44,263,319 |

## InterProScan on all proteins of representative and typical strains

To be filled by other projects

```shell
mkdir -p ~/data/Bacteria/STRAINS

```

