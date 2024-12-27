# Build alignments across a eukaryotic taxonomy rank

Genus *Trichoderma* as an example.

<!-- TOC -->
* [Build alignments across a eukaryotic taxonomy rank](#build-alignments-across-a-eukaryotic-taxonomy-rank)
  * [Taxon info](#taxon-info)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
  * [Download all assemblies](#download-all-assemblies)
    * [Create .assembly.tsv](#create-assemblytsv)
    * [Count before download](#count-before-download)
    * [Download and check](#download-and-check)
    * [Rsync to hpcc](#rsync-to-hpcc)
  * [BioSample](#biosample)
  * [MinHash](#minhash)
    * [Condense branches in the minhash tree](#condense-branches-in-the-minhash-tree)
  * [Count valid species and strains](#count-valid-species-and-strains)
    * [For *genomic alignments*](#for-genomic-alignments)
    * [For *protein families*](#for-protein-families)
  * [Collect proteins](#collect-proteins)
  * [Phylogenetics with fungi61](#phylogenetics-with-fungi61)
    * [Find corresponding representative proteins by `hmmsearch`](#find-corresponding-representative-proteins-by-hmmsearch)
    * [Domain related protein sequences](#domain-related-protein-sequences)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
    * [The protein tree](#the-protein-tree)
  * [Groups and targets](#groups-and-targets)
  * [Prepare sequences for `egaz`](#prepare-sequences-for-egaz)
  * [Generate alignments](#generate-alignments)
<!-- TOC -->

## Taxon info

* [Trichoderma](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=5543)
* [Entrez records](http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5543)
* [WGS](https://www.ncbi.nlm.nih.gov/Traces/wgs/?view=wgs&search=Trichoderma) is now useless and can
  be ignored.

A nice review article about Trichoderma:

Woo, S.L. et al.
Trichoderma: a multipurpose, plant-beneficial microorganism for eco-sustainable agriculture.
Nat Rev Microbiol 21, 312–326 (2023).
https://doi.org/10.1038/s41579-022-00819-5

### List all ranks

There are no noteworthy classification ranks other than species.

```shell
nwr member Trichoderma |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    rgr md stdin --num

nwr lineage Trichoderma |
    tsv-filter --str-ne 1:clade |
    tsv-filter --str-ne "1:no rank" |
    sed -n '/kingdom\tFungi/,$p' |
    sed -E "s/\b(genus)\b/*\1*/"| # Highlight genus
    (echo -e '#rank\tsci_name\ttax_id' && cat) |
    rgr md stdin

```

| rank     | count |
|----------|------:|
| genus    |     1 |
| species  |   476 |
| no rank  |     1 |
| varietas |     2 |
| strain   |    14 |
| forma    |     2 |

| #rank      | sci_name          | tax_id |
|------------|-------------------|--------|
| kingdom    | Fungi             | 4751   |
| subkingdom | Dikarya           | 451864 |
| phylum     | Ascomycota        | 4890   |
| subphylum  | Pezizomycotina    | 147538 |
| class      | Sordariomycetes   | 147550 |
| subclass   | Hypocreomycetidae | 222543 |
| order      | Hypocreales       | 5125   |
| family     | Hypocreaceae      | 5129   |
| *genus*    | Trichoderma       | 5543   |

### Species with assemblies

The family Hypocreaceae as outgroups.

```shell
mkdir -p ~/data/Trichoderma/summary
cd ~/data/Trichoderma/summary

# should have a valid name of genus
nwr member Hypocreaceae -r genus |
    grep -v " x " |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#18 genus.list

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
            AND species NOT LIKE '% x %' -- Crossbreeding of two species
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
#   9 RS1.tsv
#  57 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     9
#GB1     172

```

## Download all assemblies

### Create .assembly.tsv

This step is pretty important

* `nwr kb formats` will give the requirements for `.assembly.tsv`.

* The naming of assemblies has two aspects:
    * for program operation they are unique identifiers;
    * for researchers, they should provide taxonomic information.

If a RefSeq assembly is available, the corresponding GenBank one will not be listed

```shell
cd ~/data/Trichoderma/summary

# Reference genome
echo "
.headers ON
    SELECT
        *
    FROM ar
    WHERE 1=1
        AND species IN ('Saccharomyces cerevisiae')
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
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% x %'
        AND genome_rep IN ('Full')
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

echo "
    SELECT
        genus || ' sp. ' || infraspecific_name || ' ' || assembly_accession AS name,
        genus || ' sp.', genus, ftp_path, biosample, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species LIKE '% sp.%'
        AND species NOT LIKE '% x %'
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

cat raw.tsv |
    rgr dedup stdin |
    datamash check
#174 lines, 7 fields

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
    > Trichoderma.assembly.tsv

datamash check < Trichoderma.assembly.tsv
#173 lines, 5 fields

# find potential duplicate strains or assemblies
cat Trichoderma.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Trichoderma.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

# Edit .assembly.tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Trichoderma.assembly.tsv
#
# Save the file to another directory to prevent accidentally changing it
# cp Trichoderma.assembly.tsv ~/Scripts/genomes/assembly

# Cleaning
rm raw*.*sv

```

### Count before download

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class

```shell
cd ~/data/Trichoderma

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
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
    rgr md stdin --num

```

| item    | count |
|---------|------:|
| strain  |   172 |
| species |    44 |
| genus   |     7 |
| family  |     2 |
| order   |     2 |
| class   |     2 |

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Cladobotryum     |        2 |        3 |
| Escovopsis       |        2 |        7 |
| Hypomyces        |        2 |        2 |
| Mycogone         |        1 |        1 |
| Saccharomyces    |        1 |        1 |
| Sphaerostilbella |        1 |        1 |
| Trichoderma      |       35 |      157 |

### Download and check

* When `rsync.sh` is interrupted, run `check.sh` before restarting
* For projects that have finished downloading, but have renamed strains, you can run `reorder.sh` to
  avoid re-downloading
    * `misplaced.tsv`
    * `remove.list`
* The parameters of `n50.sh` should be determined by the distribution of the description statistics
* `collect.sh` generates a file of type `.tsv`, which is intended to be opened by spreadsheet
  software.
    * Information of assemblies are collected from *_assembly_report.txt *after* downloading
    * **Note**: `*_assembly_report.txt` have `CRLF` at the end of the line.
* `finish.sh` generates the following files
    * `omit.lst` - no annotations
    * `collect.pass.tsv` - passes the n50 check
    * `pass.lst` - passes the n50 check
    * `rep.lst` - representative or reference strains
    * `counts.tsv`

```shell
cd ~/data/Trichoderma

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --ass

# Run
bash ASSEMBLY/rsync.sh

# Check md5; create check.lst
# rm ASSEMBLY/check.lst
bash ASSEMBLY/check.sh

## Put the misplaced directory into the right place
#bash ASSEMBLY/reorder.sh
#
## This operation will delete some files in the directory, so please be careful
#cat ASSEMBLY/remove.lst |
#    parallel --no-run-if-empty --linebuffer -k -j 1 '
#        if [[ -e "ASSEMBLY/{}" ]]; then
#            echo Remove {}
#            rm -fr "ASSEMBLY/{}"
#        fi
#    '

# N50 C S; create n50.tsv and n50.pass.tsv
# LEN_N50   N_CONTIG    LEN_SUM
bash ASSEMBLY/n50.sh 100000 1000 1000000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50" --max "C" --min "S"
#N50_min C_max   S_min
#579860  533     33215161

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "N50:0.1,0.5" --quantile "C:0.5,0.9" --quantile "S:0.1,0.5" |
    datamash transpose
#N50_pct10       103504.6
#N50_pct50       1412965.5
#C_pct50 157
#C_pct90 1044.8
#S_pct10 32206756.6
#S_pct50 37298829.5

# After the above steps are completed, run the following commands.

# Collect; create collect.tsv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    rgr md stdin --fmt

```

| #item            | fields | lines |
|------------------|-------:|------:|
| url.tsv          |      3 |   172 |
| check.lst        |      1 |   172 |
| collect.tsv      |     20 |   173 |
| n50.tsv          |      4 |   173 |
| n50.pass.tsv     |      4 |   151 |
| collect.pass.tsv |     23 |   151 |
| pass.lst         |      1 |   150 |
| omit.lst         |      1 |   128 |
| rep.lst          |      1 |    48 |
| sp.lst           |      1 |    15 |

### Rsync to hpcc

```shell
rsync -avP \
    ~/data/Trichoderma/ \
    wangq@202.119.37.251:data/Trichoderma

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Trichoderma/ \
    wangq@58.213.64.36:data/Trichoderma

# rsync -avP wangq@202.119.37.251:data/Trichoderma/ ~/data/Trichoderma

# rsync -avP -e 'ssh -p 8804' wangq@58.213.64.36:data/Trichoderma/ ~/data/Trichoderma

```

## BioSample

ENA's BioSample missed many strains, so NCBI's was used.

```shell
cd ~/data/Trichoderma

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --bs

bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 10

datamash check < BioSample/biosample.tsv
#170 lines, 39 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## MinHash

Estimate nucleotide divergences among strains.

* Abnormal strains
    * This [paper](https://doi.org/10.1038/s41467-018-07641-9) showed that >95% intra-species and
      and <83% inter-species ANI values.
    * If the maximum value of ANI between strains within a species is greater than *0.05*, the
      median and maximum value will be reported. Strains that cannot be linked by the median
      ANI, e.g., have no similar strains in the species, will be considered as abnormal strains.
    * It may consist of two scenarios:
        1. Wrong species identification
        2. Poor assembly quality

* Non-redundant strains
    * If the ANI value between two strains within a species is less than *0.005*, the two strains
      are considered to be redundant.
    * Need these files:  representative.lst and omit.lst

* MinHash tree
    * A rough tree is generated by k-mean clustering.

* These abnormal strains should be manually checked to determine whether to include them in the
  subsequent steps.

```shell
cd ~/data/Trichoderma

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --mh \
    --parallel 8 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
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
#  80 summary/NR.lst
#  51 summary/redundant.lst

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#13

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Trichoderma/

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --mh \
    --parallel 8 \
    --not-in summary/redundant.lst \
    --height 0.4

bash MinHash/dist.sh

```

### Condense branches in the minhash tree

* This phylo-tree is not really formal/correct, and shouldn't be used to interpret phylogenetic
  relationships
* It is just used to find more abnormal strains

```shell
mkdir -p ~/data/Trichoderma/tree
cd ~/data/Trichoderma/tree

nw_reroot ../MinHash/tree.nwk Sa_cer_S288C |
    nwr order stdin --nd --an \
    > minhash.reroot.newick

nwr pl-condense --map -r species \
    minhash.reroot.newick ../MinHash/species.tsv |
    nwr comment stdin -r "S=" |
    nwr comment stdin -r "member=" |
    nwr comment stdin -r "^\d+$" |
    nwr order stdin --nd --an \
    > minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

nwr tex minhash.condensed.newick --bl -o Trichoderma.minhash.tex

tectonic Trichoderma.minhash.tex

```

## Count valid species and strains

### For *genomic alignments*

```shell
cd ~/data/Trichoderma/

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    rgr md stdin --num

# Can accept N_COUNT
bash Count/lineage.sh 1

cat Count/lineage.count.tsv |
    rgr md stdin --num

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| item    | count |
|---------|------:|
| strain  |   137 |
| species |    38 |
| genus   |     5 |
| family  |     2 |
| order   |     2 |
| class   |     2 |

| genus         | #species | #strains |
|---------------|---------:|---------:|
| Cladobotryum  |        2 |        3 |
| Escovopsis    |        2 |        7 |
| Hypomyces     |        2 |        2 |
| Saccharomyces |        1 |        1 |
| Trichoderma   |       31 |      124 |

| #family            | genus         | species                     | count |
|--------------------|---------------|-----------------------------|------:|
| Hypocreaceae       | Cladobotryum  | Cladobotryum mycophilum     |     2 |
|                    |               | Cladobotryum protrusum      |     1 |
|                    | Escovopsis    | Escovopsis sp.              |     5 |
|                    |               | Escovopsis weberi           |     2 |
|                    | Hypomyces     | Hypomyces perniciosus       |     1 |
|                    |               | Hypomyces rosellus          |     1 |
|                    | Trichoderma   | Trichoderma afroharzianum   |     5 |
|                    |               | Trichoderma aggressivum     |     1 |
|                    |               | Trichoderma arundinaceum    |     4 |
|                    |               | Trichoderma asperelloides   |     3 |
|                    |               | Trichoderma asperellum      |    18 |
|                    |               | Trichoderma atrobrunneum    |     1 |
|                    |               | Trichoderma atroviride      |    10 |
|                    |               | Trichoderma breve           |     1 |
|                    |               | Trichoderma brevicrassum    |     1 |
|                    |               | Trichoderma citrinoviride   |     4 |
|                    |               | Trichoderma cornu-damae     |     1 |
|                    |               | Trichoderma endophyticum    |     4 |
|                    |               | Trichoderma erinaceum       |     2 |
|                    |               | Trichoderma gamsii          |     2 |
|                    |               | Trichoderma gracile         |     1 |
|                    |               | Trichoderma guizhouense     |     1 |
|                    |               | Trichoderma hamatum         |     1 |
|                    |               | Trichoderma harzianum       |    10 |
|                    |               | Trichoderma koningii        |     1 |
|                    |               | Trichoderma koningiopsis    |     5 |
|                    |               | Trichoderma lentiforme      |     1 |
|                    |               | Trichoderma lixii           |     1 |
|                    |               | Trichoderma longibrachiatum |     7 |
|                    |               | Trichoderma orchidacearum   |     1 |
|                    |               | Trichoderma polysporum      |     1 |
|                    |               | Trichoderma reesei          |    19 |
|                    |               | Trichoderma semiorbis       |     1 |
|                    |               | Trichoderma simmonsii       |     1 |
|                    |               | Trichoderma sp.             |     4 |
|                    |               | Trichoderma virens          |     9 |
|                    |               | Trichoderma viride          |     3 |
| Saccharomycetaceae | Saccharomyces | Saccharomyces cerevisiae    |     1 |

### For *protein families*

```shell
cd ~/data/Trichoderma/

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    rgr md stdin --num

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv

```

| item    | count |
|---------|------:|
| strain  |    35 |
| species |    20 |
| genus   |     4 |
| family  |     2 |
| order   |     2 |
| class   |     2 |

| genus         | #species | #strains |
|---------------|---------:|---------:|
| Cladobotryum  |        1 |        1 |
| Escovopsis    |        1 |        1 |
| Saccharomyces |        1 |        1 |
| Trichoderma   |       17 |       32 |

## Collect proteins

```shell
cd ~/data/Trichoderma/

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --pro \
    --in ASSEMBLY/pass.lst \
    --not-in ASSEMBLY/omit.lst

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
    tsv-summarize -H --count --sum 2-7 |
    sed 's/^count/species/' |
    datamash transpose |
    (echo -e "#item\tcount" && cat) |
    rgr md stdin --fmt

```

| #item      |   count |
|------------|--------:|
| species    |      21 |
| strain_sum |      39 |
| total_sum  | 363,402 |
| dedup_sum  | 363,402 |
| rep_sum    | 284,914 |
| fam88_sum  | 258,055 |
| fam38_sum  | 220,021 |

## Phylogenetics with fungi61

```shell
cd ~/data/Trichoderma/

mkdir -p HMM

# The Fungi HMM set
tar xvfz ~/data/HMM/fungi61/fungi61.tar.gz --directory=HMM
cp HMM/fungi61.lst HMM/marker.lst

```

### Find corresponding representative proteins by `hmmsearch`

* 61 fungal marker genes
    * Ref.: https://doi.org/10.1093/nar/gkac894

* The `E_VALUE` was manually adjusted to 1e-20 to reach a balance between sensitivity and
  speciality.

```shell
cd ~/data/Trichoderma

E_VALUE=1e-20

cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 \
    > Protein/species-f.tsv

cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ -s Protein/"${SPECIES}"/fungi61.tsv ]]; then
        continue
    fi
    if [[ ! -f Protein/"${SPECIES}"/rep_seq.fa.gz ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    cat HMM/marker.lst |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 "
            gzip -dcf Protein/${SPECIES}/rep_seq.fa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/hmm/{}.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), q({}), \$1; '
        " \
        > Protein/${SPECIES}/fungi61.tsv
done

fd --full-path "Protein/.+/fungi61.tsv" -X cat |
    tsv-summarize --group-by 1 --count |
    tsv-summarize --quantile 2:0.25,0.5,0.75
#23      25      55

# There are 461 species and 616 strains
fd --full-path "Protein/.+/fungi61.tsv" -X cat |
    tsv-summarize --group-by 1 --count |
    tsv-filter --invert --ge 2:20 --le 2:30 |
    cut -f 1 \
    > Protein/marker.omit.lst

wc -l HMM/marker.lst Protein/marker.omit.lst
# 61 HMM/marker.lst
# 27 Protein/marker.omit.lst

cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ ! -s Protein/"${SPECIES}"/fungi61.tsv ]]; then
        continue
    fi
    if [[ ! -f Protein/"${SPECIES}"/seq.sqlite ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    # single copy
    cat Protein/"${SPECIES}"/fungi61.tsv |
        grep -v -Fw -f Protein/marker.omit.lst \
        > Protein/"${SPECIES}"/fungi61.sc.tsv

    nwr seqdb -d Protein/${SPECIES} --rep f3=Protein/${SPECIES}/fungi61.sc.tsv

done

```

### Domain related protein sequences

```shell
cd ~/data/Trichoderma

mkdir -p Domain

# each assembly
cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
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
    hnsm dedup stdin |
    hnsm gz stdin -o Domain/fungi61.fa

fd --full-path "Protein/.+/seq_asm_f3.tsv" -X cat \
    > Domain/seq_asm_f3.tsv

cat Domain/seq_asm_f3.tsv |
    tsv-join -e -d 2 -f summary/redundant.lst -k 1 \
    > Domain/seq_asm_f3.NR.tsv

```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Trichoderma

# Extract proteins
cat HMM/marker.lst |
    grep -v -Fw -f Protein/marker.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

        mkdir -p Domain/{}

        hnsm some Domain/fungi61.fa.gz <(
            cat Domain/seq_asm_f3.tsv |
                tsv-filter --str-eq "3:{}" |
                tsv-select -f 1 |
                rgr dedup stdin
            ) \
            > Domain/{}/{}.pro.fa
    '

# Align each marker
cat HMM/marker.lst |
    grep -v -Fw -f Protein/marker.omit.lst |
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
    grep -v -Fw -f Protein/marker.omit.lst |
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
        tsv-filter --str-eq "3:${marker}" |
        tsv-select -f 1-2 |
        hnsm replace -s Domain/${marker}/${marker}.aln.fa stdin \
        > Domain/${marker}/${marker}.replace.fa
done

# Concat marker genes
cat HMM/marker.lst |
    grep -v -Fw -f Protein/marker.omit.lst |
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
    > Domain/fungi61.aln.fas

cat Domain/seq_asm_f3.NR.tsv |
    cut -f 2 |
    rgr dedup stdin |
    sort |
    fasops concat Domain/fungi61.aln.fas stdin -o Domain/fungi61.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Domain/fungi61.aln.fa -out Domain/fungi61.trim.fa -automated1

hnsm size Domain/fungi61.*.fa |
    rgr dedup stdin -f 2 |
    cut -f 2
#31562
#20070

# To make it faster
FastTree -fastest -noml Domain/fungi61.trim.fa > Domain/fungi61.trim.newick

```

### The protein tree

```shell
cd ~/data/Trichoderma/tree

nw_reroot  ../Domain/fungi61.trim.newick Sa_cer_S288C |
    nwr order stdin --nd --an \
    > fungi61.reroot.newick

nwr pl-condense --map -r species \
    fungi61.reroot.newick ../Count/species.tsv |
    nwr order stdin --nd --an \
    > fungi61.condensed.newick

mv condensed.tsv fungi61.condense.tsv

nwr topo --bl fungi61.condensed.newick | # remove comments
    nwr tex stdin --bl -o Trichoderma.fungi61.tex

tectonic Trichoderma.fungi61.tex

```

## Groups and targets

Grouping criteria:

* The mash tree and the marker protein tree
* `MinHash/groups.tsv`

Target selecting criteria:

* `ASSEMBLY/collect.pass.tsv`
* Prefer Sanger sequenced assemblies
* RefSeq_category with `Representative Genome`
* Assembly_level with `Complete Genome` or `Chromosome`

Create a Bash `ARRAY` manually with a format of `group::target`.

```shell
mkdir -p ~/data/Trichoderma/taxon
cd ~/data/Trichoderma/taxon

cat ../ASSEMBLY/collect.pass.tsv |
    tsv-filter -H --str-eq annotations:Yes --le C:100 |
    tsv-select -H -f name,Assembly_level,Genome_coverage,Sequencing_technology,N50,C \
    > potential-target.tsv

cat ../ASSEMBLY/collect.pass.tsv |
    tsv-filter -H --or \
        --str-eq Assembly_level:"Complete Genome" \
        --str-eq Assembly_level:"Chromosome" \
        --le C:50 |
    tsv-select -H -f name,Assembly_level,Genome_coverage,Sequencing_technology,N50,C \
    > complete-genome.tsv

echo -e "#Serial\tGroup\tTarget\tCount" > group_target.tsv

# groups according `groups.tsv`
ARRAY=(
    'C_E_H::E_web_GCA_001278495_1' # 1
    'T_afr_har::T_har_CGMCC_20739_GCA_019097725_1' # 3
    'T_asperello_asperellum::T_asperellum_FT101_GCA_020647865_1' # 5
    'T_atrov_koningio::T_atrov_IMI_206040_GCF_000171015_1' # 6
    'T_cit_lon_ree::T_ree_QM6a_GCF_000167675_1' # 7
    'T_vire::T_vire_Gv29_8_GCF_000170995_1' # 8
)

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"

    SERIAL=$(
        cat ../MinHash/groups.tsv |
            tsv-filter --str-eq 2:${TARGET_NAME} |
            tsv-select -f 1
    )

    cat ../MinHash/groups.tsv |
        tsv-filter --str-eq 1:${SERIAL} |
        tsv-select -f 2 |
        tsv-join -f ../ASSEMBLY/url.tsv -k 1 -a 3 \
        > ${GROUP_NAME}

    COUNT=$(cat ${GROUP_NAME} | wc -l )

    echo -e "${SERIAL}\t${GROUP_NAME}\t${TARGET_NAME}\t${COUNT}" >> group_target.tsv

done

# Custom groups
ARRAY=(
    'Trichoderma::T_ree_QM6a_GCF_000167675_1'
    'Trichoderma_reesei::T_ree_QM6a_GCF_000167675_1'
    'Trichoderma_asperellum::T_asperellum_FT101_GCA_020647865_1'
    'Trichoderma_harzianum::T_har_CGMCC_20739_GCA_019097725_1'
    'Trichoderma_atroviride::T_atrov_IMI_206040_GCF_000171015_1'
    'Trichoderma_virens::T_vire_Gv29_8_GCF_000170995_1'
)

SERIAL=100
for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"

    SERIAL=$((SERIAL + 1))
    GROUP_NAME_2=$(echo $GROUP_NAME | tr "_" " ")

    if [ "$GROUP_NAME" = "Trichoderma" ]; then
        cat ../ASSEMBLY/collect.pass.tsv |
            tsv-filter -H --not-blank RefSeq_category |
            sed '1d' |
            tsv-select -f 1 \
            > T.tmp
        echo "C_pro_CCMJ2080_GCA_004303015_1" >> T.tmp
        echo "E_web_EWB_GCA_003055145_1" >> T.tmp
        echo "E_web_GCA_001278495_1" >> T.tmp
        echo "H_perniciosus_HP10_GCA_008477525_1" >> T.tmp
        echo "H_ros_CCMJ2808_GCA_011799845_1" >> T.tmp
        cat T.tmp |
            rgr dedup stdin |
            tsv-join -f ../ASSEMBLY/url.tsv -k 1 -a 3 \
            > ${GROUP_NAME}

    else
        cat ../ASSEMBLY/collect.pass.tsv |
            tsv-select -f 1,2 |
            grep "${GROUP_NAME_2}" |
            tsv-select -f 1 |
            tsv-join -f ../ASSEMBLY/url.tsv -k 1 -a 3 \
            > ${GROUP_NAME}
    fi

    COUNT=$(cat ${GROUP_NAME} | wc -l )

    echo -e "${SERIAL}\t${GROUP_NAME}\t${TARGET_NAME}\t${COUNT}" >> group_target.tsv

done

cat group_target.tsv |
    rar md stdin --right 4

```

| #Serial | Group                  | Target                             | Count |
|---------|------------------------|------------------------------------|------:|
| 1       | C_E_H                  | E_web_GCA_001278495_1              |     5 |
| 3       | T_afr_har              | T_har_CGMCC_20739_GCA_019097725_1  |    20 |
| 5       | T_asperello_asperellum | T_asperellum_FT101_GCA_020647865_1 |    16 |
| 6       | T_atrov_koningio       | T_atrov_IMI_206040_GCF_000171015_1 |    15 |
| 7       | T_cit_lon_ree          | T_ree_QM6a_GCF_000167675_1         |    25 |
| 8       | T_vire                 | T_vire_Gv29_8_GCF_000170995_1      |     9 |
| 101     | Trichoderma            | T_ree_QM6a_GCF_000167675_1         |    32 |
| 102     | Trichoderma_reesei     | T_ree_QM6a_GCF_000167675_1         |    13 |
| 103     | Trichoderma_asperellum | T_asperellum_FT101_GCA_020647865_1 |    13 |
| 104     | Trichoderma_harzianum  | T_har_CGMCC_20739_GCA_019097725_1  |    10 |
| 105     | Trichoderma_atroviride | T_atrov_IMI_206040_GCF_000171015_1 |     7 |
| 106     | Trichoderma_virens     | T_vire_Gv29_8_GCF_000170995_1      |     8 |

## Prepare sequences for `egaz`

* `--perseq` for Chromosome-level assemblies and targets
    * means split fasta by names, targets or good assembles should set it

```shell
cd ~/data/Trichoderma

# /share/home/wangq/homebrew/Cellar/repeatmasker@4.1.1/4.1.1/libexec/famdb.py \
#   -i /share/home/wangq/homebrew/Cellar/repeatmasker@4.1.1/4.1.1/libexec/Libraries/RepeatMaskerLib.h5 \
#   lineage Fungi

# prep
egaz template \
    ASSEMBLY \
    --prep -o Genome \
    $( cat taxon/group_target.tsv |
        sed -e '1d' | cut -f 3 |
        parallel -j 1 echo " --perseq {} "
    ) \
    $( cat taxon/complete-genome.tsv |
        sed '1d' | cut -f 1 |
        parallel -j 1 echo " --perseq {} "
    ) \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--parallel 16"

bash Genome/0_prep.sh

# gff
for n in \
    $(cat taxon/group_target.tsv | sed -e '1d' | cut -f 3 ) \
    $( cat taxon/potential-target.tsv | sed -e '1d' | cut -f 1 ) \
    ; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"

    gzip -dc ${FILE_GFF} > Genome/${n}/chr.gff
done

```

## Generate alignments

```shell
cd ~/data/Trichoderma

cat taxon/group_target.tsv |
    sed -e '1d' |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{3}]\n"

        egaz template \
            Genome/{3} \
            $(cat taxon/{2} | cut -f 1 | grep -v -x "{3}" | xargs -I[] echo "Genome/[]") \
            --multi -o groups/{2}/ \
            --tree MinHash/tree.nwk \
            --parallel 16 -v

        bash groups/{2}/1_pair.sh
        bash groups/{2}/3_multi.sh
    '

# clean
find groups -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find groups -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type f -name "output.*" | parallel -r rm

```
