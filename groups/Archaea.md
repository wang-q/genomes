# Archaea

All genomes of *Archaea*

<!-- TOC -->
* [Archaea](#archaea)
  * [Strain info](#strain-info)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
  * [Download all assemblies](#download-all-assemblies)
    * [Create assembly.tsv](#create-assemblytsv)
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
  * [Phylogenetics with ar53](#phylogenetics-with-ar53)
    * [Find corresponding representative proteins by `hmmsearch`](#find-corresponding-representative-proteins-by-hmmsearch)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
    * [Condense branches in the protein tree](#condense-branches-in-the-protein-tree)
<!-- TOC -->

## Strain info

* [Archaea](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2157)

### List all ranks

```shell
nwr member Archaea |
    grep -v " sp." |
    grep -v " x " |
    tsv-summarize -H -g rank --count |
    rgr md stdin --num

```

| rank          | count |
|---------------|------:|
| superkingdom  |     1 |
| kingdom       |     4 |
| phylum        |    46 |
| class         |    47 |
| no rank       |   405 |
| species       |  3639 |
| order         |    79 |
| genus         |   311 |
| family        |   111 |
| clade         |    40 |
| strain        |   352 |
| species group |     2 |
| isolate       |     6 |

### Species with assemblies

* 'RefSeq'
    * RS1 - genome_rep: 'Full'
* 'Genbank'
    * GB1 - genome_rep: 'Full'

```shell
mkdir -p ~/data/Archaea/summary
cd ~/data/Archaea/summary

# should have a valid name of genus
nwr member Archaea -r genus |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#311 genus.list

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
#  726 RS1.tsv
#  890 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     1508
#GB1     2274

```

## Download all assemblies

### Create assembly.tsv

If a refseq assembly is available, the corresponding genbank one is not downloaded

```shell
cd ~/data/Archaea/summary

echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus IN ('Saccharomyces')
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
#2280 lines, 7 fields

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
    > Archaea.assembly.tsv

datamash check < Archaea.assembly.tsv
#2277 lines, 5 fields

# find potential duplicate strains or assemblies
cat Archaea.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Archaea.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

# Edit .assembly.tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Archaea.assembly.tsv
#
# Save the file to another directory to prevent accidentally changing it
# cp Archaea.assembly.tsv ~/Scripts/genomes/assembly

# Cleaning
rm raw*.*sv

```

### Count before download

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class

```shell
cd ~/data/Archaea

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
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
| strain  |  2276 |
| species |   894 |
| genus   |   293 |
| family  |   108 |
| order   |    74 |
| class   |    43 |

| genus                | #species | #strains |
|----------------------|---------:|---------:|
| Ferroplasma          |        2 |       20 |
| Haloarcula           |       32 |       59 |
| Halobacterium        |        9 |       46 |
| Haloferax            |       17 |       43 |
| Halorubrum           |       44 |      112 |
| Halorussus           |       17 |       21 |
| Metallosphaera       |        8 |       21 |
| Methanobacterium     |       14 |       38 |
| Methanobrevibacter   |       14 |      206 |
| Methanococcus        |        4 |       38 |
| Methanoculleus       |       16 |       47 |
| Methanomethylophilus |        1 |       25 |
| Methanosarcina       |       11 |      120 |
| Methanothermobacter  |        5 |       23 |
| Methanothrix         |        3 |       69 |
| Natrinema            |       23 |       37 |
| Nitrosopumilus       |       11 |       24 |
| Saccharolobus        |        3 |       20 |
| Sulfolobus           |        3 |       85 |
| Thermococcus         |       35 |       50 |

### Download and check

```shell
cd ~/data/Archaea

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
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

| #item            | fields | lines |
|------------------|-------:|------:|
| url.tsv          |      3 | 2,276 |
| check.lst        |      1 | 2,276 |
| collect.tsv      |     20 | 2,277 |
| n50.tsv          |      4 | 2,277 |
| n50.pass.tsv     |      4 | 1,407 |
| collect.pass.tsv |     23 | 1,407 |
| pass.lst         |      1 | 1,406 |
| omit.lst         |      1 |   343 |
| rep.lst          |      1 |   660 |
| sp.lst           |      0 |     0 |

### Rsync to hpcc

```shell
rsync -avP \
    ~/data/Archaea/ \
    wangq@202.119.37.251:data/Archaea

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Archaea/ \
    wangq@58.213.64.36:data/Archaea

# rsync -avP wangq@202.119.37.251:data/Archaea/ ~/data/Archaea

# rsync -avP -e 'ssh -p 8804' wangq@58.213.64.36:data/Archaea/ ~/data/Archaea

```

## BioSample

```shell
cd ~/data/Archaea

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --bs

# Run this script twice and it will re-download the failed files
bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 50

datamash check < BioSample/biosample.tsv
#2271 lines, 56 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## MinHash

```shell
cd ~/data/Archaea

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --mh \
    --parallel 8 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
    --ani-nr 0.005

# Compute assembly sketches
bash MinHash/compute.sh

# Distances within species
bash MinHash/species.sh

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#68

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
#  519 summary/NR.lst
#  365 summary/redundant.lst

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Archaea/

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --mh \
    --parallel 8 \
    --not-in MinHash/abnormal.lst \
    --not-in summary/redundant.lst \
    --height 0.4

bash MinHash/dist.sh

```

### Condense branches in the minhash tree

```shell
mkdir -p ~/data/Archaea/tree
cd ~/data/Archaea/tree

nw_reroot ../MinHash/tree.nwk Saccharom_cere_S288C Saccharom_eub |
    nwr order stdin --nd --an \
    > minhash.reroot.newick

nwr pl-condense --map -r order -r family -r genus \
    minhash.reroot.newick ../MinHash/species.tsv |
    nwr order stdin --nd --an \
    -o minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 minhash.condensed.newick |
    rsvg-convert -o Archaea.minhash.png

```

## Count valid species and strains

### For *genomic alignments*

```shell
cd ~/data/Archaea/

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --rank order --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/order.count.tsv |
    tsv-filter -H --ge "3:20" |
    rgr md stdin --num

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:20" |
    rgr md stdin --num

# Can accept N_COUNT
bash Count/lineage.sh 5

cat Count/lineage.count.tsv |
    rgr md stdin --num

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| item    | count |
|:--------|------:|
| strain  |  1338 |
| species |   739 |
| genus   |   235 |
| family  |    78 |
| order   |    52 |
| class   |    28 |

| order              | #species | #strains |
|:-------------------|---------:|---------:|
| Desulfurococcales  |       20 |       23 |
| Halobacteriales    |      363 |      594 |
| Methanobacteriales |       30 |      117 |
| Methanococcales    |       15 |       45 |
| Methanomicrobiales |       48 |       55 |
| Methanosarcinales  |       45 |       86 |
| Sulfolobales       |       28 |      136 |
| Thermococcales     |       41 |       61 |

| genus              | #species | #strains |
|:-------------------|---------:|---------:|
| Haloarcula         |       32 |       56 |
| Halobacterium      |        8 |       28 |
| Haloferax          |       16 |       41 |
| Halorubrum         |       32 |       83 |
| Halorussus         |       16 |       20 |
| Methanobacterium   |       14 |       25 |
| Methanobrevibacter |        8 |       73 |
| Methanococcus      |        4 |       33 |
| Methanoculleus     |       16 |       21 |
| Methanosarcina     |       10 |       37 |
| Natrinema          |       21 |       29 |
| Sulfolobus         |        2 |       76 |
| Thermococcus       |       34 |       41 |

| #family                         | genus                         | species                                       | count |
|:--------------------------------|:------------------------------|:----------------------------------------------|------:|
| Candidatus Haiyanarchaeaceae    | Candidatus Haiyanarchaeum     | Candidatus Haiyanarchaeum thermophilum        |     5 |
| Candidatus Methanomethylicaceae | Candidatus Methanosuratincola | Candidatus Methanosuratincola petrocarbonis   |    10 |
| Candidatus Wolframiiraptoraceae | Candidatus Terraquivivens     | Candidatus Terraquivivens tengchongensis      |     7 |
| Haloarculaceae                  | Haloarcula                    | Haloarcula argentinensis                      |     5 |
|                                 |                               | Haloarcula hispanica                          |     8 |
|                                 | Halocatena                    | Halocatena marina                             |     6 |
| Halobacteriaceae                | Halobacterium                 | Halobacterium hubeiense                       |     8 |
|                                 |                               | Halobacterium salinarum                       |    12 |
| Haloferacaceae                  | Haloferax                     | Haloferax mediterranei                        |     5 |
|                                 |                               | Haloferax volcanii                            |    15 |
|                                 | Halogeometricum               | Halogeometricum borinquense                   |     6 |
|                                 | Halorubrum                    | Halorubrum ezzemoulense                       |    38 |
|                                 |                               | Halorubrum lacusprofundi                      |     5 |
| Methanobacteriaceae             | Methanobacterium              | Methanobacterium formicicum                   |     6 |
|                                 | Methanobrevibacter            | Methanobrevibacter smithii                    |    54 |
|                                 |                               | Methanobrevibacter woesei                     |     7 |
|                                 | Methanothermobacter           | Methanothermobacter wolfeii                   |     5 |
| Methanococcaceae                | Methanococcus                 | Methanococcus maripaludis                     |    27 |
| Methanomassiliicoccaceae        | Methanomassiliicoccus         | Candidatus Methanomassiliicoccus intestinalis |     5 |
| Methanomethylophilaceae         | Methanomethylophilus          | Methanomethylophilus alvi                     |     7 |
| Methanomicrobiaceae             | Methanoculleus                | Methanoculleus bourgensis                     |     5 |
| Methanosarcinaceae              | Methanosarcina                | Methanosarcina mazei                          |    13 |
|                                 |                               | Methanosarcina thermophila                    |    12 |
| Sulfolobaceae                   | Metallosphaera                | Metallosphaera sedula                         |    11 |
|                                 | Saccharolobus                 | Saccharolobus solfataricus                    |    14 |
|                                 | Sulfolobus                    | Sulfolobus acidocaldarius                     |    55 |
|                                 |                               | Sulfolobus islandicus                         |    21 |
| Thermococcaceae                 | Pyrococcus                    | Pyrococcus kukulkanii                         |    11 |

### For *protein families*

```shell
cd ~/data/Archaea/

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
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
    tsv-filter -H --ge "3:10" |
    rgr md stdin --num

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv

```

| item    | count |
|:--------|------:|
| strain  |  1279 |
| species |   720 |
| genus   |   223 |
| family  |    74 |
| order   |    48 |
| class   |    26 |

| genus               | #species | #strains |
|:--------------------|---------:|---------:|
| Haladaptatus        |        7 |       11 |
| Haloarcula          |       32 |       55 |
| Halobacterium       |        8 |       28 |
| Halobaculum         |       10 |       16 |
| Halobellus          |        8 |       13 |
| Haloferax           |       16 |       41 |
| Halomarina          |        7 |       12 |
| Haloplanus          |        9 |       14 |
| Halorientalis       |        7 |       10 |
| Halorubrum          |       32 |       83 |
| Halorussus          |       16 |       20 |
| Metallosphaera      |        7 |       18 |
| Methanobacterium    |       14 |       24 |
| Methanobrevibacter  |        8 |       66 |
| Methanococcus       |        4 |       31 |
| Methanoculleus      |       16 |       21 |
| Methanohalophilus   |        6 |       14 |
| Methanosarcina      |       10 |       35 |
| Methanothermobacter |        4 |       11 |
| Natrinema           |       21 |       29 |
| Natronorubrum       |        9 |       13 |
| Nitrosopumilus      |       10 |       12 |
| Pyrococcus          |        5 |       18 |
| Saccharolobus       |        3 |       19 |
| Sulfolobus          |        2 |       75 |
| Thermococcus        |       34 |       41 |

## Collect proteins

```shell
cd ~/data/Archaea/

nwr template ~/Scripts/genomes/assembly/Archaea.assembly.tsv \
    --pro \
    --parallel 8 \
    --in ASSEMBLY/pass.lst \
    --not-in ASSEMBLY/omit.lst

# collect proteins and
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
    perl -nla -F"\t" -MNumber::Format -e '
        printf qq(%s\t%s\n), $F[0], Number::Format::format_number($F[1], 0,);
        ' |
    (echo -e "#item\tcount" && cat) |
    rgr md stdin -r 2

```

| #item      |     count |
|------------|----------:|
| species    |       735 |
| strain_sum |     1,342 |
| total_sum  | 3,857,410 |
| dedup_sum  | 2,733,465 |
| rep_sum    | 2,391,856 |
| fam88_sum  | 2,297,232 |
| fam38_sum  | 1,934,238 |

## Phylogenetics with ar53

### Find corresponding representative proteins by `hmmsearch`

```shell
cd ~/data/Archaea

# The Archaea61 HMM set
nwr kb ar53 -o HMM
cp HMM/ar53.lst HMM/marker.lst

mkdir -p Domain

cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 \
    > Protein/species-f.tsv

cat Protein/species-f.tsv |
    tsv-select -f 2 |
    tsv-uniq |
while read SPECIES; do
    if [[ -s Protein/"${SPECIES}"/ar53.tsv ]]; then
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
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\t%s\n), q({}), \$1; '
        " \
        > Protein/${SPECIES}/ar53.tsv
done

cat Protein/species-f.tsv |
    tsv-select -f 2 |
    tsv-uniq |
while read SPECIES; do
    if [[ ! -s Protein/"${SPECIES}"/ar53.tsv ]]; then
        continue
    fi
    if [[ ! -f Protein/"${SPECIES}"/seq.sqlite ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    nwr seqdb -d Protein/${SPECIES} --rep f3=Protein/${SPECIES}/ar53.tsv

done

# each assembly
cat Protein/species-f.tsv |
    tsv-select -f 2 |
    tsv-uniq |
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
                tsv-uniq
        )
done |
    hnsm dedup stdin |
    hnsm gz stdin -o Domain/ar53.fa

fd --full-path "Protein/.*/seq_asm_f3.tsv" -X cat \
    > Domain/seq_asm_f3.tsv

cat Domain/seq_asm_f3.tsv |
    tsv-join -e -d 2 -f summary/NR.lst -k 1 \
    > Domain/seq_asm_f3.NR.tsv

```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Archaea

# Extract proteins
cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

        mkdir -p Domain/{}

        hnsm some Domain/ar53.fa.gz <(
            cat Domain/seq_asm_f3.tsv |
                tsv-filter --str-eq "3:{}" |
                tsv-select -f 1 |
                tsv-uniq
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
        tsv-filter --str-eq "3:${marker}" |
        tsv-select -f 1-2 |
        hnsm replace -s Domain/${marker}/${marker}.aln.fa stdin \
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
    > Domain/ar53.aln.fas

cat Domain/seq_asm_f3.NR.tsv |
    cut -f 2 |
    tsv-uniq |
    sort |
    fasops concat Domain/ar53.aln.fas stdin -o Domain/ar53.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Domain/ar53.aln.fa -out Domain/ar53.trim.fa -automated1

hnsm size Domain/ar53.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#45929
#12258

# To make it faster
FastTree -fastest -noml Domain/ar53.trim.fa > Domain/ar53.trim.newick

```

### Condense branches in the protein tree

```shell
cd ~/data/Archaea/tree

nw_reroot  ../Domain/ar53.trim.newick Saccharom_cere_S288C Saccharom_eub |
    nwr order stdin --nd --an \
    > ar53.reroot.newick

nwr pl-condense --map -r order -r family -r genus \
    ar53.reroot.newick ../Count/species.tsv |
    nwr order stdin --nd --an \
    -o ar53.condensed.newick

mv condensed.tsv ar53.condense.tsv

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 ar53.condensed.newick |
    rsvg-convert -o Archaea.ar53.png

```

