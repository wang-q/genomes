# Fungi

Download all genomes and analyze representative strains.

[TOC levels=2-4]: #

- [Taxon info](#taxon-info)
  - [List all ranks](#list-all-ranks)
  - [Species with assemblies](#species-with-assemblies)
  - [Model organisms](#model-organisms)
- [Download all assemblies](#download-all-assemblies)
  - [Create assembly.tsv](#create-assemblytsv)
  - [Count before download](#count-before-download)
  - [Download and check](#download-and-check)
  - [Rsync to hpcc](#rsync-to-hpcc)
- [BioSample](#biosample)
- [Divergence of Fungi](#divergence-of-fungi)
  - [ReRoot](#reroot)
- [MinHash](#minhash)
  - [Condense branches in the minhash tree](#condense-branches-in-the-minhash-tree)
- [Count valid species and strains](#count-valid-species-and-strains)
  - [For genomic alignments](#for-genomic-alignments)
  - [For protein families](#for-protein-families)
- [Collect proteins](#collect-proteins)
- [Phylogenetics with fungi61](#phylogenetics-with-fungi61)
  - [Find corresponding proteins by ](#find-corresponding-proteins-by-)
  - [Find corresponding representative proteins by ](#find-corresponding-representative-proteins-by-)
  - [Domain related protein sequences](#domain-related-protein-sequences)
  - [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
  - [Condense branches in the protein tree](#condense-branches-in-the-protein-tree)
- [InterProScan on all proteins of representative and typical strains](#interproscan-on-all-proteins-of-representative-and-typical-strains)
- [Groups and targets](#groups-and-targets)

## Taxon info

- [Fungi](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4751)

Open nomenclature (开放命名法)

- Sp. (pl. spp.; short for "species") indicates potentially new species
- Sp. aff. or aff. (short for "species affinis") indicates a potentially new and undescribed
  species has an affinity to, but is not identical to, the named species
- Cf. (short for the Latin: confer, "compare with") or a question mark (?, also inc., species
  incerta) signify varying degrees or types of uncertainty
- Candidatus, a candidate taxon proposed from metagenomics or other incomplete information
- Incertae sedis, a taxon of uncertain position in a classification
- Nomen invalidum, a name that is not validly published

### List all ranks

```shell
nwr member Fungi |
    grep -v " sp." |
    grep -v " aff." |
    grep -v " cf." |
    grep -v " x " |
    tva stats -H -g rank --count |
    tva to md --num

nwr member Fungi | grep " sp\."
nwr member Fungi | grep " aff\."
nwr member Fungi | grep " cf\."
nwr member Fungi -r species | grep "\." | grep -v " sp\." | grep -v " aff\." | grep -v " cf\."

nwr member Fungi | grep " x "
```

| rank          | count |
| ------------- | ----: |
| clade         |    27 |
| class         |    71 |
| family        |   985 |
| forma         |   255 |
| genus         |  8095 |
| isolate       |     5 |
| kingdom       |     1 |
| morph         |     2 |
| no rank       |  5768 |
| order         |   273 |
| phylum        |    10 |
| section       |    39 |
| species       | 66685 |
| species group |     1 |
| strain        |  2228 |
| subclass      |    18 |
| subfamily     |    19 |
| subgenus      |    10 |
| subkingdom    |     1 |
| suborder      |    25 |
| subphylum     |    14 |
| subspecies    |   143 |
| superfamily   |     1 |
| tribe         |     3 |
| varietas      |  1115 |

### Species with assemblies

In the vast majority of fungal species, only one genome was selected for refseq.

- 'RefSeq'
- 'Genbank'

```shell
mkdir -p ~/data/Fungi/summary
cd ~/data/Fungi/summary

# should have a valid name of genus
nwr member Fungi -r genus |
    grep -v " x " |
    sed '1d' |
    tva sort -n -k 1 \
    > genus.list.tsv

wc -l genus.list.tsv
#8095 genus.list

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
#    661 RS1.tsv
#   5448 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tva stats --sum 3
        fi
    done
done
# RS1     705
# GB1     22016
```

- The names of some genera are abnormal

```shell
cd ~/data/Fungi/summary

cat RS*.tsv GB*.tsv |
    tva select -f 1,2 |
    tva uniq |
    grep "\[" |
    nwr append stdin -r genus |
    head
# 561895  [Candida] subhashii     Spathaspora
# 566037  [Ashbya] aceris (nom. inval.)   Eremothecium
# 45518   [Candida] aaseri        Yamadazyma
# 1171601 [Candida] adriatica     Cyberlindnera
# 456249  [Candida] andamanensis (nom. inval.)    Yamadazyma
# 148631  [Candida] anglica       Kurtzmaniella
# 130810  [Candida] arabinofermentans     Ogataea
# 391823  [Candida] ascalaphidarum        Yamadazyma
# 224031  [Candida] aurita        Kurtzmaniella
# 551772  [Candida] awuae Pichia
```

### Model organisms

There is one true model organism in Fungi, Saccharomyces cerevisiae S288C

```shell
cd ~/data/Fungi/summary

GENUS=$(
    cat genus.list.tsv |
        tva select -f 1 |
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
        AND genus_id IN ($GENUS)
    GROUP BY refseq_category
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
# refseq_category count
# na      5
# reference genome        700
```

## Download all assemblies

### Create assembly.tsv

If a refseq assembly is available, the corresponding genbank one is not listed

Group all questionable strains (sp., cf., aff.) under `Genus sp.`.

```shell
cd ~/data/Fungi/summary

# RS
SPECIES=$(
    cat RS1.tsv |
        tva select -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
.header ON
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND genome_rep IN ('Full')
        AND species NOT LIKE '% x %'
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > raw.tsv

# Preference for refseq
cat raw.tsv |
    tva select -H -f "assembly_accession" \
    > rs.acc.tsv

# GB
SPECIES=$(
    cat GB1.tsv |
        tva select -f 1 |
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
        AND genome_rep IN ('Full')
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% aff.%'
        AND species NOT LIKE '% cf.%'
        AND species NOT LIKE '% x %'
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tva join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

echo "
    SELECT
        genus || ' sp. ' || infraspecific_name || ' ' || assembly_accession AS name,
        genus || ' sp.', genus, ftp_path, biosample, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND genome_rep IN ('Full')
        AND (
            species LIKE '% sp.%'
            OR species LIKE '% aff.%'
            OR species LIKE '% cf.%'
        )
        AND species NOT LIKE '% x %'
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tva join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

cat raw.tsv |
    tva uniq |
    tva check
# 21941 lines, 7 fields

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
    > Fungi.assembly.tsv

tva check < Fungi.assembly.tsv
# 21933 lines, 5 fields

# find potential duplicate strains or assemblies
cat Fungi.assembly.tsv |
    tva uniq -f 1 --repeated

cat Fungi.assembly.tsv |
    tva filter --str-not-in-fld 2:ftp

# Edit .assembly.tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Fungi.assembly.tsv
#
# Save the file to another directory to prevent accidentally changing it
# cp Fungi.assembly.tsv ~/Scripts/genomes/assembly

# Cleaning
rm raw*.*sv
```

### Count before download

- `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class

```shell
cd ~/data/Fungi

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
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

| item    |  count |
| ------- | -----: |
| strain  | 21,925 |
| species |  5,396 |
| genus   |  1,514 |
| family  |    503 |
| order   |    184 |
| class   |     61 |

| genus             | #species | #strains |
| ----------------- | -------: | -------: |
| Alternaria        |       38 |      247 |
| Aspergillus       |      183 |    1,440 |
| Aureobasidium     |       14 |      177 |
| Beauveria         |       15 |      375 |
| Botryosphaeria    |        3 |      142 |
| Candida           |       34 |      384 |
| Candidozyma       |       10 |      705 |
| Clavispora        |        9 |      103 |
| Colletotrichum    |       84 |      328 |
| Cryphonectria     |        7 |      115 |
| Cryptococcus      |       12 |      233 |
| Exophiala         |       13 |      102 |
| Fusarium          |      204 |    1,895 |
| Komagataella      |        7 |      193 |
| Lachancea         |       11 |      131 |
| Macrophomina      |        4 |      491 |
| Malassezia        |       18 |      102 |
| Metschnikowia     |       64 |      151 |
| Neurospora        |       11 |      124 |
| Ogataea           |       56 |      149 |
| Parastagonospora  |        2 |      190 |
| Penicillium       |      117 |      614 |
| Pichia            |       32 |      176 |
| Psilocybe         |       15 |      144 |
| Pyricularia       |        5 |      606 |
| Rhizopus          |        5 |      107 |
| Rhodotorula       |       13 |      206 |
| Saccharomyces     |       10 |    1,924 |
| Scheffersomyces   |       20 |      142 |
| Torulaspora       |       20 |      153 |
| Trichoderma       |       53 |      208 |
| Trichophyton      |       15 |      184 |
| Venturia          |        8 |      107 |
| Verticillium      |       10 |      126 |
| Yamadazyma        |       55 |      104 |
| Zygosaccharomyces |       14 |      159 |

### Download and check

```shell
cd ~/data/Fungi

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
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

find ASSEMBLY/ -name "*_genomic.fna.gz" |
    grep -v "_from_" |
    wc -l
#21930

# N50 C S; create n50.tsv and n50.pass.tsv
bash ASSEMBLY/n50.sh 100000 2000 1000000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tva filter -H --str-in-fld "name:_GCF_" |
    tva stats -H --min "N50" --max "C" --min "S"
#N50_min C_max   S_min
#32179   2016    2187595

cat ASSEMBLY/n50.tsv |
    tva stats -H --quantile "N50:0.1,0.5" --quantile "C:0.5,0.9" --quantile "S:0.1,0.5" |
    tva transpose
# N50_quantile_0.1        24176.5
# N50_quantile_0.5        321544.5
# C_quantile_0.5  350
# C_quantile_0.9  5084.5
# S_quantile_0.1  11665362.4
# S_quantile_0.5  32559599.5

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

| #item            | fields |  lines |
| ---------------- | -----: | -----: |
| url.tsv          |      3 | 21,932 |
| check.lst        |      1 | 21,932 |
| collect.tsv      |     20 | 21,933 |
| n50.pass.tsv     |      4 | 14,971 |
| collect.pass.tsv |     23 | 14,971 |
| pass.lst         |      1 | 14,970 |
| omit.lst         |      1 | 16,573 |
| rep.lst          |      1 |  3,976 |
| sp.lst           |      1 |    216 |

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Fungi/ \
    wangq@202.119.37.251:data/Fungi

rsync -avP \
    /Volumes/data-1/Fungi/ \
    ~/data/Fungi


# rsync -avP wangq@202.119.37.251:data/Fungi/ ~/data/Fungi
```

## BioSample

```shell
cd ~/data/Fungi

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --bs

# Run this script twice and it will re-download the failed files
bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 10

tva check < BioSample/biosample.tsv
#21555 lines, 239 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/
```

## Divergence of Fungi

Ref.:

1. A genome-scale phylogeny of the kingdom Fungi. Cur Biology, 2021. 
   doi.org/10.1016/j.cub.2021.01.074

![1-s2.0-S0960982221001391-fx1_lrg.jpg](images/1-s2.0-S0960982221001391-fx1_lrg.jpg)

### ReRoot

- 小孢子虫 (Microsporidia) 可以当作真菌的基部类群, 也可以当作真菌的姊妹类群
- 芽枝霉门 (Blastocladiomycota) 是真菌的基部类群
- For the latter steps, use the following two as the outgroups
    - Enc_hellem_ATCC_50504_GCF_000277815_2
    - Nematoc_ausu_ERTm6_GCF_000738915_1

```shell
cd ~/data/Fungi

cat summary/collect.pass.tsv |
    tva select -f 1,3 |
    sed '1d' |
    grep -v -Fw -f summary/omit.lst |
    nwr append stdin -c 2 -r species -r phylum -r subkingdom |
    tva filter --str-ne "5:Dikarya" | # 双核亚界
    tva filter --str-ne "4:Mucoromycota" | # 毛霉门
    tva filter --str-ne "4:Zoopagomycota" | #
    tva filter --str-ne "4:Chytridiomycota" | # 壶菌门
    tva select -f 1,3,4 |
    tva sort -k 3,1 |
    tva stats -g 3,2 --count
# Blastocladiomycota      Allomyces arbusculus    1
# Blastocladiomycota      Allomyces javanicus     1
# Blastocladiomycota      Allomyces macrogynus    1
# Blastocladiomycota      Blastocladiella emersonii       1
# Blastocladiomycota      Catenaria anguillulae   1
# Blastocladiomycota      Paraphysoderma sedebokerense    1
# Blastocladiomycota      Sorochytrium milnesiophthora    1
# Microsporidia   Binucleata daphniae     1
# Microsporidia   Ecytonucleospora hepatopenaei   1
# Microsporidia   Edhazardia aedis        1
# Microsporidia   Encephalitozoon cuniculi        3
# Microsporidia   Encephalitozoon hellem  4
# Microsporidia   Encephalitozoon intestinalis    5
# Microsporidia   Encephalitozoon romaleae        1
# Microsporidia   Glugoides intestinalis  1
# Microsporidia   Gurleya vavrai  1
# Microsporidia   Hamiltosporidium tvaerminnensis 1
# Microsporidia   Mitosporidium daphniae  1
# Microsporidia   Nematocida ausubeli     2
# Microsporidia   Nematocida displodere   1
# Microsporidia   Nematocida major        1
# Microsporidia   Nematocida parisii      2
# Microsporidia   Ordospora colligata     5
# Microsporidia   Ordospora pajunii       1
# Microsporidia   Vairimorpha bombi       1
# Microsporidia   Vairimorpha ceranae     1
# Microsporidia   Vairimorpha necatrix    1
# Microsporidia   Vittaforma corneae      1

cat summary/collect.pass.tsv |
    tva filter -H --not-blank RefSeq_category |
    tva filter -H --or \
        --str-in-fld "2:Allomyces" \
        --str-in-fld "2:Blastocladiella" \
        --str-in-fld "2:Encephalitozoon" \
        --str-in-fld "2:Nematocida" |
    grep -v -Fw -f summary/omit.lst |
    tva select -H -f "#name,Assembly_level,Assembly_method,Genome_coverage,Sequencing_technology"
```

## MinHash

```shell
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --mh \
    --parallel 8 \
    --in summary/pass.lst \
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
    # 3757 summary/NR.lst
    # 6946 summary/redundant.lst

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#  203

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --mh \
    --parallel 8 \
    --not-in summary/sp.lst \
    --not-in summary/omit.lst \
    --not-in MinHash/abnormal.lst \
    --not-in summary/redundant.lst \
    --height 0.4

bash MinHash/dist.sh
```

### Condense branches in the minhash tree

```shell
mkdir -p ~/data/Fungi/tree
cd ~/data/Fungi/tree

pgr nwk reroot ../MinHash/tree.nwk -n Enc_hellem_ATCC_50504_GCF_000277815_2 -n Nematoc_ausu_ERTm6_GCF_000738915_1 |
    pgr nwk order stdin --nd --an \
    > minhash.reroot.newick

pgr pl condense --map -t ../Count/strains.taxon.tsv -r 5 -r 4 -r 3 \
    minhash.reroot.newick |
    pgr nwk order stdin --nd --an \
    > minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

# # pdf
# pgr nwk to-tex minhash.condensed.newick --bl |
#     tectonic - &&
#     mv texput.pdf Fungi.minhash.pdf

# svg
pgr nwk topo --bl minhash.condensed.newick | # remove comments
    nw_display -s -b 'visibility:hidden' -w 1200 -v 20 - \
    > Fungi.minhash.svg

```

## Count valid species and strains

### For *genomic alignments*

```shell
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --count \
    --in summary/pass.lst \
    --rank order --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    tva to md --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/order.count.tsv |
    tva filter -H --ge "3:100" |
    tva to md --num

cat Count/genus.count.tsv |
    tva filter -H --ge "3:100" |
    tva to md --num

# Can accept N_COUNT
bash Count/lineage.sh 50

cat Count/lineage.count.tsv |
    tva to md --num

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv
cp Count/genus.count.tsv summary/genus.genome.tsv
```

| item    | count |
| ------- | ----: |
| strain  | 14968 |
| species |  3903 |
| genus   |  1134 |
| family  |   418 |
| order   |   163 |
| class   |    55 |

| order             | #species | #strains |
| ----------------- | -------: | -------: |
| Agaricales        |      221 |      389 |
| Botryosphaeriales |       31 |      257 |
| Chaetothyriales   |       59 |      170 |
| Diaporthales      |       57 |      181 |
| Dipodascales      |      145 |      268 |
| Dothideales       |       20 |      110 |
| Eurotiales        |      327 |     1991 |
| Glomerellales     |       87 |      325 |
| Helotiales        |      101 |      170 |
| Hypocreales       |      373 |     1912 |
| Magnaporthales    |       13 |      363 |
| Mucorales         |       61 |      112 |
| Mycosphaerellales |       75 |      214 |
| Onygenales        |       46 |      267 |
| Phaffomycetales   |       81 |      139 |
| Pichiales         |      134 |      505 |
| Pleosporales      |      174 |      723 |
| Polyporales       |       81 |      130 |
| Saccharomycetales |      129 |     2239 |
| Serinales         |      339 |     1492 |
| Sordariales       |       76 |      160 |
| Sporidiobolales   |       28 |      221 |
| Tremellales       |       47 |      335 |
| Ustilaginales     |       38 |      102 |
| Xylariales        |       82 |      160 |

| genus            | #species | #strains |
| ---------------- | -------: | -------: |
| Alternaria       |       38 |      233 |
| Aspergillus      |      165 |     1240 |
| Aureobasidium    |       14 |      101 |
| Candida          |       24 |      220 |
| Candidozyma      |       10 |      528 |
| Clavispora       |        9 |      100 |
| Colletotrichum   |       70 |      207 |
| Cryptococcus     |       12 |      226 |
| Fusarium         |      136 |     1233 |
| Komagataella     |        7 |      190 |
| Lachancea        |       10 |      128 |
| Ogataea          |       55 |      136 |
| Parastagonospora |        2 |      169 |
| Penicillium      |      114 |      574 |
| Pyricularia      |        4 |      346 |
| Rhodotorula      |       11 |      197 |
| Saccharomyces    |       10 |     1661 |
| Scheffersomyces  |       19 |      115 |
| Torulaspora      |       20 |      140 |
| Trichoderma      |       48 |      191 |
| Trichophyton     |       14 |      111 |
| Verticillium     |       10 |      107 |
| Yamadazyma       |       55 |      103 |

| #family                  | genus               | species                   | count |
| ------------------------ | ------------------- | ------------------------- | ----: |
| Aspergillaceae           | Aspergillus         | Aspergillus flavus        |   287 |
|                          |                     | Aspergillus fumigatus     |   270 |
|                          |                     | Aspergillus niger         |   117 |
|                          |                     | Aspergillus oryzae        |   120 |
|                          | Penicillium         | Penicillium chrysogenum   |    88 |
|                          |                     | Penicillium commune       |    65 |
| Botryosphaeriaceae       | Botryosphaeria      | Botryosphaeria dothidea   |    77 |
|                          | Macrophomina        | Macrophomina phaseolina   |    64 |
| Cryphonectriaceae        | Cryphonectria       | Cryphonectria parasitica  |    73 |
| Cryptococcaceae          | Cryptococcus        | Cryptococcus neoformans   |   189 |
| Debaryomycetaceae        | Candida             | Candida albicans          |    70 |
|                          |                     | Candida parapsilosis      |    80 |
|                          | Scheffersomyces     | Scheffersomyces spartinae |    85 |
| Metschnikowiaceae        | Candidozyma         | Candidozyma auris         |   507 |
|                          | Clavispora          | Clavispora lusitaniae     |    88 |
| Nectriaceae              | Fusarium            | Fusarium asiaticum        |   252 |
|                          |                     | Fusarium graminearum      |   130 |
|                          |                     | Fusarium oxysporum        |   358 |
|                          |                     | Fusarium verticillioides  |    61 |
| Onygenaceae              | Ophidiomyces        | Ophidiomyces ophidiicola  |    73 |
| Phaeosphaeriaceae        | Parastagonospora    | Parastagonospora nodorum  |   167 |
| Pichiaceae               | Komagataella        | Komagataella phaffii      |   135 |
| Pleosporaceae            | Alternaria          | Alternaria alternata      |   132 |
| Pyriculariaceae          | Pyricularia         | Pyricularia oryzae        |   338 |
| Saccharomycetaceae       | Lachancea           | Lachancea thermotolerans  |   108 |
|                          | Nakaseomyces        | Nakaseomyces glabratus    |    76 |
|                          | Saccharomyces       | Saccharomyces cerevisiae  |  1540 |
|                          | Torulaspora         | Torulaspora delbrueckii   |    67 |
| Saccotheciaceae          | Aureobasidium       | Aureobasidium pullulans   |    50 |
| Schizosaccharomycetaceae | Schizosaccharomyces | Schizosaccharomyces pombe |    89 |
| Sporidiobolaceae         | Rhodotorula         | Rhodotorula mucilaginosa  |   120 |
| Trimorphomycetaceae      | Saitozyma           | Saitozyma podzolica       |    50 |

### For *protein families*

```shell
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --count \
    --in summary/pass.lst \
    --not-in MinHash/abnormal.lst \
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

| item    | count |
| ------- | ----: |
| strain  |  4067 |
| species |  1770 |
| genus   |   736 |
| family  |   337 |
| order   |   145 |
| class   |    49 |

| genus            | #species | #strains |
| ---------------- | -------: | -------: |
| Alternaria       |       21 |       65 |
| Aspergillus      |      116 |      422 |
| Candida          |       15 |       58 |
| Colletotrichum   |       45 |       76 |
| Cryptococcus     |       12 |       63 |
| Exophiala        |        9 |       57 |
| Fusarium         |       52 |      189 |
| Ogataea          |        7 |       55 |
| Ophidiomyces     |        1 |       69 |
| Parastagonospora |        1 |      153 |
| Penicillium      |       85 |      224 |
| Saccharomyces    |        8 |      421 |
| Trichoderma      |       32 |       54 |

## Collect proteins

```shell
cd ~/data/Fungi/

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
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

| #item      |      count |
| ---------- | ---------: |
| species    |      1,776 |
| strain_sum |      3,854 |
| total_sum  | 41,268,126 |
| dedup_sum  | 41,268,126 |
| rep_sum    | 23,464,695 |
| fam88_sum  | 21,630,459 |
| fam38_sum  | 18,205,855 |

## Phylogenetics with fungi61

### Find corresponding proteins by `hmmsearch`

```shell
cd ~/data/Fungi

mkdir -p HMM

# The Fungi HMM set
tar xvfz ~/data/HMM/fungi61/fungi61.tar.gz --directory=HMM
cp HMM/fungi61.lst HMM/marker.lst
```

### Find corresponding representative proteins by `hmmsearch`

```shell
cd ~/data/Fungi

# The profiles do not have NC fields
E_VALUE=1e-20

cat Protein/species.tsv |
    tva join -f summary/pass.lst -k 1 |
    tva join -f summary/rep.lst -k 1 |
    tva join -f summary/NR.lst -k 1 |
    tva join -e -f MinHash/abnormal.lst -k 1 |
    tva join -e -f summary/omit.lst -k 1 \
    > Protein/species-f.tsv

cat Protein/species-f.tsv |
    tva select -f 2 |
    tva uniq |
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
    tva stats --group-by 1 --count |
    tva stats --quantile 2:0.25,0.5,0.75
#627	662	1639

# There are 461 species and 616 strains
fd --full-path "Protein/.+/fungi61.tsv" -X cat |
    tva stats --group-by 1 --count |
    tva filter --invert --ge 2:600 --le 2:750 |
    cut -f 1 \
    > Protein/marker.omit.lst

wc -l HMM/marker.lst Protein/marker.omit.lst
# 61 HMM/marker.lst
# 32 Protein/marker.omit.lst

cat Protein/species-f.tsv |
    tva select -f 2 |
    tva uniq |
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
cd ~/data/Fungi

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

    pgr fa some Protein/"${SPECIES}"/pro.fa.gz <(
            tva select -f 1 Protein/"${SPECIES}"/seq_asm_f3.tsv |
            tva uniq
        )
done |
    pgr fa dedup stdin |
    pgr fa gz stdin -o Domain/fungi61.fa.gz

fd --full-path "Protein/.+/seq_asm_f3.tsv" -X cat \
    > Domain/seq_asm_f3.tsv

cat Domain/seq_asm_f3.tsv |
    tva join -e -d 2 -f summary/redundant.lst -k 1 \
    > Domain/seq_asm_f3.NR.tsv
```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Fungi

# Extract proteins
cat HMM/marker.lst |
    grep -v -Fw -f Protein/marker.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

        mkdir -p Domain/{}

        pgr fa some Domain/fungi61.fa.gz <(
            cat Domain/seq_asm_f3.tsv |
                tva filter --str-eq "3:{}" |
                tva select -f 1 |
                tva uniq
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
        tva filter --str-eq "3:${marker}" |
        tva select -f 1-2 |
        pgr fa replace -s Domain/${marker}/${marker}.aln.fa stdin \
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
    tva uniq |
    sort |
    fasops concat Domain/fungi61.aln.fas stdin -o Domain/fungi61.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Domain/fungi61.aln.fa -out Domain/fungi61.trim.fa -automated1

pgr fa size Domain/fungi61.*.fa |
    tva uniq -f 2 |
    cut -f 2
# 194230
# 11241

# To make it faster
FastTree -fastest -noml Domain/fungi61.trim.fa > Domain/fungi61.trim.newick
```

### Condense branches in the protein tree

```shell
cd ~/data/Fungi/tree

pgr nwk reroot ../Domain/fungi61.trim.newick -n Enc_hellem_ATCC_50604_GCA_024399255_1 -n Nematoc_paris_ERTm3_GCA_000190615_1 |
    pgr nwk order stdin --nd --an \
    > fungi61.reroot.newick

pgr pl condense --map -t ../Count/strains.taxon.tsv -r 5 -r 4 -r 3 -r 2 \
    fungi61.reroot.newick |
    pgr nwk order stdin --nd --an \
    > fungi61.condensed.newick

mv condensed.tsv fungi61.condense.tsv

# pdf
pgr nwk topo --bl fungi61.condensed.newick | # remove comments
    pgr nwk to-tex stdin --bl -o Fungi.fungi61.tex

tectonic Fungi.fungi61.tex
```

## InterProScan on all proteins of representative and typical strains

To be filled by other projects

```shell
mkdir -p ~/data/Fungi/STRAINS
```

## Groups and targets

Review `summary/collect.pass.tsv` and `tree/groups.tsv`

