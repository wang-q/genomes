# Fungi

Download all genomes and analyze representative strains.

<!-- TOC -->

* [Fungi](#fungi)
    * [Taxon info](#taxon-info)
        * [List all ranks](#list-all-ranks)
        * [Species with assemblies](#species-with-assemblies)
        * [Model organisms](#model-organisms)
    * [Download all assemblies](#download-all-assemblies)
        * [Create assembly.tsv](#create-assemblytsv)
        * [Count before download](#count-before-download)
        * [Download and check](#download-and-check)
        * [Rsync to hpcc](#rsync-to-hpcc)
    * [BioSample](#biosample)
    * [Divergence of Fungi](#divergence-of-fungi)
        * [ReRoot](#reroot)
    * [MinHash](#minhash)
        * [Condense branches in the minhash tree](#condense-branches-in-the-minhash-tree)
    * [Count valid species and strains](#count-valid-species-and-strains)
        * [For *genomic alignments*](#for-genomic-alignments)
        * [For *protein families*](#for-protein-families)
    * [Collect proteins](#collect-proteins)
    * [Phylogenetics with fungi61](#phylogenetics-with-fungi61)
        * [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
        * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
        * [Condense branches in the protein tree](#condense-branches-in-the-protein-tree)
    * [InterProScan on all proteins of representative and typical strains](#interproscan-on-all-proteins-of-representative-and-typical-strains)
    * [Groups and targets](#groups-and-targets)

<!-- TOC -->

## Taxon info

* [Fungi](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4751)

Open nomenclature (开放命名法)

* Sp. (pl. spp.; short for "species") indicates potentially new species
* Sp. aff. or aff. (short for "species affinis") indicates a potentially new and undescribed species
  has an affinity to, but is not identical to, the named species
* Cf. (short for the Latin: confer, "compare with") or a question mark (?, also inc., species
  incerta) signify varying degrees or types of uncertainty
* Candidatus, a candidate taxon proposed from metagenomics or other incomplete information
* Incertae sedis, a taxon of uncertain position in a classification
* Nomen invalidum, a name that is not validly published

### List all ranks

```shell
nwr member Fungi |
    grep -v " sp." |
    grep -v " aff." |
    grep -v " cf." |
    grep -v " x " |
    tsv-summarize -H -g rank --count |
    rgr md stdin --num

nwr member Fungi | grep " sp\."
nwr member Fungi | grep " aff\."
nwr member Fungi | grep " cf\."
nwr member Fungi -r species | grep "\." | grep -v " sp\." | grep -v " aff\." | grep -v " cf\."

nwr member Fungi | grep " x "

```

| rank          | count |
|---------------|------:|
| kingdom       |     1 |
| no rank       |  5306 |
| species       | 61884 |
| subkingdom    |     1 |
| class         |    71 |
| order         |   256 |
| family        |   926 |
| genus         |  7817 |
| phylum        |    10 |
| subphylum     |    14 |
| strain        |  2258 |
| varietas      |  1044 |
| forma         |   222 |
| isolate       |     5 |
| subspecies    |   144 |
| subclass      |    19 |
| suborder      |    24 |
| subfamily     |    19 |
| subgenus      |    10 |
| clade         |    28 |
| section       |    38 |
| species group |     1 |
| tribe         |     3 |
| superfamily   |     1 |
| morph         |     2 |

### Species with assemblies

In the vast majority of fungal species, only one genome was selected for refseq.

* 'RefSeq'
* 'Genbank'

```shell
mkdir -p ~/data/Fungi/summary
cd ~/data/Fungi/summary

# should have a valid name of genus
nwr member Fungi -r genus |
    grep -v " x " |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#7817 genus.list

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
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB1.tsv

wc -l RS*.tsv GB*.tsv
#   627 RS1.tsv
#  4919 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     632
#GB1     17810

```

* The names of some genera are abnormal

```shell
cd ~/data/Fungi/summary

cat RS*.tsv GB*.tsv |
    cut -f 1,2 |
    tsv-uniq |
    grep "\[" |
    nwr append stdin -r genus |
    head
#561895  [Candida] subhashii     Spathaspora
#566037  [Ashbya] aceris (nom. inval.)   Eremothecium
#45518   [Candida] aaseri        Yamadazyma
#1171601 [Candida] adriatica     Cyberlindnera
#456249  [Candida] andamanensis (nom. inval.)    Yamadazyma
#148631  [Candida] anglica       Kurtzmaniella
#130810  [Candida] arabinofermentans     Ogataea
#391823  [Candida] ascalaphidarum        Yamadazyma
#45521   [Candida] atlantica     Yamadazyma
#45522   [Candida] atmosphaerica Yamadazyma

```

### Model organisms

There is one true model organism in Fungi, Saccharomyces cerevisiae S288C

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
        refseq_category,
        COUNT(*) AS count
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
    GROUP BY refseq_category
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
#refseq_category count
#na      5
#reference genome        627

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
        cut -f 1 |
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
        AND genome_rep IN ('Full')
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% aff.%'
        AND species NOT LIKE '% cf.%'
        AND species NOT LIKE '% x %'
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
        AND genome_rep IN ('Full')
        AND (
            species LIKE '% sp.%'
            OR species LIKE '% aff.%'
            OR species LIKE '% cf.%'
        )
        AND species NOT LIKE '% x %'
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

cat raw.tsv |
    tsv-uniq |
    datamash check
#17737 lines, 7 fields

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
    > Fungi.assembly.tsv

datamash check < Fungi.assembly.tsv
#17425 lines, 5 fields

# find potential duplicate strains or assemblies
cat Fungi.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Fungi.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

# Edit .assembly.tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Fungi.assembly.tsv
#
# Save the file to another directory to prevent accidentally changing it
# cp Fungi.assembly.tsv ~/Scripts/genomes/assembly

# Cleaning
rm raw*.*sv

```

### Count before download

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class

```shell
cd ~/data/Fungi

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --count \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --num

# .lst and .count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    keep-header -- tsv-sort -k1,1 |
    tsv-filter -H --ge 3:100 |
    rgr md stdin --num

```

| item    | count |
|---------|------:|
| strain  | 17419 |
| species |  4838 |
| genus   |  1389 |
| family  |   485 |
| order   |   176 |
| class   |    61 |

| genus             | #species | #strains |
|-------------------|---------:|---------:|
| Alternaria        |       35 |      180 |
| Aspergillus       |      143 |     1280 |
| Aureobasidium     |       12 |      167 |
| Beauveria         |       11 |      330 |
| Botryosphaeria    |        3 |      141 |
| Candida           |       34 |      270 |
| Candidozyma       |        9 |      221 |
| Colletotrichum    |       81 |      277 |
| Cryphonectria     |        7 |      110 |
| Cryptococcus      |       12 |      215 |
| Exophiala         |       13 |      102 |
| Fusarium          |      189 |     1587 |
| Komagataella      |        7 |      189 |
| Metschnikowia     |       63 |      145 |
| Ogataea           |       55 |      146 |
| Parastagonospora  |        2 |      189 |
| Penicillium       |      110 |      482 |
| Pichia            |       30 |      134 |
| Psilocybe         |       12 |      141 |
| Pyricularia       |        5 |      410 |
| Rhodotorula       |       10 |      178 |
| Saccharomyces     |       11 |     1829 |
| Torulaspora       |        8 |      108 |
| Trichoderma       |       34 |      138 |
| Zygosaccharomyces |       12 |      146 |

### Download and check

```shell
cd ~/data/Fungi

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
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
#17423

# N50 C S; create n50.tsv and n50.pass.tsv
bash ASSEMBLY/n50.sh 100000 2000 1000000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50" --max "C" --min "S"
#N50_min C_max   S_min
#32179   2016    2187595

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "N50:0.1,0.5" --quantile "C:0.5,0.9" --quantile "S:0.1,0.5"
#N50_pct10       N50_pct50       C_pct50 C_pct90 S_pct10 S_pct50
#20869.6 327499  348     4815.2  11659596        32316518

# After the above steps are completed, run the following commands.

# Collect; create collect.tsv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    rgr md stdin --right 2-3

```

| #item            | fields |  lines |
|------------------|-------:|-------:|
| url.tsv          |      3 | 17,423 |
| check.lst        |      1 | 17,423 |
| collect.tsv      |     20 | 17,424 |
| n50.tsv          |      4 | 17,424 |
| n50.pass.tsv     |      4 | 12,035 |
| collect.pass.tsv |     23 | 12,035 |
| pass.lst         |      1 | 12,034 |
| omit.lst         |      1 | 12,704 |
| rep.lst          |      1 |  3,385 |
| sp.lst           |      1 |    157 |

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Fungi/ \
    wangq@202.119.37.251:data/Fungi

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

datamash check < BioSample/biosample.tsv
#17097 lines, 208 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## Divergence of Fungi

Ref.:

1. A genome-scale phylogeny of the kingdom Fungi. Cur Biology, 2021.
   https://doi.org/10.1016/j.cub.2021.01.074

![1-s2.0-S0960982221001391-fx1_lrg.jpg](images/1-s2.0-S0960982221001391-fx1_lrg.jpg)

### ReRoot

* 小孢子虫 (Microsporidia) 可以当作真菌的基部类群, 也可以当作真菌的姊妹类群

* 芽枝霉门 (Blastocladiomycota) 是真菌的基部类群

* For the latter steps, use the following two as the outgroups
    * Enc_hellem_ATCC_50504_GCF_000277815_2
    * Nematoc_ausu_ERTm6_GCF_000738915_1

```shell
cd ~/data/Fungi

cat summary/collect.pass.tsv |
    tsv-select -f 1,3 |
    sed '1d' |
    grep -v -Fw -f ASSEMBLY/omit.lst |
    nwr append stdin -c 2 -r species -r phylum -r subkingdom |
    tsv-filter --str-ne "5:Dikarya" | # 双核亚界
    tsv-filter --str-ne "4:Mucoromycota" | # 毛霉门
    tsv-filter --str-ne "4:Zoopagomycota" | #
    tsv-filter --str-ne "4:Chytridiomycota" | # 壶菌门
    tsv-select -f 1,3,4 |
    tsv-sort -k3,3 -k1,1 |
    tsv-summarize -g 3,2 --count
#Blastocladiomycota      Allomyces arbusculus    1
#Blastocladiomycota      Allomyces javanicus     1
#Blastocladiomycota      Allomyces macrogynus    1
#Blastocladiomycota      Blastocladiella emersonii       1
#Blastocladiomycota      Catenaria anguillulae   1
#Blastocladiomycota      Paraphysoderma sedebokerense    1
#Microsporidia   Edhazardia aedis        1
#Microsporidia   Encephalitozoon cuniculi        3
#Microsporidia   Encephalitozoon hellem  4
#Microsporidia   Encephalitozoon intestinalis    2
#Microsporidia   Encephalitozoon romaleae        1
#Microsporidia   Enterocytozoon hepatopenaei     1
#Microsporidia   Hamiltosporidium tvaerminnensis 1
#Microsporidia   Nematocida ausubeli     2
#Microsporidia   Nematocida displodere   1
#Microsporidia   Nematocida major        1
#Microsporidia   Nematocida parisii      2
#Microsporidia   Ordospora colligata     4
#Microsporidia   Ordospora pajunii       1
#Microsporidia   Vairimorpha bombi       1
#Microsporidia   Vairimorpha ceranae     1
#Microsporidia   Vairimorpha necatrix    1
#Microsporidia   Vittaforma corneae      1

cat summary/collect.pass.tsv |
    tsv-filter -H --not-blank RefSeq_category |
    tsv-filter -H --or \
        --str-in-fld "2:Allomyces" \
        --str-in-fld "2:Blastocladiella" \
        --str-in-fld "2:Encephalitozoon" \
        --str-in-fld "2:Nematocida" |
    grep -v -Fw -f ASSEMBLY/omit.lst |
    tsv-select -H -f "#name,Assembly_level,Assembly_method,Genome_coverage,Sequencing_technology"

```

## MinHash

```shell
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
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
#  3222 summary/NR.lst
#  6573 summary/redundant.lst

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#188

find MinHash -name "redundant.lst" -size +0 | wc -l
#778

find MinHash -name "redundant.lst" -empty | wc -l
#394

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --mh \
    --parallel 8 \
    --not-in ASSEMBLY/sp.lst \
    --not-in ASSEMBLY/omit.lst \
    --not-in MinHash/abnormal.lst \
    --not-in summary/redundant.lst \
    --height 0.4

bash MinHash/dist.sh

```

### Condense branches in the minhash tree

```shell
mkdir -p ~/data/Fungi/tree
cd ~/data/Fungi/tree

nw_reroot ../MinHash/tree.nwk Enc_hellem_ATCC_50504_GCF_000277815_2 Nematoc_ausu_ERTm6_GCF_000738915_1 |
    nwr order stdin --nd --an \
    > minhash.reroot.newick

nwr pl-condense --map -r order -r family -r genus \
    minhash.reroot.newick ../MinHash/species.tsv |
    nwr order stdin --nd --an \
    > minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

# png
nwr topo --bl minhash.condensed.newick | # remove comments
    nw_display -s -b 'visibility:hidden' -w 1200 -v 20 - |
    rsvg-convert -o Fungi.minhash.png

```

## Count valid species and strains

### For *genomic alignments*

```shell
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --rank order --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/order.count.tsv |
    tsv-filter -H --ge "3:100" |
    rgr md stdin --num

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:100" |
    rgr md stdin --num

# Can accept N_COUNT
bash Count/lineage.sh 50

cat Count/lineage.count.tsv |
    rgr md stdin --num

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| item    | count |
|---------|------:|
| strain  | 12033 |
| species |  3410 |
| genus   |  1021 |
| family  |   392 |
| order   |   154 |
| class   |    54 |

| order             | #species | #strains |
|-------------------|---------:|---------:|
| Agaricales        |      178 |      286 |
| Botryosphaeriales |       27 |      186 |
| Chaetothyriales   |       48 |      153 |
| Diaporthales      |       46 |      149 |
| Dipodascales      |      138 |      245 |
| Eurotiales        |      276 |     1682 |
| Glomerellales     |       82 |      240 |
| Helotiales        |       92 |      142 |
| Hypocreales       |      304 |     1537 |
| Magnaporthales    |       10 |      192 |
| Mycosphaerellales |       65 |      195 |
| Onygenales        |       42 |      183 |
| Phaffomycetales   |       77 |      108 |
| Pichiales         |      131 |      479 |
| Pleosporales      |      148 |      592 |
| Polyporales       |       71 |      106 |
| Saccharomycetales |      114 |     1925 |
| Serinales         |      333 |      927 |
| Sordariales       |       73 |      133 |
| Sporidiobolales   |       16 |      178 |
| Tremellales       |       42 |      305 |
| Xylariales        |      108 |      196 |

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Alternaria       |       35 |      176 |
| Aspergillus      |      127 |     1092 |
| Candida          |       23 |      168 |
| Candidozyma      |        9 |      204 |
| Colletotrichum   |       67 |      176 |
| Cryptococcus     |       12 |      208 |
| Fusarium         |      120 |     1050 |
| Komagataella     |        7 |      186 |
| Ogataea          |       54 |      133 |
| Parastagonospora |        2 |      168 |
| Penicillium      |      106 |      452 |
| Pyricularia      |        4 |      185 |
| Rhodotorula      |       10 |      171 |
| Saccharomyces    |       11 |     1569 |
| Trichoderma      |       30 |      124 |

| #family                  | genus               | species                   | count |
|--------------------------|---------------------|---------------------------|------:|
| Aspergillaceae           | Aspergillus         | Aspergillus flavus        |   265 |
|                          |                     | Aspergillus fumigatus     |   260 |
|                          |                     | Aspergillus niger         |   112 |
|                          |                     | Aspergillus oryzae        |   116 |
|                          | Penicillium         | Penicillium chrysogenum   |    83 |
| Botryosphaeriaceae       | Botryosphaeria      | Botryosphaeria dothidea   |    76 |
| Cryphonectriaceae        | Cryphonectria       | Cryphonectria parasitica  |    68 |
| Cryptococcaceae          | Cryptococcus        | Cryptococcus neoformans   |   175 |
| Debaryomycetaceae        | Candida             | Candida albicans          |    68 |
| Metschnikowiaceae        | Candidozyma         | Candidozyma auris         |   188 |
|                          | Clavispora          | Clavispora lusitaniae     |    73 |
| Nectriaceae              | Fusarium            | Fusarium asiaticum        |   250 |
|                          |                     | Fusarium graminearum      |   125 |
|                          |                     | Fusarium oxysporum        |   281 |
|                          |                     | Fusarium verticillioides  |    57 |
| Onygenaceae              | Ophidiomyces        | Ophidiomyces ophidiicola  |    73 |
| Phaeosphaeriaceae        | Parastagonospora    | Parastagonospora nodorum  |   166 |
| Pichiaceae               | Komagataella        | Komagataella phaffii      |   134 |
| Pleosporaceae            | Alternaria          | Alternaria alternata      |    87 |
| Pyriculariaceae          | Pyricularia         | Pyricularia oryzae        |   177 |
| Saccharomycetaceae       | Nakaseomyces        | Nakaseomyces glabratus    |    56 |
|                          | Saccharomyces       | Saccharomyces cerevisiae  |  1453 |
|                          | Torulaspora         | Torulaspora delbrueckii   |    72 |
| Schizosaccharomycetaceae | Schizosaccharomyces | Schizosaccharomyces pombe |    84 |
| Sporidiobolaceae         | Rhodotorula         | Rhodotorula mucilaginosa  |   117 |
| Trimorphomycetaceae      | Saitozyma           | Saitozyma podzolica       |    50 |

### For *protein families*

```shell
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
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
    tsv-filter -H --ge "3:50" |
    rgr md stdin --num

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv

```

| item    | count |
|---------|------:|
| strain  |  3545 |
| species |  1563 |
| genus   |   692 |
| family  |   317 |
| order   |   135 |
| class   |    48 |

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Aspergillus      |       80 |      371 |
| Candida          |       14 |       50 |
| Colletotrichum   |       44 |       73 |
| Cryptococcus     |       11 |       62 |
| Exophiala        |        9 |       57 |
| Fusarium         |       46 |      154 |
| Ogataea          |        7 |       55 |
| Ophidiomyces     |        1 |       69 |
| Parastagonospora |        1 |      153 |
| Penicillium      |       78 |      209 |
| Saccharomyces    |        9 |      411 |

## Collect proteins

```shell
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
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

| #item      |      count |
|------------|-----------:|
| species    |      1,568 |
| strain_sum |      3,602 |
| total_sum  | 38,555,598 |
| dedup_sum  | 38,555,598 |
| rep_sum    | 21,106,780 |
| fam88_sum  | 19,379,724 |
| fam38_sum  | 16,301,299 |

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

E_VALUE=1e-20

cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -f ASSEMBLY/rep.lst -k 1 |
    tsv-join -f summary/NR.lst -k 1 |
    tsv-join -e -f MinHash/abnormal.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 \
    > Protein/species-f.tsv

cat Protein/species-f.tsv |
    tsv-select -f 2 |
    tsv-uniq |
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
#628     664     1659

# There are 461 species and 616 strains
fd --full-path "Protein/.+/fungi61.tsv" -X cat |
    tsv-summarize --group-by 1 --count |
    tsv-filter --invert --ge 2:600 --le 2:750 |
    cut -f 1 \
    > Protein/marker.omit.lst

wc -l HMM/marker.lst Protein/marker.omit.lst
# 61 HMM/marker.lst
# 33 Protein/marker.omit.lst

cat Protein/species-f.tsv |
    tsv-select -f 2 |
    tsv-uniq |
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
    hnsm gz stdin -o Domain/fungi61.fa

fd --full-path "Protein/.+/seq_asm_f3.tsv" -X cat \
    > Domain/seq_asm_f3.tsv

cat Domain/seq_asm_f3.tsv |
    tsv-join -e -d 2 -f summary/NR.lst -k 1 \
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

        hnsm some Domain/fungi61.fa.gz <(
            cat Domain/seq_asm_f3.tsv |
                tsv-filter --str-eq "3:{}" |
                tsv-select -f 1 |
                tsv-uniq
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
    tsv-uniq |
    sort |
    fasops concat Domain/fungi61.aln.fas stdin -o Domain/fungi61.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Domain/fungi61.aln.fa -out Domain/fungi61.trim.fa -automated1

hnsm size Domain/fungi61.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#192590
#24940

# To make it faster
FastTree -fastest -noml Domain/fungi61.trim.fa > Domain/fungi61.trim.newick

```

### Condense branches in the protein tree

```shell
cd ~/data/Fungi/tree

nw_reroot ../Domain/fungi61.trim.newick Enc_hellem_ATCC_50604_GCA_024399255_1 Nematoc_paris_ERTm3_GCA_000190615_1 |
    nwr order stdin --nd --an \
    > fungi61.reroot.newick

nwr pl-condense -r order -r family -r genus -r species \
    fungi61.reroot.newick ../Count/species.tsv --map \
    -o fungi61.condensed.newick

mv condensed.tsv fungi61.condense.tsv

# pdf
nwr tex fungi61.condensed.newick --bl -o Fungi.fungi61.tex

tectonic Fungi.fungi61.tex

```

## InterProScan on all proteins of representative and typical strains

To be filled by other projects

```shell
mkdir -p ~/data/Fungi/STRAINS

```

## Groups and targets

Review `summary/collect.pass.tsv` and `tree/groups.tsv`
