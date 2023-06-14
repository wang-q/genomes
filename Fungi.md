# Fungi

Download all genomes and analyze representative strains.

<!-- toc -->

- [Taxon info](#taxon-info)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
    * [Model organisms](#model-organisms)
- [Download all assemblies](#download-all-assemblies)
    * [Create assembly.tsv](#create-assemblytsv)
    * [Rsync and check](#rsync-and-check)
    * [Check N50 of assemblies](#check-n50-of-assemblies)
    * [Rsync to hpcc](#rsync-to-hpcc)
- [BioSample](#biosample)
- [Count species and strains](#count-species-and-strains)
    * [Order](#order)
    * [Genus](#genus)
    * [ReRoot](#reroot)
- [MinHash](#minhash)
    * [Compute MinHash](#compute-minhash)
    * [Raw phylo-tree of representative assemblies](#raw-phylo-tree-of-representative-assemblies)
    * [Tweak the mash tree](#tweak-the-mash-tree)
- [Non-redundant strains within species](#non-redundant-strains-within-species)
    * [Genus](#genus-1)
- [Groups and targets](#groups-and-targets)
- [Collect proteins](#collect-proteins)
    * [`all.pro.fa`](#allprofa)
    * [`all.replace.fa`](#allreplacefa)
    * [`all.info.tsv`](#allinfotsv)
- [Phylogenetics with fungi61](#phylogenetics-with-fungi61)
    * [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
    * [Tweak the concat tree](#tweak-the-concat-tree)
- [InterProScan on all proteins of representative and typical strains](#interproscan-on-all-proteins-of-representative-and-typical-strains)
- [Divergence of Fungi](#divergence-of-fungi)

<!-- tocstop -->

## Taxon info

* [Fungi](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4751)

### List all ranks

```shell
nwr member Fungi |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

```

| rank          | count |
|---------------|------:|
| kingdom       |     1 |
| no rank       |  4623 |
| species       | 65297 |
| subkingdom    |     1 |
| class         |    65 |
| order         |   242 |
| family        |   904 |
| genus         |  7466 |
| phylum        |    10 |
| subphylum     |    14 |
| strain        |  2236 |
| varietas      |  1060 |
| subspecies    |   162 |
| forma         |   210 |
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
    * RS1 - assembly_level: 'Complete Genome', 'Chromosome'
    * RS2 - genome_rep: 'Full'
* 'Genbank'
    * '>= 20 genomes'
        * GB1 - With strain ID; assembly_level: 'Complete Genome', 'Chromosome'
        * GB2 - assembly_level: 'Complete Genome', 'Chromosome'
        * GB3 - genome_rep: 'Full'
    * '>= 1 genomes'
        * GB4 - assembly_level: 'Complete Genome', 'Chromosome'
        * GB5 - genome_rep: 'Full'

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
#7466 genus.list

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
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
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
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
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND assembly_level IN ('Complete Genome', 'Chromosome')
            AND genome_rep IN ('Full')
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
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND assembly_level IN ('Complete Genome', 'Chromosome')
            AND genome_rep IN ('Full')
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
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
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
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND assembly_level IN ('Complete Genome', 'Chromosome')
        GROUP BY species_id
        HAVING count >= 1
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
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB5.tsv

wc -l RS*.tsv GB*.tsv
#    91 RS1.tsv
#   499 RS2.tsv
#     1 GB1.tsv
#     4 GB2.tsv
#    85 GB3.tsv
#   281 GB4.tsv
#  3395 GB5.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     94
#RS2     504
#GB1     109
#GB2     770
#GB3     7260
#GB4     1325
#GB5     13371

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
#498019  [Candida] auris Candida/Metschnikowiaceae
#1231522 [Candida] duobushaemulonis      Candida/Metschnikowiaceae
#45357   [Candida] haemuloni     Candida/Metschnikowiaceae
#418784  [Candida] pseudohaemulonii      Candida/Metschnikowiaceae
#561895  [Candida] subhashii     Spathaspora
#566037  [Ashbya] aceris (nom. inval.)   Eremothecium
#312227  [Candida] hispaniensis  Candida/Saccharomycetales
#45354   [Candida] intermedia    Candida/Metschnikowiaceae
#2093215 [Candida] vulturna (nom. inval.)        Clavispora
#45518   [Candida] aaseri        Yamadazyma

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
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > reference.tsv

cp reference.tsv ~/Scripts/genomes/assembly/Fungi.reference.tsv

```

## Download all assemblies

### Create assembly.tsv

Two levels:

* 'RefSeq'
    * RS2 - genome_rep: 'Full'
* 'Genbank'
    * GB5 - genome_rep: 'Full'

If a refseq assembly is available, the corresponding genbank one is not downloaded

```shell
cd ~/data/Fungi/summary

cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession \
    > raw.tsv

# RS2
SPECIES=$(
    cat RS2.tsv |
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

# Preference for refseq
cat raw.tsv |
    tsv-select -H -f "assembly_accession" \
    > rs.acc.tsv

# GB5
SPECIES=$(
    cat GB5.tsv |
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
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

cat raw.tsv |
    tsv-uniq |
    datamash check
#13376 lines, 7 fields

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
#13109 lines, 5 fields

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

# strains.taxon.tsv
bash Count/strains.sh

# genus.lst and genus.count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    keep-header -- tsv-sort -k1,1 |
    tsv-filter -H --ge 3:100 |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```

| genus                     | #species | #strains |
|---------------------------|---------:|---------:|
| Alternaria                |       31 |      157 |
| Aspergillus               |      133 |     1072 |
| Aureobasidium             |       11 |      154 |
| Botryosphaeria            |        3 |      141 |
| Candida                   |       18 |      186 |
| Candida/Metschnikowiaceae |        9 |      158 |
| Colletotrichum            |       60 |      225 |
| Cryphonectria             |        7 |      110 |
| Cryptococcus              |       11 |      149 |
| Fusarium                  |      188 |     1216 |
| Komagataella              |        7 |      181 |
| Parastagonospora          |        2 |      189 |
| Penicillium               |      103 |      459 |
| Pyricularia               |        3 |      381 |
| Rhodotorula               |       10 |      158 |
| Saccharomyces             |       11 |     1765 |
| Torulaspora               |        8 |      101 |
| Trichoderma               |       31 |      105 |
| Zygosaccharomyces         |        9 |      136 |

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
#13108

# N50 C S; create n50.tsv and n50.pass.tsv
bash ASSEMBLY/n50.sh 100000 2000 1000000

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

```

| #item            | fields |  lines |
|------------------|-------:|-------:|
| url.tsv          |      3 | 13,108 |
| check.lst        |      1 | 13,108 |
| collect.tsv      |     20 | 13,109 |
| n50.tsv          |      4 | 13,109 |
| n50.pass.tsv     |      4 |  9,345 |
| collect.pass.tsv |     23 |  9,345 |
| pass.lst         |      1 |  9,344 |
| omit.lst         |      1 |  8,926 |
| rep.lst          |      1 |  2,543 |

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Fungi/ \
    wangq@202.119.37.251:data/Fungi

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Fungi/ \
    wangq@58.213.64.36:data/Fungi

# scp -P 8804 ./nwr wangq@58.213.64.36:bin

# rsync -avP wangq@202.119.37.251:data/Fungi/ ~/data/Fungi

# rsync -avP -e 'ssh -p 8804' wangq@58.213.64.36:data/Fungi/ ~/data/Fungi

```

For such huge collections, we can rsync files in parallel.

```shell
# Copy Directory Structure
rsync -avP \
    -e 'ssh -p 8804' \
    -f"+ */" -f"- *" \
    wangq@58.213.64.36:data/Fungi/ \
    ~/data/Fungi

# Transfer species directories in parallel
cat ~/data/Fungi/ASSEMBLY/url.tsv |
    tsv-select -f 3 |
    tsv-uniq |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo -e "\n==> {}"
        rsync -avP \
            -e "ssh -p 8804" \
            wangq@58.213.64.36:data/Fungi/ASSEMBLY/{}/ \
            ~/data/Fungi/ASSEMBLY/{}
    '

# rsync other files
rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Fungi/ \
    wangq@58.213.64.36:data/Fungi

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
#12865 lines, 169 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## Divergence of Fungi

Ref.:

1. A genome-scale phylogeny of the kingdom Fungi. Cur Biology, 2021.
   https://doi.org/10.1016/j.cub.2021.01.074

![1-s2.0-S0960982221001391-fx1_lrg.jpg](images%2F1-s2.0-S0960982221001391-fx1_lrg.jpg)

### ReRoot

* 小孢子虫 (Microsporidia) 可以当作真菌的基部类群, 也可以当作真菌的姊妹类群

* 芽枝霉门 (Blastocladiomycota) 是真菌的基部类群

* For the latter steps, use the following two as the outgroups
    * Enc_hell_ATCC_50504_GCF_000277815_2
    * Nemat_ausu_ERTm6_GCF_000738915_1

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
#Microsporidia   Nematocida ausubeli     2
#Microsporidia   Nematocida displodere   1
#Microsporidia   Nematocida major        1
#Microsporidia   Nematocida parisii      2
#Microsporidia   Nosema ceranae  1
#Microsporidia   Ordospora colligata     4
#Microsporidia   Vittaforma corneae      1

cat summary/collect.pass.tsv |
    tsv-filter -H --not-blank RefSeq_category |
    tsv-filter -H --or \
        --str-in-fld "2:Allomyces" \
        --str-in-fld "2:Blastocladiella" \
        --str-in-fld "2:Encephalitozoon" \
        --str-in-fld "2:Nematocida" |
    grep -v -Fw -f ASSEMBLY/omit.lst |
    tsv-select -H -f "name,Assembly_level,Assembly_method,Genome_coverage,Sequencing_technology"

```

## MinHash

```shell
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
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
#359

# Non-redundant strains within species
bash MinHash/nr.sh

find MinHash -name "mash.dist.tsv" -size +0 | wc -l
#808

find MinHash -name "redundant.lst" -size +0 | wc -l
#503

find MinHash -name "redundant.lst" -empty | wc -l
#304

find MinHash -name "NR.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/NR.lst
wc -l summary/NR.lst
#2403

## All representative should be in NR
#cat ASSEMBLY/rep.lst |
#    grep -v -F -f summary/NR.lst

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
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
mkdir -p ~/data/Fungi/tree
cd ~/data/Fungi/tree

nw_reroot ../MinHash/tree.nwk Enc_hell_ATCC_50504_GCF_000277815_2 Nemat_ausu_ERTm6_GCF_000738915_1 |
    nw_order -c n - \
    > minhash.reroot.newick

# Avoid "Candida/Metschnikowiaceae"
cat ../Count/strains.taxon.tsv |
    sed "s/\//_/g" \
    > strains.taxon.tsv

# rank::col
ARRAY=(
    'order::5'
    'family::4'
    'genus::3'
#    'species::2'
)

rm minhash.condensed.map
CUR_TREE=minhash.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/genomes/bin/condense_tree.sh ${CUR_TREE} strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick minhash.${GROUP_NAME}.newick
    cat condense.map >> minhash.condensed.map

    CUR_TREE=minhash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 minhash.genus.newick |
    rsvg-convert -o Fungi.minhash.png

```

## Count valid species and strains

### For *genomic alignments*

```shell
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --rank order --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv
bash Count/strains.sh

# .lst and .count.tsv
bash Count/rank.sh

cat Count/order.count.tsv |
    tsv-filter -H --ge "3:50" |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:50" |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

# Can accept N_COUNT
bash Count/lineage.sh 50

cat Count/lineage.count.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| order             | #species | #strains |
|-------------------|---------:|---------:|
| Agaricales        |      138 |      192 |
| Botryosphaeriales |       19 |      134 |
| Chaetothyriales   |       41 |      107 |
| Diaporthales      |       41 |      135 |
| Dothideales       |       13 |       65 |
| Eurotiales        |      246 |     1408 |
| Glomerellales     |       65 |      194 |
| Helotiales        |       83 |      107 |
| Hypocreales       |      271 |     1003 |
| Magnaporthales    |        8 |      177 |
| Mucorales         |       52 |       74 |
| Mycosphaerellales |       55 |      148 |
| Onygenales        |       39 |      177 |
| Ophiostomatales   |       54 |       82 |
| Pleosporales      |      129 |      523 |
| Polyporales       |       60 |       82 |
| Russulales        |       32 |       51 |
| Saccharomycetales |      436 |     2791 |
| Sordariales       |       21 |       62 |
| Sporidiobolales   |       10 |      141 |
| Tremellales       |       32 |      204 |
| Trichosporonales  |       29 |       51 |
| Ustilaginales     |       34 |       85 |
| Xylariales        |       86 |      152 |

| genus                     | #species | #strains |
|---------------------------|---------:|---------:|
| Alternaria                |       31 |      150 |
| Aspergillus               |      116 |      896 |
| Aureobasidium             |       11 |       62 |
| Botryosphaeria            |        1 |       70 |
| Calonectria               |       15 |       54 |
| Candida                   |       17 |      118 |
| Candida/Metschnikowiaceae |        9 |      141 |
| Colletotrichum            |       50 |      133 |
| Cryphonectria             |        6 |       82 |
| Cryptococcus              |       11 |      128 |
| Exophiala                 |       11 |       50 |
| Fusarium                  |      115 |      666 |
| Komagataella              |        7 |      167 |
| Metschnikowia             |       30 |       50 |
| Nakaseomyces              |        6 |       52 |
| Ogataea                   |       23 |       79 |
| Ophidiomyces              |        1 |       73 |
| Parastagonospora          |        1 |      166 |
| Penicillium               |       96 |      405 |
| Pichia                    |       12 |       51 |
| Pyricularia               |        3 |      172 |
| Rhodotorula               |        9 |      140 |
| Saccharomyces             |       11 |     1497 |
| Torulaspora               |        7 |       76 |
| Trichoderma               |       26 |       88 |
| Trichophyton              |       11 |       52 |
| Ustilago                  |       14 |       56 |
| Verticillium              |       10 |       53 |
| Zymoseptoria              |        3 |       50 |

| #family            | genus                     | species                  | count |
|--------------------|---------------------------|--------------------------|------:|
| Aspergillaceae     | Aspergillus               | Aspergillus flavus       |   190 |
|                    |                           | Aspergillus fumigatus    |   246 |
|                    |                           | Aspergillus niger        |   106 |
|                    |                           | Aspergillus oryzae       |   102 |
|                    | Penicillium               | Penicillium chrysogenum  |    83 |
| Botryosphaeriaceae | Botryosphaeria            | Botryosphaeria dothidea  |    70 |
| Cryphonectriaceae  | Cryphonectria             | Cryphonectria parasitica |    68 |
| Cryptococcaceae    | Cryptococcus              | Cryptococcus neoformans  |    98 |
| Debaryomycetaceae  | Candida                   | Candida albicans         |    60 |
| Metschnikowiaceae  | Candida/Metschnikowiaceae | [Candida] auris          |   128 |
| Nectriaceae        | Fusarium                  | Fusarium graminearum     |   120 |
|                    |                           | Fusarium oxysporum       |   253 |
| Onygenaceae        | Ophidiomyces              | Ophidiomyces ophidiicola |    73 |
| Phaeosphaeriaceae  | Parastagonospora          | Parastagonospora nodorum |   166 |
| Phaffomycetaceae   | Komagataella              | Komagataella phaffii     |   133 |
| Pleosporaceae      | Alternaria                | Alternaria alternata     |    84 |
| Pyriculariaceae    | Pyricularia               | Pyricularia oryzae       |   167 |
| Saccharomycetaceae | Saccharomyces             | Saccharomyces cerevisiae |  1398 |
|                    | Torulaspora               | Torulaspora delbrueckii  |    60 |
| Sporidiobolaceae   | Rhodotorula               | Rhodotorula mucilaginosa |   106 |

### For *protein families*

```shell
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst \
    --rank genus

# strains.taxon.tsv
bash Count/strains.sh

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:50" |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv

```

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Aspergillus      |       79 |      355 |
| Colletotrichum   |       29 |       54 |
| Cryptococcus     |       10 |       53 |
| Fusarium         |       43 |      140 |
| Ogataea          |        7 |       55 |
| Ophidiomyces     |        1 |       69 |
| Parastagonospora |        1 |      153 |
| Penicillium      |       71 |      193 |
| Saccharomyces    |        9 |      410 |

## Collect proteins

```shell
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --pro \
    --in ASSEMBLY/pass.lst \
    --in ASSEMBLY/rep.lst \
    --in summary/NR.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst

# collect proteins
bash Protein/collect.sh

cat Protein/counts.tsv |
    mlr --itsv --omd cat

```

| #item                          | count     |
|--------------------------------|-----------|
| Proteins                       | 3,761,901 |
| Unique headers and annotations | 3,761,901 |
| Unique proteins                | 3,761,901 |
| all.replace.fa                 | 3,761,901 |
| all.annotation.tsv             | 3,761,902 |
| all.info.tsv                   | 3,761,902 |

## Phylogenetics with fungi61

### Find corresponding proteins by `hmmsearch`

```shell
cd ~/data/Fungi

# The fungi61 HMM set
nwr kb fungi61 -o HMM
cp HMM/fungi61.lst HMM/marker.lst

E_VALUE=1e-20

# Find all genes
for marker in $(cat HMM/marker.lst); do
    echo >&2 "==> marker [${marker}]"

    mkdir -p Protein/${marker}

    cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -f ASSEMBLY/rep.lst -k 1 |
        tsv-join -f summary/NR.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/hmm/${marker}.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/${marker}/replace.tsv

    echo >&2
done

```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Fungi

cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat Protein/{}/replace.tsv |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
#367     382     930

cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        cat Protein/{}/replace.tsv |
            wc -l
    ' |
    paste - - |
    tsv-filter --invert --ge 2:300 --le 2:500 |
    cut -f 1 \
    > Protein/marker.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process
cat HMM/marker.lst |
    grep -v -Fx -f Protein/marker.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

        cat Protein/{}/replace.tsv \
            > Protein/{}/{}.replace.tsv

        faops some Protein/all.uniq.fa.gz <(
            cat Protein/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > Protein/{}/{}.pro.fa
    '

# Align each markers with muscle
cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo >&2 "==> marker [{}]"
        if [ ! -s Protein/{}/{}.pro.fa ]; then
            exit
        fi
        if [ -s Protein/{}/{}.aln.fa ]; then
            exit
        fi

        muscle -quiet -in Protein/{}/{}.pro.fa -out Protein/{}/{}.aln.fa
    '

for marker in $(cat HMM/marker.lst); do
    echo >&2 "==> marker [${marker}]"
    if [ ! -s Protein/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sometimes `muscle` can not produce alignments
    if [ ! -s Protein/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # 1 name to many names
    cat Protein/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s Protein/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > Protein/${marker}/${marker}.replace.fa
done

# Concat marker genes
for marker in $(cat HMM/marker.lst); do
    if [ ! -s Protein/${marker}/${marker}.pro.fa ]; then
        continue
    fi
    if [ ! -s Protein/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # sequences in one line
    faops filter -l 0 Protein/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > Protein/fungi61.aln.fas

cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -f ASSEMBLY/rep.lst -k 1 |
    tsv-join -f summary/NR.lst -k 1 |
    tsv-join -e -f MinHash/abnormal.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
    cut -f 1 |
    fasops concat Protein/fungi61.aln.fas stdin -o Protein/fungi61.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Protein/fungi61.aln.fa -out Protein/fungi61.trim.fa -automated1

faops size Protein/fungi61.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#91210
#13874

# To make it faster
FastTree -fastest -noml Protein/fungi61.trim.fa > Protein/fungi61.trim.newick

```

### Condense branches in the protein tree

```shell
cd ~/data/Fungi/tree

nw_reroot ../Protein/fungi61.trim.newick Mit_dap_GCF_000760515_2 No_cera_GCF_000988165_1 |
    nw_order -c n - \
    > fungi61.reroot.newick

# Avoid "Candida/Metschnikowiaceae"
cat ../Count/strains.taxon.tsv |
    sed "s/\//_/g" \
    > strains.taxon.tsv

# rank::col
ARRAY=(
    'order::5'
    'family::4'
    'genus::3'
#    'species::2'
)

rm fungi61.condensed.map
CUR_TREE=fungi61.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/genomes/bin/condense_tree.sh ${CUR_TREE} strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick fungi61.${GROUP_NAME}.newick
    cat condense.map >> fungi61.condensed.map

    CUR_TREE=fungi61.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 fungi61.genus.newick |
    rsvg-convert -o Fungi.fungi61.png

```


```shell
cd ~/data/Fungi/tree

nw_reroot ../MProtein/tree.nwk Enc_hell_ATCC_50504_GCF_000277815_2 Nemat_major_JUm2507_GCF_021653875_1 |
    nw_order -c n - \
    > minhash.reroot.newick

rm minhash.condensed.map
CUR_TREE=minhash.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/genomes/bin/condense_tree.sh ${CUR_TREE} strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick minhash.${GROUP_NAME}.newick
    cat condense.map >> minhash.condensed.map

    CUR_TREE=minhash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 minhash.genus.newick |
    rsvg-convert -o Fungi.minhash.png

```

## InterProScan on all proteins of representative and typical strains

To be filled by other projects

```shell
mkdir -p ~/data/Fungi/STRAINS

```

## Groups and targets

Review `summary/collect.pass.tsv` and `tree/groups.tsv`
