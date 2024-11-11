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
#11602116.3      32228284.5      21666.1 315609  375     5088.5

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
| url.tsv          |      3 | 15,812 |
| check.lst        |      1 | 15,811 |
| collect.tsv      |     20 | 15,812 |
| n50.tsv          |      4 | 15,813 |
| n50.pass.tsv     |      4 | 10,884 |
| collect.pass.tsv |     23 | 10,884 |
| pass.lst         |      1 | 10,883 |
| omit.lst         |      1 | 11,262 |
| rep.lst          |      1 |  3,294 |
| sp.lst           |      1 |    156 |

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Fungi/ \
    wangq@202.119.37.251:data/Fungi

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Fungi/ \
    wangq@58.213.64.36:data/Fungi

# rsync -avP wangq@202.119.37.251:data/Fungi/ ~/data/Fungi

# rsync -avP -e 'ssh -p 8804' wangq@58.213.64.36:data/Fungi/ ~/data/Fungi

```

For such huge collections, we can rsync files inside ASSEMBLY/ in parallel.

```shell
# Copy Directory Structure
rsync -avP \
    -e 'ssh -p 8804' \
    -f"+ */" -f"- *" \
    wangq@58.213.64.36:data/Fungi/ASSEMBLY/ \
    ~/data/Fungi/ASSEMBLY

# Transfer species directories in parallel
cat ~/data/Fungi/ASSEMBLY/url.tsv |
    tsv-select -f 3 |
    tsv-uniq |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo -e "\n==> {}"
        rsync -avP \
            -e "ssh -p 8804" \
            wangq@58.213.64.36:data/Fungi/ASSEMBLY/{}/ \
            ~/data/Fungi/ASSEMBLY/{}
    '

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
#15539 lines, 196 fields

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
#Microsporidia   Ecytonucleospora hepatopenaei   1
#Microsporidia   Edhazardia aedis        1
#Microsporidia   Encephalitozoon cuniculi        3
#Microsporidia   Encephalitozoon hellem  4
#Microsporidia   Encephalitozoon intestinalis    2
#Microsporidia   Encephalitozoon romaleae        1
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
#467

# Non-redundant strains within species
bash MinHash/nr.sh

find MinHash -name "mash.dist.tsv" -size +0 | wc -l
#1083

find MinHash -name "redundant.lst" -size +0 | wc -l
#722

find MinHash -name "redundant.lst" -empty | wc -l
#361

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
#  3055 summary/NR.lst
#  5754 summary/redundant.lst

## All representative should be in NR
#cat ASSEMBLY/rep.lst |
#    grep -v -F -f summary/NR.lst

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --mh \
    --parallel 16 \
    --in ASSEMBLY/rep.lst \
    --not-in ASSEMBLY/sp.lst \
    --not-in ASSEMBLY/omit.lst \
    --not-in MinHash/abnormal.lst \
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
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

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

| item    | count |
|---------|------:|
| strain  | 10882 |
| species |  3205 |
| genus   |   958 |
| family  |   380 |
| order   |   140 |
| class   |    48 |

| order             | #species | #strains |
|-------------------|---------:|---------:|
| Agaricales        |      171 |      259 |
| Boletales         |       46 |       54 |
| Botryosphaeriales |       27 |      179 |
| Chaetothyriales   |       48 |      148 |
| Diaporthales      |       44 |      142 |
| Dothideales       |       13 |       89 |
| Eurotiales        |      263 |     1557 |
| Glomerellales     |       82 |      236 |
| Helotiales        |       87 |      126 |
| Hypocreales       |      294 |     1193 |
| Magnaporthales    |        9 |      188 |
| Malasseziales     |       18 |       70 |
| Mucorales         |       54 |       83 |
| Mycosphaerellales |       65 |      191 |
| Onygenales        |       41 |      182 |
| Ophiostomatales   |       54 |       83 |
| Orbiliales        |       24 |       55 |
| Pleosporales      |      142 |      570 |
| Polyporales       |       70 |      103 |
| Russulales        |       35 |       56 |
| Saccharomycetales |      828 |     3512 |
| Sordariales       |       64 |      121 |
| Sporidiobolales   |       15 |      172 |
| Tremellales       |       40 |      304 |
| Trichosporonales  |       31 |       67 |
| Ustilaginales     |       37 |       96 |
| Xylariales        |      100 |      186 |

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Alternaria       |       34 |      164 |
| Aspergillus      |      124 |      986 |
| Aureobasidium    |       11 |       85 |
| Beauveria        |        5 |       58 |
| Botryosphaeria   |        2 |       78 |
| Calonectria      |       15 |       55 |
| Candida          |       23 |      162 |
| Clavispora       |        4 |       80 |
| Colletotrichum   |       67 |      173 |
| Cryphonectria    |        6 |       82 |
| Cryptococcus     |       12 |      206 |
| Exophiala        |       13 |       80 |
| Fusarium         |      115 |      733 |
| Kazachstania     |       30 |       53 |
| Komagataella     |        7 |      183 |
| Malassezia       |       18 |       70 |
| Metschnikowia    |       34 |       74 |
| Nakaseomyces     |        7 |       69 |
| Neurospora       |        8 |       51 |
| Ogataea          |       54 |      131 |
| Ophidiomyces     |        1 |       73 |
| Parastagonospora |        2 |      168 |
| Penicillium      |      103 |      444 |
| Pichia           |       23 |       81 |
| Pyricularia      |        3 |      181 |
| Rhodotorula      |       10 |      166 |
| Saccharomyces    |       11 |     1524 |
| Saitozyma        |        1 |       50 |
| Talaromyces      |       22 |       59 |
| Torulaspora      |        8 |       94 |
| Trichoderma      |       29 |      113 |
| Trichophyton     |       12 |       56 |
| Ustilago         |       15 |       60 |
| Verticillium     |       10 |       55 |
| Wickerhamiella   |       29 |       51 |
| Yamadazyma       |       50 |       71 |
| Yarrowia         |       13 |       58 |
| Zymoseptoria     |        3 |       50 |

| #family             | genus            | species                  | count |
|---------------------|------------------|--------------------------|------:|
| Aspergillaceae      | Aspergillus      | Aspergillus flavus       |   228 |
|                     |                  | Aspergillus fumigatus    |   248 |
|                     |                  | Aspergillus niger        |   112 |
|                     |                  | Aspergillus oryzae       |   110 |
|                     | Penicillium      | Penicillium chrysogenum  |    83 |
| Botryosphaeriaceae  | Botryosphaeria   | Botryosphaeria dothidea  |    76 |
| Cryphonectriaceae   | Cryphonectria    | Cryphonectria parasitica |    68 |
| Cryptococcaceae     | Cryptococcus     | Cryptococcus neoformans  |   174 |
| Debaryomycetaceae   | Candida          | Candida albicans         |    63 |
| Metschnikowiaceae   | Clavispora       | Clavispora lusitaniae    |    73 |
| Nectriaceae         | Fusarium         | Fusarium graminearum     |   124 |
|                     |                  | Fusarium oxysporum       |   273 |
| Onygenaceae         | Ophidiomyces     | Ophidiomyces ophidiicola |    73 |
| Phaeosphaeriaceae   | Parastagonospora | Parastagonospora nodorum |   166 |
| Phaffomycetaceae    | Komagataella     | Komagataella phaffii     |   134 |
| Pleosporaceae       | Alternaria       | Alternaria alternata     |    86 |
| Pyriculariaceae     | Pyricularia      | Pyricularia oryzae       |   174 |
| Saccharomycetaceae  | Nakaseomyces     | Nakaseomyces glabratus   |    56 |
|                     | Saccharomyces    | Saccharomyces cerevisiae |  1420 |
|                     | Torulaspora      | Torulaspora delbrueckii  |    71 |
| Sporidiobolaceae    | Rhodotorula      | Rhodotorula mucilaginosa |   114 |
| Trimorphomycetaceae | Saitozyma        | Saitozyma podzolica      |    50 |

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
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:50" |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv

```

| item    | count |
|---------|------:|
| strain  |  3310 |
| species |  1441 |
| genus   |   642 |
| family  |   305 |
| order   |   120 |
| class   |    43 |

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Aspergillus      |       74 |      352 |
| Colletotrichum   |       44 |       71 |
| Cryptococcus     |       11 |       54 |
| Exophiala        |        9 |       51 |
| Fusarium         |       43 |      149 |
| Ogataea          |        7 |       55 |
| Ophidiomyces     |        1 |       69 |
| Parastagonospora |        1 |      153 |
| Penicillium      |       70 |      192 |
| Saccharomyces    |        9 |      410 |

## Collect proteins

```shell
cd ~/data/Fungi/

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --pro \
    --in ASSEMBLY/pass.lst \
    --not-in ASSEMBLY/omit.lst \
    --clust-id 0.95 \
    --clust-cov 0.95

# collect proteins
bash Protein/collect.sh

# clustering
bash Protein/compute.sh

# counts
bash Protein/count.sh

cat Protein/counts.tsv |
    tsv-summarize -H --count --sum 2-5 |
    sed 's/^count/species/' |
    datamash transpose |
    perl -nla -F"\t" -MNumber::Format -e '
        printf qq(%s\t%s\n), $F[0], Number::Format::format_number($F[1], 0,);
        ' |
    (echo -e "#item\tcount" && cat) |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

```

| #item      |      count |
|------------|-----------:|
| species    |      1,498 |
| strain_sum |      3,455 |
| total_sum  | 37,177,773 |
| dedup_sum  | 37,155,186 |
| rep_sum    | 20,265,558 |

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

nwr reroot ../Protein/fungi61.trim.newick -n Enc_hell_ATCC_50504_GCF_000277815_2 -n Nemat_ausu_ERTm6_GCF_000738915_1 |
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
