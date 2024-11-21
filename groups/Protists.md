# Various genera from Protists

Eukaryota other than Plants, Metazoa, and Fungi

<!-- TOC -->
* [Various genera from Protists](#various-genera-from-protists)
  * [Taxon info](#taxon-info)
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
  * [Phylogenetics with BUSCO](#phylogenetics-with-busco)
    * [Find corresponding representative proteins by `hmmsearch`](#find-corresponding-representative-proteins-by-hmmsearch)
    * [Domain related protein sequences](#domain-related-protein-sequences)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
    * [Condense branches in the protein tree](#condense-branches-in-the-protein-tree)
  * [Groups and targets](#groups-and-targets)
<!-- TOC -->

## Taxon info

* Amoebozoa 变形虫
    * Discosea
    * Evosea
        * Eumycetozoa
            * Dictyostelium 网柄菌属
            * Entamoeba 内变形虫属
* Apusozoa
    * Apusomonadida
* Breviatea
    * Lenisia
* Cryptophyceae (cryptomonads) 隐藻纲
    * Guillardia
* Discoba
    * Euglenozoa 眼虫门
        * Kinetoplastea 动基体
            * Leishmania 利什曼虫属
            * Trypanosoma 锥虫属
            * Angomonas
            * Porcisia
            * Crithidia 短膜虫属
            * Endotrypanum 内锥虫属
    * Heterolobosea 异叶足纲
        * Naegleria 耐格里属
        * Tetramitus 四鞭毛虫属
* Haptista
    * Haptophyta
* Metamonada 后滴门
    * Fornicata
        * Diplomonadida 双滴虫目
            * Giardia 贾第虫属
    * Parabasalia (parabasalids) 副基体门
        * Trichomonas 毛滴虫属
    * Preaxostyla
* Opisthokonta
    * Choanoflagellata
    * Filasterea
    * Fungi
    * Ichthyosporea
    * Metazoa (metazoans)
    * Rotosphaerida
    * Opisthokonta incertae sedis
* Rhodophyta (red algae) 红藻门
    * Bangiophyceae 红毛菜纲
        * Cyanidioschyzon
        * Neoporphyra
    * Florideophyceae
* Sar
    * Alveolata (alveolates) 囊泡虫类
        * Apicomplexa 顶复门
        * Aconoidasida 孢子虫纲
            * Babesia 巴倍虫属
            * Plasmodium 疟原虫属
            * Theileria 泰勒虫属
        * Conoidasida 类锥体纲
            * Cryptosporidium 隐孢子虫属
            * Eimeria 艾美球虫
            * Neospora 新孢子虫属
            * Toxoplasma 弓形虫属
        * Ciliophora 纤毛虫门
            * Oligohymenophorea 寡膜纲
                * Tetrahymena 四膜虫
    * Rhizaria 有孔虫类
        * Cercozoa
            * Chlorarachniophyceae
                * Bigelowiella 比奇洛藻属
    * Stramenopiles (heterokonts) 不等鞭毛类
        * Bacillariophyta 硅藻门
            * Bacillariophyceae 硅藻纲
                * Fragilariopsis 拟脆杆藻属
                * Phaeodactylum 褐指藻属
                * Epithemia 窗纹藻属
            * Coscinodiscophyceae 圆筛藻纲
                * Thalassiosira 海链藻属
        * Eustigmatophyceae 眼点藻纲
            * Nannochloropsis 微拟球藻
        * Phaeophyceae 褐藻纲
            * Ectocarpus 外子藻属
            * Undaria 裙带菜属
        * Oomycota 卵菌门
            * Aphanomyces 丝囊霉属
            * Elongisporangium
            * Globisporangium
            * Hyaloperonospora
            * Plasmopara 单轴霉属
            * Phytophthora 疫霉属
            * Pythium 腐霉属
            * Albugo
            * Peronospora 霜霉属
            * Phytopythium
* Viridiplantae 绿色植物
    * Chlorophyta (green algae)
    * Streptophyta 链形植物

### Species with assemblies

```shell
mkdir -p ~/data/Protists/summary
cd ~/data/Protists/summary

# should have a valid name of genus
nwr member Eukaryota -r genus |
    nwr restrict -e Viridiplantae |
    nwr restrict -e Metazoa |
    nwr restrict -e Fungi |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#4099 genus.list

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
#  122 RS1.tsv
#  905 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     123
#GB1     2059

```

## Download all assemblies

### Create .assembly.tsv

```shell
cd ~/data/Protists/summary

# Reference genome
echo "
.headers ON
    SELECT
        *
    FROM ar
    WHERE 1=1
        AND species IN ('Caenorhabditis elegans', 'Drosophila melanogaster', 'Arabidopsis thaliana', 'Saccharomyces cerevisiae')
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
#2065 lines, 7 fields

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
    > Protists.assembly.tsv

datamash check < Protists.assembly.tsv
#2056 lines, 5 fields

# find potential duplicate strains or assemblies
cat Protists.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Protists.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

# Edit .assembly.tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Protists.assembly.tsv
#
# Save the file to another directory to prevent accidentally changing it
# cp Protists.assembly.tsv ~/Scripts/genomes/assembly

# Cleaning
rm raw*.*sv

```

### Count before download

```shell
cd ~/data/Protists

nwr template ~/Scripts/genomes/assembly/Protists.assembly.tsv \
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
    tsv-filter -H --ge 3:10 |
    rgr md stdin --num

```

| item    | count |
|---------|------:|
| strain  |  2041 |
| species |   844 |
| genus   |   350 |
| family  |   183 |
| order   |   115 |
| class   |    50 |

| genus           | #species | #strains |
|-----------------|---------:|---------:|
| Acanthamoeba    |       17 |       28 |
| Albugo          |        2 |       12 |
| Aphanomyces     |        5 |       27 |
| Babesia         |       10 |       30 |
| Blastocystis    |        2 |       74 |
| Crithidia       |        8 |       11 |
| Cryptosporidium |       16 |       77 |
| Cyclospora      |        1 |       41 |
| Cyclotella      |        7 |       10 |
| Dictyostelium   |        9 |       11 |
| Ectocarpus      |        3 |       12 |
| Eimeria         |       11 |       19 |
| Entamoeba       |        5 |       16 |
| Galdieria       |        4 |       15 |
| Giardia         |        2 |       38 |
| Globisporangium |       48 |       82 |
| Isotricha       |        3 |       13 |
| Leishmania      |       26 |      123 |
| Naegleria       |        3 |       14 |
| Nannochloropsis |        6 |       17 |
| Ochromonas      |        2 |       11 |
| Ophryoscolex    |        1 |       10 |
| Paramecium      |        9 |       18 |
| Peronospora     |        6 |       18 |
| Phytophthora    |       77 |      249 |
| Phytopythium    |       14 |       18 |
| Plasmodiophora  |        1 |       51 |
| Plasmodium      |       21 |      163 |
| Plasmopara      |        4 |       11 |
| Pythium         |       45 |       69 |
| Skeletonema     |        8 |       12 |
| Symbiodinium    |        6 |       13 |
| Thalassiosira   |       15 |       20 |
| Theileria       |        6 |       16 |
| Toxoplasma      |        1 |       29 |
| Trypanosoma     |       22 |       75 |

### Download and check

```shell
cd ~/data/Protists

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Protists.assembly.tsv \
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

find ASSEMBLY/ -name "*_genomic.fna.gz" |
    grep -v "_from_" |
    wc -l
#2054

# N50 C S; create n50.tsv and n50.pass.tsv
# LEN_N50   N_CONTIG    LEN_SUM
bash ASSEMBLY/n50.sh 5000 30000 5000000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50" --max "C" --min "S"
#N50_min C_max   S_min
#7707    29495   6434485

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "N50:0.1,0.5" --quantile "C:0.5,0.9" --quantile "S:0.1,0.5" |
    datamash transpose
#N50_pct10       3293.6
#N50_pct50       64426.5
#C_pct50 1790
#C_pct90 32978.6
#S_pct10 11412417.1
#S_pct50 39091805.5

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
| url.tsv          |      3 | 2,055 |
| check.lst        |      1 | 2,054 |
| collect.tsv      |     20 | 2,055 |
| n50.tsv          |      4 | 2,055 |
| n50.pass.tsv     |      4 | 1,702 |
| collect.pass.tsv |     23 | 1,702 |
| pass.lst         |      1 | 1,701 |
| omit.lst         |      1 | 1,516 |
| rep.lst          |      1 |   691 |
| sp.lst           |      1 |   155 |

### Rsync to hpcc

```shell
rsync -avP \
    ~/data/Protists/ \
    wangq@202.119.37.251:data/Protists

# rsync -avP wangq@202.119.37.251:data/Protists/ ~/data/Protists

```

## BioSample

ENA's BioSample missed many strains, so NCBI's was used.

```shell
cd ~/data/Protists

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Protists.assembly.tsv \
    --bs

# Run this script twice and it will re-download the failed files
bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 10

datamash check < BioSample/biosample.tsv
#2023 lines, 121 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## MinHash

```shell
cd ~/data/Protists

nwr template ~/Scripts/genomes/assembly/Protists.assembly.tsv \
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
#  621 summary/NR.lst
#  623 summary/redundant.lst

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#33

# Distances between all selected sketches, then hierarchical clustering
nwr template ~/Scripts/genomes/assembly/Protists.assembly.tsv \
    --mh \
    --parallel 8 \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in summary/redundant.lst \
    --height 0.4

bash MinHash/dist.sh

```

### Condense branches in the minhash tree

```shell
mkdir -p ~/data/Protists/tree
cd ~/data/Protists/tree

nwr order ../MinHash/tree.nwk --nd --an \
    > minhash.order.newick

nwr pl-condense --map -r order -r family -r genus -r species \
    minhash.order.newick ../MinHash/species.tsv |
    nwr order stdin --nd --an \
    > minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

# svg
nwr topo --bl minhash.condensed.newick | # remove comments
    nw_display -s -b 'visibility:hidden' -w 1200 -v 20 - \
    > Protists.minhash.svg

```

## Count valid species and strains

### For *genomic alignments*

```shell
cd ~/data/Protists/

nwr template ~/Scripts/genomes/assembly/Protists.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --rank order --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv
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
bash Count/lineage.sh 10

cat Count/lineage.count.tsv |
    rgr md stdin --num

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| item    | count |
|---------|------:|
| strain  |  1659 |
| species |   672 |
| genus   |   258 |
| family  |   140 |
| order   |    96 |
| class   |    49 |

| order            | #species | #strains |
|------------------|---------:|---------:|
| Diplomonadida    |        4 |       36 |
| Eucoccidiorida   |       32 |      170 |
| Haemosporida     |       23 |      151 |
| Opalinata        |        3 |       53 |
| Peronosporales   |       85 |      271 |
| Piroplasmida     |       17 |       46 |
| Plasmodiophorida |        3 |       53 |
| Pythiales        |      111 |      156 |
| Saprolegniales   |        9 |       27 |
| Trypanosomatida  |       86 |      238 |

| genus           | #species | #strains |
|-----------------|---------:|---------:|
| Aphanomyces     |        5 |       23 |
| Babesia         |       10 |       29 |
| Blastocystis    |        2 |       52 |
| Cryptosporidium |       16 |       76 |
| Cyclospora      |        1 |       41 |
| Giardia         |        2 |       34 |
| Globisporangium |       48 |       70 |
| Leishmania      |       26 |      117 |
| Phytophthora    |       66 |      232 |
| Plasmodiophora  |        1 |       51 |
| Plasmodium      |       21 |      149 |
| Pythium         |       41 |       59 |
| Toxoplasma      |        1 |       28 |
| Trypanosoma     |       22 |       64 |

| #family           | genus           | species                    | count |
|-------------------|-----------------|----------------------------|------:|
| Albuginaceae      | Albugo          | Albugo candida             |    11 |
| Babesiidae        | Babesia         | Babesia microti            |    10 |
| Blastocystidae    | Blastocystis    | Blastocystis sp.           |    48 |
| Cryptosporidiidae | Cryptosporidium | Cryptosporidium hominis    |    15 |
|                   |                 | Cryptosporidium parvum     |    26 |
| Eimeriidae        | Cyclospora      | Cyclospora cayetanensis    |    41 |
| Entamoebidae      | Entamoeba       | Entamoeba histolytica      |    12 |
| Hexamitidae       | Giardia         | Giardia intestinalis       |    33 |
| Parameciidae      | Paramecium      | Paramecium bursaria        |    10 |
| Peronosporaceae   | Phytophthora    | Phytophthora cactorum      |    28 |
|                   |                 | Phytophthora capsici       |    16 |
|                   |                 | Phytophthora fragariae     |    13 |
|                   |                 | Phytophthora kernoviae     |    12 |
|                   |                 | Phytophthora nicotianae    |    13 |
|                   |                 | Phytophthora ramorum       |    28 |
|                   |                 | Phytophthora sojae         |    11 |
| Plasmodiidae      | Plasmodium      | Plasmodium falciparum      |    59 |
|                   |                 | Plasmodium vinckei         |    10 |
|                   |                 | Plasmodium vivax           |    19 |
|                   |                 | Plasmodium yoelii          |    14 |
| Plasmodiophoridae | Plasmodiophora  | Plasmodiophora brassicae   |    51 |
| Pythiaceae        | Globisporangium | Globisporangium irregulare |    15 |
|                   | Pythium         | Pythium insidiosum         |    11 |
| Saprolegniaceae   | Aphanomyces     | Aphanomyces astaci         |    10 |
| Sarcocystidae     | Toxoplasma      | Toxoplasma gondii          |    28 |
| Trypanosomatidae  | Leishmania      | Leishmania donovani        |    10 |
|                   |                 | Leishmania infantum        |    39 |
|                   | Trypanosoma     | Trypanosoma cruzi          |    34 |

### For *protein families*

```shell
cd ~/data/Protists/

nwr template ~/Scripts/genomes/assembly/Protists.assembly.tsv \
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
|---------|------:|
| strain  |   493 |
| species |   272 |
| genus   |   142 |
| family  |    89 |
| order   |    69 |
| class   |    35 |

| genus           | #species | #strains |
|-----------------|---------:|---------:|
| Aphanomyces     |        5 |       21 |
| Babesia         |       10 |       10 |
| Cryptosporidium |       13 |       23 |
| Leishmania      |       15 |       29 |
| Peronospora     |        5 |       10 |
| Phytophthora    |       19 |       72 |
| Plasmodium      |       20 |       78 |
| Toxoplasma      |        1 |       14 |
| Trypanosoma     |       10 |       22 |

## Collect proteins

```shell
cd ~/data/Protists/

nwr template ~/Scripts/genomes/assembly/Protists.assembly.tsv \
    --pro \
    --parallel 8 \
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

| #item      |     count |
|------------|----------:|
| species    |       273 |
| strain_sum |       503 |
| total_sum  | 6,904,114 |
| dedup_sum  | 6,904,114 |
| rep_sum    | 4,526,186 |
| fam88_sum  | 3,960,800 |
| fam38_sum  | 3,333,428 |

## Phylogenetics with BUSCO

```shell
cd ~/data/Protists/

rm -fr BUSCO

curl -L https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2024-01-08.tar.gz |
    tar xvz

mv eukaryota_odb10/ BUSCO

```

### Find corresponding representative proteins by `hmmsearch`

```shell
cd ~/data/Protists

cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 \
    > Protein/species-f.tsv

cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ -s Protein/"${SPECIES}"/busco.tsv ]]; then
        continue
    fi
    if [[ ! -f Protein/"${SPECIES}"/rep_seq.fa.gz ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    cat BUSCO/scores_cutoff |
        parallel --colsep '\s+' --no-run-if-empty --linebuffer -k -j 4 "
            gzip -dcf Protein/${SPECIES}/rep_seq.fa.gz |
                hmmsearch -T {2} --domT {2} --noali --notextw BUSCO/hmms/{1}.hmm - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), q({1}), \$1; '
        " \
        > Protein/${SPECIES}/busco.tsv
done

fd --full-path "Protein/.+/busco.tsv" -X cat |
    tsv-summarize --group-by 1 --count |
    tsv-summarize --quantile 2:0.25,0.5,0.75
#207     267     335.5

# There are 273 species and 503 strains
fd --full-path "Protein/.+/busco.tsv" -X cat |
    tsv-summarize --group-by 1 --count |
    tsv-filter --invert --ge 2:260 --le 2:520 |
    cut -f 1 \
    > Protein/marker.omit.lst

cat BUSCO/scores_cutoff |
    parallel --colsep '\s+' --no-run-if-empty --linebuffer -k -j 1 "
        echo {1}
    " \
    > Protein/marker.lst

wc -l Protein/marker.lst Protein/marker.omit.lst
# 255 Protein/marker.lst
# 146 Protein/marker.omit.lst

cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ ! -s Protein/"${SPECIES}"/busco.tsv ]]; then
        continue
    fi
    if [[ ! -f Protein/"${SPECIES}"/seq.sqlite ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    # single copy
    cat Protein/"${SPECIES}"/busco.tsv |
        grep -v -Fw -f Protein/marker.omit.lst \
        > Protein/"${SPECIES}"/busco.sc.tsv

    nwr seqdb -d Protein/${SPECIES} --rep f3=Protein/${SPECIES}/busco.sc.tsv

done

```

### Domain related protein sequences

```shell
cd ~/data/Protists

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
    hnsm gz stdin -o Domain/busco.fa

fd --full-path "Protein/.+/seq_asm_f3.tsv" -X cat \
    > Domain/seq_asm_f3.tsv

cat Domain/seq_asm_f3.tsv |
    tsv-join -e -d 2 -f summary/redundant.lst -k 1 \
    > Domain/seq_asm_f3.NR.tsv

```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Protists

# Extract proteins
cat Protein/marker.lst |
    grep -v -Fw -f Protein/marker.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

        mkdir -p Domain/{}

        hnsm some Domain/busco.fa.gz <(
            cat Domain/seq_asm_f3.tsv |
                tsv-filter --str-eq "3:{}" |
                tsv-select -f 1 |
                rgr dedup stdin
            ) \
            > Domain/{}/{}.pro.fa
    '

# Align each marker
cat Protein/marker.lst |
    grep -v -Fw -f Protein/marker.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
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

cat Protein/marker.lst |
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
cat Protein/marker.lst |
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
    > Domain/busco.aln.fas

cat Domain/seq_asm_f3.NR.tsv |
    cut -f 2 |
    rgr dedup stdin |
    sort |
    fasops concat Domain/busco.aln.fas stdin -o Domain/busco.aln.fa

trimal -in Domain/busco.aln.fa -out Domain/busco.trim.fa -automated1

hnsm size Domain/busco.*.fa |
    rgr dedup stdin -f 2 |
    cut -f 2
#787011
#170

# To make it faster
FastTree -fastest -noml Domain/busco.trim.fa > Domain/busco.trim.newick

```

### Condense branches in the protein tree

```shell
cd ~/data/Protists/tree

nwr order ../Domain/busco.trim.newick --nd --an \
    > busco.order.newick

nwr pl-condense --map -r order -r family -r genus \
    busco.order.newick ../Count/species.tsv |
    nwr order stdin --nd --an \
    > busco.condensed.newick

mv condensed.tsv busco.condense.tsv

# svg
nwr topo --bl busco.condensed.newick | # remove comments
    nw_display -s -b 'visibility:hidden' -w 1200 -v 20 - \
    > Protists.busco.svg

```

## Groups and targets

Review `ASSEMBLY/Protists.assembly.pass.csv` and `mash/groups.tsv`

| #Serial | Group              | Count | Target                 | Sequencing     |
|:--------|:-------------------|:------|:-----------------------|:---------------|
| 1       | A_Ei               | 4     | A_cas                  |                |
| 2       | Babesia            | 7     | Ba_bov_T2Bo            | 8X Sanger, WGS |
| 4       | Blastocystis       | 8     | Bl_hom                 |                |
| 8       | Cryptosporidium    | 11    | Cry_parvum_Iowa_II     |                |
| 7       | Dictyostelium      | 11    | D_disc_AX4             |                |
| 10      | Giardia            | 6     | G_intes                |                |
| 11      | Leishmania         | 24    | L_maj_Friedlin         |                |
| 12      | Nannochloropsis    | 8     | N_oce                  |                |
| 13      | Phytophthora       | 10    | Ph_soj                 |                |
| 15      | Pl_ber_cha_vin_yoe | 10    | Pl_yoe                 |                |
| 16      | Pl_kno_viv         | 12    | Pl_viv                 |                |
| 17      | Pl_falcip          | 31    | Pl_falcip_3D7          |                |
| 18      | Pythium            | 4     | Py_gui                 |                |
| 19      | Theileria          | 6     | Th_parva_Muguga        |                |
| 20      | Toxoplasma         | 16    | To_gondii_ME49         |                |
| 21      | Trypanosoma        | 5     | Tr_bruc_brucei_TREU927 |                |

```shell
$(brew --prefix repeatmasker)/libexec/util/queryRepeatDatabase.pl \
    -species Eukaryota -stat
#174 ancestral and ubiquitous sequence(s) with a total length of 45804 bp
#0 eukaryota specific repeats with a total length of 0 bp
#31314 lineage specific sequence(s) with a total length of 81842705 bp

$(brew --prefix repeatmasker)/libexec/util/queryRepeatDatabase.pl \
    -species Alveolata -stat
#174 ancestral and ubiquitous sequence(s) with a total length of 45804 bp
#3 alveolata specific repeats with a total length of 7681 bp
#30 lineage specific sequence(s) with a total length of 81127 bp

$(brew --prefix repeatmasker)/libexec/util/queryRepeatDatabase.pl \
    -species Kinetoplastida -stat
#176 ancestral and ubiquitous sequence(s) with a total length of 48931 bp
#0 kinetoplastida specific repeats with a total length of 0 bp
#44 lineage specific sequence(s) with a total length of 67956 bp

$(brew --prefix repeatmasker)/libexec/util/queryRepeatDatabase.pl \
    -species Stramenopiles -stat
#174 ancestral and ubiquitous sequence(s) with a total length of 45804 bp
#3 stramenopiles specific repeats with a total length of 6250 bp
#620 lineage specific sequence(s) with a total length of 2090320 bp

```

