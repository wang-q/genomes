# Various genera from Protists

Eukaryota other than Plants, Metazoa, and Fungi

<!-- toc -->

* [Species with assemblies](#species-with-assemblies)

- [Download all assemblies](#download-all-assemblies)
    * [Create assembly.tsv](#create-assemblytsv)
    * [rsync and check](#rsync-and-check)
- [Count strains](#count-strains)
- [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash)
- [Groups and targets](#groups-and-targets)
- [Protists: prepare](#protists-prepare)
- [plasmodium: run](#plasmodium-run)
- [Protists: run](#protists-run)

<!-- tocstop -->

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
#4058 genus.list

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
#   97 RS1.tsv
#  831 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     98
#GB1     1930

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
        AND genus IN ('Caenorhabditis', 'Drosophila', 'Arabidopsis', 'Saccharomyces')
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
    tsv-uniq |
    datamash check
#1936 lines, 7 fields

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
    > Protists.assembly.tsv

datamash check < Protists.assembly.tsv
#1927 lines, 5 fields

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
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# genus.lst and genus.count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    keep-header -- tsv-sort -k1,1 |
    tsv-filter -H --ge 3:10 |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```

| item    | count |
|---------|------:|
| strain  |  1912 |
| species |   781 |
| genus   |   322 |
| family  |   165 |
| order   |   108 |
| class   |    48 |

| genus           | #species | #strains |
|-----------------|---------:|---------:|
| Acanthamoeba    |       17 |       28 |
| Albugo          |        2 |       11 |
| Aphanomyces     |        5 |       27 |
| Babesia         |       10 |       25 |
| Blastocystis    |        2 |       74 |
| Crithidia       |        8 |       11 |
| Cryptosporidium |       16 |       76 |
| Cyclospora      |        1 |       41 |
| Cyclotella      |        7 |       10 |
| Dictyostelium   |        9 |       11 |
| Eimeria         |       11 |       18 |
| Entamoeba       |        5 |       16 |
| Galdieria       |        4 |       15 |
| Giardia         |        2 |       38 |
| Globisporangium |       47 |       81 |
| Isotricha       |        3 |       13 |
| Leishmania      |       21 |      105 |
| Naegleria       |        3 |       14 |
| Nannochloropsis |        6 |       16 |
| Ochromonas      |        2 |       11 |
| Ophryoscolex    |        1 |       10 |
| Paramecium      |        9 |       16 |
| Peronospora     |        6 |       18 |
| Phytophthora    |       62 |      225 |
| Phytopythium    |       14 |       18 |
| Plasmodiophora  |        1 |       51 |
| Plasmodium      |       21 |      158 |
| Plasmopara      |        4 |       11 |
| Pythium         |       44 |       68 |
| Symbiodinium    |        6 |       11 |
| Thalassiosira   |       15 |       18 |
| Theileria       |        6 |       16 |
| Toxoplasma      |        1 |       29 |
| Trypanosoma     |       19 |       69 |

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
#1926

# N50 C S; create n50.tsv and n50.pass.tsv
# LEN_N50   N_CONTIG    LEN_SUM
bash ASSEMBLY/n50.sh 5000 30000 5000000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50,S" --max "C"
#N50_min S_min   C_max
#7707    6434485 29495

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "S:0.1,0.5" --quantile "N50:0.1,0.5"  --quantile "C:0.5,0.9"
#S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
#11412690.5      38915944        3313.5  63464.5 1790    29762

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

| #item            | fields | lines |
|------------------|-------:|------:|
| url.tsv          |      3 | 1,926 |
| check.lst        |      1 | 1,926 |
| collect.tsv      |     20 | 1,927 |
| n50.tsv          |      4 | 1,927 |
| n50.pass.tsv     |      4 | 1,611 |
| collect.pass.tsv |     23 | 1,611 |
| pass.lst         |      1 | 1,610 |
| omit.lst         |      1 | 1,424 |
| rep.lst          |      1 |   686 |
| sp.lst           |      1 |   140 |

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Protists/ \
    wangq@202.119.37.251:data/Protists

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Protists/ \
    wangq@58.213.64.36:data/Protists

# rsync -avP wangq@202.119.37.251:data/Protists/ ~/data/Protists

# rsync -avP -e 'ssh -p 8804' wangq@58.213.64.36:data/Protists/ ~/data/Protists

```

## BioSample

ENA's BioSample missed many strains, so NCBI's was used.

```shell
cd ~/data/Protists

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Protists.assembly.tsv \
    --bs

bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 10

datamash check < BioSample/biosample.tsv
#1897 lines, 92 fields

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

# Distances within species
bash MinHash/species.sh

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#112

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
#  592 summary/NR.lst
#  588 summary/redundant.lst

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
    -o minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

# png
nwr topo --bl minhash.condensed.newick | # remove comments
    nw_display -s -b 'visibility:hidden' -w 1200 -v 20 - |
    rsvg-convert -o Protists.minhash.png

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

# .lst and .count.tsv
bash Count/rank.sh

cat Count/order.count.tsv |
    tsv-filter -H --ge "3:10" |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:10" |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

# Can accept N_COUNT
bash Count/lineage.sh 10

cat Count/lineage.count.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| order            | #species | #strains |
|------------------|---------:|---------:|
| Albuginales      |        2 |       11 |
| Bacillariales    |       10 |       11 |
| Dictyosteliales  |       10 |       11 |
| Diplomonadida    |        4 |       34 |
| Eucoccidiorida   |       32 |      166 |
| Eustigmatales    |        6 |       18 |
| Galdieriales     |        4 |       11 |
| Haemosporida     |       23 |      145 |
| Longamoebia      |        7 |       14 |
| Mastigamoebida   |        6 |       17 |
| Opalinata        |        3 |       43 |
| Peniculida       |        9 |       16 |
| Peronosporales   |       76 |      255 |
| Piroplasmida     |       16 |       37 |
| Plasmodiophorida |        3 |       53 |
| Pythiales        |      106 |      148 |
| Saprolegniales   |        9 |       27 |
| Stephanodiscales |       10 |       12 |
| Suessiales       |        8 |       10 |
| Thalassiosirales |       13 |       15 |
| Trypanosomatida  |       76 |      200 |
| Vestibuliferida  |        3 |       10 |

| genus           | #species | #strains |
|-----------------|---------:|---------:|
| Acanthamoeba    |        6 |       11 |
| Albugo          |        2 |       11 |
| Aphanomyces     |        5 |       23 |
| Babesia         |       10 |       24 |
| Blastocystis    |        2 |       42 |
| Crithidia       |        8 |       11 |
| Cryptosporidium |       16 |       73 |
| Cyclospora      |        1 |       41 |
| Eimeria         |        9 |       15 |
| Entamoeba       |        5 |       16 |
| Galdieria       |        4 |       11 |
| Giardia         |        2 |       32 |
| Globisporangium |       46 |       67 |
| Isotricha       |        3 |       10 |
| Leishmania      |       21 |       91 |
| Naegleria       |        3 |       14 |
| Nannochloropsis |        5 |       15 |
| Paramecium      |        9 |       16 |
| Peronospora     |        6 |       18 |
| Phytophthora    |       58 |      217 |
| Phytopythium    |       14 |       17 |
| Plasmodiophora  |        1 |       51 |
| Plasmodium      |       21 |      143 |
| Pythium         |       38 |       54 |
| Theileria       |        5 |       12 |
| Toxoplasma      |        1 |       28 |
| Trypanosoma     |       19 |       58 |

| #family           | genus           | species                    | count |
|-------------------|-----------------|----------------------------|------:|
| Albuginaceae      | Albugo          | Albugo candida             |    10 |
| Blastocystidae    | Blastocystis    | Blastocystis sp.           |    40 |
| Cryptosporidiidae | Cryptosporidium | Cryptosporidium hominis    |    15 |
|                   |                 | Cryptosporidium parvum     |    25 |
| Eimeriidae        | Cyclospora      | Cyclospora cayetanensis    |    41 |
| Entamoebidae      | Entamoeba       | Entamoeba histolytica      |    12 |
| Hexamitidae       | Giardia         | Giardia intestinalis       |    31 |
| Peronosporaceae   | Phytophthora    | Phytophthora cactorum      |    28 |
|                   |                 | Phytophthora capsici       |    16 |
|                   |                 | Phytophthora fragariae     |    13 |
|                   |                 | Phytophthora kernoviae     |    12 |
|                   |                 | Phytophthora nicotianae    |    13 |
|                   |                 | Phytophthora ramorum       |    28 |
|                   |                 | Phytophthora sojae         |    11 |
| Plasmodiidae      | Plasmodium      | Plasmodium falciparum      |    57 |
|                   |                 | Plasmodium vinckei         |    10 |
|                   |                 | Plasmodium vivax           |    18 |
|                   |                 | Plasmodium yoelii          |    11 |
| Plasmodiophoridae | Plasmodiophora  | Plasmodiophora brassicae   |    51 |
| Pythiaceae        | Globisporangium | Globisporangium irregulare |    15 |
|                   | Pythium         | Pythium insidiosum         |    11 |
| Saprolegniaceae   | Aphanomyces     | Aphanomyces astaci         |    10 |
| Sarcocystidae     | Toxoplasma      | Toxoplasma gondii          |    28 |
| Trypanosomatidae  | Leishmania      | Leishmania donovani        |    10 |
|                   |                 | Leishmania infantum        |    31 |
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

# strains.taxon.tsv
bash Count/strains.sh

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:10" |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv

```

| genus           | #species | #strains |
|-----------------|---------:|---------:|
| Aphanomyces     |        5 |       21 |
| Cryptosporidium |       13 |       22 |
| Leishmania      |       11 |       19 |
| Peronospora     |        5 |       10 |
| Phytophthora    |       19 |       72 |
| Plasmodium      |       20 |       78 |
| Toxoplasma      |        1 |       14 |
| Trypanosoma     |       10 |       22 |

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

