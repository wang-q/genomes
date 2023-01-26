# Various genera from Protists

Eukaryota without Plants, Metazoa, and Fungi

[TOC levels=1-3]: # ""

- [Aligning various genera from Protists](#aligning-various-genera-from-protists)
    - [Strain info](#strain-info)
    - [NCBI Assembly](#ncbi-assembly)
    - [Count strains](#count-strains)
    - [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash)
    - [Groups and targets](#groups-and-targets)
    - [Protists: prepare](#protists-prepare)
    - [plasmodium: run](#plasmodium-run)
    - [Protists: run](#protists-run)

## Strain info

### Species with assemblies

```shell
mkdir -p ~/data/Protists/summary
cd ~/data/Protists/summary

# should have a valid name of genus
nwr member Eukaryota -r genus |
    grep -v -i "Candidatus " |
    grep -v -i "candidate " |
    nwr restrict -e Viridiplantae |
    nwr restrict -e Metazoa |
    nwr restrict -e Fungi |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#3916 genus.list

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
    GB=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND genus_id = ${RANK_ID}
            " |
            sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
    )

    FULL=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND genus_id = ${RANK_ID}
                AND genome_rep IN ('Full')
            " |
            sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
    )

    if [[ ${GB} -gt 0 ]]; then
        echo -e "${RANK_ID}\t${GB}\t${FULL}"
    fi
done |
    nwr append stdin |
    tsv-select -f 1,4,2-3 |
    tsv-sort -k3,3nr -k4,4nr -k2,2 |
    (echo -e '#tax_id\tgenus\tGB\tFULL' && cat) \
    > genus.count.tsv

wc -l genus.count.tsv
#255 genus.count.tsv

cat genus.count.tsv |
    tsv-filter -H --ge GB:5 --ge FULL:1 |
    mlr --itsv --omd cat

cat genus.count.tsv |
    tsv-filter -H --ge GB:5 --ge FULL:1 |
    nwr append stdin -r phylum -r class |
    keep-header -- sort -k5,5 -k6,6

```

| #tax_id | genus            | GB  | FULL |
|---------|------------------|-----|------|
| 4783    | Phytophthora     | 189 | 189  |
| 5820    | Plasmodium       | 152 | 151  |
| 5658    | Leishmania       | 72  | 66   |
| 4797    | Pythium          | 68  | 68   |
| 5806    | Cryptosporidium  | 65  | 64   |
| 5690    | Trypanosoma      | 57  | 55   |
| 1448052 | Globisporangium  | 53  | 53   |
| 37359   | Plasmodiophora   | 50  | 49   |
| 44417   | Cyclospora       | 40  | 40   |
| 5810    | Toxoplasma       | 28  | 28   |
| 100860  | Aphanomyces      | 27  | 27   |
| 5754    | Acanthamoeba     | 25  | 25   |
| 12967   | Blastocystis     | 23  | 23   |
| 5740    | Giardia          | 22  | 22   |
| 5864    | Babesia          | 21  | 21   |
| 795339  | Phytopythium     | 18  | 18   |
| 5758    | Entamoeba        | 16  | 16   |
| 5884    | Paramecium       | 16  | 16   |
| 83373   | Galdieria        | 15  | 15   |
| 5748    | Nannochloropsis  | 15  | 15   |
| 5873    | Theileria        | 15  | 14   |
| 5761    | Naegleria        | 14  | 14   |
| 5800    | Eimeria          | 14  | 13   |
| 5986    | Isotricha        | 13  | 13   |
| 2985    | Ochromonas       | 11  | 11   |
| 70742   | Peronospora      | 11  | 11   |
| 4780    | Plasmopara       | 11  | 11   |
| 2949    | Symbiodinium     | 11  | 11   |
| 65356   | Albugo           | 10  | 10   |
| 5782    | Dictyostelium    | 10  | 10   |
| 47895   | Ophryoscolex     | 10  | 10   |
| 5655    | Crithidia        | 8   | 8    |
| 40635   | Entodinium       | 8   | 8    |
| 5935    | Euplotes         | 8   | 8    |
| 28000   | Perkinsus        | 8   | 8    |
| 47893   | Polyplastron     | 8   | 8    |
| 1003337 | Angomonas        | 6   | 6    |
| 47888   | Diplodinium      | 6   | 6    |
| 40637   | Epidinium        | 6   | 6    |
| 5890    | Tetrahymena      | 6   | 6    |
| 5721    | Trichomonas      | 6   | 6    |
| 33652   | Cafeteria        | 5   | 5    |
| 40805   | Dasytricha       | 5   | 5    |
| 1448050 | Elongisporangium | 5   | 5    |
| 184462  | Hyaloperonospora | 5   | 5    |
| 35127   | Thalassiosira    | 5   | 5    |

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
          * Aphanomyces   丝囊霉属
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

## Download all assemblies

### Create assembly.tsv


```shell
cd ~/data/Protists/summary

GENUS=$(
    cat genus.count.tsv |
        grep -v "^#" |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-filter --regex '2:^[A-Za-z ]+$' \
    > raw.tsv

# Create abbr.
cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[5]}++; # abbr_name
        $seen{$F[5]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $F[5], $F[3], $F[1], $F[4];
        ' |
    keep-header -- sort -k3,3 -k1,1 \
    > Protists.assembly.tsv

datamash check < Protists.assembly.tsv
#1408 lines, 4 fields

# find potential duplicate strains or assemblies
cat Protists.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Protists.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Protists.assembly.tsv
# cp Protists.assembly.tsv ~/Scripts/genomes/assembly

# Comment out unneeded strains

# Cleaning
rm raw*.*sv

```

### rsync and check

```shell
cd ~/data/Protists

cat ~/Scripts/genomes/assembly/Protists.assembly.tsv |
    tsv-filter -v --str-in-fld 2:http

nwr assembly ~/Scripts/genomes/assembly/Protists.assembly.tsv \
    -o ASSEMBLY

# Remove dirs not in the list
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    tr "/" "\t" |
    cut -f 2 |
    tsv-join --exclude -k 1 -f ASSEMBLY/url.tsv -d 1 |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        echo Remove {}
        rm -fr ASSEMBLY/{}
    '

# Run
proxychains4 bash ASSEMBLY/rsync.sh

# md5
# rm ASSEMBLY/check.list
bash ASSEMBLY/check.sh

# collect
bash ASSEMBLY/collect.sh

```

## Count strains

```bash
cd ~/data/alignment/Protists

for dir in $(find ASSEMBLY -maxdepth 1 -mindepth 1 -type d | sort); do
    1>&2 echo "==> ${dir}"
    name=$(basename ${dir})

    find ${dir} -type f -name "*_genomic.fna.gz" |
        grep -v "_from_" | # exclude CDS and rna
        xargs cat |
        faops n50 -C -S stdin |
        (echo -e "name\t${name}" && cat) |
        datamash transpose
done |
    tsv-uniq |
    tee ASSEMBLY/n50.tsv

cat ASSEMBLY/n50.tsv |
    tsv-filter \
        -H --or \
        --le 4:3000 \
        --ge 2:100000 |
    tsv-filter -H --ge 3:1000000 |
    tr "\t" "," \
    > ASSEMBLY/n50.pass.csv

wc -l ASSEMBLY/n50*
#  204 ASSEMBLY/n50.pass.csv
#  296 ASSEMBLY/n50.tsv

tsv-join \
    ASSEMBLY/Protists.assembly.collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv \
    > ASSEMBLY/Protists.assembly.pass.csv

wc -l ASSEMBLY/Protists.assembly*csv
#   293 ASSEMBLY/Protists.assembly.collect.csv
#   201 ASSEMBLY/Protists.assembly.pass.csv

# find potential duplicated strains names
cat ASSEMBLY/Protists.assembly.pass.csv |
    cut -d, -f 7 |
    sort |
    uniq -c |
    sort -nr

```

```bash
cd ~/data/alignment/Protists

parallel --no-run-if-empty --linebuffer -k -j 4 '
    n_species=$(cat ASSEMBLY/Protists.assembly.pass.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        uniq |
        wc -l)

    n_strains=$(cat ASSEMBLY/Protists.assembly.pass.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        wc -l)

    printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' ::: $(
        cat ASSEMBLY/Protists.assembly.pass.csv |
            sed -e '1d' |
            cut -d"," -f 2 |
            grep -v "Candidatus" |
            cut -d" " -f 1 |
            sort |
            uniq
    )

#Acanthamoeba    1       2
#Babesia 6       7
#Blastocystis    2       8
#Crithidia       2       2
#Cryptosporidium 13      16
#Dictyostelium   4       4
#Eimeria 2       2
#Entamoeba       2       3
#Giardia 3       6
#Leishmania      17      25
#Nannochloropsis 4       8
#Phytophthora    11      14
#Plasmodium      20      60
#Pythium 4       4
#Theileria       4       6
#Toxoplasma      1       16
#Trypanosoma     8       13

```

## Raw phylogenetic tree by MinHash

```bash
mkdir -p ~/data/alignment/Protists/mash
cd ~/data/alignment/Protists/mash

for name in $(cat ../ASSEMBLY/Protists.assembly.pass.csv | sed -e '1d' | cut -d"," -f 1 ); do
    2>&1 echo "==> ${name}"

    if [[ -e ${name}.msh ]]; then
        continue
    fi

    find ../ASSEMBLY/${name} -name "*.fsa_nt.gz" -or -name "*_genomic.fna.gz" |
        grep -v "_from_" |
        xargs cat |
        mash sketch -k 21 -s 100000 -p 8 - -I "${name}" -o ${name}
done

mash triangle -E -p 8 -l <(
    cat ../ASSEMBLY/Protists.assembly.pass.csv | sed -e '1d' | cut -d"," -f 1 | parallel echo "{}.msh"
    ) \
    > dist.tsv

# fill matrix with lower triangle
tsv-select -f 1-3 dist.tsv |
    (tsv-select -f 2,1,3 dist.tsv && cat) |
    (
        cut -f 1 dist.tsv |
            tsv-uniq |
            parallel -j 1 --keep-order 'echo -e "{}\t{}\t0"' &&
        cat
    ) \
    > dist_full.tsv

cat dist_full.tsv |
    Rscript -e '
        library(readr);
        library(tidyr);
        library(ape);
        pair_dist <- read_tsv(file("stdin"), col_names=F);
        tmp <- pair_dist %>%
            pivot_wider( names_from = X2, values_from = X3, values_fill = list(X3 = 1.0) )
        tmp <- as.matrix(tmp)
        mat <- tmp[,-1]
        rownames(mat) <- tmp[,1]

        dist_mat <- as.dist(mat)
        clusters <- hclust(dist_mat, method = "ward.D2")
        tree <- as.phylo(clusters)
        write.tree(phy=tree, file="tree.nwk")

        group <- cutree(clusters, h=0.5) # k=3
        groups <- as.data.frame(group)
        groups$ids <- rownames(groups)
        rownames(groups) <- NULL
        groups <- groups[order(groups$group), ]
        write_tsv(groups, "groups.tsv")
    '

nw_display -s -b 'visibility:hidden' -w 600 -v 30 tree.nwk |
    rsvg-convert -o ~/Scripts/withncbi/image/Protists.png

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

```bash
mkdir -p ~/data/alignment/Protists/taxon
cd ~/data/alignment/Protists/taxon

cp ../mash/tree.nwk .

# manually combine Ba_mic
# manually remove bad assemblies
cat ../mash/groups.tsv |
    grep -v "Ba_mic" |
    grep -v "Cry_bai" |
    grep -v "En_inv_IP1" |
    grep -v "L_sp_A" \
    > groups.tsv
echo -e "2\tBa_mic" >> groups.tsv
echo -e "2\tBa_mic_RI" >> groups.tsv
#echo -e "11\tCr_win_CBS_7118" >> groups.tsv
#echo -e "11\tCr_dep_CBS_7841" >> groups.tsv
#echo -e "11\tCr_dep_CBS_7855" >> groups.tsv

ARRAY=(
    'A_Ei::A_cas'
    'Babesia::Ba_bov_T2Bo'
    'Blastocystis::Bl_hom'
    'Cryptosporidium::Cry_parvum_Iowa_II'
    'Dictyostelium::D_disc_AX4'
    'Giardia::G_intes'
    'Leishmania::L_maj_Friedlin'
    'Nannochloropsis::N_oce'
    'Phytophthora::Ph_soj'
    'Pl_ber_cha_vin_yoe::Pl_yoe'
    'Pl_kno_viv::Pl_viv'
    'Pl_falcip::Pl_falcip_3D7'
    'Pythium::Py_gui'
    'Theileria::Th_parva_Muguga'
    'Toxoplasma::To_gondii_ME49'
    'Trypanosoma::Tr_bruc_brucei_TREU927'
)

echo -e "#Serial\tGroup\tCount\tTarget\tSequencing" > group_target.tsv

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"

    SERIAL=$(
        cat ../mash/groups.tsv |
            tsv-filter --str-eq 2:${TARGET_NAME} |
            tsv-select -f 1
    )

    cat groups.tsv |
        tsv-filter --str-eq 1:${SERIAL} |
        tsv-select -f 2 \
        > ${GROUP_NAME}

    COUNT=$(cat ${GROUP_NAME} | wc -l )

    echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}\t" >> group_target.tsv

done

mlr --itsv --omd cat group_target.tsv

cat <<'EOF' > chr-level.list
Ba_big
Ba_bov_T2Bo
Ba_mic_RI
Cry_parvum_Iowa_II
D_disc_AX4
L_braz_MHOM_BR_75_M2904
L_don
L_infa_JPCM5
L_maj_Friedlin
L_mex_MHOM_GT_2001_U1103
Pl_falcip_3D7
Pl_kno_H
Pl_viv
Th_ann
Th_equi_WA
Th_ori_Shintoku
Th_parva
Th_parva_Muguga
Tr_bruc_brucei_TREU927
Tr_bruc_gambiense_DAL972
EOF

```

## Protists: prepare

* Rsync to hpcc

```bash
rsync -avP \
    ~/data/alignment/Protists/ \
    wangq@202.119.37.251:data/alignment/Protists

# rsync -avP wangq@202.119.37.251:data/alignment/Protists/ ~/data/alignment/Protists

```

* `--perseq` for Chromosome-level assemblies and targets
* Use `Stramenopiles` as `Eukaryota` takes too long to compute

```bash
cd ~/data/alignment/Protists/

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

# prep
egaz template \
    ASSEMBLY \
    --prep -o GENOMES \
    $( cat taxon/group_target.tsv | sed -e '1d' | cut -f 4 | parallel -j 1 echo " --perseq {} " ) \
    $( cat taxon/chr-level.list | parallel -j 1 echo " --perseq {} " ) \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--species Stramenopiles --parallel 24"

bsub -q mpi -n 24 -J "Protists-0_prep" "bash GENOMES/0_prep.sh"

ls -t output.* | head -n 1 | xargs tail -f | grep "==>"

# gff
for n in $(cat taxon/group_target.tsv | sed -e '1d' | cut -f 4 ) \
    $( cat taxon/chr-level.list ) \
    ; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"

    gzip -d -c ${FILE_GFF} > GENOMES/${n}/chr.gff
done

```

## plasmodium: run

Plasmodium distributed on many branches.

```bash
cd ~/data/alignment/Protists/

# sanger
egaz template \
    GENOMES/Pl_falcip_3D7 \
    GENOMES/Pl_yoe \
    GENOMES/Pl_viv \
    GENOMES/Pl_ber_ANKA \
    GENOMES/Pl_cha_chabaudi \
    GENOMES/Pl_cyn_B \
    GENOMES/Pl_kno_H \
    GENOMES/Pl_ovale \
    GENOMES/Pl_rel \
    --multi -o groups/plasmodium/ \
    --multiname sanger \
    --tree taxon/tree.nwk \
    --parallel 24 -v

bsub -q mpi -n 24 -J "plasmodium-1_pair" "bash groups/plasmodium/1_pair.sh"
bsub  -w "ended(plasmodium-1_pair)" \
    -q mpi -n 24 -J "plasmodium-3_multi" "bash groups/plasmodium/3_multi.sh"

```

## Protists: run

```bash
cd ~/data/alignment/Protists/

cat taxon/group_target.tsv |
    sed -e '1d' |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            --multi -o groups/{2}/ \
            --tree taxon/tree.nwk \
            --parallel 24 -v

        bsub -q mpi -n 24 -J "{2}-1_pair" "bash groups/{2}/1_pair.sh"
        bsub -w "ended({2}-1_pair)" \
            -q mpi -n 24 -J "{2}-3_multi" "bash groups/{2}/3_multi.sh"
    '

# clean
find groups -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find groups -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type f -name "output.*" | parallel -r rm

```

