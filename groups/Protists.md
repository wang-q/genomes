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
    grep -v -i "Candidatus " |
    grep -v -i "candidate " |
    nwr restrict -e Viridiplantae |
    nwr restrict -e Metazoa |
    nwr restrict -e Fungi |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#3966 genus.list

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
#   94 RS1.tsv
#  591 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     95
#GB1     1437

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
        AND genus IN ('Caenorhabditis', 'Arabidopsis', 'Saccharomyces')
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
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

cat raw.tsv |
    tsv-uniq |
    datamash check
#1442 lines, 7 fields

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
#1439 lines, 5 fields

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

# genus.lst and genus.count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    keep-header -- tsv-sort -k1,1 |
    tsv-filter -H --ge 3:5 |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Acanthamoeba     |       16 |       25 |
| Albugo           |        2 |       10 |
| Angomonas        |        3 |        6 |
| Aphanomyces      |        5 |       27 |
| Babesia          |        8 |       21 |
| Cafeteria        |        2 |        5 |
| Crithidia        |        5 |        8 |
| Cryptosporidium  |       14 |       63 |
| Cyclospora       |        1 |       40 |
| Dasytricha       |        1 |        5 |
| Dictyostelium    |        9 |       10 |
| Diplodinium      |        2 |        6 |
| Eimeria          |        9 |       10 |
| Elongisporangium |        5 |        5 |
| Entamoeba        |        5 |       16 |
| Entodinium       |        3 |        8 |
| Epidinium        |        2 |        6 |
| Euplotes         |        5 |        7 |
| Galdieria        |        3 |       15 |
| Giardia          |        2 |       32 |
| Globisporangium  |       47 |       53 |
| Hyaloperonospora |        3 |        7 |
| Isotricha        |        2 |        6 |
| Leishmania       |       20 |       62 |
| Naegleria        |        3 |       14 |
| Nannochloropsis  |        5 |       14 |
| Ophryoscolex     |        1 |       10 |
| Paramecium       |        9 |       16 |
| Perkinsus        |        3 |        7 |
| Peronospora      |        6 |       18 |
| Phytophthora     |       53 |      191 |
| Phytopythium     |       14 |       18 |
| Plasmodiophora   |        1 |       49 |
| Plasmodium       |       20 |      148 |
| Plasmopara       |        4 |       11 |
| Polyplastron     |        1 |        8 |
| Pythium          |       43 |       66 |
| Symbiodinium     |        5 |        7 |
| Tetrahymena      |        4 |        6 |
| Theileria        |        5 |       15 |
| Toxoplasma       |        1 |       28 |
| Trichomonas      |        4 |        7 |
| Trypanosoma      |       11 |       56 |

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

