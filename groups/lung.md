# Fungi infecting the lungs

Genomic alignment of fungi infecting the lungs was performed to obtain consensus sequences.

<!-- toc -->

- [Strain info](#strain-info)
    * [Symlink](#symlink)
    * [Infections](#infections)
    * [List strains within families of target genara](#list-strains-within-families-of-target-genara)
    * [List all ranks within the genus of interest](#list-all-ranks-within-the-genus-of-interest)
- [All assemblies](#all-assemblies)
    * [List strains of the target genus and remove abnormal strains](#list-strains-of-the-target-genus-and-remove-abnormal-strains)

<!-- tocstop -->

## Strain info

* [Aspergillus](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=5052)
* [Candida](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=5475)
* [Cryptococcus](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=5206)
* [Pneumocystis](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4753)
* [Candida auris](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=498019)

### Symlink

```shell
mkdir -p ~/data/lung
cd ~/data/lung

rm -fr ASSEMBLY
rm -fr NR
rm -fr STRAINS

ln -s ../Fungi/ASSEMBLY ASSEMBLY
ln -s ../Fungi/NR NR
ln -s ../Fungi/STRAINS STRAINS

```

### Infections

* https://patient.info/doctor/fungal-lung-infections

* https://www.cdc.gov/fungal/diseases/index.html

```shell
cd ~/data/lung/
mkdir -p summary

cat <<EOF | grep -v '^#' > summary/infections.csv
## Ascomycetes
#Saccharomyces,酵母菌属
Aspergillus,曲霉菌属
#Blastomyces,芽生菌属 (lung)
Candida,念珠菌属
Candida/Metschnikowiaceae, ([Candida] auris 耳念珠菌)
Coccidioides,球孢子菌属 (lung)
#Colletotrichum,刺盘孢属 (炭疽)
#Epichloe,
#Fusarium,镰刀菌
Hanseniaspora,有孢汉逊酵母
#Histoplasma,组织胞浆菌属 (lung)
#Kazachstania,
Metschnikowia,梅奇酵母属
Nakaseomyces, (Candida glabrata 光滑念珠菌)
#Ogataea,
Paracoccidioides,副球孢子菌属 (lung)
#Penicillium,青霉菌属
#Pichia,毕赤酵母属
Pneumocystis,肺孢子菌属 (lung)
#Pyricularia,梨孢属 (稻瘟病)
#Sporothrix,孢子丝菌属 (skin)
#Talaromyces,踝节菌属 (lung)
#Trichoderma,木霉属
#Trichophyton,毛癣菌属
#Verticillium,轮枝菌属 (黄萎病)
Yarrowia,耶氏酵母
#Zymoseptoria,
## Basidiomycetes
Cryptococcus,隐球菌属 (脑膜炎)
#Malassezia,马拉色菌属 (skin)
#Puccinia,柄锈菌属 (条锈病)
#Rhodotorula,红酵母属
#Ustilago,黑粉菌属
## Other
#Mucor,毛霉菌属
## Top
Candida albicans,白色念珠菌
Nakaseomyces glabratus,光滑念珠菌
Candida tropicalis,热带念珠菌
Candida parapsilosis,近平滑念珠菌
[Candida] auris,耳念珠菌
EOF

cat summary/infections.csv |
    tr ',' '\t' |
    tsv-join -f ~/data/Fungi/summary/genus.count.tsv -k 2 -d 1 --append-fields 1,3,4 |
    tsv-select -f 3,1,4,5,2 |
    (echo -e '#tax_id\tgenus\t#species\t#strains\tcomment' && cat) |
    mlr --itsv --omd cat

cat ~/data/Fungi/summary/strains.taxon.tsv |
    tsv-summarize -g 3 --count |
    tsv-join -k 1 -f <(
        cat summary/infections.csv |
            tr ',' '\t'
        ) \
        -a 2 |
    (echo -e '#species\t#strains\tcomment' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | genus                     | #species | #strains | comment                  |
|---------|---------------------------|----------|----------|--------------------------|
| 5052    | Aspergillus               | 70       | 468      | 曲霉菌属                     |
| 5475    | Candida                   | 55       | 150      | 念珠菌属                     |
| 2964429 | Candida/Metschnikowiaceae | 6        | 55       | ([Candida] auris 耳念珠菌)   |
| 5500    | Coccidioides              | 11       | 12       | 球孢子菌属 (lung)             |
| 27320   | Metschnikowia             | 10       | 61       | 梅奇酵母属                    |
| 374468  | Nakaseomyces              | 2        | 42       | (Candida glabrata 光滑念珠菌) |
| 38946   | Paracoccidioides          | 2        | 2        | 副球孢子菌属 (lung)            |
| 4753    | Pneumocystis              | 8        | 11       | 肺孢子菌属 (lung)             |
| 4951    | Yarrowia                  | 4        | 23       | 耶氏酵母                     |
| 5206    | Cryptococcus              | 27       | 55       | 隐球菌属 (脑膜炎)               |

| #species               | #strains | comment |
|------------------------|----------|---------|
| Candida albicans       | 53       | 白色念珠菌   |
| Candida parapsilosis   | 21       | 近平滑念珠菌  |
| Candida tropicalis     | 9        | 热带念珠菌   |
| [Candida] auris        | 46       | 耳念珠菌    |
| Nakaseomyces glabratus | 42       | 光滑念珠菌   |

### List strains within families of target genara

```shell
cd ~/data/lung/

cat ../Fungi/summary/collect.pass.csv |
    sed -e '1d' |
    tr "," "\t" |
    tsv-select -f 1,3 |
    nwr append stdin -c 2 -r species -r genus -r family -r order |
    grep -Fw -f <( # Families of target genara
        cat summary/infections.csv |
            tr ',' '\t' |
            cut -f 1 |
            nwr append stdin -r family |
            cut -f 2 |
            grep -v "NA" |
            tsv-uniq
        ) \
    > summary/strains.taxon.tsv

cat summary/strains.taxon.tsv |
    tsv-select -f 5,4,3 |
    tsv-sort -k1,1 -k2,2 -k3,3 |
    tsv-summarize -g 1,2,3 --count |
    tsv-filter --gt 4:1 |
    perl -nla -F'\t' -e '
            BEGIN { our $family = q(); our $genus = q(); }

            # record the current family
            if ($F[0] eq $family) {
                printf qq(\t);
            } else {
                $family = $F[0];
                printf qq($family\t);
            }
            # record the current genus
            if ($F[1] eq $genus) {
                printf qq(\t);
            } else {
                $genus = $F[1];
                printf qq($genus\t);
            }

            print join qq(\t), ($F[2], $F[3]);
        ' |
    (echo -e '#family\tgenus\tspecies\tcount' && cat) |
    mlr --itsv --omd cat

```

| #family            | genus                     | species                               | count |
|--------------------|---------------------------|---------------------------------------|-------|
| Aspergillaceae     | Aspergillus               | Aspergillus chevalieri                | 4     |
|                    |                           | Aspergillus felis                     | 3     |
|                    |                           | Aspergillus flavus                    | 136   |
|                    |                           | Aspergillus fumigatus                 | 74    |
|                    |                           | Aspergillus luchuensis                | 12    |
|                    |                           | Aspergillus nidulans                  | 3     |
|                    |                           | Aspergillus niger                     | 94    |
|                    |                           | Aspergillus oryzae                    | 90    |
|                    |                           | Aspergillus puulaauensis              | 2     |
|                    |                           | Aspergillus sojae                     | 5     |
|                    | Monascus                  | Monascus purpureus                    | 6     |
|                    | Penicillium               | Penicillium chrysogenum               | 77    |
|                    |                           | Penicillium citrinum                  | 23    |
|                    |                           | Penicillium digitatum                 | 5     |
|                    |                           | Penicillium nalgiovense               | 29    |
|                    |                           | Penicillium oxalicum                  | 10    |
|                    |                           | Penicillium polonicum                 | 8     |
|                    |                           | Penicillium roqueforti                | 36    |
|                    |                           | Penicillium salamii                   | 21    |
| Cryptococcaceae    | Cryptococcus              | Cryptococcus floricola                | 2     |
|                    |                           | Cryptococcus gattii VGI               | 7     |
|                    |                           | Cryptococcus gattii VGII              | 10    |
|                    |                           | Cryptococcus neoformans               | 32    |
|                    |                           | Cryptococcus wingfieldii              | 2     |
| Debaryomycetaceae  | Candida                   | Candida albicans                      | 53    |
|                    |                           | Candida orthopsilosis                 | 4     |
|                    |                           | Candida parapsilosis                  | 21    |
|                    |                           | Candida tropicalis                    | 9     |
|                    | Debaryomyces              | Debaryomyces hansenii                 | 11    |
|                    | Hyphopichia               | Hyphopichia burtonii                  | 3     |
|                    | Scheffersomyces           | Scheffersomyces stipitis              | 2     |
| Dipodascaceae      | Geotrichum                | Geotrichum candidum                   | 17    |
|                    | Yarrowia                  | Yarrowia lipolytica                   | 23    |
| Metschnikowiaceae  | Candida/Metschnikowiaceae | [Candida] auris                       | 46    |
|                    |                           | [Candida] haemuloni                   | 4     |
|                    |                           | [Candida] intermedia                  | 3     |
|                    | Clavispora                | Clavispora lusitaniae                 | 32    |
|                    | Metschnikowia             | Metschnikowia reukaufii               | 2     |
|                    |                           | Metschnikowia zobellii                | 2     |
| Onygenaceae        | Coccidioides              | Coccidioides posadasii                | 11    |
|                    | Ophidiomyces              | Ophidiomyces ophidiicola              | 63    |
| Pneumocystidaceae  | Pneumocystis              | Pneumocystis canis                    | 3     |
|                    |                           | Pneumocystis jirovecii                | 3     |
| Saccharomycetaceae | Eremothecium              | Eremothecium gossypii                 | 2     |
|                    | Kluyveromyces             | Kluyveromyces lactis                  | 3     |
|                    |                           | Kluyveromyces marxianus               | 15    |
|                    | Lachancea                 | Lachancea fermentati                  | 2     |
|                    |                           | Lachancea thermotolerans              | 4     |
|                    | Nakaseomyces              | Nakaseomyces glabratus                | 42    |
|                    | Naumovozyma               | Naumovozyma castellii                 | 7     |
|                    | Saccharomyces             | Saccharomyces arboricola              | 2     |
|                    |                           | Saccharomyces bayanus                 | 5     |
|                    |                           | Saccharomyces boulardii (nom. inval.) | 7     |
|                    |                           | Saccharomyces cerevisiae              | 111   |
|                    |                           | Saccharomyces eubayanus               | 13    |
|                    |                           | Saccharomyces kudriavzevii            | 9     |
|                    |                           | Saccharomyces mikatae                 | 3     |
|                    |                           | Saccharomyces paradoxus               | 25    |
|                    |                           | Saccharomyces pastorianus             | 18    |
|                    |                           | Saccharomyces uvarum                  | 20    |
|                    | Torulaspora               | Torulaspora delbrueckii               | 23    |
|                    | Zygosaccharomyces         | Zygosaccharomyces rouxii              | 4     |
|                    | Zygotorulaspora           | Zygotorulaspora mrakii                | 2     |
| Saccharomycodaceae | Saccharomycodes           | Saccharomycodes ludwigii              | 3     |

### List all ranks within the genus of interest

下面 6 个属是这次研究的主要目标

There are no noteworthy classification ranks other than species.

```shell

nwr member Aspergillus Candida "Candida/Metschnikowiaceae" Nakaseomyces Pneumocystis Cryptococcus |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

for N in Aspergillus Candida "Candida/Metschnikowiaceae" Nakaseomyces Pneumocystis Cryptococcus; do
    nwr lineage "${N}" |
        tsv-filter --str-ne 1:clade |
        tsv-filter --str-ne "1:no rank" |
        sed -n '/kingdom\tFungi/,$p'
done |
    tsv-uniq |
    sed -E "s/\b(genus)\b/*\1*/"| # Highlight genus
    (echo -e '#rank\tsci_name\ttax_id' && cat) |
    mlr --itsv --omd cat

```

| rank       | count |
|------------|-------|
| genus      | 6     |
| species    | 738   |
| subgenus   | 6     |
| strain     | 730   |
| no rank    | 7     |
| varietas   | 43    |
| subspecies | 1     |
| isolate    | 4     |

| #rank      | sci_name                  | tax_id  |
|------------|---------------------------|---------|
| kingdom    | Fungi                     | 4751    |
| subkingdom | Dikarya                   | 451864  |
| phylum     | Ascomycota                | 4890    |
| subphylum  | Pezizomycotina            | 147538  |
| class      | Eurotiomycetes            | 147545  |
| subclass   | Eurotiomycetidae          | 451871  |
| order      | Eurotiales                | 5042    |
| family     | Aspergillaceae            | 1131492 |
| *genus*    | Aspergillus               | 5052    |
| subphylum  | Saccharomycotina          | 147537  |
| class      | Saccharomycetes           | 4891    |
| order      | Saccharomycetales         | 4892    |
| family     | Debaryomycetaceae         | 766764  |
| *genus*    | Candida                   | 5475    |
| family     | Metschnikowiaceae         | 27319   |
| *genus*    | Candida/Metschnikowiaceae | 2964429 |
| family     | Saccharomycetaceae        | 4893    |
| *genus*    | Nakaseomyces              | 374468  |
| subphylum  | Taphrinomycotina          | 451866  |
| class      | Pneumocystidomycetes      | 147553  |
| order      | Pneumocystidales          | 37987   |
| family     | Pneumocystidaceae         | 44281   |
| *genus*    | Pneumocystis              | 4753    |
| phylum     | Basidiomycota             | 5204    |
| subphylum  | Agaricomycotina           | 5302    |
| class      | Tremellomycetes           | 155616  |
| order      | Tremellales               | 5234    |
| family     | Cryptococcaceae           | 1884633 |
| *genus*    | Cryptococcus              | 5206    |

## All assemblies

### List strains of the target families and remove abnormal strains

```shell
cd ~/data/lung

mkdir -p summary

# Target families
FAMILY=(
    # Ascomycota
    # Pezizomycotina
    Aspergillaceae

    # Saccharomycotina
    # All families inside Saccharomycetales...

    # Taphrinomycotina
    Pneumocystidaceae

    # Basidiomycota
    Cryptococcaceae
)
FAMILY+=( $(nwr member Saccharomycetales -r family | sed 1d | cut -f 2) )
#FAMILY=$(IFS=, ; echo "${FAMILY[*]}")

FAMILY_ID=$(
    for G in "${FAMILY[@]}"; do echo $G; done |
        nwr append stdin --id -r family |
        cut -f 3 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        DISTINCT assembly_accession
    FROM ar
    WHERE 1=1
        AND family_id IN ($FAMILY_ID)
        AND species NOT LIKE '% sp.%'
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    grep -v -i "symbiont " \
    > tmp.lst

echo "
    SELECT
        DISTINCT assembly_accession
    FROM ar
    WHERE 1=1
        AND family_id IN ($FAMILY_ID)
        AND species NOT LIKE '% sp.%'
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    grep -v -i "symbiont " \
    >> tmp.lst

# all.taxon.tsv
cat ../Fungi/summary/collect.pass.csv |
    grep -Fw -f tmp.lst |
    tsv-select -d, -f 1,3 |
    tr "," "\t" |
    sort |
    tsv-uniq |
    nwr append stdin -c 2 -r species -r genus -r family -r order \
    > summary/all.taxon.tsv

```

### Remove abnormal strains

If the maximum value of ANI between strains within a species was greater than 0.05, the median and
maximum ANI would be reported. Strains that cannot be linked by the median ANI will be considered as
abnormal strains.

```shell
cd ~/data/lung

ANI_VALUE_THRESHOLD=0.05

# Abnormal strains
cat summary/all.taxon.tsv | tsv-select -f 3 | tsv-uniq | #head -n 1 |
while read SPECIES; do
    SPECIES_=$(
        echo "${SPECIES}" |
            tr " " "_"
    )

    if [[ ! -e NR/${SPECIES_}/assembly.lst ]]; then
        continue
    fi

    # Number of assemblies >= 3
    N_ASM=$(
        cat NR/${SPECIES_}/assembly.lst | wc -l
    )
    if [[ $N_ASM -lt 3 ]]; then
        continue
    fi

    D_MAX=$(
        cat NR/${SPECIES_}/mash.dist.tsv |
            tsv-summarize --max 3
    )
    if (( $(echo "$D_MAX < $ANI_VALUE_THRESHOLD" | bc -l) )); then
        continue
    fi

    # Link assemblies with median ANI
    D_MEDIAN=$(
        cat NR/${SPECIES_}/mash.dist.tsv |
            tsv-filter --lt "3:$ANI_VALUE_THRESHOLD" |
            tsv-summarize --median 3
    )
    cat "NR/${SPECIES_}/mash.dist.tsv" |
        tsv-filter --ff-str-ne 1:2 --le "3:$D_MEDIAN" |
        perl -nla -F"\t" -MGraph::Undirected -e '
            BEGIN {
                our $g = Graph::Undirected->new;
            }

            $g->add_edge($F[0], $F[1]);

            END {
                for my $cc ( $g->connected_components ) {
                    print join qq{\n}, sort @{$cc};
                }
            }
        ' \
        > "NR/${SPECIES_}/median.cc.lst"

    1>&2 echo -e "==> ${SPECIES_}\t${D_MEDIAN}\t${D_MAX}"
    cat NR/${SPECIES_}/assembly.lst |
        grep -v -Fw -f "NR/${SPECIES_}/median.cc.lst"
done |
    tee summary/abnormal.lst

cat summary/abnormal.lst
#Asp_flav_GCA_023653635_1
#Asp_nig_GCA_027923785_1
#Candida_haemuloni_GCA_019332025_1
#Crypt_gattii_VGI_GCA_000786445_1
#Crypt_gattii_VGI_GCA_014964225_1
#De_han_GCA_016097515_1
#De_han_GCA_016097625_1
#Geot_candidum_GCA_013305805_1
#Mal_globosa_GCA_001264795_1
#Mal_globosa_GCA_001264815_1
#Penicillium_polo_GCA_003344595_1
#Penicillium_polo_GCA_015585905_1
#Penicillium_sal_GCA_911639275_1
#Pn_canis_GCA_017911135_1
#Saccharomyces_uva_GCA_013180255_1
#Saccharomyces_uva_GCA_918221285_1
#Saccharomyces_uva_GCA_918280685_1
#Yar_lip_GCA_003571375_1
#Zygos_rou_GCA_020521415_1

# Recreate summary/strains.taxon.tsv to avoid these assemblies
cat summary/all.taxon.tsv |
    cut -f 1,2 |
    grep -v -F -w -f summary/abnormal.lst |
    nwr append stdin -c 2 -r species -r genus -r family -r order \
    > summary/strains.taxon.tsv

wc -l summary/*.taxon.tsv
#  1557 summary/all.taxon.tsv
#  1540 summary/strains.taxon.tsv

# other lists
cat summary/strains.taxon.tsv | cut -f 1 | sort | uniq \
    > summary/strains.lst

cat summary/strains.taxon.tsv | tsv-select -f 4 | sort | uniq \
    > summary/genus.lst

cat ../Fungi/summary/collect.pass.csv |
    grep -F -w -f <(
        cat ~/Scripts/genomes/assembly/Fungi.reference.tsv |
            tsv-select -H -f assembly_accession |
            sed '1d'
    ) |
    cut -d, -f 1 \
    > summary/reference.lst

```

### Extract from `../Fungi`

```shell
cd ~/data/lung

head -n 1 ../Fungi/summary/collect.pass.csv \
    > summary/collect.pass.csv

cat ../Fungi/summary/collect.pass.csv |
    grep -Fw -f <(cat summary/strains.lst summary/reference.lst) \
    >> summary/collect.pass.csv

# biosample.tsv
cp ../Fungi/summary/attributes.lst summary/

head -n 1 ../Fungi/summary/biosample.tsv \
    > summary/biosample.tsv

cat ../Fungi/summary/biosample.tsv |
    grep -Fw -f <(
        cat summary/collect.pass.csv |
            tsv-select -H -d, -f BioSample |
            sort | uniq |
            grep "^SAM"
        ) \
    >> summary/biosample.tsv

# NR.lst and representative.lst
cat ../Fungi/summary/NR.lst |
    grep -Fw -f <(cat summary/strains.lst) \
    > summary/NR.lst

cat ../Fungi/summary/representative.lst |
    grep -Fw -f <(cat summary/strains.lst) \
    > summary/representative.lst

# All representative should be in NR
cat summary/representative.lst |
    grep -v -F -f summary/NR.lst

wc -l \
    summary/strains.taxon.tsv \
    summary/collect.pass.csv \
    summary/biosample.tsv \
    summary/NR.lst \
    summary/representative.lst
#   4182 summary/strains.taxon.tsv
#   4196 summary/collect.pass.csv
#   4203 summary/biosample.tsv
#   1100 summary/NR.lst
#    131 summary/representative.lst

```

### Count strains - Genus

```shell
cd ~/data/lung

cat summary/genus.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(
            cat summary/collect.pass.csv |
                sed "1d" |
                tsv-select -d, -f 3 |
                nwr append stdin -r genus -r species |
                grep -w {} |
                tsv-select -f 1,3 |
                tsv-uniq |
                wc -l
        )

        n_strains=$(
            cat summary/collect.pass.csv |
                sed "1d" |
                tsv-select -d, -f 3 |
                nwr append stdin -r genus |
                grep -w {} |
                wc -l
        )

        n_nr=$(
            cat summary/collect.pass.csv |
                grep -F -w -f summary/NR.lst |
                tsv-select -d, -f 3 |
                nwr append stdin -r genus |
                grep -w {} |
                wc -l
        )

        printf "%s\t%d\t%d\t%d\n" {} ${n_species} ${n_strains} ${n_nr}
    ' |
    nwr append stdin --id |
    tsv-select -f 6,5,2,3,4 |
    tsv-sort -k2,2 |
    tsv-filter --ge 4:5 |
    (echo -e '#tax_id\tgenus\t#species\t#strains\t#NR' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | genus                     | #species | #strains | #NR |
|---------|---------------------------|----------|----------|-----|
| 5052    | Aspergillus               | 70       | 466      | 56  |
| 13366   | Brettanomyces             | 3        | 7        | 5   |
| 5475    | Candida                   | 54       | 147      | 31  |
| 2964429 | Candida/Metschnikowiaceae | 6        | 54       | 8   |
| 36910   | Clavispora                | 3        | 33       | 7   |
| 5206    | Cryptococcus              | 26       | 53       | 10  |
| 604195  | Cyberlindnera             | 5        | 7        | 5   |
| 4958    | Debaryomyces              | 3        | 10       | 6   |
| 33170   | Eremothecium              | 5        | 5        | 1   |
| 43987   | Geotrichum                | 1        | 16       | 3   |
| 4910    | Kluyveromyces             | 4        | 18       | 8   |
| 460517  | Komagataella              | 4        | 128      | 2   |
| 490731  | Kwoniella                 | 5        | 5        | 0   |
| 300275  | Lachancea                 | 9        | 12       | 4   |
| 27320   | Metschnikowia             | 4        | 6        | 2   |
| 5097    | Monascus                  | 1        | 6        | 1   |
| 374468  | Nakaseomyces              | 2        | 42       | 6   |
| 278028  | Naumovozyma               | 4        | 8        | 1   |
| 461281  | Ogataea                   | 6        | 23       | 3   |
| 5073    | Penicillium               | 23       | 211      | 20  |
| 4919    | Pichia                    | 2        | 23       | 2   |
| 4753    | Pneumocystis              | 7        | 9        | 3   |
| 4930    | Saccharomyces             | 128      | 211      | 38  |
| 4943    | Saccharomycopsis          | 2        | 5        | 4   |
| 4948    | Torulaspora               | 3        | 25       | 5   |
| 599737  | Wickerhamomyces           | 3        | 7        | 6   |
| 4951    | Yarrowia                  | 4        | 22       | 1   |

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/lung/ \
    wangq@202.119.37.251:data/lung

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/lung/ \
    wangq@58.213.64.36:data/lung

```

## Collect proteins

### `all.pro.fa`

```shell script
cd ~/data/lung

mkdir -p PROTEINS

cat summary/strains.lst | grep -v -Fw -f ASSEMBLY/omit.lst |
while read STRAIN; do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz
done |
    pigz -p4 \
    > PROTEINS/all.pro.fa.gz

gzip -dcf PROTEINS/all.pro.fa.gz |
    perl -nl -e '
        BEGIN { our %seen; our $h; }

        if (/^>/) {
            $h = (split(" ", $_))[0];
            $seen{$h}++;
            $_ = $h;
        }
        print if $seen{$h} == 1;
    ' |
    pigz -p4 \
    > PROTEINS/all.uniq.fa.gz

# counting proteins
gzip -dcf PROTEINS/all.pro.fa.gz |
    grep "^>" |
    wc -l |
    numfmt --to=si
#4.4M

gzip -dcf PROTEINS/all.pro.fa.gz |
    grep "^>" |
    tsv-uniq |
    wc -l |
    numfmt --to=si
#4.4M

# annotations may be different
gzip -dcf PROTEINS/all.uniq.fa.gz |
    grep "^>" |
    wc -l |
    numfmt --to=si
#4.4M

```

### `all.replace.fa`

```shell
cd ~/data/lung

rm PROTEINS/all.strain.tsv PROTEINS/all.replace.fa.gz

cat summary/strains.lst | grep -v -Fw -f ASSEMBLY/omit.lst |
while read STRAIN; do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
        grep "^>" |
        cut -d" " -f 1 |
        sed "s/^>//" |
        STRAIN=${STRAIN} perl -nl -e '
            $n = $_;
            $s = $n;
            $s =~ s/\.\d+//;
            printf qq{%s\t%s_%s\t%s\n}, $n, $ENV{STRAIN}, $s, $ENV{STRAIN};
        ' \
    > PROTEINS/${STRAIN}.replace.tsv

    cut -f 2,3 PROTEINS/${STRAIN}.replace.tsv >> PROTEINS/all.strain.tsv

    faops replace -s \
        ASSEMBLY/${STRAIN}/*_protein.faa.gz \
        <(cut -f 1,2 PROTEINS/${STRAIN}.replace.tsv) \
        stdout |
        pigz -p4 \
        >> PROTEINS/all.replace.fa.gz

    rm PROTEINS/${STRAIN}.replace.tsv
done

gzip -dcf PROTEINS/all.replace.fa.gz |
    grep "^>" |
    wc -l |
    numfmt --to=si
#4.4M

(echo -e "#name\tstrain" && cat PROTEINS/all.strain.tsv)  \
    > temp &&
    mv temp PROTEINS/all.strain.tsv

faops size PROTEINS/all.replace.fa.gz > PROTEINS/all.replace.sizes

(echo -e "#name\tsize" && cat PROTEINS/all.replace.sizes) > PROTEINS/all.size.tsv

rm PROTEINS/all.replace.sizes

```

### `all.info.tsv`

```shell
cd ~/data/lung

cat summary/strains.lst | grep -v -Fw -f ASSEMBLY/omit.lst |
while read STRAIN; do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
        grep "^>" |
        sed "s/^>//" |
        perl -nl -e '/\[.+\[/ and s/\[/\(/; print' |
        perl -nl -e '/\].+\]/ and s/\]/\)/; print' |
        perl -nl -e 's/\s+\[.+?\]$//g; print' |
        perl -nl -e 's/MULTISPECIES: //g; print' |
        STRAIN=${STRAIN} perl -nl -e '
            /^(\w+)\.\d+\s+(.+)$/ or next;
            printf qq{%s_%s\t%s\n}, $ENV{STRAIN}, $1, $2;
        '
done \
    > PROTEINS/all.annotation.tsv

cat PROTEINS/all.annotation.tsv |
    wc -l |
    numfmt --to=si
#4.4M

(echo -e "#name\tannotation" && cat PROTEINS/all.annotation.tsv) \
    > temp &&
    mv temp PROTEINS/all.annotation.tsv

# check differences
cat PROTEINS/all.size.tsv |
    grep -v -F -f <(cut -f 1 PROTEINS/all.annotation.tsv)

tsv-join \
    PROTEINS/all.strain.tsv \
    --data-fields 1 \
    -f PROTEINS/all.size.tsv \
    --key-fields 1 \
    --append-fields 2 \
    > PROTEINS/all.strain_size.tsv

tsv-join \
    PROTEINS/all.strain_size.tsv \
    --data-fields 1 \
    -f PROTEINS/all.annotation.tsv \
    --key-fields 1 \
    --append-fields 2 \
    > PROTEINS/all.info.tsv

cat PROTEINS/all.info.tsv |
    wc -l |
    numfmt --to=si
#21M

```
