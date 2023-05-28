# Fungi infecting the lungs

Genomic alignment of fungi infecting the lungs was performed to obtain consensus sequences.

<!-- toc -->

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
Blastomyces,芽生菌属 (lung)
Candida,念珠菌属
Candida/Metschnikowiaceae, ([Candida] auris 耳念珠菌)
Coccidioides,球孢子菌属 (lung)
#Colletotrichum,刺盘孢属 (炭疽)
#Epichloe,
#Fusarium,镰刀菌
Hanseniaspora,有孢汉逊酵母
Histoplasma,组织胞浆菌属 (lung)
#Kazachstania,
Metschnikowia,梅奇酵母属
Nakaseomyces, (Candida glabrata 光滑念珠菌)
#Ogataea,
Paracoccidioides,副球孢子菌属 (lung)
#Penicillium,青霉菌属
#Pichia,毕赤酵母属
Pneumocystis,肺孢子菌属 (lung)
#Pyricularia,梨孢属 (稻瘟病)
Sporothrix,孢子丝菌属 (skin)
Talaromyces,踝节菌属 (lung)
#Trichoderma,木霉属
Trichophyton,毛癣菌属
#Verticillium,轮枝菌属 (黄萎病)
Yarrowia,耶氏酵母
#Zymoseptoria,
## Basidiomycetes
Cryptococcus,隐球菌属 (脑膜炎)
Malassezia,马拉色菌属 (skin)
#Puccinia,柄锈菌属 (条锈病)
#Rhodotorula,红酵母属
#Ustilago,黑粉菌属
## Other
Mucor,毛霉菌属
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
| 5052    | Aspergillus               | 66       | 454      | 曲霉菌属                     |
| 229219  | Blastomyces               | 2        | 2        | 芽生菌属 (lung)              |
| 5475    | Candida                   | 51       | 137      | 念珠菌属                     |
| 2964429 | Candida/Metschnikowiaceae | 5        | 52       | ([Candida] auris 耳念珠菌)   |
| 5500    | Coccidioides              | 2        | 2        | 球孢子菌属 (lung)             |
| 5036    | Histoplasma               | 5        | 8        | 组织胞浆菌属 (lung)            |
| 27320   | Metschnikowia             | 6        | 53       | 梅奇酵母属                    |
| 374468  | Nakaseomyces              | 2        | 42       | (Candida glabrata 光滑念珠菌) |
| 38946   | Paracoccidioides          | 2        | 2        | 副球孢子菌属 (lung)            |
| 4753    | Pneumocystis              | 8        | 11       | 肺孢子菌属 (lung)             |
| 29907   | Sporothrix                | 2        | 2        | 孢子丝菌属 (skin)             |
| 5094    | Talaromyces               | 8        | 10       | 踝节菌属 (lung)              |
| 5550    | Trichophyton              | 12       | 25       | 毛癣菌属                     |
| 4951    | Yarrowia                  | 4        | 23       | 耶氏酵母                     |
| 5206    | Cryptococcus              | 9        | 35       | 隐球菌属 (脑膜炎)               |
| 55193   | Malassezia                | 7        | 16       | 马拉色菌属 (skin)             |
| 4830    | Mucor                     | 1        | 1        | 毛霉菌属                     |

| #species               | #strains | comment |
|------------------------|----------|---------|
| Candida albicans       | 53       | 白色念珠菌   |
| Candida parapsilosis   | 21       | 近平滑念珠菌  |
| Candida tropicalis     | 1        | 热带念珠菌   |
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
| Ajellomycetaceae   | Histoplasma               | Histoplasma capsulatum                | 7     |
| Arthrodermataceae  | Trichophyton              | Trichophyton rubrum                   | 23    |
| Aspergillaceae     | Aspergillus               | Aspergillus flavus                    | 136   |
|                    |                           | Aspergillus fumigatus                 | 74    |
|                    |                           | Aspergillus luchuensis                | 12    |
|                    |                           | Aspergillus niger                     | 94    |
|                    |                           | Aspergillus oryzae                    | 90    |
|                    | Monascus                  | Monascus purpureus                    | 6     |
|                    | Penicillium               | Penicillium chrysogenum               | 77    |
|                    |                           | Penicillium citrinum                  | 23    |
|                    |                           | Penicillium nalgiovense               | 29    |
|                    |                           | Penicillium oxalicum                  | 10    |
|                    |                           | Penicillium roqueforti                | 36    |
|                    |                           | Penicillium salamii                   | 21    |
| Cryptococcaceae    | Cryptococcus              | Cryptococcus neoformans               | 32    |
| Debaryomycetaceae  | Candida                   | Candida albicans                      | 53    |
|                    |                           | Candida orthopsilosis                 | 4     |
|                    |                           | Candida parapsilosis                  | 21    |
|                    | Debaryomyces              | Debaryomyces hansenii                 | 11    |
|                    | Scheffersomyces           | Scheffersomyces stipitis              | 2     |
| Dipodascaceae      | Geotrichum                | Geotrichum candidum                   | 17    |
|                    | Yarrowia                  | Yarrowia lipolytica                   | 23    |
| Malasseziaceae     | Malassezia                | Malassezia furfur                     | 2     |
|                    |                           | Malassezia restricta                  | 4     |
|                    |                           | Malassezia sympodialis                | 8     |
| Metschnikowiaceae  | Candida/Metschnikowiaceae | [Candida] auris                       | 46    |
|                    |                           | [Candida] intermedia                  | 3     |
|                    | Clavispora                | Clavispora lusitaniae                 | 32    |
| Onygenaceae        | Ophidiomyces              | Ophidiomyces ophidiicola              | 63    |
| Pneumocystidaceae  | Pneumocystis              | Pneumocystis canis                    | 3     |
|                    |                           | Pneumocystis jirovecii                | 3     |
| Saccharomycetaceae | Eremothecium              | Eremothecium gossypii                 | 2     |
|                    | Kluyveromyces             | Kluyveromyces lactis                  | 3     |
|                    |                           | Kluyveromyces marxianus               | 15    |
|                    | Nakaseomyces              | Nakaseomyces glabratus                | 42    |
|                    | Saccharomyces             | Saccharomyces boulardii (nom. inval.) | 7     |
|                    |                           | Saccharomyces cerevisiae              | 111   |
|                    |                           | Saccharomyces eubayanus               | 13    |
|                    |                           | Saccharomyces kudriavzevii            | 9     |
|                    |                           | Saccharomyces paradoxus               | 25    |
|                    |                           | Saccharomyces pastorianus             | 18    |
|                    |                           | Saccharomyces uvarum                  | 20    |
|                    | Torulaspora               | Torulaspora delbrueckii               | 23    |
| Saccharomycodaceae | Saccharomycodes           | Saccharomycodes ludwigii              | 3     |
| Trichocomaceae     | Talaromyces               | Talaromyces marneffei                 | 5     |

### List all ranks within the genus of interest

下面 5 个属是这次研究的主要目标

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

### List strains of the target genus and remove abnormal strains

If the maximum value of ANI between strains within a species was greater than 0.05, the median and
maximum ANI would be reported. Strains that cannot be linked by the median ANI will be considered as
abnormal strains.

```shell
cd ~/data/lung

ANI_VALUE_THRESHOLD=0.05

# Abnormal strains
cat summary/strains.taxon.tsv | tsv-select -f 3 | tsv-uniq | #head -n 1 |
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
        grep -v -F -w -f "NR/${SPECIES_}/median.cc.lst"
done |
    tee summary/abnormal.lst

cat summary/abnormal.lst
#Asp_flav_GCA_023653635_1
#Asp_nig_GCA_027923785_1
#De_han_GCA_016097515_1
#De_han_GCA_016097625_1
#Geot_candidum_GCA_013305805_1
#Penicillium_sal_GCA_911639275_1
#Yar_lip_GCA_003571375_1

# Recreate summary/strains.taxon.tsv to avoid these assemblies
cat ../Fungi/summary/collect.pass.csv |
    grep -F -f tmp.lst |
    tsv-select -d, -f 1,3 |
    tr "," "\t" |
    sort |
    tsv-uniq |
    grep -v -F -w -f summary/abnormal.lst |
    nwr append stdin -c 2 -r species -r genus -r family -r order \
    > summary/strains.taxon.tsv

# other lists
cat summary/strains.taxon.tsv | cut -f 1 | sort | uniq \
    > strains.lst

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
