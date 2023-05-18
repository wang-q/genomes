# Fungi infecting the lungs

Genomic alignment of fungi infecting the lungs was performed to obtain consensus sequences.

<!-- toc -->

## Strain info

* [Aspergillus](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=5052)
* [Candida](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=5475)
* [Cryptococcus](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=5206)
* [Pneumocystis](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4753)
* [Candida auris](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=498019)

### Infections

* https://patient.info/doctor/fungal-lung-infections

* https://www.cdc.gov/fungal/diseases/index.html

```shell
cd ~/data/Fungi

cat <<EOF | grep -v '^#' > summary/infections.csv
# Ascomycetes
Saccharomyces,酵母菌属
Aspergillus,曲霉菌属
Blastomyces,芽生菌属 (lung)
Candida,念珠菌属
Candida/Metschnikowiaceae,
Coccidioides,球孢子菌属 (lung)
Colletotrichum,炭疽菌属
Epichloe,
Fusarium,镰刀菌
Hanseniaspora,有孢汉逊酵母
Histoplasma,组织胞浆菌属 (lung)
Kazachstania,
Metschnikowia,梅奇酵母属
Ogataea,
Paracoccidioides,副球孢子菌属 (lung)
Penicillium,青霉菌属
Pichia,毕赤酵母属
Pneumocystis,肺孢子菌属 (lung)
Pyricularia,梨孢属
Sporothrix,孢子丝菌属 (skin)
Talaromyces,踝节菌属 (lung)
Trichoderma,木霉属
Trichophyton,毛癣菌属
Verticillium,轮枝菌属
Yarrowia,耶氏酵母
Zymoseptoria,
# Basidiomycetes
Cryptococcus,隐球菌属 (脑膜炎)
Malassezia,马拉色菌属
Puccinia,柄锈菌属
Rhodotorula,红酵母属
Ustilago,黑粉菌属
# Other
Mucor,毛霉菌属
EOF

cat summary/infections.csv |
    tr ',' '\t' |
    tsv-join -f summary/genus.count.tsv -k 2 -d 1 --append-fields 1,3,4 |
    tsv-select -f 3,1,4,5,2 |
    (echo -e '#tax_id\tgenus\t#species\t#strains\tcomment' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | genus                     | #species | #strains | comment       |
|---------|---------------------------|----------|----------|---------------|
| 4930    | Saccharomyces             | 121      | 203      | 酵母菌属          |
| 5052    | Aspergillus               | 66       | 454      | 曲霉菌属          |
| 229219  | Blastomyces               | 2        | 2        | 芽生菌属 (lung)   |
| 5475    | Candida                   | 51       | 137      | 念珠菌属          |
| 2964429 | Candida/Metschnikowiaceae | 5        | 52       |               |
| 5500    | Coccidioides              | 2        | 2        | 球孢子菌属 (lung)  |
| 5455    | Colletotrichum            | 14       | 14       | 炭疽菌属          |
| 5112    | Epichloe                  | 4        | 4        |               |
| 5506    | Fusarium                  | 40       | 224      | 镰刀菌           |
| 5036    | Histoplasma               | 5        | 8        | 组织胞浆菌属 (lung) |
| 71245   | Kazachstania              | 3        | 3        |               |
| 27320   | Metschnikowia             | 6        | 53       | 梅奇酵母属         |
| 461281  | Ogataea                   | 5        | 20       |               |
| 38946   | Paracoccidioides          | 2        | 2        | 副球孢子菌属 (lung) |
| 5073    | Penicillium               | 19       | 202      | 青霉菌属          |
| 4919    | Pichia                    | 2        | 23       | 毕赤酵母属         |
| 4753    | Pneumocystis              | 3        | 3        | 肺孢子菌属 (lung)  |
| 48558   | Pyricularia               | 8        | 172      | 梨孢属           |
| 29907   | Sporothrix                | 2        | 2        | 孢子丝菌属 (skin)  |
| 5094    | Talaromyces               | 8        | 10       | 踝节菌属 (lung)   |
| 5543    | Trichoderma               | 12       | 31       | 木霉属           |
| 5550    | Trichophyton              | 12       | 25       | 毛癣菌属          |
| 1036719 | Verticillium              | 7        | 20       | 轮枝菌属          |
| 4951    | Yarrowia                  | 4        | 23       | 耶氏酵母          |
| 1047167 | Zymoseptoria              | 7        | 21       |               |
| 5206    | Cryptococcus              | 9        | 35       | 隐球菌属 (脑膜炎)    |
| 55193   | Malassezia                | 7        | 16       | 马拉色菌属         |
| 5296    | Puccinia                  | 7        | 14       | 柄锈菌属          |
| 5533    | Rhodotorula               | 3        | 67       | 红酵母属          |
| 5269    | Ustilago                  | 3        | 36       | 黑粉菌属          |
| 4830    | Mucor                     | 1        | 1        | 毛霉菌属          |

### List all ranks

There are no noteworthy classification ranks other than species.

```shell
mkdir -p ~/data/lung
cd ~/data/lung

nwr member Aspergillus Candida "Candida/Metschnikowiaceae" Pneumocystis Cryptococcus |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

for N in Aspergillus Candida "Candida/Metschnikowiaceae" Pneumocystis Cryptococcus; do
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
| genus      | 5     |
| species    | 728   |
| subgenus   | 6     |
| strain     | 715   |
| no rank    | 6     |
| varietas   | 43    |
| subspecies | 1     |
| isolate    | 3     |

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
