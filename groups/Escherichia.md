# *Escherichia*

<!-- TOC -->
* [*Escherichia*](#escherichia)
  * [Strain info](#strain-info)
    * [Symlink](#symlink)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
  * [All assemblies](#all-assemblies)
    * [Extract from `../Bacteria` and create assembly.tsv](#extract-from-bacteria-and-create-assemblytsv)
    * [Count `assembly.tsv`](#count-assemblytsv)
  * [MinHash](#minhash)
  * [Count valid species and strains](#count-valid-species-and-strains)
    * [For *genomic alignments* and *protein families*](#for-genomic-alignments-and-protein-families)
    * [Count strains - Genus](#count-strains---genus)
  * [Collect proteins](#collect-proteins)
<!-- TOC -->

## Strain info

* [Escherichia](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=561)
* [Shigella](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=620)
* [Salmonella](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=590)
* [Klebsiella](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=570)

### Symlink

```bash
mkdir -p ~/data/Escherichia
cd ~/data/Escherichia

rm -fr ASSEMBLY
rm -fr STRAINS

ln -s ../Bacteria/ASSEMBLY ASSEMBLY
ln -s ../Bacteria/STRAINS STRAINS

```

### List all ranks

```shell
mkdir -p ~/data/Escherichia
cd ~/data/Escherichia

nwr member Enterobacteriaceae -r genus |
    keep-header -- tsv-sort -k2,2 |
    rgr md stdin --num

nwr member \
    Enterobacteriaceae |
    tsv-summarize -H -g 3 --count |
    rgr md stdin --fmt

nwr member \
    Enterobacteriaceae \
    -r "species group" -r "species subgroup" |
    tsv-select -f 1-3 |
    keep-header -- tsv-sort -k3,3 -k2,2 |
    rgr md stdin --num

```

| #tax_id | sci_name                      | rank  | division |
|--------:|-------------------------------|-------|----------|
| 3142780 | Apirhabdus                    | genus | Bacteria |
|   58337 | Aranicola                     | genus | Bacteria |
| 1903434 | Atlantibacter                 | genus | Bacteria |
|  347014 | Averyella                     | genus | Bacteria |
|   82976 | Buttiauxella                  | genus | Bacteria |
| 1302410 | Candidatus Annandia           | genus | Bacteria |
| 2592573 | Candidatus Arocatia           | genus | Bacteria |
| 1485663 | Candidatus Aschnera           | genus | Bacteria |
|  203804 | Candidatus Blochmanniella     | genus | Bacteria |
|  690591 | Candidatus Curculioniphilus   | genus | Bacteria |
|  862793 | Candidatus Cuticobacterium    | genus | Bacteria |
| 3121261 | Candidatus Dasytiphilus       | genus | Bacteria |
| 1906657 | Candidatus Doolittlea         | genus | Bacteria |
| 1177217 | Candidatus Gillettellia       | genus | Bacteria |
| 1906661 | Candidatus Gullanella         | genus | Bacteria |
|  568987 | Candidatus Hamiltonella       | genus | Bacteria |
| 1449912 | Candidatus Hartigia           | genus | Bacteria |
| 1906659 | Candidatus Hoaglandella       | genus | Bacteria |
| 2592572 | Candidatus Ischnodemia        | genus | Bacteria |
|  409304 | Candidatus Ishikawella        | genus | Bacteria |
|  862775 | Candidatus Kleidoceria        | genus | Bacteria |
| 2137422 | Candidatus Kotejella          | genus | Bacteria |
| 2974993 | Candidatus Lightella          | genus | Bacteria |
|  400524 | Candidatus Macropleicola      | genus | Bacteria |
| 3343722 | Candidatus Malihini           | genus | Bacteria |
| 1906660 | Candidatus Mikella            | genus | Bacteria |
| 1048757 | Candidatus Moranella          | genus | Bacteria |
| 1410382 | Candidatus Parastrichiiphilus | genus | Bacteria |
|   57962 | Candidatus Phlomobacter       | genus | Bacteria |
| 1177214 | Candidatus Profftia           | genus | Bacteria |
| 1408812 | Candidatus Puchtella          | genus | Bacteria |
|  472825 | Candidatus Purcelliella       | genus | Bacteria |
|  568988 | Candidatus Regiella           | genus | Bacteria |
|  401618 | Candidatus Riesia             | genus | Bacteria |
| 1008944 | Candidatus Rohrkolberia       | genus | Bacteria |
| 1414835 | Candidatus Rosenkranzia       | genus | Bacteria |
| 1081630 | Candidatus Schneideria        | genus | Bacteria |
| 2608261 | Candidatus Stammera           | genus | Bacteria |
|  442155 | Candidatus Stammerula         | genus | Bacteria |
| 1682492 | Candidatus Tachikawaea        | genus | Bacteria |
| 1699619 | Candidatus Westeberhardia     | genus | Bacteria |
|  158483 | Cedecea                       | genus | Bacteria |
|     544 | Citrobacter                   | genus | Bacteria |
|  413496 | Cronobacter                   | genus | Bacteria |
| 3021684 | Dryocola                      | genus | Bacteria |
| 2678248 | Edaphovirga                   | genus | Bacteria |
|     547 | Enterobacter                  | genus | Bacteria |
| 2943317 | Entomohabitans                | genus | Bacteria |
|     561 | Escherichia                   | genus | Bacteria |
| 1649295 | Franconibacter                | genus | Bacteria |
| 2976153 | Huaxiibacter                  | genus | Bacteria |
| 2899543 | Intestinirhabdus              | genus | Bacteria |
| 2815296 | Jejubacter                    | genus | Bacteria |
|     570 | Klebsiella                    | genus | Bacteria |
|     579 | Kluyvera                      | genus | Bacteria |
| 1330547 | Kosakonia                     | genus | Bacteria |
|   83654 | Leclercia                     | genus | Bacteria |
| 1330545 | Lelliottia                    | genus | Bacteria |
|  451512 | Mangrovibacter                | genus | Bacteria |
|  447792 | Phytobacter                   | genus | Bacteria |
|     702 | Plesiomonas                   | genus | Bacteria |
| 1330546 | Pluralibacter                 | genus | Bacteria |
| 2994443 | Pseudenterobacter             | genus | Bacteria |
| 2055880 | Pseudescherichia              | genus | Bacteria |
| 1504576 | Pseudocitrobacter             | genus | Bacteria |
|  160674 | Raoultella                    | genus | Bacteria |
|     590 | Salmonella                    | genus | Bacteria |
| 2726810 | Scandinavium                  | genus | Bacteria |
|     620 | Shigella                      | genus | Bacteria |
| 1335483 | Shimwellia                    | genus | Bacteria |
| 1649298 | Siccibacter                   | genus | Bacteria |
| 3016631 | Silvania                      | genus | Bacteria |
| 2303321 | Superficieibacter             | genus | Bacteria |
|     801 | Symbiopectobacterium          | genus | Bacteria |
| 2918836 | Tenebrionibacter              | genus | Bacteria |
| 2943312 | Tenebrionicola                | genus | Bacteria |
|  158851 | Trabulsiella                  | genus | Bacteria |
|  158876 | Yokenella                     | genus | Bacteria |

| rank             |  count |
|------------------|-------:|
| family           |      1 |
| genus            |     78 |
| species          | 17,002 |
| no rank          |  1,748 |
| strain           |  6,212 |
| subspecies       |     27 |
| isolate          |     27 |
| forma specialis  |     23 |
| clade            |      5 |
| serotype         |    196 |
| varietas         |      1 |
| serogroup        |    110 |
| species group    |      2 |
| species subgroup |      6 |

| #tax_id | sci_name                             | rank             |
|--------:|--------------------------------------|------------------|
| 1344959 | Citrobacter freundii complex         | species group    |
|  354276 | Enterobacter cloacae complex         | species group    |
| 2302426 | Enterobacter cloacae complex clade K | species subgroup |
| 2302427 | Enterobacter cloacae complex clade L | species subgroup |
| 2302428 | Enterobacter cloacae complex clade N | species subgroup |
| 2302429 | Enterobacter cloacae complex clade O | species subgroup |
| 2302430 | Enterobacter cloacae complex clade P | species subgroup |
| 2302431 | Enterobacter cloacae complex clade S | species subgroup |

### Species with assemblies

```bash
mkdir -p ~/data/Escherichia/summary
cd ~/data/Escherichia/summary

SPECIES=$(
    nwr member \
        Enterobacteriaceae \
        -r species |
        sed '1d' |
        cut -f 1 |
        sort |
        uniq
)

for S in $SPECIES; do
    RS=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    CHR=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
                AND assembly_level IN ('Complete Genome', 'Chromosome')
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    if [[ ${RS} -gt 0 ]]; then
        echo -e "$S\t$RS\t$CHR"
    fi
done |
    nwr append stdin |
    tsv-select -f 1,4,2-3 |
    tsv-sort -k3,3nr -k4,4nr -k2,2 |
    (echo -e 'species_id\tspecies\tRS\tCHR' && cat) \
    > species.count.tsv

cat species.count.tsv |
    tsv-filter -H --ge CHR:20 |
    rgr md stdin --num

```

| species_id | species                    |    RS |  CHR |
|-----------:|----------------------------|------:|-----:|
|        562 | Escherichia coli           | 40138 | 3964 |
|        573 | Klebsiella pneumoniae      | 21886 | 2500 |
|      28901 | Salmonella enterica        | 14204 | 1923 |
|     158836 | Enterobacter hormaechei    |  3763 |  399 |
|        624 | Shigella sonnei            |  1594 |   54 |
|    1463165 | Klebsiella quasipneumoniae |  1283 |  175 |
|        546 | Citrobacter freundii       |  1107 |  188 |
|     244366 | Klebsiella variicola       |  1017 |  172 |
|        623 | Shigella flexneri          |   859 |  129 |
|    1134687 | Klebsiella michiganensis   |   633 |   83 |
|        548 | Klebsiella aerogenes       |   626 |   78 |
|      28141 | Cronobacter sakazakii      |   546 |   25 |
|    1812935 | Enterobacter roggenkampii  |   511 |   53 |
|        550 | Enterobacter cloacae       |   511 |   47 |
|      61645 | Enterobacter asburiae      |   471 |   54 |
|        571 | Klebsiella oxytoca         |   460 |   50 |
|     208224 | Enterobacter kobei         |   375 |   29 |
|     208962 | Escherichia albertii       |   276 |   77 |
|    1639133 | Citrobacter portucalensis  |   270 |   46 |
|     881260 | Enterobacter bugandensis   |   221 |   22 |
|    2058152 | Klebsiella grimontii       |   214 |   28 |
|      57706 | Citrobacter braakii        |   190 |   46 |
|        545 | Citrobacter koseri         |   188 |   32 |
|        564 | Escherichia fergusonii     |   187 |   74 |
|     299767 | Enterobacter ludwigii      |   186 |   39 |
|      54291 | Raoultella ornithinolytica |   184 |   46 |
|      83655 | Leclercia adecarboxylata   |    93 |   22 |
|        622 | Shigella dysenteriae       |    69 |   28 |

## All assemblies

### Extract from `../Bacteria` and create assembly.tsv

```bash
cd ~/data/Escherichia

mkdir -p summary

cat ../Bacteria/summary/collect.pass.tsv |
    tsv-join -H -d RefSeq_assembly_accession -f ~/Scripts/genomes/assembly/Bacteria.reference.tsv -k assembly_accession \
    > summary/collect.pass.tsv

cat ../Bacteria/summary/collect.pass.tsv | # 316617
    nwr restrict Enterobacteriaceae -f stdin -c 3 | # restrict to Enterobacteriaceae 91393
    sed '1d' |
    tsv-join -e -f ../Bacteria/ASSEMBLY/omit.lst -k 1 | # 91392
    tsv-join -e -f ../Bacteria/MinHash/abnormal.lst -k 1 | # 91334
    sort \
    >> summary/collect.pass.tsv

cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    tsv-join -H -f summary/collect.pass.tsv -k 1 \
    > summary/assembly.tsv

# biosample.tsv
cp ../Bacteria/summary/attributes.lst summary/

cat ../Bacteria/summary/biosample.tsv |
    grep -Fw -f <(cat summary/collect.pass.tsv | tsv-select -H -f BioSample | sort | uniq) \
    > summary/biosample.tsv

```

### Count `assembly.tsv`

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class

```bash
cd ~/data/Escherichia

nwr template summary/assembly.tsv \
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
    rgr md stdin --num

```

| item    | count |
|---------|------:|
| strain  | 91356 |
| species |   170 |
| genus   |    57 |
| family  |    18 |
| order   |    13 |
| class   |     6 |

| genus                     | #species | #strains |
|---------------------------|---------:|---------:|
| Acinetobacter             |        1 |        1 |
| Atlantibacter             |        2 |      116 |
| Bacillus                  |        3 |        3 |
| Bacteroides               |        1 |        1 |
| Bifidobacterium           |        1 |        1 |
| Brucella                  |        1 |        1 |
| Buttiauxella              |        5 |       16 |
| Campylobacter             |        1 |        1 |
| Candidatus Blochmanniella |        0 |        0 |
| Candidatus Hamiltonella   |        0 |        0 |
| Candidatus Moranella      |        0 |        0 |
| Candidatus Purcelliella   |        0 |        0 |
| Candidatus Regiella       |        0 |        0 |
| Candidatus Riesia         |        0 |        0 |
| Candidatus Schneideria    |        0 |        0 |
| Cedecea                   |        3 |       18 |
| Citrobacter               |       18 |     2032 |
| Cronobacter               |        7 |      744 |
| Dryocola                  |        1 |        3 |
| Enterobacter              |       27 |     6168 |
| Enterococcus              |        1 |        1 |
| Escherichia               |        6 |    39441 |
| Francisella               |        1 |        1 |
| Franconibacter            |        3 |       14 |
| Klebsiella                |       14 |    25695 |
| Kluyvera                  |        6 |       53 |
| Kosakonia                 |        7 |       73 |
| Leclercia                 |        3 |       88 |
| Lelliottia                |        4 |       41 |
| Listeria                  |        1 |        1 |
| Listeria                  |        1 |        1 |
| Listeria                  |        1 |        1 |
| Mangrovibacter            |        1 |        3 |
| Mycobacterium             |        1 |        1 |
| Mycobacteroides           |        1 |        1 |
| Phocaeicola               |        1 |        1 |
| Phytobacter               |        3 |       23 |
| Plesiomonas               |        1 |       48 |
| Pluralibacter             |        2 |       30 |
| Proteus                   |        1 |        1 |
| Pseudenterobacter         |        1 |        4 |
| Pseudescherichia          |        1 |       11 |
| Pseudocitrobacter         |        1 |        2 |
| Escherichia               |        1 |        1 |
| Raoultella                |        4 |      311 |
| Salmonella                |        2 |    14009 |
| Scandinavium              |        4 |        9 |
| Shigella                  |        4 |     2291 |
| Shimwellia                |        1 |        4 |
| Siccibacter               |        2 |       10 |
| Staphylococcus            |        1 |        1 |
| Streptococcus             |        2 |        2 |
| Superficieibacter         |        1 |        2 |
| Tenebrionicola            |        1 |        2 |
| Trabulsiella              |        1 |        8 |
| Vibrio                    |        1 |        1 |
| Xanthomonas               |        1 |        1 |
| Yersinia                  |        1 |        1 |
| Yokenella                 |        1 |       14 |

## MinHash

```bash
cd ~/data/Escherichia/

# relaxed thresholds
nwr template summary/assembly.tsv \
    --mh \
    --parallel 8 \
    --ani-ab 0.12 \
    --ani-nr 0.01

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
#   2217 summary/NR.lst
#  89115 summary/redundant.lst

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#2

```

## Count valid species and strains

### For *genomic alignments* and *protein families*

```bash
cd ~/data/Escherichia/

nwr template summary/assembly.tsv \
    --count \
    --not-in MinHash/abnormal.lst \
    --rank order --rank genus \
    --lineage genus

# strains.taxon.tsv
bash Count/strains.sh

# .lst and .count.tsv
bash Count/rank.sh

cat Count/order.count.tsv |
    tsv-filter -H --ge "3:2" |
    rgr md stdin --fmt

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:2" |
    rgr md stdin --fmt

# Can accept N_COUNT
bash Count/lineage.sh 50

cat Count/lineage.count.tsv |
    rgr md stdin --fmt

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| order            | #species | #strains |
|------------------|---------:|---------:|
| Bacillales       |        5 |        5 |
| Bacteroidales    |        2 |        2 |
| Enterobacterales |      150 |   91,334 |
| Lactobacillales  |        3 |        3 |
| Mycobacteriales  |        2 |        2 |

| genus             | #species | #strains |
|-------------------|---------:|---------:|
| Atlantibacter     |        2 |      116 |
| Bacillus          |        3 |        3 |
| Buttiauxella      |        5 |       16 |
| Cedecea           |        3 |       18 |
| Citrobacter       |       18 |    2,031 |
| Cronobacter       |        7 |      744 |
| Cronobacter       |        7 |      744 |
| Dryocola          |        1 |        3 |
| Enterobacter      |       27 |    6,167 |
| Escherichia       |        6 |   39,441 |
| Franconibacter    |        3 |       14 |
| Klebsiella        |       14 |   25,695 |
| Kluyvera          |        6 |       53 |
| Kosakonia         |        7 |       73 |
| Leclercia         |        3 |       88 |
| Lelliottia        |        4 |       41 |
| Mangrovibacter    |        1 |        3 |
| Phytobacter       |        3 |       23 |
| Plesiomonas       |        1 |       48 |
| Pluralibacter     |        2 |       30 |
| Pseudenterobacter |        1 |        4 |
| Pseudescherichia  |        1 |       11 |
| Pseudocitrobacter |        1 |        2 |
| Raoultella        |        4 |      311 |
| Salmonella        |        2 |   14,009 |
| Scandinavium      |        4 |        9 |
| Shigella          |        4 |    2,291 |
| Shimwellia        |        1 |        4 |
| Siccibacter       |        2 |       10 |
| Streptococcus     |        2 |        2 |
| Superficieibacter |        1 |        2 |
| Tenebrionicola    |        1 |        2 |
| Trabulsiella      |        1 |        8 |
| Yokenella         |        1 |       14 |

| #genus        | species                    |  count |
|---------------|----------------------------|-------:|
| Atlantibacter | Atlantibacter hermannii    |    104 |
| Citrobacter   | Citrobacter amalonaticus   |     77 |
|               | Citrobacter braakii        |    189 |
|               | Citrobacter freundii       |  1,083 |
|               | Citrobacter koseri         |    186 |
|               | Citrobacter portucalensis  |    268 |
|               | Citrobacter werkmanii      |     64 |
| Cronobacter   | Cronobacter dublinensis    |     66 |
|               | Cronobacter malonaticus    |     87 |
|               | Cronobacter sakazakii      |    538 |
| Enterobacter  | Enterobacter asburiae      |    468 |
|               | Enterobacter bugandensis   |    220 |
|               | Enterobacter cloacae       |    479 |
|               | Enterobacter hormaechei    |  3,697 |
|               | Enterobacter kobei         |    369 |
|               | Enterobacter ludwigii      |    180 |
|               | Enterobacter roggenkampii  |    506 |
| Escherichia   | Escherichia albertii       |    275 |
|               | Escherichia coli           | 38,841 |
|               | Escherichia fergusonii     |    186 |
|               | Escherichia marmotae       |    114 |
| Klebsiella    | Klebsiella aerogenes       |    604 |
|               | Klebsiella grimontii       |    210 |
|               | Klebsiella michiganensis   |    622 |
|               | Klebsiella oxytoca         |    447 |
|               | Klebsiella pasteurii       |     53 |
|               | Klebsiella pneumoniae      | 21,442 |
|               | Klebsiella quasipneumoniae |  1,262 |
|               | Klebsiella variicola       |    998 |
| Leclercia     | Leclercia adecarboxylata   |     80 |
| Raoultella    | Raoultella ornithinolytica |    179 |
|               | Raoultella planticola      |     79 |
|               | Raoultella terrigena       |     50 |
| Salmonella    | Salmonella enterica        | 13,978 |
| Shigella      | Shigella boydii            |     76 |
|               | Shigella flexneri          |    722 |
|               | Shigella sonnei            |  1,449 |

### Count strains - Genus

```bash
cd ~/data/Escherichia

cat Count/genus.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(
            cat summary/collect.pass.tsv |
                sed "1d" |
                tsv-select -f 3 |
                nwr append stdin -r genus -r species |
                grep -w {} |
                tsv-select -f 1,3 |
                rgr dedup stdin |
                wc -l
        )

        n_strains=$(
            cat summary/collect.pass.tsv |
                sed "1d" |
                tsv-select -f 3 |
                nwr append stdin -r genus |
                grep -w {} |
                wc -l
        )

        n_nr=$(
            cat summary/collect.pass.tsv |
                grep -Fw -f summary/NR.lst |
                tsv-select -f 3 |
                nwr append stdin -r genus |
                grep -w {} |
                wc -l
        )

        printf "%s\t%d\t%d\t%d\n" {} ${n_species} ${n_strains} ${n_nr}
    ' |
    nwr append stdin --id |
    tsv-select -f 6,5,2,3,4 |
    tsv-sort -k2,2 |
    tsv-filter --ge 4:2 |
    (echo -e '#tax_id\tgenus\t#species\t#strains\t#NR' && cat) |
    rgr md stdin

```

| #tax_id | genus                     | #species | #strains | #NR |
|---------|---------------------------|----------|----------|-----|
| 1903434 | Atlantibacter             | 3        | 116      | 28  |
| 1386    | Bacillus                  | 3        | 3        | 0   |
| 82976   | Buttiauxella              | 10       | 16       | 12  |
| 203804  | Candidatus Blochmanniella | 5        | 11       | 7   |
| 568987  | Candidatus Hamiltonella   | 4        | 20       | 7   |
| 1048757 | Candidatus Moranella      | 2        | 2        | 1   |
| 472825  | Candidatus Purcelliella   | 1        | 5        | 2   |
| 568988  | Candidatus Regiella       | 2        | 4        | 2   |
| 401618  | Candidatus Riesia         | 4        | 6        | 3   |
| 1081630 | Candidatus Schneideria    | 1        | 3        | 1   |
| 158483  | Cedecea                   | 5        | 18       | 14  |
| 544     | Citrobacter               | 30       | 2032     | 313 |
| 413496  | Cronobacter               | 25       | 744      | 141 |
| 3021684 | Dryocola                  | 1        | 3        | 1   |
| 547     | Enterobacter              | 60       | 6168     | 643 |
| 561     | Escherichia               | 2216     | 39443    | 392 |
| 1649295 | Franconibacter            | 9        | 14       | 6   |
| 570     | Klebsiella                | 397      | 25697    | 181 |
| 579     | Kluyvera                  | 10       | 53       | 30  |
| 1330547 | Kosakonia                 | 12       | 73       | 45  |
| 83654   | Leclercia                 | 4        | 88       | 39  |
| 1330545 | Lelliottia                | 5        | 41       | 26  |
| 451512  | Mangrovibacter            | 1        | 3        | 1   |
| 447792  | Phytobacter               | 4        | 23       | 14  |
| 702     | Plesiomonas               | 2        | 48       | 38  |
| 1330546 | Pluralibacter             | 4        | 30       | 16  |
| 2994443 | Pseudenterobacter         | 1        | 4        | 3   |
| 2055880 | Pseudescherichia          | 2        | 11       | 9   |
| 1504576 | Pseudocitrobacter         | 1        | 2        | 2   |
| 160674  | Raoultella                | 11       | 311      | 15  |
| 590     | Salmonella                | 1389     | 14010    | 146 |
| 2726810 | Scandinavium              | 4        | 9        | 7   |
| 620     | Shigella                  | 123      | 2292     | 59  |
| 1335483 | Shimwellia                | 2        | 4        | 1   |
| 1649298 | Siccibacter               | 3        | 10       | 8   |
| 1301    | Streptococcus             | 2        | 2        | 0   |
| 2303321 | Superficieibacter         | 1        | 2        | 1   |
| 2943312 | Tenebrionicola            | 1        | 2        | 1   |
| 158851  | Trabulsiella              | 1        | 8        | 4   |
| 158876  | Yokenella                 | 3        | 14       | 4   |


## Collect proteins

```bash
cd ~/data/Escherichia/

ulimit -n `ulimit -Hn`

nwr template summary/assembly.tsv \
    --pro \
    --parallel 8

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
