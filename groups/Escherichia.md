# *Escherichia*

## Strain info

* [Escherichia](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=561)
* [Shigella](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=620)
* [Salmonella](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=590)
* [Klebsiella](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=570)

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
