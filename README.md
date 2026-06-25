# Genomes

## Software

```shell
cbp install nwr
cbp install sqlite3 # 3.34 or above

cbp install hmmer easel
cbp install mafft muscle trimal
cbp install samtools
cbp install fasttree iqtree2
cbp install mash mmseqs

cbp install fd pigz
cbp install jq tectonic
cbp install tva 
cbp install pgr
```

### Other Packages

- Pangenome
    - `PPanGGOLiN` is used in this project. Installation steps can be found
      [here](https://github.com/wang-q/dotfiles/blob/master/others.sh).

## Strain info

- [Bacteria](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2)
- [Archaea](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2157)
- [Fungi](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4751)
- [Oomycota](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4762)

## New Prokaryote Systematics Reference

Göker and Oren (2024) validly published names for
**two domains** and **seven kingdoms**of prokaryotes under
the International Code of Nomenclature of Prokaryotes (ICNP)(
[IJSEM](https://www.microbiologyresearch.org/content/journal/ijsem/10.1099/ijsem.0.006242))
. This established a formal higher-rank classification system for Bacteria and
Archaea.

| Domain   | Kingdom           | Former/Informal Name | Key Phyla                                  |
|----------|-------------------|----------------------|--------------------------------------------|
| Bacteria | Bacillati         | Terrabacteria        | Bacillota, Actinomycetota, Cyanobacteriota |
|          | Pseudomonadati    | Hydrobacteria        | Pseudomonadota, Bacteroidota               |
|          | Fusobacteriati    | —                    | Fusobacteriota                             |
|          | Thermotogati      | —                    | Thermotogota, Aquificota                   |
| Archaea  | Methanobacteriati | Euryarchaeota        | Methanobacteriota, Thermoplasmatota        |
|          | Nanobdellati      | DPANN group          | Nanobdellota, Micrarchaeota                |
|          | Thermoproteati    | TACK group           | Thermoproteota, Asgardarchaeota            |

This classification was adopted by NCBI Taxonomy in October 2024, introducing
the rank of kingdom for prokaryotes for the first time.

### Count qualified assemblies of Prokaryote groups

```bash
PHYLUM=$(
    nwr member Bacteria Archaea -r phylum |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        sed '1d' |
        tva select -f 2
)

for P in $PHYLUM; do
    GENUS=$(
        nwr member ${P} -r genus |
            grep -v -i "Candidatus " |
            grep -v -i "candidate " |
            sed '1d' |
            tva select -f 1 |
            tr "\n" "," |
            sed 's/,$/\)/' |
            sed 's/^/\(/'
    )

    if [[ ${#GENUS} -lt 3 ]]; then
        >&2 echo $P has no genera
        continue
    fi

    printf "$P\t"

    echo "
        SELECT
            COUNT(*)
        FROM ar
        WHERE 1=1
            AND genus_id IN $GENUS
            AND assembly_level IN ('Complete Genome', 'Chromosome')
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
        tr "\n" "\t"
        
    for L in 'Scaffold' 'Contig'; do
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND genus_id IN $GENUS
                AND assembly_level IN ('$L')
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    done |
        tr "\n" "\t" |
        sed 's/\t$//'

    echo
done \
    > groups.tsv

nwr append groups.tsv -r domain -r kingdom |
    tva select -f 5-6,1-4 |
    tva sort -k 2,3 |
    tva sort -r -k 1 |
    tva blank -f 1 -f 2 |
    (echo -e 'Domain\tKingdom\tPhylum\tComplete\tScaffold\tContig' && cat ) |
    tva to md --num
```

| Domain   | Kingdom           | Phylum                  | Complete | Scaffold | Contig |
|----------|-------------------|-------------------------|---------:|---------:|-------:|
| Bacteria | Bacillati         | Actinomycetota          |     7026 |    26124 |  20668 |
|          |                   | Armatimonadota          |        7 |        4 |      8 |
|          |                   | Bacillota               |    15717 |    42966 |  69650 |
|          |                   | Chloroflexota           |       55 |       63 |    109 |
|          |                   | Cyanobacteriota         |      460 |      803 |   1331 |
|          |                   | Mycoplasmatota          |     1024 |      382 |   1135 |
|          |                   | Vulcanimicrobiota       |        1 |        0 |      0 |
|          | Fusobacteriati    | Fusobacteriota          |      271 |      211 |    472 |
|          | NA                | Minisyncoccota          |        1 |        0 |      0 |
|          | Pseudomonadati    | Abditibacteriota        |        1 |        0 |      1 |
|          |                   | Acidobacteriota         |       58 |       38 |     67 |
|          |                   | Aquificota              |       27 |       26 |     67 |
|          |                   | Atribacterota           |        3 |        1 |      2 |
|          |                   | Bacteroidota            |     2281 |     7928 |  10745 |
|          |                   | Balneolota              |        4 |       15 |     39 |
|          |                   | Bdellovibrionota        |       59 |       48 |     44 |
|          |                   | Caldisericota           |        1 |        9 |      2 |
|          |                   | Calditrichota           |        2 |        0 |      3 |
|          |                   | Campylobacterota        |     1598 |     2584 |   8148 |
|          |                   | Chlamydiota             |      393 |       54 |    193 |
|          |                   | Chlorobiota             |       17 |        9 |     36 |
|          |                   | Chrysiogenota           |        3 |        5 |      0 |
|          |                   | Coprothermobacterota    |        1 |        1 |      2 |
|          |                   | Deferribacterota        |        9 |        9 |     22 |
|          |                   | Dictyoglomota           |        7 |        6 |      1 |
|          |                   | Elusimicrobiota         |        4 |        0 |      1 |
|          |                   | Fibrobacterota          |        2 |       23 |     60 |
|          |                   | Fidelibacterota         |        1 |        0 |      0 |
|          |                   | Gemmatimonadota         |       11 |        9 |     48 |
|          |                   | Ignavibacteriota        |        3 |        5 |     12 |
|          |                   | Kiritimatiellota        |        2 |        0 |      6 |
|          |                   | Lentisphaerota          |        2 |        1 |     23 |
|          |                   | Myxococcota             |      140 |       37 |    148 |
|          |                   | Nitrospinota            |        1 |        1 |     10 |
|          |                   | Nitrospirota            |       24 |       19 |     24 |
|          |                   | Planctomycetota         |      116 |       61 |    117 |
|          |                   | Pseudomonadota          |    36901 |    71037 | 157832 |
|          |                   | Rhodothermota           |       22 |       41 |     99 |
|          |                   | Spirochaetota           |      751 |      373 |   1411 |
|          |                   | Thermodesulfobacteriota |      198 |      279 |    487 |
|          |                   | Thermodesulfobiota      |        2 |        0 |      2 |
|          |                   | Thermomicrobiota        |        2 |        3 |      9 |
|          |                   | Thermosulfidibacterota  |        1 |        0 |      0 |
|          |                   | Verrucomicrobiota       |      158 |      237 |    272 |
|          |                   | Zhurongbacterota        |        1 |        0 |      0 |
|          | Thermotogati      | Deinococcota            |      118 |      142 |    234 |
|          |                   | Synergistota            |       16 |       49 |    110 |
|          |                   | Thermotogota            |       62 |      105 |     99 |
| Archaea  | Methanobacteriati | Methanobacteriota       |      544 |      547 |   1147 |
|          |                   | Thermoplasmatota        |       16 |        9 |     75 |
|          | Nanobdellati      | Microcaldota            |        0 |        0 |      0 |
|          |                   | Nanobdellota            |        1 |        0 |      0 |
|          | Promethearchaeati | Promethearchaeota       |        1 |        0 |      0 |
|          | Thermoproteati    | Nitrososphaerota        |       24 |       12 |     26 |
|          |                   | Thermoproteota          |      139 |      117 |    127 |

```bash
echo "
    SELECT
        family,
        COUNT(*) AS complete
    FROM ar
    WHERE 1=1
        AND family != 'NA'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    GROUP BY family
    HAVING complete >= 100
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    nwr restrict Bacteria |
    nwr append stdin -r phylum -r class -r order |
    tva select -f 3-5,1-2 |
    tva sort -k 1-4 |
    tva blank -f 1 -f 2 -f 3 -f 4 |
    (echo -e 'Phylum\tClass\tOrder\tFamily\tComplete' && cat ) |
    tva to md --num
```

| Phylum            | Class                 | Order                | Family                 | Complete |
|-------------------|-----------------------|----------------------|------------------------|---------:|
| Actinomycetota    | Actinomycetes         | Actinomycetales      | Actinomycetaceae       |      157 |
|                   |                       | Bifidobacteriales    | Bifidobacteriaceae     |      498 |
|                   |                       | Kitasatosporales     | Streptomycetaceae      |     1599 |
|                   |                       | Micrococcales        | Microbacteriaceae      |      630 |
|                   |                       |                      | Micrococcaceae         |      374 |
|                   |                       | Micromonosporales    | Micromonosporaceae     |      201 |
|                   |                       | Mycobacteriales      | Corynebacteriaceae     |      651 |
|                   |                       |                      | Mycobacteriaceae       |     1554 |
|                   |                       |                      | Nocardiaceae           |      322 |
|                   |                       | Propionibacteriales  | Propionibacteriaceae   |      163 |
|                   |                       | Pseudonocardiales    | Pseudonocardiaceae     |      157 |
| Bacillota         | Bacilli               | Bacillales           | Bacillaceae            |     3321 |
|                   |                       |                      | Caryophanaceae         |      102 |
|                   |                       |                      | Listeriaceae           |      871 |
|                   |                       |                      | Paenibacillaceae       |      481 |
|                   |                       |                      | Staphylococcaceae      |     3318 |
|                   |                       | Lactobacillales      | Enterococcaceae        |     1221 |
|                   |                       |                      | Lactobacillaceae       |     2034 |
|                   |                       |                      | Streptococcaceae       |     2418 |
|                   | Clostridia            | Eubacteriales        | Clostridiaceae         |      465 |
|                   |                       |                      | Oscillospiraceae       |      124 |
|                   |                       | Lachnospirales       | Lachnospiraceae        |      250 |
|                   |                       | Peptostreptococcales | Peptostreptococcaceae  |      325 |
| Bacteroidota      | Bacteroidia           | Bacteroidales        | Bacteroidaceae         |      233 |
|                   |                       |                      | Prevotellaceae         |      135 |
|                   | Flavobacteriia        | Flavobacteriales     | Flavobacteriaceae      |      801 |
|                   |                       |                      | Weeksellaceae          |      415 |
|                   | Sphingobacteriia      | Sphingobacteriales   | Sphingobacteriaceae    |      163 |
| Campylobacterota  | Epsilonproteobacteria | Campylobacterales    | Arcobacteraceae        |      104 |
|                   |                       |                      | Campylobacteraceae     |      897 |
|                   |                       |                      | Helicobacteraceae      |      517 |
| Chlamydiota       | Chlamydiia            | Chlamydiales         | Chlamydiaceae          |      388 |
| Fusobacteriota    | Fusobacteriia         | Fusobacteriales      | Fusobacteriaceae       |      240 |
| Mycoplasmatota    | Mollicutes            | Entomoplasmatales    | Spiroplasmataceae      |      121 |
|                   |                       | Mycoplasmatales      | Mycoplasmataceae       |      191 |
|                   | NA                    | Mycoplasmoidales     | Metamycoplasmataceae   |      483 |
|                   |                       |                      | Mycoplasmoidaceae      |      188 |
| Pseudomonadota    | Alphaproteobacteria   | Acetobacterales      | Acetobacteraceae       |      132 |
|                   |                       | Caulobacterales      | Caulobacteraceae       |      135 |
|                   |                       | Hyphomicrobiales     | Bartonellaceae         |      127 |
|                   |                       |                      | Brucellaceae           |      379 |
|                   |                       |                      | Nitrobacteraceae       |      425 |
|                   |                       |                      | Phyllobacteriaceae     |      139 |
|                   |                       |                      | Rhizobiaceae           |      617 |
|                   |                       | Rhodobacterales      | Paracoccaceae          |      202 |
|                   |                       |                      | Roseobacteraceae       |      307 |
|                   |                       | Rickettsiales        | Anaplasmataceae        |      452 |
|                   |                       |                      | Rickettsiaceae         |      203 |
|                   |                       | Sphingomonadales     | Erythrobacteraceae     |      103 |
|                   |                       |                      | Sphingomonadaceae      |      182 |
|                   | Betaproteobacteria    | Burkholderiales      | Alcaligenaceae         |     1183 |
|                   |                       |                      | Burkholderiaceae       |     1344 |
|                   |                       |                      | Comamonadaceae         |      333 |
|                   |                       |                      | Oxalobacteraceae       |      141 |
|                   |                       | Neisseriales         | Neisseriaceae          |      585 |
|                   | Gammaproteobacteria   | Aeromonadales        | Aeromonadaceae         |      469 |
|                   |                       | Alteromonadales      | Alteromonadaceae       |      143 |
|                   |                       |                      | Pseudoalteromonadaceae |      123 |
|                   |                       |                      | Shewanellaceae         |      160 |
|                   |                       | Enterobacterales     | Enterobacteriaceae     |    14630 |
|                   |                       |                      | Erwiniaceae            |      521 |
|                   |                       |                      | Hafniaceae             |      243 |
|                   |                       |                      | Morganellaceae         |      559 |
|                   |                       |                      | Pectobacteriaceae      |      269 |
|                   |                       |                      | Yersiniaceae           |      887 |
|                   |                       | Legionellales        | Coxiellaceae           |      100 |
|                   |                       |                      | Legionellaceae         |      276 |
|                   |                       | Lysobacterales       | Lysobacteraceae        |     1530 |
|                   |                       | Moraxellales         | Moraxellaceae          |     1709 |
|                   |                       | Oceanospirillales    | Halomonadaceae         |      239 |
|                   |                       | Pasteurellales       | Pasteurellaceae        |     1255 |
|                   |                       | Pseudomonadales      | Pseudomonadaceae       |     3757 |
|                   |                       | Thiotrichales        | Francisellaceae        |      141 |
|                   |                       |                      | Piscirickettsiaceae    |      120 |
|                   |                       | Vibrionales          | Vibrionaceae           |     1151 |
| Spirochaetota     | Leptospiria           | Leptospirales        | Leptospiraceae         |      159 |
|                   | Spirochaetia          | Spirochaetales       | Borreliaceae           |      286 |
|                   |                       |                      | Treponemataceae        |      254 |
| Verrucomicrobiota | Verrucomicrobiia      | Verrucomicrobiales   | Akkermansiaceae        |      116 |

## Download all valid Bacteria and Archaea genomes

- [Trichoderma](./groups/Trichoderma.md) is a good example of familiarizing
  yourself with the processing steps.
- [Bacteria](./groups/Bacteria.md): All genomes of **Bacteria**, species by
  species
- [Archaea](./groups/Archaea.md): All genomes of **Archaea**, species by species
- [Fungi](./groups/Fungi.md): All genomes of **Fungi**, species by species

## Prokaryote groups

- [Pseudomonas](groups/Pseudomonas.md)
- MTBC
- Tenericutes

## Eukaryote groups

- [Trichoderma](groups/Trichoderma.md)
- [Protistis](groups/Protists.md)

