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

Göker and Oren (2024) validly published names for **two domains** and **seven kingdoms**
of prokaryotes under the International Code of Nomenclature of Prokaryotes (ICNP)
([IJSEM](https://www.microbiologyresearch.org/content/journal/ijsem/10.1099/ijsem.0.006242)). This
established a formal higher-rank classification system for Bacteria and Archaea.

| Domain   | Kingdom           | Former/Informal Name | Key Phyla                                                  |
|----------|-------------------|----------------------|------------------------------------------------------------|
| Bacteria | Bacillati         | Terrabacteria        | Bacillota (Firmicutes), Actinomycetota, Cyanobacteriota    |
|          | Pseudomonadati    | Hydrobacteria        | Pseudomonadota (Proteobacteria), Bacteroidota, Chlorobiota |
|          | Fusobacteriati    | —                    | Fusobacteriota                                             |
|          | Thermotogati      | —                    | Thermotogota, Aquificota                                   |
| Archaea  | Methanobacteriati | Euryarchaeota        | Methanobacteriota, Thermoplasmatota, Halobacteriota        |
|          | Nanobdellati      | DPANN group          | Nanobdellota, Micrarchaeota                                |
|          | Thermoproteati    | TACK group           | Thermoproteota, Nitrososphaerota, Asgardarchaeota          |

This classification was adopted by NCBI Taxonomy in October 2024, introducing the rank of kingdom
for prokaryotes for the first time.

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

    for L in 'Complete Genome' 'Chromosome' 'Scaffold' 'Contig'; do
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
    tva select -f 6-7,1-5 |
    tva sort -k 2,3 |
    tva sort -r -k 1 |
    tva blank -f 1 -f 2 |
    (echo -e 'Domain\tKingdom\tPhylum\tComplete Genome\tChromosome\tScaffold\tContig' && cat ) |
    tva to md --num

```

| Domain   | Kingdom           | Phylum                  | Complete Genome | Chromosome | Scaffold | Contig |
| -------- | ----------------- | ----------------------- | --------------: | ---------: | -------: | -----: |
| Bacteria | Bacillati         | Actinomycetota          |            6050 |        976 |    26124 |  20668 |
|          |                   | Armatimonadota          |               3 |          4 |        4 |      8 |
|          |                   | Bacillota               |           14115 |       1602 |    42966 |  69650 |
|          |                   | Chloroflexota           |              54 |          1 |       63 |    109 |
|          |                   | Cyanobacteriota         |             416 |         44 |      803 |   1331 |
|          |                   | Mycoplasmatota          |             953 |         71 |      382 |   1135 |
|          |                   | Vulcanimicrobiota       |               1 |          0 |        0 |      0 |
|          | Fusobacteriati    | Fusobacteriota          |             262 |          9 |      211 |    472 |
|          | NA                | Minisyncoccota          |               1 |          0 |        0 |      0 |
|          | Pseudomonadati    | Abditibacteriota        |               1 |          0 |        0 |      1 |
|          |                   | Acidobacteriota         |              47 |         11 |       38 |     67 |
|          |                   | Aquificota              |              25 |          2 |       26 |     67 |
|          |                   | Atribacterota           |               3 |          0 |        1 |      2 |
|          |                   | Bacteroidota            |            1997 |        284 |     7928 |  10745 |
|          |                   | Balneolota              |               3 |          1 |       15 |     39 |
|          |                   | Bdellovibrionota        |              49 |         10 |       48 |     44 |
|          |                   | Caldisericota           |               1 |          0 |        9 |      2 |
|          |                   | Calditrichota           |               1 |          1 |        0 |      3 |
|          |                   | Campylobacterota        |            1482 |        116 |     2584 |   8148 |
|          |                   | Chlamydiota             |             303 |         90 |       54 |    193 |
|          |                   | Chlorobiota             |              16 |          1 |        9 |     36 |
|          |                   | Chrysiogenota           |               3 |          0 |        5 |      0 |
|          |                   | Coprothermobacterota    |               1 |          0 |        1 |      2 |
|          |                   | Deferribacterota        |               9 |          0 |        9 |     22 |
|          |                   | Dictyoglomota           |               7 |          0 |        6 |      1 |
|          |                   | Elusimicrobiota         |               4 |          0 |        0 |      1 |
|          |                   | Fibrobacterota          |               2 |          0 |       23 |     60 |
|          |                   | Fidelibacterota         |               1 |          0 |        0 |      0 |
|          |                   | Gemmatimonadota         |              10 |          1 |        9 |     48 |
|          |                   | Ignavibacteriota        |               3 |          0 |        5 |     12 |
|          |                   | Kiritimatiellota        |               2 |          0 |        0 |      6 |
|          |                   | Lentisphaerota          |               2 |          0 |        1 |     23 |
|          |                   | Myxococcota             |             131 |          9 |       37 |    148 |
|          |                   | Nitrospinota            |               1 |          0 |        1 |     10 |
|          |                   | Nitrospirota            |              24 |          0 |       19 |     24 |
|          |                   | Planctomycetota         |              86 |         30 |       61 |    117 |
|          |                   | Pseudomonadota          |           33304 |       3597 |    71037 | 157832 |
|          |                   | Rhodothermota           |              19 |          3 |       41 |     99 |
|          |                   | Spirochaetota           |             467 |        284 |      373 |   1411 |
|          |                   | Thermodesulfobacteriota |             186 |         12 |      279 |    487 |
|          |                   | Thermodesulfobiota      |               2 |          0 |        0 |      2 |
|          |                   | Thermomicrobiota        |               2 |          0 |        3 |      9 |
|          |                   | Thermosulfidibacterota  |               1 |          0 |        0 |      0 |
|          |                   | Verrucomicrobiota       |             149 |          9 |      237 |    272 |
|          |                   | Zhurongbacterota        |               1 |          0 |        0 |      0 |
|          | Thermotogati      | Deinococcota            |             113 |          5 |      142 |    234 |
|          |                   | Synergistota            |              12 |          4 |       49 |    110 |
|          |                   | Thermotogota            |              61 |          1 |      105 |     99 |
| Archaea  | Methanobacteriati | Methanobacteriota       |             523 |         21 |      547 |   1147 |
|          |                   | Thermoplasmatota        |              16 |          0 |        9 |     75 |
|          | Nanobdellati      | Microcaldota            |               0 |          0 |        0 |      0 |
|          |                   | Nanobdellota            |               1 |          0 |        0 |      0 |
|          | Promethearchaeati | Promethearchaeota       |               1 |          0 |        0 |      0 |
|          | Thermoproteati    | Nitrososphaerota        |              21 |          3 |       12 |     26 |
|          |                   | Thermoproteota          |             133 |          6 |      117 |    127 |

Table: refseq - Prokaryotes

## Download all valid Bacteria and Archaea genomes

- [Trichoderma](./groups/Trichoderma.md) is a good example of familiarizing yourself with the
  processing steps.
- [Bacteria](./groups/Bacteria.md): All genomes of **Bacteria** and **Archaea**, species by species
- [Fungi](./groups/Fungi.md): All genomes of **Fungi**, species by species

## Prokaryote groups

- [Pseudomonas](groups/Pseudomonas.md)
- MTBC
- Tenericutes

## Eukaryote groups

- [Trichoderma](groups/Trichoderma.md)
- [Protistis](groups/Protists.md)

