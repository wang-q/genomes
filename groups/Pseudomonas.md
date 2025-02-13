# *Pseudomonas* HGT

<!-- TOC -->
* [*Pseudomonas* HGT](#pseudomonas-hgt)
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
    * [Typical strains](#typical-strains)
    * [Rsync to hpcc](#rsync-to-hpcc)
  * [Collect proteins](#collect-proteins)
  * [Phylogenetics with bac120](#phylogenetics-with-bac120)
    * [Find corresponding representative proteins by `hmmsearch`](#find-corresponding-representative-proteins-by-hmmsearch)
    * [Domain related protein sequences](#domain-related-protein-sequences)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
    * [Condense branches in the protein tree](#condense-branches-in-the-protein-tree)
  * [InterProScan on all proteins of representative and typical strains](#interproscan-on-all-proteins-of-representative-and-typical-strains)
  * [Count members of protein families](#count-members-of-protein-families)
    * [Pseudomonas](#pseudomonas)
    * [P. aeruginosa](#p-aeruginosa)
    * [P. putida](#p-putida)
    * [Acinetobacter](#acinetobacter)
    * [A. baumannii](#a-baumannii)
  * [Protein families](#protein-families)
    * [IPR007416 - YggL 50S ribosome-binding protein](#ipr007416---yggl-50s-ribosome-binding-protein)
<!-- TOC -->

## Strain info

* [Pseudomonas](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=286)
* [Acinetobacter](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=469)
* [Stenotrophomonas](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=40323)

According to a 2020 [study](https://doi.org/10.1128/mSystems.00543-20), there are some order-level
changes in Gammaproteobacteria. We include both old and new orders.

* Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
* New ones: Moraxellales, Kangiellales, and Pseudomonadales

![msystems.00543-20-f0003.jpeg](../images/msystems.00543-20-f0003.jpeg)

Another [paper](https://doi.org/10.1099/ijsem.0.005542) proposed that an outlier clade of P.
aeruginosa (containing PA7) split into a new species, Pseudomonas paraeruginosa.

Similarly, a recent [publication](https://doi.org/10.1016/j.syapm.2021.126289) suggested that some
species of Pseudomonas should be transferred to Halopseudomonas and Stutzerimonas.

Additionally, Azotobacter is actually inside other Pseudomonas.

In 2024, [Rudra and Gupta](https://doi.org/10.3389/fmicb.2023.1273665) proposed splitting
Pseudomonas into several genera, a change that has been accepted by NCBI.

### Symlink

```shell
mkdir -p ~/data/Pseudomonas
cd ~/data/Pseudomonas

rm -fr ASSEMBLY
rm -fr STRAINS

ln -s ../Bacteria/ASSEMBLY ASSEMBLY
ln -s ../Bacteria/STRAINS STRAINS

```

### List all ranks

```shell
mkdir -p ~/data/Pseudomonas
cd ~/data/Pseudomonas

nwr member Pseudomonas Acinetobacter |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    rgr md stdin --num

nwr member Pseudomonas Acinetobacter -r "species group" -r "species subgroup" |
    tsv-select -f 1-3 |
    keep-header -- tsv-sort -k3,3 -k2,2 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    rgr md stdin

```

| rank             | count |
|------------------|------:|
| genus            |     2 |
| species          |   568 |
| strain           |  1813 |
| subspecies       |    13 |
| no rank          |   124 |
| species group    |     7 |
| species subgroup |     7 |
| isolate          |     3 |

| #tax_id | sci_name                                 | rank             |
|---------|------------------------------------------|------------------|
| 909768  | A. calcoaceticus/baumannii complex       | species group    |
| 2839056 | A. Taxon 24                              | species group    |
| 136841  | P. aeruginosa group                      | species group    |
| 136842  | P. chlororaphis group                    | species group    |
| 136843  | P. fluorescens group                     | species group    |
| 136845  | P. putida group                          | species group    |
| 136849  | P. syringae group                        | species group    |
| 2839060 | A. Taxon 24C                             | species subgroup |
| 2839057 | A. Taxon 24D                             | species subgroup |
| 2839061 | A. Taxon 24E                             | species subgroup |
| 627141  | P. nitroreducens/multiresinivorans group | species subgroup |
| 1232139 | P. oleovorans/pseudoalcaligenes group    | species subgroup |
| 251695  | P. syringae group genomosp. 1            | species subgroup |
| 251698  | P. syringae group genomosp. 2            | species subgroup |

### Species with assemblies

* Order Pseudomonadales and Moraxellales
    * Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
    * New ones: Moraxellales, Kangiellales, and Pseudomonadales
    * Pseudomonas aeruginosa
    * Acinetobacter baumannii

```shell
mkdir -p ~/data/Pseudomonas/summary
cd ~/data/Pseudomonas/summary

SPECIES=$(
    nwr member \
        Cellvibrionales Oceanospirillales Alteromonadales \
        Moraxellales Kangiellales Pseudomonadales \
        -r species |
        grep -v " sp." |
        grep -v -E "\bbacterium\b" |
        grep -v -E "symbiont\b" |
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
    tsv-filter -H --invert --str-in-fld species:Pseudomonas --lt RS:50 |
    tsv-filter -H --invert --str-in-fld species:Acinetobacter --lt RS:50 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    rgr md stdin --num

```

| species_id | species                |    RS |  CHR |
|-----------:|------------------------|------:|-----:|
|        287 | P. aeruginosa          | 10053 | 1010 |
|        470 | A. baumannii           |  9699 |  808 |
|      48296 | A. pittii              |   698 |   66 |
|        317 | P. syringae            |   666 |   59 |
|        303 | P. putida              |   306 |   92 |
|      38313 | Shewanella algae       |   214 |   26 |
|        294 | P. fluorescens         |   197 |   31 |
|        316 | Stutzerimonas stutzeri |   169 |   33 |
|     587753 | P. chlororaphis        |   151 |   93 |
|     756892 | A. indicus             |   130 |   23 |
|     380021 | P. protegens           |   105 |   32 |
|      40214 | A. johnsonii           |    92 |   22 |
|    1530123 | A. seifertii           |    84 |   34 |
|        476 | Moraxella bovis        |    39 |   37 |

## All assemblies

### Extract from `../Bacteria` and create assembly.tsv

```shell
cd ~/data/Pseudomonas

mkdir -p summary

# Target genus
GENUS=(
    # Pseudomonadales - All of Pseudomonadaceae
    $(
        nwr member Pseudomonadaceae -r genus |
            sed '1d' |
            cut -f 2
    )

    # Alteromonadales
    Shewanella

    # Moraxellales
    Acinetobacter
    Moraxella
)
#GENUS=$(IFS=, ; echo "${GENUS[*]}")

cat ../Bacteria/summary/collect.pass.tsv |
    tsv-join -H -d RefSeq_assembly_accession -f ~/Scripts/genomes/assembly/Bacteria.reference.tsv -k assembly_accession \
    > summary/collect.pass.tsv

cat ../Bacteria/summary/collect.pass.tsv | # 316617
    nwr restrict ${GENUS[*]} -f stdin -c 3 | # restrict to these genera 28409
    tsv-filter -H --le "C:20" --ge "N50:500000" | # more stringent parameters 4486
    sed '1d' |
    tsv-join -e -f ../Bacteria/ASSEMBLY/omit.lst -k 1 | # 4485
    tsv-join -e -f ../Bacteria/MinHash/abnormal.lst -k 1 | # 4411
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

```shell
cd ~/data/Pseudomonas

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
| strain  |  4437 |
| species |   341 |
| genus   |    36 |
| family  |    19 |
| order   |    14 |
| class   |     6 |

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Acinetobacter    |       51 |     1506 |
| Aquipseudomonas  |        1 |        8 |
| Azorhizophilus   |        1 |        2 |
| Azotobacter      |        3 |       12 |
| Bacillus         |        3 |        3 |
| Bacteroides      |        1 |        1 |
| Bifidobacterium  |        1 |        1 |
| Brucella         |        1 |        1 |
| Campylobacter    |        1 |        1 |
| Denitrificimonas |        1 |        1 |
| Ectopseudomonas  |        7 |       28 |
| Enterococcus     |        1 |        1 |
| Entomomonas      |        1 |        1 |
| Escherichia      |        1 |        2 |
| Francisella      |        1 |        1 |
| Halopseudomonas  |        6 |        8 |
| Klebsiella       |        2 |        2 |
| Listeria         |        1 |        1 |
| Metapseudomonas  |        4 |       13 |
| Moraxella        |       14 |       97 |
| Mycobacterium    |        1 |        1 |
| Mycobacteroides  |        1 |        1 |
| Permianibacter   |        1 |        1 |
| Phocaeicola      |        1 |        1 |
| Proteus          |        1 |        1 |
| Pseudomonas      |      187 |     2546 |
| Salmonella       |        1 |        1 |
| Shewanella       |       29 |      105 |
| Shigella         |        1 |        1 |
| Staphylococcus   |        1 |        1 |
| Streptococcus    |        2 |        2 |
| Stutzerimonas    |        9 |       75 |
| Thiopseudomonas  |        1 |        8 |
| Vibrio           |        1 |        1 |
| Xanthomonas      |        1 |        1 |
| Yersinia         |        1 |        1 |

## MinHash

```shell
cd ~/data/Pseudomonas/

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
#  1161 summary/NR.lst
#  3153 summary/redundant.lst

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#16

```

## Count valid species and strains

### For *genomic alignments* and *protein families*

```shell
cd ~/data/Pseudomonas/

nwr template summary/assembly.tsv \
    --count \
    --not-in MinHash/abnormal.lst \
    --rank order --rank genus \
    --lineage family --lineage genus

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
bash Count/lineage.sh 10

cat Count/lineage.count.tsv |
    rgr md stdin --fmt

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| order            | #species | #strains |
|------------------|---------:|---------:|
| Alteromonadales  |       29 |      105 |
| Bacillales       |        5 |        5 |
| Bacteroidales    |        2 |        2 |
| Enterobacterales |        7 |        8 |
| Lactobacillales  |        3 |        3 |
| Moraxellales     |       65 |    1,603 |
| Mycobacteriales  |        2 |        2 |
| Pseudomonadales  |      222 |    2,687 |

| genus           | #species | #strains |
|-----------------|---------:|---------:|
| Acinetobacter   |       51 |    1,506 |
| Aquipseudomonas |        1 |        8 |
| Azorhizophilus  |        1 |        2 |
| Azotobacter     |        3 |       12 |
| Bacillus        |        3 |        3 |
| Ectopseudomonas |        7 |       28 |
| Escherichia     |        1 |        2 |
| Halopseudomonas |        6 |        8 |
| Klebsiella      |        2 |        2 |
| Metapseudomonas |        4 |       13 |
| Moraxella       |       14 |       97 |
| Pseudomonas     |      187 |    2,536 |
| Shewanella      |       29 |      105 |
| Streptococcus   |        2 |        2 |
| Stutzerimonas   |        9 |       69 |
| Thiopseudomonas |        1 |        8 |

| #family          | genus         | species                                | count |
|------------------|---------------|----------------------------------------|------:|
| Moraxellaceae    | Acinetobacter | Acinetobacter baumannii                | 1,038 |
|                  |               | Acinetobacter bereziniae               |    12 |
|                  |               | Acinetobacter calcoaceticus            |    13 |
|                  |               | Acinetobacter haemolyticus             |    20 |
|                  |               | Acinetobacter indicus                  |    32 |
|                  |               | Acinetobacter johnsonii                |    51 |
|                  |               | Acinetobacter junii                    |    21 |
|                  |               | Acinetobacter lwoffii                  |    21 |
|                  |               | Acinetobacter nosocomialis             |    25 |
|                  |               | Acinetobacter pittii                   |    96 |
|                  |               | Acinetobacter radioresistens           |    11 |
|                  |               | Acinetobacter schindleri               |    10 |
|                  |               | Acinetobacter seifertii                |    35 |
|                  |               | Acinetobacter towneri                  |    13 |
|                  |               | Acinetobacter ursingii                 |    11 |
|                  | Moraxella     | Moraxella bovis                        |    38 |
|                  |               | Moraxella catarrhalis                  |    21 |
|                  |               | Moraxella osloensis                    |    14 |
| Pseudomonadaceae | Pseudomonas   | Pseudomonas aeruginosa                 | 1,506 |
|                  |               | Pseudomonas amygdali                   |    16 |
|                  |               | Pseudomonas asiatica                   |    21 |
|                  |               | Pseudomonas atacamensis                |    11 |
|                  |               | Pseudomonas brassicacearum             |    13 |
|                  |               | Pseudomonas chlororaphis               |   112 |
|                  |               | Pseudomonas fluorescens                |    28 |
|                  |               | Pseudomonas fragi                      |    16 |
|                  |               | Pseudomonas fulva                      |    11 |
|                  |               | Pseudomonas monteilii                  |    15 |
|                  |               | Pseudomonas oryzihabitans              |    14 |
|                  |               | Pseudomonas poae                       |    22 |
|                  |               | Pseudomonas protegens                  |    54 |
|                  |               | Pseudomonas putida                     |   102 |
|                  |               | Pseudomonas rhodesiae                  |    13 |
|                  |               | Pseudomonas simiae                     |    11 |
|                  |               | Pseudomonas synxantha                  |    12 |
|                  |               | Pseudomonas syringae                   |    82 |
|                  |               | Pseudomonas syringae group genomosp. 3 |    10 |
|                  |               | Pseudomonas trivialis                  |    24 |
|                  |               | Pseudomonas viridiflava                |    21 |
|                  | Stutzerimonas | Stutzerimonas balearica                |    10 |
|                  |               | Stutzerimonas frequens                 |    15 |
|                  |               | Stutzerimonas stutzeri                 |    29 |
| Shewanellaceae   | Shewanella    | Shewanella algae                       |    30 |
|                  |               | Shewanella baltica                     |    15 |
|                  |               | Shewanella xiamenensis                 |    13 |

### Count strains - Genus

```shell
cd ~/data/Pseudomonas

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

| #tax_id | genus           | #species | #strains | #NR |
|---------|-----------------|----------|----------|-----|
| 469     | Acinetobacter   | 153      | 1507     | 411 |
| 3236652 | Aquipseudomonas | 1        | 8        | 8   |
| 157913  | Azorhizophilus  | 1        | 2        | 1   |
| 352     | Azotobacter     | 7        | 12       | 5   |
| 1386    | Bacillus        | 3        | 3        | 0   |
| 3236654 | Ectopseudomonas | 11       | 28       | 25  |
| 561     | Escherichia     | 2        | 2        | 2   |
| 2901189 | Halopseudomonas | 6        | 8        | 3   |
| 570     | Klebsiella      | 2        | 2        | 0   |
| 3236656 | Metapseudomonas | 4        | 13       | 11  |
| 475     | Moraxella       | 15       | 97       | 31  |
| 286     | Pseudomonas     | 430      | 2547     | 520 |
| 22      | Shewanella      | 41       | 105      | 80  |
| 1301    | Streptococcus   | 2        | 2        | 0   |
| 2901164 | Stutzerimonas   | 14       | 75       | 61  |
| 1654787 | Thiopseudomonas | 1        | 8        | 5   |

### Typical strains

* Typical (Popular) strains in [pseudomonas.com](https://www.pseudomonas.com/strain/list)
    * 107 - Pseudomonas aeruginosa PAO1
    * 109 - Pseudomonas aeruginosa UCBPP-PA14
    * 119 - Pseudomonas aeruginosa PA7
    * 6441 - Pseudomonas aeruginosa PAK
    * 125 - Pseudomonas aeruginosa LESB58
    * 7360 - Pseudomonas putida KT2440
    * 118 - Pseudomonas putida F1
    * 479 - Pseudomonas chlororaphis subsp. aureofaciens 30-84
    * 116 - Pseudomonas fluorescens SBW25
    * 111 - Pseudomonas syringae pv. tomato DC3000
    * 112 - Pseudomonas syringae pv. syringae B728a
    * 114 - Pseudomonas savastanoi pv. phaseolicola 1448A
    * 123 - Pseudomonas stutzeri A1501
    * 113 - Pseudomonas protegens Pf-5
    * 117 - Pseudomonas entomophila L48

    * Pseudomonas fluorescens SBW25 was removed from refseq

    * Pseudomonas aeruginosa CF39S is Pseudom_aeruginosa_GCF_011466835_1

* Typical strains of Acin_bau was listed
  in [here](https://www.biorxiv.org/content/10.1101/2022.02.27.482139v3.full)
    * ATCC17978-VUB
    * ATCC19606-VUB
    * DSM30011-VUB https://www.ncbi.nlm.nih.gov/assembly/GCF_001936675.2
        * Assembly level: Scaffold
    * AB5075-VUB https://www.ncbi.nlm.nih.gov/assembly/GCF_016919505.2

```shell
cd ~/data/Pseudomonas

# Typical strains
for S in \
    Pseudom_aeruginosa_PAO1 \
    Pseudom_aeruginosa_UCBPP_PA14_GCF_000014625_1 \
    Pseudom_aeruginosa_PA7_GCF_000017205_1 \
    Pseudom_aeruginosa_PAK_GCF_000408865_1 \
    Pseudom_aeruginosa_CF39S_GCF_011466835_1 \
    Pseudom_aeruginosa_LESB58_GCF_000026645_1 \
    Pseudom_viridifl_CFBP_1590_GCF_900184295_1 \
    Pseudom_DC3000_GCF_000007805_1 \
    Pseudom_syringae_B728a_GCF_000012245_1 \
    Pseudom_fluore_ATCC_13525_GCF_900215245_1 \
    Pseudom_putida_KT2440_GCF_000007565_2 \
    Pseudom_putida_NBRC_14164_GCF_000412675_1 \
    Stu_stutz_A1501_GCF_000013785_1 \
    Pseudom_chloror_30_84_GCF_000281915_1 \
    Pseudom_amygdali_ATCC_11528_GCF_000145945_2 \
    Pseudom_proteg_Pf_5_GCF_000012265_1 \
    Pseudom_proteg_CHA0_GCF_900560965_1 \
    Pseudom_entomophila_L48_GCF_000026105_1 \
    Acin_baum_K09_14_GCF_008632635_1 \
    Acin_baum_ATCC_17978_Lab_WT_GCF_004794235_2 \
    Acin_baum_ATCC_19606_GCF_019331655_1 \
    Acin_baum_DSM_30011_GCF_001936675_2 \
    Acin_baum_AB5075_VUB_GCF_016919505_2 \
    ; do
    echo ${S}
done \
    > summary/typical.manual.lst


```

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Pseudomonas/ \
    wangq@202.119.37.251:data/Pseudomonas

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Pseudomonas/ \
    wangq@58.213.64.36:data/Pseudomonas

```

## Collect proteins

```shell
cd ~/data/Pseudomonas/

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

| #item      |      count |
|------------|-----------:|
| species    |        341 |
| strain_sum |      4,437 |
| total_sum  | 21,652,155 |
| dedup_sum  |  4,967,821 |
| rep_sum    |  2,370,952 |
| fam88_sum  |  2,026,703 |
| fam38_sum  |  1,606,432 |

## Phylogenetics with bac120

```shell
cd ~/data/Pseudomonas/

# The Bacteria HMM set
nwr kb bac120 -o HMM
cp HMM/bac120.lst HMM/marker.lst

```

### Find corresponding representative proteins by `hmmsearch`

```shell
cd ~/data/Pseudomonas

cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -e -f ASSEMBLY/sp.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 \
    > Protein/species-f.tsv

cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ -s Protein/"${SPECIES}"/bac120.tsv ]]; then
        continue
    fi
    if [[ ! -f Protein/"${SPECIES}"/rep_seq.fa.gz ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    cat HMM/marker.lst |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 "
            gzip -dcf Protein/${SPECIES}/rep_seq.fa.gz |
                hmmsearch --cut_nc --noali --notextw HMM/hmm/{}.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), q({}), \$1; '
        " \
        > Protein/${SPECIES}/bac120.tsv
done

cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ ! -s Protein/"${SPECIES}"/bac120.tsv ]]; then
        continue
    fi
    if [[ ! -f Protein/"${SPECIES}"/seq.sqlite ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    nwr seqdb -d Protein/${SPECIES} --rep f3=Protein/${SPECIES}/bac120.tsv

done

```

### Domain related protein sequences

```shell
cd ~/data/Pseudomonas

mkdir -p Domain

# each assembly
cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ ! -f Protein/"${SPECIES}"/seq.sqlite ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    echo "
        SELECT
            seq.name,
            asm.name,
            rep.f3
        FROM asm_seq
        JOIN rep_seq ON asm_seq.seq_id = rep_seq.seq_id
        JOIN seq ON asm_seq.seq_id = seq.id
        JOIN rep ON rep_seq.rep_id = rep.id
        JOIN asm ON asm_seq.asm_id = asm.id
        WHERE 1=1
            AND rep.f3 IS NOT NULL
        ORDER BY
            asm.name,
            rep.f3
        " |
        sqlite3 -tabs Protein/${SPECIES}/seq.sqlite \
        > Protein/${SPECIES}/seq_asm_f3.tsv

    hnsm some Protein/"${SPECIES}"/pro.fa.gz <(
            tsv-select -f 1 Protein/"${SPECIES}"/seq_asm_f3.tsv |
            rgr dedup stdin
        )
done |
    hnsm dedup stdin |
    hnsm gz stdin -o Domain/bac120.fa

fd --full-path "Protein/.+/seq_asm_f3.tsv" -X cat \
    > Domain/seq_asm_f3.tsv

cat Domain/seq_asm_f3.tsv |
    tsv-join -e -d 2 -f summary/redundant.lst -k 1 |
    tsv-join -e -d 2 -f ASSEMBLY/sp.lst -k 1 \
    > Domain/seq_asm_f3.NR.tsv

```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Pseudomonas

# Extract proteins
cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

        mkdir -p Domain/{}

        hnsm some Domain/bac120.fa.gz <(
            cat Domain/seq_asm_f3.tsv |
                tsv-filter --str-eq "3:{}" |
                tsv-select -f 1 |
                rgr dedup stdin
            ) \
            > Domain/{}/{}.pro.fa
    '

# Align each marker
cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo >&2 "==> marker [{}]"
        if [ ! -s Domain/{}/{}.pro.fa ]; then
            exit
        fi
        if [ -s Domain/{}/{}.aln.fa ]; then
            exit
        fi

#        muscle -quiet -in Domain/{}/{}.pro.fa -out Domain/{}/{}.aln.fa
        mafft --auto Domain/{}/{}.pro.fa > Domain/{}/{}.aln.fa
    '

cat HMM/marker.lst |
while read marker; do
    echo >&2 "==> marker [${marker}]"
    if [ ! -s Domain/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sometimes `muscle` can not produce alignments
    if [ ! -s Domain/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # Only NR strains
    # 1 name to many names
    cat Domain/seq_asm_f3.NR.tsv |
        tsv-filter --str-eq "3:${marker}" |
        tsv-select -f 1-2 |
        hnsm replace -s Domain/${marker}/${marker}.aln.fa stdin \
        > Domain/${marker}/${marker}.replace.fa
done

# Concat marker genes
cat HMM/marker.lst |
while read marker; do
    if [ ! -s Domain/${marker}/${marker}.pro.fa ]; then
        continue
    fi
    if [ ! -s Domain/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    cat Domain/${marker}/${marker}.replace.fa

    # empty line for .fas
    echo
done \
    > Domain/bac120.aln.fas

cat Domain/seq_asm_f3.NR.tsv |
    cut -f 2 |
    rgr dedup stdin |
    sort |
    fasops concat Domain/bac120.aln.fas stdin -o Domain/bac120.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Domain/bac120.aln.fa -out Domain/bac120.trim.fa -automated1

hnsm size Domain/bac120.*.fa |
    rgr dedup stdin -f 2 |
    cut -f 2
#73147
#45083

# To make it faster
FastTree -fastest -noml Domain/bac120.trim.fa > Domain/bac120.trim.newick

```

### Condense branches in the protein tree

```shell
mkdir -p ~/data/Pseudomonas/tree
cd ~/data/Pseudomonas/tree

cat ../Domain/bac120.trim.newick |
    nwr order stdin --nd --an \
    > bac120.order.newick

nwr pl-condense --map -r species \
    bac120.order.newick ../Count/species.tsv |
    nwr order stdin --nd --an \
    > bac120.condensed.newick

mv condensed.tsv bac120.condense.tsv

# svg
nwr topo --bl bac120.condensed.newick | # remove comments
    nw_display -s -b 'visibility:hidden' -w 1200 -v 20 - \
    > Pseudomonas.bac120.svg

```

## InterProScan on all proteins of representative and typical strains

```shell
cd ~/data/Pseudomonas

cat summary/NR.lst summary/reference.lst summary/typical.manual.lst |
    sort |
    uniq \
    > summary/ips.lst

wc -l summary/ips.lst
#1119

for S in $(cat summary/ips.lst); do
    if [ -s STRAINS/${S}/family.tsv ]; then
        continue
    fi

    >&2 echo "==> ${S}"

    mkdir -p STRAINS/${S}

    # Too many sequences would make ips crash
    faops split-about ASSEMBLY/${S}/*_protein.faa.gz 500000 STRAINS/${S}/
done

# max job number is 200
N_JOBS=$(bjobs -w | wc -l)
COUNT=
cat summary/ips.lst | #sort -r |
    grep -v -F -w -f <(bjobs -w | tr -s " " | cut -d " " -f 7 | cut -d "/" -f 2) | # Job exists
while read S; do
    for f in $(find STRAINS/${S}/ -maxdepth 1 -type f -name "[0-9]*.fa" | sort); do
        if [ -e ${f}.tsv ]; then
            continue
        fi

        COUNT=$((COUNT + 1))

        if ((  COUNT + N_JOBS > 195 )); then
            >&2 echo MAX JOBS
            break 2 # In a nested loop
        fi

        >&2 echo "==> ${f}"

        bsub -q mpi -n 24 -J "${f}" "
            interproscan.sh --cpu 24 -dp -f tsv,json -i ${f} --output-file-base ${f}
        "
    done
done

find STRAINS/ -type f -name "*.json" | sort |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        if [ $(({#} % 10)) -eq "0" ]; then
            >&2 printf "."
        fi
        pigz -p 4 {}
    '

find STRAINS/ -type f -name "*.json.gz" -size +1M | wc -l

# same protein may have multiple families
for S in $(cat summary/ips.lst); do
    if [ -s STRAINS/${S}/family.tsv ]; then
        continue
    fi

    for f in $(find STRAINS/${S}/ -maxdepth 1 -type f -name "[0-9]*.json.gz" | sort); do
        >&2 echo "==> ${f}"
        gzip -dcf ${f} |
            jq .results |
            jq -r -c '
                .[] |
                .xref as $name |
                .matches[] |
                .signature.entry |
                select(.type == "FAMILY") |
                [$name[0].name, .accession, .description] |
                @tsv
            ' |
            rgr dedup stdin
    done \
        > STRAINS/${S}/family.tsv
done

find STRAINS/ -type f -name "family.tsv" -empty

find STRAINS/ -type f -name "family.tsv" | wc -l
#1121

find STRAINS/ -type d |
    grep -v -Fw -f summary/ips.lst

find STRAINS/ -type f -name "family.tsv" |
    xargs wc -l |
    sed 's/^\s*//g' |
    grep -v "total$" |
    cut -d" " -f 1 |
    tsv-summarize --quantile 1:0,0.25,0.5,0.75,1
#735     2644    3471    3823    9168

find STRAINS/ -type f -name "family.tsv" |
    xargs wc -l |
    sed 's/^\s*//g' |
    grep -v "total$" |
    tsv-filter -d ' ' --lt 1:1500
#1253 STRAINS/Cox_burn_RSA_493/family.tsv
#1269 STRAINS/Cam_jej_jejuni_NCTC_11168_ATCC_700819/family.tsv
#735 STRAINS/Chlamydia_tracho_D_UW_3_CX/family.tsv

# seems OK

# Count protein family per strains
for S in $(cat summary/ips.lst); do
    if [ ! -s STRAINS/${S}/family.tsv ]; then
        continue
    fi
    if [ -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi

    cat STRAINS/${S}/family.tsv |
        tsv-summarize -g 2,3 --count \
        > STRAINS/${S}/family-count.tsv

done

```

## Count members of protein families

### Pseudomonas

```shell
cd ~/data/Pseudomonas

mkdir -p count

cat summary/strains.taxon.tsv |
    tsv-join -k 1 -f summary/ips.lst |
    tsv-filter --or --str-eq 4:Pseudomonas --str-eq 4:Halopseudomonas --str-eq 4:Stutzerimonas --str-eq 4:Azotobacter \
    > count/Pseudomonas.taxon.tsv

N_SPECIES=$(
    cat count/Pseudomonas.taxon.tsv |
        tsv-select -f 3 |
        rgr dedup stdin |
        wc -l
)
echo ${N_SPECIES}
# 91

cat count/Pseudomonas.taxon.tsv |
    tsv-select -f 3 |
    rgr dedup stdin |
    sort | # head |
while read SPECIES; do
    N_STRAINS=$(
        cat count/Pseudomonas.taxon.tsv |
            tsv-filter --str-eq "3:${SPECIES}" |
            wc -l
    )

    cat count/Pseudomonas.taxon.tsv |
        tsv-filter --str-eq "3:${SPECIES}" |
        tsv-select -f 1 |
        SPECIES=${SPECIES} N_STRAINS=${N_STRAINS} perl -nl -MPath::Tiny -e '
            BEGIN { our %count_of; }

            my @lines = path(qq{STRAINS/$_/family-count.tsv})->lines({chomp => 1});
            for my $l (@lines) {
                my ($f, $d, $c) = split /\t/, $l;
                $count_of{qq($f\t$d)} += $c;
            }

            #print $ENV{SPECIES}, $ENV{N_STRAINS};

            END {
                for my $f (sort keys %count_of) {
                    print join qq(\t), $f, $ENV{SPECIES}, sprintf(q{%.1f}, $count_of{$f} / $ENV{N_STRAINS} ) ;
                }
            }
        '
done |
    tsv-filter --istr-not-in-fld 2:"probable" |
    tsv-filter --istr-not-in-fld 2:"putative" |
    tsv-filter --istr-not-in-fld 2:"Uncharacterised" |
    tsv-filter --istr-not-in-fld 2:" DUF" \
    > count/Pseudomonas.family.tsv

# Families occurring in at least 81/91 species
cat count/Pseudomonas.family.tsv |
    tsv-summarize -g 1,2 --count |
    tsv-filter --gt "3:$((N_SPECIES - 10))" \
    > count/Pseudomonas.universal.tsv

```

### P. aeruginosa

```shell
cd ~/data/Pseudomonas

cat count/Pseudomonas.taxon.tsv |
    tsv-filter --str-eq "3:Pseudomonas aeruginosa" |
    wc -l
#12

# All other species should have only 1 family member
# The number is shrunk to be more tolerable, such as Pseudomonas paraeruginosa
cat count/Pseudomonas.family.tsv |
    tsv-filter --str-ne "3:Pseudomonas aeruginosa" |
    tsv-filter --le 4:1.5 |
    tsv-summarize -g 1,2 --count |
    tsv-filter --gt "3:$((N_SPECIES - 15))" \
    > count/Paer.1.tsv

# All strains of target species should have multiple family members
cat count/Pseudomonas.family.tsv |
    tsv-filter --str-eq "3:Pseudomonas aeruginosa" |
    tsv-filter --ge 4:1.8 |
    tsv-filter --lt 4:3 \
    > count/Paer.n.tsv

wc -l \
    count/Pseudomonas.universal.tsv\
    count/Paer.1.tsv \
    count/Paer.n.tsv
#  1640 count/Pseudomonas.universal.tsv
#  1150 count/Paer.1.tsv
#   305 count/Paer.n.tsv

cat count/Paer.n.tsv |
    tsv-join -k 1 -f count/Pseudomonas.universal.tsv |
    tsv-join -k 1 -f count/Paer.1.tsv |
    tsv-select -f 1,2,4 |
    tsv-sort |
    (echo -e "#family\tdesc\tcount" && cat) |
    mlr --itsv --omd cat

```

| #family   | desc                                                          | count |
|-----------|---------------------------------------------------------------|-------|
| IPR000440 | NADH:ubiquinone/plastoquinone oxidoreductase, chain 3         | 2.0   |
| IPR000473 | Ribosomal protein L36                                         | 2.0   |
| IPR000813 | 7Fe ferredoxin                                                | 2.1   |
| IPR001353 | Proteasome, subunit alpha/beta                                | 2.2   |
| IPR001404 | Heat shock protein Hsp90 family                               | 2.1   |
| IPR002307 | Tyrosine-tRNA ligase                                          | 2.2   |
| IPR002480 | DAHP synthetase, class II                                     | 2.8   |
| IPR002830 | UbiD decarboxylyase family                                    | 2.2   |
| IPR003672 | CobN/magnesium chelatase                                      | 2.2   |
| IPR004633 | Na/Pi-cotransporter II-related/YqeW-like protein              | 1.8   |
| IPR004685 | Branched-chain amino acid transport system II carrier protein | 2.2   |
| IPR004769 | Adenylosuccinate lyase                                        | 2.1   |
| IPR005955 | Glutathione S-transferases, class Zeta                        | 2.3   |
| IPR005999 | Glycerol kinase                                               | 2.1   |
| IPR006261 | dNTP triphosphohydrolase                                      | 2.2   |
| IPR006684 | Acyl-CoA thioester hydrolase YbgC/YbaW family                 | 2.2   |
| IPR007416 | YggL 50S ribosome-binding protein                             | 2.1   |
| IPR007692 | DNA helicase, DnaB type                                       | 1.8   |
| IPR011548 | 3-hydroxyisobutyrate dehydrogenase                            | 2.2   |
| IPR011757 | Lytic transglycosylase MltB                                   | 2.2   |
| IPR014311 | Guanine deaminase                                             | 2.2   |
| IPR017039 | Virulence factor BrkB                                         | 2.2   |
| IPR023695 | Thiosulfate sulfurtransferase, bacterial                      | 2.3   |
| IPR024088 | Tyrosine-tRNA ligase, bacterial-type                          | 2.2   |
| IPR026968 | 3-oxoadipate enol-lactonase 1/2                               | 2.0   |
| IPR028883 | tRNA-specific adenosine deaminase                             | 2.0   |

### P. putida

```shell
cd ~/data/Pseudomonas

cat count/Pseudomonas.taxon.tsv |
    tsv-filter --str-eq "3:Pseudomonas putida" |
    wc -l
#55

# All other species should have only 1 family member
cat count/Pseudomonas.family.tsv |
    tsv-filter --str-ne "3:Pseudomonas putida" |
    tsv-filter --le 4:1.5 |
    tsv-summarize -g 1,2 --count |
    tsv-filter --gt "3:$((N_SPECIES - 15))" \
    > count/Pput.1.tsv

# All strains of target species should have multiple family members
cat count/Pseudomonas.family.tsv |
    tsv-filter --str-eq "3:Pseudomonas putida" |
    tsv-filter --ge 4:1.8 |
    tsv-filter --lt 4:3 \
    > count/Pput.n.tsv

wc -l \
    count/Pseudomonas.universal.tsv\
    count/Pput.1.tsv \
    count/Pput.n.tsv
#  1640 count/Pseudomonas.universal.tsv
#  1147 count/Pput.1.tsv
#   239 count/Pput.n.tsv

cat count/Pput.n.tsv |
    tsv-join -k 1 -f count/Pseudomonas.universal.tsv |
    tsv-join -k 1 -f count/Pput.1.tsv |
    tsv-select -f 1,2,4 |
    tsv-sort |
    (echo -e "#family\tdesc\tcount" && cat) |
    mlr --itsv --omd cat

```

| #family   | desc                                                       | count |
|-----------|------------------------------------------------------------|-------|
| IPR001783 | Lumazine-binding protein                                   | 1.9   |
| IPR002033 | Sec-independent periplasmic protein translocase TatC       | 1.9   |
| IPR002446 | Lipocalin, bacterial                                       | 1.9   |
| IPR003171 | Methylenetetrahydrofolate reductase-like                   | 2.0   |
| IPR006312 | Sec-independent protein translocase protein TatA/E         | 1.9   |
| IPR012794 | Beta-ketoadipate transcriptional regulator, PcaR/PcaU/PobR | 1.8   |
| IPR018448 | Sec-independent protein translocase protein TatB           | 1.9   |
| IPR018550 | Lipid A 3-O-deacylase-related                              | 2.3   |

### Acinetobacter

```shell
cd ~/data/Pseudomonas

mkdir -p count

cat summary/strains.taxon.tsv |
    tsv-join -k 1 -f summary/ips.lst |
    tsv-filter --or --str-eq 4:Acinetobacter \
    > count/Acinetobacter.taxon.tsv

N_SPECIES=$(
    cat count/Acinetobacter.taxon.tsv |
        tsv-select -f 3 |
        rgr dedup stdin |
        wc -l
)
echo ${N_SPECIES}
# 22

cat count/Acinetobacter.taxon.tsv |
    tsv-select -f 3 |
    rgr dedup stdin |
    sort | # head |
while read SPECIES; do
    N_STRAINS=$(
        cat count/Acinetobacter.taxon.tsv |
            tsv-filter --str-eq "3:${SPECIES}" |
            wc -l
    )

    cat count/Acinetobacter.taxon.tsv |
        tsv-filter --str-eq "3:${SPECIES}" |
        tsv-select -f 1 |
        SPECIES=${SPECIES} N_STRAINS=${N_STRAINS} perl -nl -MPath::Tiny -e '
            BEGIN { our %count_of; }

            my @lines = path(qq{STRAINS/$_/family-count.tsv})->lines({chomp => 1});
            for my $l (@lines) {
                my ($f, $d, $c) = split /\t/, $l;
                $count_of{qq($f\t$d)} += $c;
            }

            #print $ENV{SPECIES}, $ENV{N_STRAINS};

            END {
                for my $f (sort keys %count_of) {
                    print join qq(\t), $f, $ENV{SPECIES}, sprintf(q{%.1f}, $count_of{$f} / $ENV{N_STRAINS} ) ;
                }
            }
        '
done |
    tsv-filter --istr-not-in-fld 2:"probable" |
    tsv-filter --istr-not-in-fld 2:"putative" |
    tsv-filter --istr-not-in-fld 2:"Uncharacterised" |
    tsv-filter --istr-not-in-fld 2:" DUF" \
    > count/Acinetobacter.family.tsv

# Families occurring in at least 18/22 species
cat count/Acinetobacter.family.tsv |
    tsv-summarize -g 1,2 --count |
    tsv-filter --gt "3:$((N_SPECIES - 4))" \
    > count/Acinetobacter.universal.tsv

```

### A. baumannii

```shell
cd ~/data/Pseudomonas

cat count/Acinetobacter.taxon.tsv |
    tsv-filter --str-eq "3:Acinetobacter baumannii" |
    wc -l
#73

# All other species should have only 1 family member
cat count/Acinetobacter.family.tsv |
    tsv-filter --str-ne "3:Pseudomonas baumannii" |
    tsv-filter --le 4:1.5 |
    tsv-summarize -g 1,2 --count |
    tsv-filter --gt "3:$((N_SPECIES - 6))" \
    > count/Abau.1.tsv

# All strains of target species should have multiple family members
cat count/Acinetobacter.family.tsv |
    tsv-filter --str-eq "3:Acinetobacter baumannii" |
    tsv-filter --ge 4:1.8 |
    tsv-filter --lt 4:3 \
    > count/Abau.n.tsv

wc -l \
    count/Acinetobacter.universal.tsv\
    count/Abau.1.tsv \
    count/Abau.n.tsv
#  1613 count/Acinetobacter.universal.tsv
#  1132 count/Abau.1.tsv
#   177 count/Abau.n.tsv

cat count/Abau.n.tsv |
    tsv-join -k 1 -f count/Acinetobacter.universal.tsv |
    tsv-join -k 1 -f count/Abau.1.tsv |
    tsv-select -f 1,2,4 |
    tsv-sort |
    (echo -e "#family\tdesc\tcount" && cat) |
    mlr --itsv --omd cat

```

| #family   | desc                                                          | count |
|-----------|---------------------------------------------------------------|-------|
| IPR001930 | Peptidase M1, alanine aminopeptidase/leukotriene A4 hydrolase | 1.9   |
| IPR002129 | Pyridoxal phosphate-dependent decarboxylase                   | 2.0   |
| IPR006424 | Glyceraldehyde-3-phosphate dehydrogenase, type I              | 1.8   |
| IPR030048 | Survival protein SurE                                         | 1.8   |
| IPR030664 | FAD-dependent oxidoreductase SdhA/FrdA/AprA                   | 2.0   |
| IPR047109 | Cinnamyl alcohol dehydrogenase-like                           | 1.9   |

## Protein families

### IPR007416 - YggL 50S ribosome-binding protein

<https://www.ebi.ac.uk/interpro/entry/InterPro/IPR007416/>

* Pfam:     PF04320 - YggL_50S_bp
* PANTHER:  PTHR38778 - CYTOPLASMIC PROTEIN-RELATED (PTHR38778)

```shell
cd ~/data/Pseudomonas

cat STRAINS/Pseudom_aeruginosa_PAO1/*.tsv |
    grep "IPR007416"

mkdir -p Targets/YggL/HMM

curl -L https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF04320?annotation=hmm |
    gzip -dc \
    > Targets/YggL/HMM/YggL_50S_bp.hmm
curl -L www.pantherdb.org/panther/exportHmm.jsp?acc=PTHR38778 \
    > Targets/YggL/HMM/PTHR38778.hmm

# Ribosomal protein L10 and S8
curl -L https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00466?annotation=hmm |
    gzip -dc \
    > Targets/YggL/HMM/Ribosomal_L10.hmm
curl -L https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00410?annotation=hmm |
    gzip -dc \
    > Targets/YggL/HMM/Ribosomal_S8.hmm

E_VALUE=1e-20
for domain in YggL_50S_bp PTHR38778 Ribosomal_L10 Ribosomal_S8 ; do
    >&2 echo "==> domain [${domain}]"

    if [ -e Targets/YggL/${domain}.replace.tsv ]; then
        continue;
    fi

    for ORDER in $(cat summary/order.lst); do
        >&2 echo "==> ORDER [${ORDER}]"

        cat taxon/${ORDER} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw Targets/YggL/HMM/${domain}.hmm - |
                    grep '>>' |
                    perl -nl -e '
                        m{>>\s+(\S+)} or next;
                        \$n = \$1;
                        \$s = \$n;
                        \$s =~ s/\.\d+//;
                        printf qq{%s\t%s_%s\n}, \$n, {}, \$s;
                    '
            "
    done \
        > Targets/YggL/${domain}.replace.tsv

    >&2 echo
done

tsv-join Targets/YggL/YggL_50S_bp.replace.tsv \
    -f Targets/YggL/PTHR38778.replace.tsv \
    > Targets/YggL/YggL.replace.tsv

wc -l Targets/YggL/*.tsv
#  2397 Targets/YggL/PTHR38778.replace.tsv
#  2721 Targets/YggL/Ribosomal_L10.replace.tsv
#  2779 Targets/YggL/Ribosomal_S8.replace.tsv
#  2396 Targets/YggL/YggL.replace.tsv
#  2401 Targets/YggL/YggL_50S_bp.replace.tsv

faops some PROTEINS/all.replace.fa.gz <(tsv-select -f 2 Targets/YggL/YggL.replace.tsv) Targets/YggL/YggL.fa

muscle -in Targets/YggL/YggL.fa -out Targets/YggL/YggL.aln.fa

FastTree Targets/YggL/YggL.aln.fa > Targets/YggL/YggL.aln.newick

nw_reroot Targets/YggL/YggL.aln.newick $(nw_labels Targets/YggL/YggL.aln.newick | grep -E "Baci_subti|Sta_aure") |
    nw_order -c n - \
    > Targets/YggL/YggL.reoot.newick

```
