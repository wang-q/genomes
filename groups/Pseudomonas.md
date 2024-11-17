# *Pseudomonas* HGT

<!-- toc -->

- [Strain info](#strain-info)
    * [Symlink](#symlink)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
- [All assemblies](#all-assemblies)
    * [List strains of the target genus and remove abnormal strains](#list-strains-of-the-target-genus-and-remove-abnormal-strains)
    * [Extract from `../Bacteria`](#extract-from-bacteria)
    * [Count strains - Genus](#count-strains---genus)
    * [Typical strains](#typical-strains)
    * [Rsync to hpcc](#rsync-to-hpcc)
- [NCBI taxonomy](#ncbi-taxonomy)
- [Collect proteins](#collect-proteins)
    * [`all.pro.fa`](#allprofa)
    * [`all.replace.fa`](#allreplacefa)
    * [`all.info.tsv`](#allinfotsv)
- [Phylogenetics with bac120](#phylogenetics-with-bac120)
    * [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
    * [Tweak the concat tree](#tweak-the-concat-tree)
- [InterProScan on all proteins of representative and typical strains](#interproscan-on-all-proteins-of-representative-and-typical-strains)
- [Count members of protein families](#count-members-of-protein-families)
    * [Pseudomonas](#pseudomonas)
    * [P. aeruginosa](#p-aeruginosa)
    * [P. putida](#p-putida)
    * [Acinetobacter](#acinetobacter)
    * [A. baumannii](#a-baumannii)
- [Protein families](#protein-families)
    * [IPR007416 - YggL 50S ribosome-binding protein](#ipr007416---yggl-50s-ribosome-binding-protein)

<!-- tocstop -->

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
    # Pseudomonadales
    Pseudomonas
    Halopseudomonas
    Stutzerimonas
    Azotobacter

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

cat ../Bacteria/summary/collect.pass.tsv | # 65357
    nwr restrict ${GENUS[*]} -f stdin -c 3 | # restrict to these genera 8877
    tsv-filter -H --le "C:20" --ge "N50:500000" | # more stringent parameters 4047
    sed '1d' |
    tsv-join -e -f ../Bacteria/ASSEMBLY/omit.lst -k 1 | # 4047
    tsv-join -e -f ../Bacteria/MinHash/abnormal.lst -k 1 | # 3892
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
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# .lst and .count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    keep-header -- tsv-sort -k1,1 |
    tsv-filter -H --ge 3:2 |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```

| item    | count |
|---------|------:|
| strain  |  3904 |
| species |   172 |
| genus   |    20 |
| family  |    15 |
| order   |    11 |
| class   |     7 |

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Acinetobacter    |       22 |      830 |
| Azotobacter      |        2 |        6 |
| Bordetella       |        9 |      807 |
| Burkholderia     |       26 |      605 |
| Escherichia      |        1 |        2 |
| Pseudomonas      |       81 |     1280 |
| Serratia         |       11 |      191 |
| Stenotrophomonas |        5 |      126 |
| Stutzerimonas    |        4 |       46 |

## MinHash

```shell
cd ~/data/Pseudomonas/

nwr template summary/assembly.tsv \
    --mh \
    --parallel 16 \
    --ani-ab 0.12 \
    --ani-nr 0.01

# Compute assembly sketches
bash MinHash/compute.sh

# Distances within species
bash MinHash/species.sh

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#8

# Non-redundant strains within species
bash MinHash/nr.sh

find MinHash -name "mash.dist.tsv" -size +0 | wc -l
#162

find MinHash -name "redundant.lst" -size +0 | wc -l
#97

find MinHash -name "redundant.lst" -empty | wc -l
#64

find MinHash -name "NR.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/NR.lst
wc -l summary/NR.lst
#874

# All representative should be in NR
cat summary/assembly.tsv |
    tsv-join -f ASSEMBLY/rep.lst -k 1 |
    cut -f 1 |
    grep -v -F -f summary/NR.lst

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
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

cat Count/genus.count.tsv |
    tsv-filter -H --ge "3:2" |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

# Can accept N_COUNT
bash Count/lineage.sh 10

cat Count/lineage.count.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| order            | #species | #strains |
|------------------|---------:|---------:|
| Bacillales       |        3 |        3 |
| Burkholderiales  |       35 |     1412 |
| Enterobacterales |       15 |      196 |
| Moraxellales     |       22 |      830 |
| Pseudomonadales  |       87 |     1324 |
| Xanthomonadales  |        5 |      126 |

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Acinetobacter    |       22 |      830 |
| Azotobacter      |        2 |        6 |
| Bordetella       |        9 |      807 |
| Burkholderia     |       26 |      605 |
| Escherichia      |        1 |        2 |
| Pseudomonas      |       81 |     1276 |
| Serratia         |       11 |      191 |
| Stenotrophomonas |        5 |      126 |
| Stutzerimonas    |        4 |       42 |

| #family          | genus            | species                      | count |
|------------------|------------------|------------------------------|------:|
| Alcaligenaceae   | Bordetella       | Bordetella bronchiseptica    |    25 |
|                  |                  | Bordetella hinzii            |    18 |
|                  |                  | Bordetella holmesii          |    66 |
|                  |                  | Bordetella parapertussis     |    90 |
|                  |                  | Bordetella pertussis         |   593 |
| Burkholderiaceae | Burkholderia     | Burkholderia ambifaria       |    15 |
|                  |                  | Burkholderia cenocepacia     |    87 |
|                  |                  | Burkholderia cepacia         |    19 |
|                  |                  | Burkholderia contaminans     |    16 |
|                  |                  | Burkholderia gladioli        |    22 |
|                  |                  | Burkholderia glumae          |    51 |
|                  |                  | Burkholderia mallei          |    53 |
|                  |                  | Burkholderia multivorans     |   102 |
|                  |                  | Burkholderia pseudomallei    |   152 |
|                  |                  | Burkholderia thailandensis   |    25 |
|                  |                  | Burkholderia vietnamiensis   |    18 |
| Moraxellaceae    | Acinetobacter    | Acinetobacter baumannii      |   601 |
|                  |                  | Acinetobacter haemolyticus   |    15 |
|                  |                  | Acinetobacter indicus        |    21 |
|                  |                  | Acinetobacter johnsonii      |    18 |
|                  |                  | Acinetobacter junii          |    11 |
|                  |                  | Acinetobacter nosocomialis   |    17 |
|                  |                  | Acinetobacter pittii         |    59 |
|                  |                  | Acinetobacter seifertii      |    25 |
| Pseudomonadaceae | Pseudomonas      | Pseudomonas aeruginosa       |   706 |
|                  |                  | Pseudomonas amygdali         |    14 |
|                  |                  | Pseudomonas chlororaphis     |   101 |
|                  |                  | Pseudomonas fluorescens      |    28 |
|                  |                  | Pseudomonas protegens        |    24 |
|                  |                  | Pseudomonas putida           |    77 |
|                  |                  | Pseudomonas synxantha        |    10 |
|                  |                  | Pseudomonas syringae         |    69 |
|                  |                  | Pseudomonas viridiflava      |    21 |
|                  | Stutzerimonas    | Stutzerimonas stutzeri       |    31 |
| Xanthomonadaceae | Stenotrophomonas | Stenotrophomonas maltophilia |   116 |
| Yersiniaceae     | Serratia         | Serratia marcescens          |   125 |
|                  |                  | Serratia plymuthica          |    18 |
|                  |                  | Serratia ureilytica          |    15 |

### Count strains - Genus

```shell
cd ~/data/Pseudomonas

cat summary/genus.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(
            cat summary/collect.pass.csv |
                sed "1d" |
                tsv-select -d, -f 3 |
                nwr append stdin -r genus -r species |
                grep -w {} |
                tsv-select -f 1,3 |
                rgr dedup stdin |
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
    tsv-sort -k2,2 | #    tsv-filter --ge 4:100 |
    (echo -e '#tax_id\tgenus\t#species\t#strains\t#NR' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | genus            | #species | #strains | #NR |
|---------|------------------|----------|----------|-----|
| 469     | Acinetobacter    | 51       | 791      | 240 |
| 352     | Azotobacter      | 5        | 6        | 4   |
| 517     | Bordetella       | 23       | 807      | 12  |
| 32008   | Burkholderia     | 103      | 692      | 123 |
| 286     | Pseudomonas      | 233      | 1418     | 495 |
| 613     | Serratia         | 21       | 191      | 69  |
| 40323   | Stenotrophomonas | 11       | 246      | 129 |
| 2901164 | Stutzerimonas    | 8        | 31       | 28  |

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
    Pseudom_aeruginosa_PAK_GCF_000568855_2 \
    Pseudom_aeruginosa_GCF_011466835_1 \
    Pseudom_aeruginosa_LESB58_GCF_000026645_1 \
    Pseudom_viridif_GCF_900184295_1 \
    Pseudom_syringae_pv_tomato_DC3000_GCF_000007805_1 \
    Pseudom_syringae_pv_syringae_B728a_GCF_000012245_1 \
    Pseudom_fluo_GCF_900215245_1 \
    Pseudom_putida_KT2440_GCF_000007565_2 \
    Pseudom_putida_NBRC_14164_GCF_000412675_1 \
    Stu_stut_A1501_GCF_000013785_1 \
    Pseudom_chl_aureofaciens_30_84_GCF_000281915_1 \
    Pseudom_amyg_pv_tabaci_ATCC_11528_GCF_000145945_2 \
    Pseudom_proteg_Pf_5_GCF_000012265_1 \
    Pseudom_proteg_CHA0_GCF_900560965_1 \
    Pseudom_entomophila_L48_GCF_000026105_1 \
    Acin_bau_GCF_008632635_1 \
    Acin_bau_ATCC_17978_GCF_004794235_2 \
    Acin_bau_ATCC_19606_CIP_70_34_JCM_6841_GCF_019331655_1 \
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

Call the pipeline script.

```shell
cd ~/data/Pseudomonas

bash ~/Scripts/genomes/bin/pl_collect_protein.sh

cat PROTEINS/counts.tsv |
    mlr --itsv --omd cat

```

## Phylogenetics with bac120

### Find corresponding proteins by `hmmsearch`

```shell
E_VALUE=1e-20

cd ~/data/Pseudomonas

# Find all genes
for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    >&2 echo "==> marker [${marker}]"

    mkdir -p PROTEINS/${marker}

    cat summary/NR.lst |
        parallel --no-run-if-empty --linebuffer -k -j 8 "
            gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/HMM/bac120/HMM/${marker}.HMM - |
                grep '>>' |
                perl -nl -e ' m{>>\s+(\S+)} and printf qq{%s\t%s\n}, \$1, {}; '
        " \
        > PROTEINS/${marker}/replace.tsv

    >&2 echo

done

```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Pseudomonas

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat PROTEINS/{}/replace.tsv |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
#2196    2200    3326.5

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        cat PROTEINS/{}/replace.tsv |
            wc -l
    ' |
    paste - - |
    tsv-filter --invert --ge 2:1800 --le 2:2600 |
    cut -f 1 \
    > PROTEINS/bac120.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process
cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    grep -v -Fx -f PROTEINS/bac120.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        cat PROTEINS/{}/replace.tsv \
            > PROTEINS/{}/{}.replace.tsv

        faops some PROTEINS/all.uniq.fa.gz <(
            cat PROTEINS/{}/{}.replace.tsv |
                cut -f 1 |
                rgr dedup stdin
            ) stdout \
            > PROTEINS/{}/{}.pro.fa
    '

# Align each markers with muscle
cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        >&2 echo "==> marker [{}]"
        if [ ! -s PROTEINS/{}/{}.pro.fa ]; then
            exit
        fi
        if [ -s PROTEINS/{}/{}.aln.fa ]; then
            exit
        fi

        muscle -quiet -in PROTEINS/{}/{}.pro.fa -out PROTEINS/{}/{}.aln.fa
    '

for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    >&2 echo "==> marker [${marker}]"
    if [ ! -s PROTEINS/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sometimes `muscle` can not produce alignments
    if [ ! -s PROTEINS/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # 1 name to many names
    cat PROTEINS/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s PROTEINS/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > PROTEINS/${marker}/${marker}.replace.fa
done

# Concat marker genes
for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    if [ ! -s PROTEINS/${marker}/${marker}.pro.fa ]; then
        continue
    fi
    if [ ! -s PROTEINS/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # sequences in one line
    faops filter -l 0 PROTEINS/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > PROTEINS/bac120.aln.fas

fasops concat PROTEINS/bac120.aln.fas summary/NR.lst -o PROTEINS/bac120.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in PROTEINS/bac120.aln.fa -out PROTEINS/bac120.trim.fa -automated1

faops size PROTEINS/bac120.*.fa |
    rgr dedup stdin -f 2 |
    cut -f 2
#30002
#24150

# To make it faster
FastTree -fastest -noml PROTEINS/bac120.trim.fa > PROTEINS/bac120.trim.newick

```

### Tweak the concat tree

```shell
cd ~/data/Pseudomonas/tree

nwr order --nd --an ../PROTEINS/bac120.trim.newick \
    > bac120.order.newick

nwr pl-condense -r species \
    bac120.order.newick ../Count/species.tsv \
    -o bac120.condensed.newick

mv condensed.tsv bac120.condense.tsv

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 bac120.condensed.newick |
    rsvg-convert -o Pseudomonas.bac120.png

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
