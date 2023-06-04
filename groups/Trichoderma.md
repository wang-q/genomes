# Build alignments across a eukaryotic taxonomy rank

Genus *Trichoderma* as an example.

<!-- toc -->

- [Strain info](#strain-info)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
- [Download all assemblies](#download-all-assemblies)
    * [Create assembly.tsv](#create-assemblytsv)
    * [Download and check](#download-and-check)
    * [Rsync to hpcc](#rsync-to-hpcc)
- [BioSample](#biosample)
- [Count species and strains](#count-species-and-strains)
- [MinHash](#minhash)
    * [Tweak the minhash tree](#tweak-the-minhash-tree)
- [Groups and targets](#groups-and-targets)
- [Prepare sequences for `egaz`](#prepare-sequences-for-egaz)
- [Generate alignments](#generate-alignments)

<!-- tocstop -->

## Strain info

* [Trichoderma](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=5543)
* [Entrez records](http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5543)
* [WGS](https://www.ncbi.nlm.nih.gov/Traces/wgs/?view=wgs&search=Trichoderma) is now useless and can
  be ignored.

A nice review article about Trichoderma:

Woo, S.L. et al.
Trichoderma: a multipurpose, plant-beneficial microorganism for eco-sustainable agriculture.
Nat Rev Microbiol 21, 312–326 (2023).
https://doi.org/10.1038/s41579-022-00819-5

### List all ranks

There are no noteworthy classification ranks other than species.

```shell

nwr member Trichoderma |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat

nwr lineage Trichoderma |
    tsv-filter --str-ne 1:clade |
    tsv-filter --str-ne "1:no rank" |
    sed -n '/kingdom\tFungi/,$p' |
    sed -E "s/\b(genus)\b/*\1*/"| # Highlight genus
    (echo -e '#rank\tsci_name\ttax_id' && cat) |
    mlr --itsv --omd cat

```

| rank     | count |
|----------|-------|
| genus    | 1     |
| species  | 443   |
| no rank  | 1     |
| varietas | 2     |
| strain   | 14    |
| forma    | 2     |

| #rank      | sci_name          | tax_id |
|------------|-------------------|--------|
| kingdom    | Fungi             | 4751   |
| subkingdom | Dikarya           | 451864 |
| phylum     | Ascomycota        | 4890   |
| subphylum  | Pezizomycotina    | 147538 |
| class      | Sordariomycetes   | 147550 |
| subclass   | Hypocreomycetidae | 222543 |
| order      | Hypocreales       | 5125   |
| family     | Hypocreaceae      | 5129   |
| *genus*    | Trichoderma       | 5543   |

### Species with assemblies

Check also the family Hypocreaceae for outgroups.

```shell
mkdir -p ~/data/Trichoderma/summary
cd ~/data/Trichoderma/summary

# should have a valid name of genus
nwr member Hypocreaceae -r genus |
    grep -v -i "Candidatus " |
    grep -v -i "candidate " |
    grep -v " sp." |
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#18 genus.list

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
#   8 RS1.tsv
#  37 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     8
#GB1     113

```

## Download all assemblies

### Create assembly.tsv

If a refseq assembly is available, the corresponding genbank one is not downloaded

```shell
cd ~/data/Trichoderma/summary

echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus IN ('Saccharomyces')
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
#115 lines, 6 fields

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
    > Trichoderma.assembly.tsv

datamash check < Trichoderma.assembly.tsv
#114 lines, 5 fields

# find potential duplicate strains or assemblies
cat Trichoderma.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Trichoderma.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

# Edit .assembly.tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Trichoderma.assembly.tsv
#
# Save the file to another directory to prevent accidentally changing it
# cp Trichoderma.assembly.tsv ~/Scripts/genomes/assembly

# Cleaning
rm raw*.*sv

```

### Download and check

* When `rsync.sh` is interrupted, run `check.sh` before restarting
* For projects that have finished downloading, but have renamed strains, you can run `reorder.sh` to
  avoid re-downloading
* The parameters of `n50.sh` should be determined by the distribution of the description statistics
* `collect.sh` generates a file of type `.csv`, which is intended to be opened by spreadsheet
  software.
* `finish.sh` generates the following files
    * `omit.lst` - no annotations
    * `collect.pass.csv` - passes the n50 check
    * `rep.lst` - representative or reference strains
    * `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class
    * `counts.tsv`

```shell
cd ~/data/Trichoderma

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --ass \
    -o .

# Run
bash ASSEMBLY/rsync.sh

# Check md5; create check.lst
# rm ASSEMBLY/check.lst
bash ASSEMBLY/check.sh

# Put the misplaced directory into the right place
#bash ASSEMBLY/reorder.sh
#
# This operation will delete some files in the directory, so please be careful
#cat ASSEMBLY/remove.list |
#    parallel --no-run-if-empty --linebuffer -k -j 1 '
#        if [[ -e "ASSEMBLY/{}" ]]; then
#            echo Remove {}
#            rm -fr "ASSEMBLY/{}"
#        fi
#    '

# N50 C S; create n50.tsv and n50.pass.csv
bash ASSEMBLY/n50.sh 100000 1000 1000000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50,S" --max "C"
#N50_min S_min   C_max
#697391  33215161        533

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "S:0.1,0.5" --quantile "N50:0.1,0.5"  --quantile "C:0.5,0.9"
#S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
#32255196.8      37316984        98936.2 1332095 167     1276.4

# After the above steps are completed, run the following commands.

# Collect; create collect.csv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.csv summary/
cp ASSEMBLY/strains.taxon.tsv summary/

cat ASSEMBLY/counts.tsv |
    mlr --itsv --omd cat |
    sed 's/-\s*|$/-:|/'

```

| #item             | count |
|-------------------|------:|
| url.tsv           |   113 |
| check.lst         |   112 |
| collect.csv       |   113 |
| n50.tsv           |   114 |
| n50.pass.csv      |    97 |
| collect.pass.csv  |    96 |
| omit.lst          |    81 |
| rep.lst           |    30 |
| strains.taxon.tsv |    95 |

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Trichoderma/ \
    wangq@202.119.37.251:data/Trichoderma

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Trichoderma/ \
    wangq@58.213.64.36:data/Trichoderma

# rsync -avP wangq@202.119.37.251:data/Trichoderma/ ~/data/Trichoderma

# rsync -avP -e 'ssh -p 8804' wangq@58.213.64.36:data/Trichoderma/ ~/data/Trichoderma

```

## BioSample

```shell
cd ~/data/Trichoderma

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --bs \
    -o .

bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 10

datamash check < BioSample/biosample.tsv
#111 lines, 37 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## Count species and strains

```shell
cd ~/data/Trichoderma/

cat summary/collect.pass.csv |
    sed -e '1d' |
    tr "," "\t" |
    tsv-select -f 1,3 |
    nwr append stdin -c 2 -r species -r genus -r family -r order \
    > summary/strains.taxon.tsv

bash ~/Scripts/genomes/bin/taxon_count.sh summary/strains.taxon.tsv 1

```

| #family            | genus         | species                     | count |
|--------------------|---------------|-----------------------------|------:|
| Hypocreaceae       | Cladobotryum  | Cladobotryum protrusum      |     1 |
|                    | Escovopsis    | Escovopsis weberi           |     2 |
|                    | Hypomyces     | Hypomyces perniciosus       |     1 |
|                    |               | Hypomyces rosellus          |     1 |
|                    | Trichoderma   | Trichoderma afroharzianum   |     4 |
|                    |               | Trichoderma arundinaceum    |     4 |
|                    |               | Trichoderma asperelloides   |     2 |
|                    |               | Trichoderma asperellum      |    13 |
|                    |               | Trichoderma atrobrunneum    |     1 |
|                    |               | Trichoderma atroviride      |     7 |
|                    |               | Trichoderma breve           |     1 |
|                    |               | Trichoderma brevicrassum    |     1 |
|                    |               | Trichoderma citrinoviride   |     3 |
|                    |               | Trichoderma erinaceum       |     2 |
|                    |               | Trichoderma gamsii          |     2 |
|                    |               | Trichoderma gracile         |     1 |
|                    |               | Trichoderma guizhouense     |     1 |
|                    |               | Trichoderma hamatum         |     1 |
|                    |               | Trichoderma harzianum       |    10 |
|                    |               | Trichoderma koningii        |     1 |
|                    |               | Trichoderma koningiopsis    |     4 |
|                    |               | Trichoderma lentiforme      |     1 |
|                    |               | Trichoderma longibrachiatum |     5 |
|                    |               | Trichoderma pseudokoningii  |     1 |
|                    |               | Trichoderma reesei          |    13 |
|                    |               | Trichoderma semiorbis       |     1 |
|                    |               | Trichoderma simmonsii       |     1 |
|                    |               | Trichoderma virens          |     8 |
|                    |               | Trichoderma viride          |     1 |
| Saccharomycetaceae | Saccharomyces | Saccharomyces cerevisiae    |     1 |

## MinHash

Estimate nucleotide divergences among strains.

* Abnormal strains
    * If the maximum value of ANI between strains within a species is greater than *0.05*, the
      median and maximum value will be reported. Strains that cannot be linked by median ANI will be
      considered as abnormal strains.
    * It may consist of two scenarios:
        1. Wrong species identification
        2. Poor assembly quality

* Non-redundant strains
    * If the ANI value between two strains within a species is less than *0.005*, the two strains
      are considered to be redundant.
    * Need these files:  representative.lst and omit.lst

* MinHash tree
    * A rough tree is generated by k-mean clustering.

* These abnormal strains should be manually checked to determine whether to include them in the
  subsequent steps.

```shell
cd ~/data/Trichoderma

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --mh \
    --pass \
    --ani_ab 0.05 \
    --ani_nr 0.005 \
    --height 0.4 \
    -o .

# Compute assembly sketches
bash MinHash/compute.sh

# Distances within species
bash MinHash/species.sh

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst
#T_har_CGMCC_20739_GCA_019097725_1
#T_har_Tr1_GCA_002894145_1
#T_har_ZL_811_GCA_021186515_1

# Non-redundant strains within species
bash MinHash/nr.sh

# Distances between all selected sketches, then hierarchical clustering
bash MinHash/dist.sh

```

### Tweak the minhash tree

```shell
mkdir -p ~/data/Trichoderma/tree
cd ~/data/Trichoderma/tree

nw_reroot ../MinHash/tree.nwk Sa_cer_S288C |
    nw_order -c n - \
    > minhash.reroot.newick

# rank::col
ARRAY=(
#    'order::6'
#    'family::5'
#    'genus::4'
    'species::3'
)

rm minhash.condensed.map
CUR_TREE=minhash.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/genomes/bin/condense_tree.sh ${CUR_TREE} ../ASSEMBLY/strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick minhash.${GROUP_NAME}.newick
    cat condense.map >> minhash.condensed.map

    CUR_TREE=minhash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 30 minhash.species.newick |
    rsvg-convert -o Trichoderma.minhash.png

```

## Groups and targets

Review `ASSEMBLY/Trichoderma.assembly.pass.csv` and `MinHash/groups.tsv`.

Create `ARRAY` manually with a format `group::target`.

Target selecting criteria:

* Prefer Sanger sequenced assemblies
* RefSeq_category with `Representative Genome`
* Assembly_level with `Complete Genome` or `Chromosome`

```shell
mkdir -p ~/data/alignment/Trichoderma/taxon
cd ~/data/alignment/Trichoderma/taxon

cp ../mash/tree.nwk .
cp ../mash/groups.tsv .

echo -e "#Serial\tGroup\tCount\tTarget" > group_target.tsv

# groups accroding `groups.tsv`
ARRAY=(
    'C_E_H::H_ros_GCA_011799845_1'
    'T_har_vire::T_vire_Gv29_8_GCA_000170995_2'
    'T_asperellum_atrov::T_atrov_IMI_206040_GCA_000171015_2'
    'T_lon_ree::T_ree_QM6a_GCA_000167675_2'
)

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

    echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}" >> group_target.tsv

done

# Custom groups
ARRAY=(
    'Trichoderma::T_ree_QM6a_GCA_000167675_2'
    'Trichoderma_reesei::T_ree_QM6a_GCA_000167675_2'
    'Trichoderma_asperellum::T_asperellum_CBS_433_97_GCA_003025105_1'
    'Trichoderma_harzianum::T_har_CBS_226_95_GCA_003025095_1'
    'Trichoderma_atroviride::T_atrov_IMI_206040_GCA_000171015_2'
    'Trichoderma_virens::T_vire_Gv29_8_GCA_000170995_2'
)

SERIAL=100
for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"

    SERIAL=$((SERIAL + 1))
    GROUP_NAME_2=$(echo $GROUP_NAME | tr "_" " ")

    if [ "$GROUP_NAME" = "Trichoderma" ]; then
        cat ../ASSEMBLY/Trichoderma.assembly.pass.csv |
            tsv-filter -H -d"," --not-blank RefSeq_category |
            sed '1d' |
            cut -d"," -f 1 \
            > ${GROUP_NAME}
        echo "C_pro_GCA_004303015_1" >> ${GROUP_NAME}
        echo "E_web_GCA_003055145_1" >> ${GROUP_NAME}
        echo "H_per_GCA_008477525_1" >> ${GROUP_NAME}
        echo "H_ros_GCA_011799845_1" >> ${GROUP_NAME}
    else
        cat ../ASSEMBLY/Trichoderma.assembly.pass.csv |
            cut -d"," -f 1,2 |
            grep "${GROUP_NAME_2}" |
            cut -d"," -f 1 \
            > ${GROUP_NAME}
    fi

    COUNT=$(cat ${GROUP_NAME} | wc -l )

    echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}" >> group_target.tsv

done

mlr --itsv --omd cat group_target.tsv

cat ../ASSEMBLY/Trichoderma.assembly.pass.csv |
    tsv-filter -H -d, --str-eq Assembly_level:"Complete Genome" |
    tsv-select -H -d, -f name \
    > complete-genome.lst


```

| #Serial | Group                  | Count | Target                                  |
|---------|------------------------|-------|-----------------------------------------|
| 1       | C_E_H                  | 6     | H_ros_GCA_011799845_1                   |
| 2       | T_har_vire             | 20    | T_vire_Gv29_8_GCA_000170995_2           |
| 3       | T_asperellum_atrov     | 23    | T_atrov_IMI_206040_GCA_000171015_2      |
| 4       | T_lon_ree              | 16    | T_ree_QM6a_GCA_000167675_2              |
| 101     | Trichoderma            | 31    | T_ree_QM6a_GCA_000167675_2              |
| 102     | Trichoderma_reesei     | 10    | T_ree_QM6a_GCA_000167675_2              |
| 103     | Trichoderma_asperellum | 10    | T_asperellum_CBS_433_97_GCA_003025105_1 |
| 104     | Trichoderma_harzianum  | 7     | T_har_CBS_226_95_GCA_003025095_1        |
| 105     | Trichoderma_atroviride | 7     | T_atrov_IMI_206040_GCA_000171015_2      |
| 106     | Trichoderma_virens     | 6     | T_vire_Gv29_8_GCA_000170995_2           |

## Prepare sequences for `egaz`

* `--perseq` for Chromosome-level assemblies and targets
    * means split fasta by names, target or good assembles should set it
* `--species Fungi` specify the species or clade of this group for RepeatMasker

```shell
cd ~/data/alignment/Trichoderma/

# prep
egaz template \
    ASSEMBLY \
    --prep -o GENOMES \
    $( cat taxon/group_target.tsv | sed -e '1d' | cut -f 4 | parallel -j 1 echo " --perseq {} " ) \
    $( cat taxon/complete-genome.lst | parallel -j 1 echo " --perseq {} " ) \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--species Fungi --parallel 24"

bash GENOMES/0_prep.sh

# gff
for n in $(cat taxon/group_target.tsv | sed -e '1d' | cut -f 4 ) \
    $( cat taxon/complete-genome.lst ) \
    ; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"

    gzip -dc ${FILE_GFF} > GENOMES/${n}/chr.gff
done

```

## Generate alignments

```shell
cd ~/data/alignment/Trichoderma/

cat taxon/group_target.tsv |
    sed -e '1d' |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            --multi -o groups/{2}/ \
            --tree taxon/tree.nwk \
            --parallel 8 -v

        bash groups/{2}/1_pair.sh
        bash groups/{2}/3_multi.sh
    '

# clean
find groups -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find groups -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type f -name "output.*" | parallel -r rm

```

Use Tatr_IMI_206040 as target

* No results for Tatr_IMI_206040vsTkon_JCM_1883

Tatr_IMI_206040;qs=Tatr_XS2015,,

```bash
cd ~/data/alignment/trichoderma

egaz template \
    GENOMES/Tatr_IMI_206040 \
    GENOMES/Tree_QM6a \
    GENOMES/Tvir_Gv29_8 \
    --multi -o multi/ \
    --multiname sanger --order \
    --parallel 8 -v

bash multi/1_pair.sh
bash multi/3_multi.sh

egaz template \
    GENOMES/Tatr_IMI_206040 \
    $(find GENOMES -maxdepth 1 -mindepth 1 -type d -path "*/????*" | grep -v "Tatr_IMI_206040"| grep -v "Tkon_JCM_1883") \
    --multi -o multi/ \
    --tree mash/tree.nwk \
    --parallel 8 -v

bash multi/1_pair.sh
bash multi/3_multi.sh

```
