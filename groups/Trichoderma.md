# Build alignments across a eukaryotic taxonomy rank

Genus *Trichoderma* as an example.

<!-- toc -->

- [Strain info](#strain-info)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
- [Download all assemblies](#download-all-assemblies)
    * [Create assembly.tsv](#create-assemblytsv)
    * [Count before download](#count-before-download)
    * [Download and check](#download-and-check)
    * [Rsync to hpcc](#rsync-to-hpcc)
- [BioSample](#biosample)
- [MinHash](#minhash)
    * [Condense branches in the minhash tree](#condense-branches-in-the-minhash-tree)
- [Count valid species and strains for *protein
  families*](#count-valid-species-and-strains-for-protein-families)
- [Collect proteins](#collect-proteins)
- [Phylogenetics with fungi61](#phylogenetics-with-fungi61)
    * [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
- [Count valid species and strains for *genomic
  alignments*](#count-valid-species-and-strains-for-genomic-alignments)
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
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

nwr lineage Trichoderma |
    tsv-filter --str-ne 1:clade |
    tsv-filter --str-ne "1:no rank" |
    sed -n '/kingdom\tFungi/,$p' |
    sed -E "s/\b(genus)\b/*\1*/"| # Highlight genus
    (echo -e '#rank\tsci_name\ttax_id' && cat) |
    mlr --itsv --omd cat

```

| rank     | count |
|----------|------:|
| genus    |     1 |
| species  |   443 |
| no rank  |     1 |
| varietas |     2 |
| strain   |    14 |
| forma    |     2 |

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

### Count before download

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class

```shell
cd ~/data/Trichoderma

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --count \
    --rank genus

# strains.taxon.tsv
bash Count/strains.sh

# genus.lst and genus.count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    mlr --itsv --omd cat

```

| genus            | #species | #strains |
|------------------|----------|----------|
| Cladobotryum     | 1        | 1        |
| Escovopsis       | 1        | 2        |
| Hypomyces        | 2        | 2        |
| Mycogone         | 1        | 1        |
| Saccharomyces    | 1        | 1        |
| Sphaerostilbella | 1        | 1        |
| Trichoderma      | 31       | 105      |

### Download and check

* When `rsync.sh` is interrupted, run `check.sh` before restarting
* For projects that have finished downloading, but have renamed strains, you can run `reorder.sh` to
  avoid re-downloading
    * `misplaced.tsv`
    * `remove.list`
* The parameters of `n50.sh` should be determined by the distribution of the description statistics
* `collect.sh` generates a file of type `.tsv`, which is intended to be opened by spreadsheet
  software.
* `finish.sh` generates the following files
    * `omit.lst` - no annotations
    * `collect.pass.tsv` - passes the n50 check
    * `pass.lst` - passes the n50 check
    * `rep.lst` - representative or reference strains
    * `counts.tsv`

```shell
cd ~/data/Trichoderma

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --ass

# Run
bash ASSEMBLY/rsync.sh

# Check md5; create check.lst
# rm ASSEMBLY/check.lst
bash ASSEMBLY/check.sh

# Put the misplaced directory into the right place
#bash ASSEMBLY/reorder.sh
#
# This operation will delete some files in the directory, so please be careful
#cat ASSEMBLY/remove.lst |
#    parallel --no-run-if-empty --linebuffer -k -j 1 '
#        if [[ -e "ASSEMBLY/{}" ]]; then
#            echo Remove {}
#            rm -fr "ASSEMBLY/{}"
#        fi
#    '

# N50 C S; create n50.tsv and n50.pass.tsv
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

# Collect; create collect.tsv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```

| #item            | fields | lines |
|------------------|-------:|------:|
| url.tsv          |      3 |   113 |
| check.lst        |      1 |   113 |
| collect.tsv      |     20 |   114 |
| n50.tsv          |      4 |   114 |
| n50.pass.tsv     |      4 |    97 |
| collect.pass.tsv |     23 |    97 |
| pass.lst         |      1 |    96 |
| omit.lst         |      1 |    81 |
| rep.lst          |      1 |    31 |

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
    --bs

bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 10

datamash check < BioSample/biosample.tsv
#111 lines, 37 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## MinHash

Estimate nucleotide divergences among strains.

* Abnormal strains
    * If the maximum value of ANI between strains within a species is greater than *0.05*, the
      median and maximum value will be reported. Strains that cannot be linked by the median
      ANI, e.g., have no similar strains in the species, will be considered as abnormal strains.
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
    --parallel 16 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
    --ani-nr 0.005 \
    --height 0.4

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

### Condense branches in the minhash tree

* This phylo-tree is not really formal/correct, and shouldn't be used to interpret phylogenetic
  relationships
* It is just used to find more abnormal strains

```shell
mkdir -p ~/data/Trichoderma/tree
cd ~/data/Trichoderma/tree

nw_reroot ../MinHash/tree.nwk Sa_cer_S288C |
    nw_order -c n - \
    > minhash.reroot.newick

# rank::col
ARRAY=(
#    'order::5'
#    'family::4'
#    'genus::3'
    'species::2'
)

rm minhash.condensed.map
CUR_TREE=minhash.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/genomes/bin/condense_tree.sh ${CUR_TREE} ../Count/strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick minhash.${GROUP_NAME}.newick
    cat condense.map >> minhash.condensed.map

    CUR_TREE=minhash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 minhash.species.newick |
    rsvg-convert -o Trichoderma.minhash.png

```

## Count valid species and strains for *protein families*

```shell
cd ~/data/Trichoderma/

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst \
    --rank genus

# strains.taxon.tsv
bash Count/strains.sh

# genus.lst genus.count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    mlr --itsv --omd cat

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv

```

| genus         | #species | #strains |
|---------------|----------|----------|
| Escovopsis    | 1        | 1        |
| Saccharomyces | 1        | 1        |
| Trichoderma   | 16       | 24       |

## Collect proteins

```shell
cd ~/data/Trichoderma/

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --pro \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst

# collect proteins
bash Protein/collect.sh

cat Protein/counts.tsv |
    mlr --itsv --omd cat

```

| #item                          | count   |
|--------------------------------|---------|
| Proteins                       | 275,985 |
| Unique headers and annotations | 275,985 |
| Unique proteins                | 275,985 |
| all.replace.fa                 | 275,985 |
| all.annotation.tsv             | 275,986 |
| all.info.tsv                   | 275,986 |

## Phylogenetics with fungi61

### Find corresponding proteins by `hmmsearch`

* 61 fungal marker genes
    * Ref.: https://doi.org/10.1093/nar/gkac894

* The `E_VALUE` was manually adjusted to 1e-20 to reach a balance between sensitivity and
  speciality.

```shell
cd ~/data/Trichoderma

# The fungi61 HMM set
nwr kb fungi61 -o HMM
cp HMM/fungi61.lst HMM/marker.lst

E_VALUE=1e-20

# Find all genes
for marker in $(cat HMM/marker.lst); do
    echo >&2 "==> marker [${marker}]"

    mkdir -p Protein/${marker}

    cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/hmm/${marker}.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/${marker}/replace.tsv

    echo >&2
done

```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Trichoderma

cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat Protein/{}/replace.tsv |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
#25      25      56

cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        cat Protein/{}/replace.tsv |
            wc -l
    ' |
    paste - - |
    tsv-filter --invert --ge 2:20 --le 2:30 |
    cut -f 1 \
    > Protein/marker.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process
cat HMM/marker.lst |
    grep -v -Fx -f Protein/marker.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

        cat Protein/{}/replace.tsv \
            > Protein/{}/{}.replace.tsv

        faops some Protein/all.uniq.fa.gz <(
            cat Protein/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > Protein/{}/{}.pro.fa
    '

# Align each markers with muscle
cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo >&2 "==> marker [{}]"
        if [ ! -s Protein/{}/{}.pro.fa ]; then
            exit
        fi
        if [ -s Protein/{}/{}.aln.fa ]; then
            exit
        fi

        muscle -quiet -in Protein/{}/{}.pro.fa -out Protein/{}/{}.aln.fa
    '

for marker in $(cat HMM/marker.lst); do
    echo >&2 "==> marker [${marker}]"
    if [ ! -s Protein/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sometimes `muscle` can not produce alignments
    if [ ! -s Protein/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # 1 name to many names
    cat Protein/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s Protein/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > Protein/${marker}/${marker}.replace.fa
done

# Concat marker genes
for marker in $(cat HMM/marker.lst); do
    if [ ! -s Protein/${marker}/${marker}.pro.fa ]; then
        continue
    fi
    if [ ! -s Protein/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # sequences in one line
    faops filter -l 0 Protein/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > Protein/fungi61.aln.fas

cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -e -f MinHash/abnormal.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
    cut -f 1 |
    fasops concat Protein/fungi61.aln.fas stdin -o Protein/fungi61.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Protein/fungi61.aln.fa -out Protein/fungi61.trim.fa -automated1

faops size Protein/fungi61.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#28706
#20432

# To make it faster
FastTree -fastest -noml Protein/fungi61.trim.fa > Protein/fungi61.trim.newick

nw_reroot Protein/fungi61.trim.newick Sa_cer_S288C |
    nw_order -c n - \
    > Protein/fungi61.reroot.newick

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 Protein/fungi61.reroot.newick |
    rsvg-convert -o tree/Trichoderma.marker.png

```

## Count valid species and strains for *genomic alignments*

```shell
cd ~/data/Trichoderma/

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv
bash Count/strains.sh

# genus.lst genus.count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    mlr --itsv --omd cat

bash Count/lineage.sh

cat Count/lineage.count.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| genus         | #species | #strains |
|---------------|----------|----------|
| Cladobotryum  | 1        | 1        |
| Escovopsis    | 1        | 2        |
| Hypomyces     | 2        | 2        |
| Saccharomyces | 1        | 1        |
| Trichoderma   | 26       | 87       |

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
|                    |               | Trichoderma harzianum       |     7 |
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

## Groups and targets

Grouping criteria:

* The mash tree and the marker protein tree
* `MinHash/groups.tsv`

Target selecting criteria:

* `ASSEMBLY/collect.pass.tsv`
* Prefer Sanger sequenced assemblies
* RefSeq_category with `Representative Genome`
* Assembly_level with `Complete Genome` or `Chromosome`

Create a Bash `ARRAY` manually with a format of `group::target`.

```shell
mkdir -p ~/data/Trichoderma/taxon
cd ~/data/Trichoderma/taxon

cat ../ASSEMBLY/collect.pass.tsv |
    tsv-filter -H --str-eq annotations:Yes --le C:100 |
    tsv-select -H -f name,Assembly_level,Genome_coverage,Sequencing_technology,N50,C \
    > potential-target.tsv

cat ../ASSEMBLY/collect.pass.tsv |
    tsv-filter -H --or \
        --str-eq Assembly_level:"Complete Genome" \
        --str-eq Assembly_level:"Chromosome" \
        --le C:50 |
    tsv-select -H -f name,Assembly_level,Genome_coverage,Sequencing_technology,N50,C \
    > complete-genome.tsv

echo -e "#Serial\tGroup\tTarget\tCount" > group_target.tsv

# groups according `groups.tsv`
ARRAY=(
    'C_E_H::E_web_GCA_001278495_1' # 1
    'T_afr_har::T_har_CGMCC_20739_GCA_019097725_1' # 3
    'T_asperello_asperellum::T_asperellum_FT101_GCA_020647865_1' # 5
    'T_atrov_koningio::T_atrov_IMI_206040_GCF_000171015_1' # 6
    'T_cit_lon_ree::T_ree_QM6a_GCF_000167675_1' # 7
    'T_vire::T_vire_Gv29_8_GCF_000170995_1' # 8
)

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"

    SERIAL=$(
        cat ../MinHash/groups.tsv |
            tsv-filter --str-eq 2:${TARGET_NAME} |
            tsv-select -f 1
    )

    cat ../MinHash/groups.tsv |
        tsv-filter --str-eq 1:${SERIAL} |
        tsv-select -f 2 |
        tsv-join -f ../ASSEMBLY/url.tsv -k 1 -a 3 \
        > ${GROUP_NAME}

    COUNT=$(cat ${GROUP_NAME} | wc -l )

    echo -e "${SERIAL}\t${GROUP_NAME}\t${TARGET_NAME}\t${COUNT}" >> group_target.tsv

done

# Custom groups
ARRAY=(
    'Trichoderma::T_ree_QM6a_GCF_000167675_1'
    'Trichoderma_reesei::T_ree_QM6a_GCF_000167675_1'
    'Trichoderma_asperellum::T_asperellum_FT101_GCA_020647865_1'
    'Trichoderma_harzianum::T_har_CGMCC_20739_GCA_019097725_1'
    'Trichoderma_atroviride::T_atrov_IMI_206040_GCF_000171015_1'
    'Trichoderma_virens::T_vire_Gv29_8_GCF_000170995_1'
)

SERIAL=100
for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"

    SERIAL=$((SERIAL + 1))
    GROUP_NAME_2=$(echo $GROUP_NAME | tr "_" " ")

    if [ "$GROUP_NAME" = "Trichoderma" ]; then
        cat ../ASSEMBLY/collect.pass.tsv |
            tsv-filter -H --not-blank RefSeq_category |
            sed '1d' |
            tsv-select -f 1 \
            > T.tmp
        echo "C_pro_CCMJ2080_GCA_004303015_1" >> T.tmp
        echo "E_web_EWB_GCA_003055145_1" >> T.tmp
        echo "E_web_GCA_001278495_1" >> T.tmp
        echo "H_perniciosus_HP10_GCA_008477525_1" >> T.tmp
        echo "H_ros_CCMJ2808_GCA_011799845_1" >> T.tmp
        cat T.tmp |
            tsv-uniq |
            tsv-join -f ../ASSEMBLY/url.tsv -k 1 -a 3 \
            > ${GROUP_NAME}

    else
        cat ../ASSEMBLY/collect.pass.tsv |
            tsv-select -f 1,2 |
            grep "${GROUP_NAME_2}" |
            tsv-select -f 1 |
            tsv-join -f ../ASSEMBLY/url.tsv -k 1 -a 3 \
            > ${GROUP_NAME}
    fi

    COUNT=$(cat ${GROUP_NAME} | wc -l )

    echo -e "${SERIAL}\t${GROUP_NAME}\t${TARGET_NAME}\t${COUNT}" >> group_target.tsv

done

cat group_target.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

```

| #Serial | Group                  | Target                             | Count |
|---------|------------------------|------------------------------------|------:|
| 1       | C_E_H                  | E_web_GCA_001278495_1              |     5 |
| 3       | T_afr_har              | T_har_CGMCC_20739_GCA_019097725_1  |    20 |
| 5       | T_asperello_asperellum | T_asperellum_FT101_GCA_020647865_1 |    16 |
| 6       | T_atrov_koningio       | T_atrov_IMI_206040_GCF_000171015_1 |    15 |
| 7       | T_cit_lon_ree          | T_ree_QM6a_GCF_000167675_1         |    25 |
| 8       | T_vire                 | T_vire_Gv29_8_GCF_000170995_1      |     9 |
| 101     | Trichoderma            | T_ree_QM6a_GCF_000167675_1         |    32 |
| 102     | Trichoderma_reesei     | T_ree_QM6a_GCF_000167675_1         |    13 |
| 103     | Trichoderma_asperellum | T_asperellum_FT101_GCA_020647865_1 |    13 |
| 104     | Trichoderma_harzianum  | T_har_CGMCC_20739_GCA_019097725_1  |    10 |
| 105     | Trichoderma_atroviride | T_atrov_IMI_206040_GCF_000171015_1 |     7 |
| 106     | Trichoderma_virens     | T_vire_Gv29_8_GCF_000170995_1      |     8 |

## Prepare sequences for `egaz`

* `--perseq` for Chromosome-level assemblies and targets
    * means split fasta by names, targets or good assembles should set it

```shell
cd ~/data/Trichoderma

# /share/home/wangq/homebrew/Cellar/repeatmasker@4.1.1/4.1.1/libexec/famdb.py \
#   -i /share/home/wangq/homebrew/Cellar/repeatmasker@4.1.1/4.1.1/libexec/Libraries/RepeatMaskerLib.h5 \
#   lineage Fungi

# prep
egaz template \
    ASSEMBLY \
    --prep -o Genome \
    $( cat taxon/group_target.tsv |
        sed -e '1d' | cut -f 3 |
        parallel -j 1 echo " --perseq {} "
    ) \
    $( cat taxon/complete-genome.tsv |
        sed '1d' | cut -f 1 |
        parallel -j 1 echo " --perseq {} "
    ) \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--parallel 16"

bash Genome/0_prep.sh

# gff
for n in \
    $(cat taxon/group_target.tsv | sed -e '1d' | cut -f 3 ) \
    $( cat taxon/potential-target.tsv | sed -e '1d' | cut -f 1 ) \
    ; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"

    gzip -dc ${FILE_GFF} > Genome/${n}/chr.gff
done

```

## Generate alignments

```shell
cd ~/data/Trichoderma

cat taxon/group_target.tsv |
    sed -e '1d' |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{3}]\n"

        egaz template \
            Genome/{3} \
            $(cat taxon/{2} | grep -v -x "{3}" | xargs -I[] echo "Genome/[]") \
            --multi -o groups/{2}/ \
            --tree MinHash/tree.nwk \
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
cd ~/data/Trichoderma/Alignment

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
