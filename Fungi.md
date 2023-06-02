# Fungi

Download all genomes and analyze representative strains.

<!-- toc -->

- [Taxon info](#taxon-info)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
    * [Model organisms](#model-organisms)
- [Download all assemblies](#download-all-assemblies)
    * [Create assembly.tsv](#create-assemblytsv)
    * [Rsync and check](#rsync-and-check)
    * [Check N50 of assemblies](#check-n50-of-assemblies)
    * [Rsync to hpcc](#rsync-to-hpcc)
- [BioSample](#biosample)
- [Count species and strains](#count-species-and-strains)
    * [Order](#order)
    * [Genus](#genus)
    * [ReRoot](#reroot)
- [MinHash](#minhash)
    * [Compute MinHash](#compute-minhash)
    * [Raw phylo-tree of representative assemblies](#raw-phylo-tree-of-representative-assemblies)
    * [Tweak the mash tree](#tweak-the-mash-tree)
- [Non-redundant strains within species](#non-redundant-strains-within-species)
    * [Genus](#genus-1)
- [Groups and targets](#groups-and-targets)
- [Collect proteins](#collect-proteins)
    * [`all.pro.fa`](#allprofa)
    * [`all.replace.fa`](#allreplacefa)
    * [`all.info.tsv`](#allinfotsv)
- [Phylogenetics with fungi61](#phylogenetics-with-fungi61)
    * [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
    * [Tweak the concat tree](#tweak-the-concat-tree)
- [InterProScan on all proteins of representative and typical strains](#interproscan-on-all-proteins-of-representative-and-typical-strains)
- [Divergence of Fungi](#divergence-of-fungi)

<!-- tocstop -->

## Taxon info

* [Fungi](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4751)

### List all ranks

```shell
nwr member Fungi |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat |
    sed 's/-\s*|$/-:|/'

```

| rank          | count |
|---------------|------:|
| kingdom       |     1 |
| no rank       |  4623 |
| species       | 65297 |
| subkingdom    |     1 |
| class         |    65 |
| order         |   242 |
| family        |   904 |
| genus         |  7466 |
| phylum        |    10 |
| subphylum     |    14 |
| strain        |  2236 |
| varietas      |  1060 |
| subspecies    |   162 |
| forma         |   210 |
| isolate       |    29 |
| subclass      |    19 |
| clade         |    26 |
| suborder      |    24 |
| subfamily     |    17 |
| subgenus      |    10 |
| section       |    37 |
| species group |     1 |
| tribe         |     3 |
| superfamily   |     1 |
| morph         |     2 |

### Species with assemblies

In the vast majority of fungal species, only one genome was selected for refseq.

* 'RefSeq'
    * RS1 - assembly_level: 'Complete Genome', 'Chromosome'
    * RS2 - genome_rep: 'Full'
* 'Genbank'
    * '>= 20 genomes'
        * GB1 - With strain ID; assembly_level: 'Complete Genome', 'Chromosome'
        * GB2 - assembly_level: 'Complete Genome', 'Chromosome'
        * GB3 - genome_rep: 'Full'
    * '>= 1 genomes'
        * GB4 - assembly_level: 'Complete Genome', 'Chromosome'
        * GB5 - genome_rep: 'Full'

```shell
mkdir -p ~/data/Fungi/summary
cd ~/data/Fungi/summary

# should have a valid name of genus
nwr member Fungi -r genus |
    grep -v -i "Candidatus " |
    grep -v -i "candidate " |
    grep -v " sp." |
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#7466 genus.list

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
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
            AND assembly_level IN ('Complete Genome', 'Chromosome')
            AND genome_rep IN ('Full') -- fully representative
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS1.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
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
    > RS2.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
    echo "
        SELECT
            species_id,
            species,
            COUNT(DISTINCT tax_id) AS count -- with strain ID
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND assembly_level IN ('Complete Genome', 'Chromosome')
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 20
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB1.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
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
            AND assembly_level IN ('Complete Genome', 'Chromosome')
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 20
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB2.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
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
        HAVING count >= 20
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB3.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
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
            AND assembly_level IN ('Complete Genome', 'Chromosome')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB4.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
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
    > GB5.tsv

wc -l RS*.tsv GB*.tsv
#    91 RS1.tsv
#   499 RS2.tsv
#     1 GB1.tsv
#     4 GB2.tsv
#    85 GB3.tsv
#   281 GB4.tsv
#  3395 GB5.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e ${C}${N}.tsv ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     94
#RS2     504
#GB1     109
#GB2     770
#GB3     7260
#GB4     1325
#GB5     13371

```

* The names of some genera are abnormal

```shell
cd ~/data/Fungi/summary

cat RS*.tsv GB*.tsv |
    cut -f 1,2 |
    tsv-uniq |
    grep "\[" |
    nwr append stdin -r genus |
    head
#498019  [Candida] auris Candida/Metschnikowiaceae
#1231522 [Candida] duobushaemulonis      Candida/Metschnikowiaceae
#45357   [Candida] haemuloni     Candida/Metschnikowiaceae
#418784  [Candida] pseudohaemulonii      Candida/Metschnikowiaceae
#561895  [Candida] subhashii     Spathaspora
#566037  [Ashbya] aceris (nom. inval.)   Eremothecium
#312227  [Candida] hispaniensis  Candida/Saccharomycetales
#45354   [Candida] intermedia    Candida/Metschnikowiaceae
#2093215 [Candida] vulturna (nom. inval.)        Clavispora
#45518   [Candida] aaseri        Yamadazyma

```

### Model organisms

There is only one model genome in Fungi, Saccharomyces cerevisiae S288C

```shell
cd ~/data/Fungi/summary

GENUS=$(
    cat genus.list.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > reference.tsv

cp reference.tsv ~/Scripts/genomes/assembly/Fungi.reference.tsv

```

## Download all assemblies

### Create assembly.tsv

Two levels:

* 'RefSeq'
    * RS2 - genome_rep: 'Full'
* 'Genbank'
    * GB5 - genome_rep: 'Full'

If a refseq assembly is available, the corresponding genbank one is not downloaded

```shell
cd ~/data/Fungi/summary

cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession \
    > raw.tsv

# RS2
SPECIES=$(
    cat RS2.tsv |
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
    cut -f 6 \
    > rs.acc.tsv

# GB5
SPECIES=$(
    cat GB5.tsv |
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
#13877 lines, 7 fields

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
    > Fungi.assembly.tsv

datamash check < Fungi.assembly.tsv
#13610 lines, 5 fields

# find potential duplicate strains or assemblies
cat Fungi.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Fungi.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

cat Fungi.assembly.tsv |
    tsv-filter --or --str-in-fld 1:genomosp --str-in-fld 1:genomovar

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Fungi.assembly.tsv
# cp Fungi.assembly.tsv ~/Scripts/genomes/assembly

# Comment out unneeded strains

# Cleaning
rm raw*.*sv

```

### Rsync and check

```shell
cd ~/data/Fungi

nwr template ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    --ass \
    -o .

# Run
bash ASSEMBLY/rsync.sh

# Check md5; create check.lst
bash ASSEMBLY/check.sh

# Put the misplaced directory in the right place
#bash ASSEMBLY/reorder.sh

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
#32322684        37330627        99165   1578216 167     1054

# Collect; create collect.csv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.csv summary/

cat ASSEMBLY/counts.tsv |
    mlr --itsv --omd cat |
    sed 's/-\s*|$/-:|/'

```

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Fungi/ \
    wangq@202.119.37.251:data/Fungi

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Fungi/ \
    wangq@58.213.64.36:data/Fungi

# rsync -avP wangq@202.119.37.251:data/Fungi/ ~/data/Fungi

# rsync -avP -e 'ssh -p 8804' wangq@58.213.64.36:data/Fungi/ ~/data/Fungi

```

## BioSample

```shell
cd ~/data/Fungi

mkdir -p biosample

ulimit -n `ulimit -Hn`

cat ASSEMBLY/collect.csv |
    tsv-select -H -d, -f BioSample |
    grep "^SAM" |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -s biosample/{}.txt ]; then
            >&2 echo {}
            curl -fsSL "https://www.ncbi.nlm.nih.gov/biosample/?term={}&report=full&format=text" -o biosample/{}.txt
#            curl -fsSL "https://www.ebi.ac.uk/biosamples/samples/{}" -o biosample/{}.json
        fi
    '

# Allowing samples not in the list
find biosample -name "SAM*.txt" | wc -l
# 4061

find biosample -name "SAM*.txt" |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat {} |
            perl -nl -e '\''
                print $1 if m{\s+\/([\w_ ]+)=};
            '\''
    ' |
    tsv-uniq --at-least 50 | # ignore rare attributes
    grep -v "^INSDC" |
    grep -v "^ENA" \
    > summary/attributes.lst

cat summary/attributes.lst |
    (echo -e "BioSample" && cat) |
    tr '\n' '\t' |
    sed 's/\t$/\n/' \
    > summary/biosample.tsv

find biosample -name "SAM*.txt" |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        >&2 echo {/.}
        cat {} |
            perl -nl -MPath::Tiny -e '\''
                BEGIN {
                    our @keys = grep {/\S/} path(q{summary/attributes.lst})->lines({chomp => 1});
                    our %stat = ();
                }

                m(\s+\/([\w_ ]+)=\"(.+)\") or next;
                my $k = $1;
                my $v = $2;
                if ( $v =~ m(\bNA|missing|Not applicable|not collected|not available|not provided|N\/A|not known|unknown\b)i ) {
                    $stat{$k} = q();
                } else {
                    $stat{$k} = $v;
                }

                END {
                    my @c;
                    for my $key ( @keys ) {
                        if (exists $stat{$key}) {
                            push @c, $stat{$key};
                        }
                        else {
                            push @c, q();
                        }
                    }
                    print join(qq{\t}, q{{/.}}, @c);
                }
            '\''
    ' \
    >> summary/biosample.tsv

```

## Count species and strains

```shell
cd ~/data/Fungi

echo "
    SELECT
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '%symbiont %'
        AND refseq_category IN ('reference genome', 'representative genome')
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > summary/assembly_accession.lst

wc -l summary/assembly_accession.lst
#5278 assembly_accession.lst

cat ASSEMBLY/collect.csv |
    tsv-select -H -d, -f Taxid |
    sed '1d' |
    tsv-uniq |
    nwr append stdin -r species |
    tsv-select -f 2 |
    tsv-uniq |
    wc -l
#659

cat summary/collect.pass.csv |
    tsv-select -H -d, -f Taxid |
    sed '1d' |
    tsv-uniq |
    nwr append stdin -r species |
    tsv-select -f 2 |
    tsv-uniq |
    wc -l
#654

cat summary/collect.pass.csv |
    tsv-filter -H -d, --or \
        --istr-eq "RefSeq_category:reference genome" --istr-eq "RefSeq_category:representative genome" |
    grep -v "symbiont " |
    tsv-select -H -d, -f name |
    sed '1d' \
    > summary/representative.lst

wc -l summary/representative.lst
#639 summary/representative.lst

cat summary/collect.pass.csv |
    sed -e '1d' |
    tr "," "\t" |
    tsv-select -f 1,3 |
    nwr append stdin -c 2 -r species -r genus -r family -r order \
    > summary/strains.taxon.tsv

cat ~/Scripts/genomes/assembly/Fungi.assembly.tsv |
    tsv-summarize -H -g 3 --count |
    tsv-filter -H --ge 2:50 |
    mlr --itsv --omd cat

cat summary/strains.taxon.tsv |
    tsv-summarize -g 3 --count |
    ( echo -e 'organism\tcount' && cat ) |
    tsv-filter -H --ge 2:50 |
    mlr --itsv --omd cat

```

| organism                 | count |
|--------------------------|-------|
| Aspergillus flavus       | 144   |
| Aspergillus fumigatus    | 80    |
| Aspergillus niger        | 96    |
| Aspergillus oryzae       | 90    |
| Botryosphaeria dothidea  | 131   |
| Candida albicans         | 64    |
| Cryphonectria parasitica | 92    |
| Fusarium graminearum     | 114   |
| Komagataella phaffii     | 129   |
| Ophidiomyces ophidiicola | 63    |
| Parastagonospora nodorum | 171   |
| Penicillium chrysogenum  | 78    |
| Pyricularia oryzae       | 251   |
| Rhodotorula mucilaginosa | 66    |
| Saccharomyces cerevisiae | 111   |
| Saitozyma podzolica      | 56    |
| Venturia inaequalis      | 85    |

| organism                 | count |
|--------------------------|-------|
| Aspergillus flavus       | 136   |
| Aspergillus fumigatus    | 74    |
| Aspergillus niger        | 94    |
| Aspergillus oryzae       | 90    |
| Botryosphaeria dothidea  | 128   |
| Candida albicans         | 53    |
| Cryphonectria parasitica | 68    |
| Fusarium graminearum     | 112   |
| Komagataella phaffii     | 127   |
| Ophidiomyces ophidiicola | 63    |
| Parastagonospora nodorum | 163   |
| Penicillium chrysogenum  | 77    |
| Pyricularia oryzae       | 170   |
| Rhodotorula mucilaginosa | 65    |
| Saccharomyces cerevisiae | 111   |

### Order

```shell
cd ~/data/Fungi

# Group by order
cat summary/collect.pass.csv |
    sed -e '1d' |
    tsv-select -d, -f 3 |
    tsv-uniq |
    nwr append stdin -r order |
    tsv-select -f 2 |
    tsv-uniq |
    grep -v "NA" |
    sort \
    > summary/order.lst

cat summary/order.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(cat summary/collect.pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r order -r species |
            grep {} |
            tsv-select -f 1,3 |
            tsv-uniq |
            wc -l)

        n_strains=$(cat summary/collect.pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r order |
            grep {} |
            wc -l)

        printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' |
    nwr append stdin --id |
    tsv-select -f 5,4,2,3 |
    tsv-sort -k2,2 |
    tsv-filter --ge 4:50 |
    (echo -e '#tax_id\torder\t#species\t#strains' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | order             | #species | #strains |
|---------|-------------------|----------|----------|
| 451869  | Botryosphaeriales | 5        | 140      |
| 34395   | Chaetothyriales   | 23       | 50       |
| 5114    | Diaporthales      | 9        | 76       |
| 5042    | Eurotiales        | 107      | 719      |
| 5125    | Hypocreales       | 137      | 440      |
| 639021  | Magnaporthales    | 10       | 176      |
| 33183   | Onygenales        | 37       | 116      |
| 92860   | Pleosporales      | 51       | 288      |
| 4892    | Saccharomycetales | 279      | 803      |
| 231213  | Sporidiobolales   | 3        | 67       |
| 5234    | Tremellales       | 36       | 112      |
| 5267    | Ustilaginales     | 13       | 52       |

### Genus

```shell
cd ~/data/Fungi

# Group by genus
cat summary/collect.pass.csv |
    sed -e '1d' |
    tsv-select -d, -f 3 |
    tsv-uniq |
    nwr append stdin -r genus |
    tsv-select -f 2 |
    tsv-uniq |
    grep -v "NA" |
    sort \
    > summary/genus.lst

cat summary/genus.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(cat summary/collect.pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r genus -r species |
            grep {} |
            tsv-select -f 1,3 |
            tsv-uniq |
            wc -l)

        n_strains=$(cat summary/collect.pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r genus |
            grep {} |
            wc -l)

        printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' |
    nwr append stdin --id |
    tsv-select -f 5,4,2,3 |
    tsv-sort -k2,2 \
    > summary/genus.count.tsv

cat summary/genus.count.tsv |
    tsv-filter --ge 4:50 |
    (echo -e '#tax_id\tgenus\t#species\t#strains' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | genus                     | #species | #strains |
|---------|---------------------------|----------|----------|
| 5598    | Alternaria                | 20       | 69       |
| 5052    | Aspergillus               | 70       | 468      |
| 45132   | Botryosphaeria            | 1        | 128      |
| 5475    | Candida                   | 55       | 150      |
| 2964429 | Candida/Metschnikowiaceae | 6        | 55       |
| 5115    | Cryphonectria             | 2        | 68       |
| 5206    | Cryptococcus              | 27       | 55       |
| 5506    | Fusarium                  | 44       | 251      |
| 460517  | Komagataella              | 4        | 128      |
| 27320   | Metschnikowia             | 10       | 61       |
| 1387562 | Ophidiomyces              | 1        | 63       |
| 1351751 | Parastagonospora          | 2        | 163      |
| 5073    | Penicillium               | 23       | 214      |
| 48558   | Pyricularia               | 9        | 175      |
| 5533    | Rhodotorula               | 3        | 67       |
| 4930    | Saccharomyces             | 128      | 214      |
| 5543    | Trichoderma               | 40       | 97       |

### ReRoot

小孢子虫 (Microsporidia) 可以当作真菌的基部类群, 也可以当作真菌的姊妹类群

```shell
cd ~/data/Fungi

cat summary/strains.taxon.tsv |
    nwr append stdin -c 2 -r phylum -r subkingdom |
    tsv-filter --str-ne "8:Dikarya" | # 双核亚界
    tsv-select -f 1,3,7 |
    tsv-sort -k3,3 -k1,1 |
    tsv-filter --or \
        --str-in-fld "2:Mitosporidium" \
        --str-in-fld "2:Nematocida" \
        --str-in-fld "2:Nosema" |
    grep -v -Fw -f ASSEMBLY/omit.lst
#Mit_dap_GCF_000760515_2 Mitosporidium daphniae  Microsporidia
#Nem_aus_GCF_000738915_1 Nematocida ausubeli     Microsporidia
#Nem_homos_GCF_024244095_1       Nematocida homosporus   Microsporidia
#Nem_maj_GCF_021653875_1 Nematocida major        Microsporidia
#Nem_minor_GCF_024244105_1       Nematocida minor        Microsporidia
#Nem_pari_ERTm1_GCF_000250985_1  Nematocida parisii      Microsporidia
#No_cera_GCF_000988165_1 Nosema ceranae  Microsporidia

```

## MinHash

### Compute MinHash

```shell
mkdir -p ~/data/Fungi/mash
cd ~/data/Fungi/mash

# Remove .msh files not in the list
find . -maxdepth 1 -mindepth 1 -type f -name "*.msh" |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        basename {} .msh
    ' |
    tsv-join --exclude -k 1 -f ../ASSEMBLY/url.tsv -d 1 |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        echo Remove {}
        rm {}.msh
    '

# Compute MinHash
cat ../ASSEMBLY/url.tsv |
    cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        if [[ -e {}.msh ]]; then
            exit
        fi

        2>&1 echo "==> {}"

        find ../ASSEMBLY/{} -name "*_genomic.fna.gz" |
            grep -v "_from_" |
            xargs cat |
            mash sketch -k 21 -s 100000 -p 2 - -I "{}" -o {}
    '

```

### Raw phylo-tree of representative assemblies

```shell
mkdir -p ~/data/Fungi/tree
cd ~/data/Fungi/tree

mash triangle -E -p 8 -l <(
    cat ../summary/representative.lst |
        parallel --no-run-if-empty --linebuffer -k -j 1 '
            if [[ -e ../mash/{}.msh ]]; then
                echo "../mash/{}.msh"
            fi
        '
    ) \
    > mash.dist.tsv

# Fill matrix with lower triangle
tsv-select -f 1-3 mash.dist.tsv |
    (tsv-select -f 2,1,3 mash.dist.tsv && cat) |
    (
        cut -f 1 mash.dist.tsv |
            tsv-uniq |
            parallel -j 1 --keep-order 'echo -e "{}\t{}\t0"' &&
        cat
    ) \
    > mash.dist_full.tsv

cat mash.dist_full.tsv |
    Rscript -e '
        library(readr);
        library(tidyr);
        library(ape);
        pair_dist <- read_tsv(file("stdin"), col_names=F);
        tmp <- pair_dist %>%
            pivot_wider( names_from = X2, values_from = X3, values_fill = list(X3 = 1.0) )
        tmp <- as.matrix(tmp)
        mat <- tmp[,-1]
        rownames(mat) <- tmp[,1]

        dist_mat <- as.dist(mat)
        clusters <- hclust(dist_mat, method = "ward.D2")
        tree <- as.phylo(clusters)
        write.tree(phy=tree, file="tree.nwk")

        group <- cutree(clusters, h=0.5) # k=3
        groups <- as.data.frame(group)
        groups$ids <- rownames(groups)
        rownames(groups) <- NULL
        groups <- groups[order(groups$group), ]
        write_tsv(groups, "groups.tsv")
    '

```

### Tweak the mash tree

```shell
cd ~/data/Fungi/tree

nw_reroot tree.nwk Mit_dap_GCF_000760515_2 No_cera_GCF_000988165_1 |
    nw_order -c n - \
    > mash.reroot.newick

# rank::col
ARRAY=(
    'order::6'
    'family::5'
#    'genus::4'
#    'species::3'
)

rm mash.condensed.map
CUR_TREE=mash.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/genomes/bin/condense_tree.sh ${CUR_TREE} ../summary/strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick mash.${GROUP_NAME}.newick
    cat condense.map >> mash.condensed.map

    CUR_TREE=mash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 mash.family.newick |
    rsvg-convert -o Fungi.mash.png

```

## Non-redundant strains within species

If the ANI value between two strains was less than 0.005, the two strains were considered as
redundant.

```shell
cd ~/data/Fungi

mkdir NR

ANI_VALUE=0.005

cat summary/strains.taxon.tsv |
    tsv-summarize -H -g 3 --count |
    tsv-filter -H --ge 2:2 \
    > NR/species.tsv

tsv-summarize NR/species.tsv -H --count --sum count
#count   count_sum
#197     3043

# each species
cat NR/species.tsv | sed '1d' | tsv-select -f 1 | #head -n 10 |
while read SPECIES; do
    1>&2 echo "==> ${SPECIES}"

    SPECIES_=$(
        echo "${SPECIES}" |
            tr " " "_"
    )

    mkdir -p NR/${SPECIES_}

    cat summary/strains.taxon.tsv |
        tsv-filter --str-eq "3:${SPECIES}" |
        tsv-select -f 1 \
        > "NR/${SPECIES_}/assembly.lst"

    1>&2 echo "    mash distances"

    if [[ ! -f "NR/${SPECIES_}/mash.dist.tsv" ]]; then
        mash triangle -E -p 8 -l <(
            cat "NR/${SPECIES_}/assembly.lst" |
                parallel --no-run-if-empty --linebuffer -k -j 1 '
                    if [[ -e mash/{}.msh ]]; then
                        echo "mash/{}.msh"
                    fi
                '
            ) \
            > "NR/${SPECIES_}/mash.dist.tsv"
    fi

    1>&2 echo "    List NR"

    cat "NR/${SPECIES_}/mash.dist.tsv" |
        tsv-filter --ff-str-ne 1:2 --le "3:${ANI_VALUE}" \
        > "NR/${SPECIES_}/redundant.dist.tsv"

    cat "NR/${SPECIES_}/redundant.dist.tsv" |
        perl -nla -F"\t" -MGraph::Undirected -e '
            BEGIN {
                our $g = Graph::Undirected->new;
            }

            $g->add_edge($F[0], $F[1]);

            END {
                for my $cc ( $g->connected_components ) {
                    print join qq{\t}, sort @{$cc};
                }
            }
        ' \
        > "NR/${SPECIES_}/connected_components.tsv"

    cat "NR/${SPECIES_}/connected_components.tsv" |
        perl -nla -MPath::Tiny -F"\t" -e '
            BEGIN {
                our %rep = map { ($_, 1) } path(q{summary/representative.lst})->lines({chomp => 1});
            }

            # Representative strains are preferred
            if ( grep { $rep{$_} } @F ) {
                @F = grep { ! $rep{$_} } @F
            }
            else {
                shift @F;
            }
            printf qq{%s\n}, $_ for @F;
            ' \
        > "NR/${SPECIES_}/redundant.lst"

    cat "NR/${SPECIES_}/assembly.lst" |
        tsv-join --exclude -f "NR/${SPECIES_}/redundant.lst" \
        > "NR/${SPECIES_}/NR.lst"

done

find NR -name "redundant.lst" -size +0 | wc -l
#152

find NR -name "redundant.lst" -empty | wc -l
#45

find NR -name "NR.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/NR.lst

wc -l summary/NR.lst
#728 summary/NR.lst

```

### Genus

```shell
cd ~/data/Fungi

cat summary/genus.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(
            cat summary/collect.pass.csv |
                sed "1d" |
                grep -F -w -f <( cat summary/NR.lst summary/representative.lst | sort | uniq ) |
                tsv-select -d, -f 3 |
                nwr append stdin -r genus -r species |
                grep {} |
                tsv-select -f 1,3 |
                tsv-uniq |
                wc -l
        )

        n_strains=$(
            cat summary/collect.pass.csv |
                sed "1d" |
                grep -F -w -f <( cat summary/NR.lst summary/representative.lst | sort | uniq ) |
                tsv-select -d, -f 3 |
                nwr append stdin -r genus |
                grep {} |
                wc -l
        )

        printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' |
    nwr append stdin --id |
    tsv-select -f 5,4,2,3 |
    tsv-sort -k2,2 |
    tsv-filter --ge 4:10 |
    (echo -e '#tax_id\tgenus\t#species\t#strains' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | genus                     | #species | #strains |
|---------|---------------------------|----------|----------|
| 5598    | Alternaria                | 20       | 44       |
| 5052    | Aspergillus               | 60       | 103      |
| 5579    | Aureobasidium             | 6        | 16       |
| 5581    | Beauveria                 | 2        | 10       |
| 45132   | Botryosphaeria            | 1        | 26       |
| 5475    | Candida                   | 30       | 41       |
| 2964429 | Candida/Metschnikowiaceae | 5        | 11       |
| 5455    | Colletotrichum            | 15       | 15       |
| 5206    | Cryptococcus              | 11       | 14       |
| 5112    | Epichloe                  | 14       | 16       |
| 5506    | Fusarium                  | 40       | 87       |
| 300275  | Lachancea                 | 9        | 10       |
| 27320   | Metschnikowia             | 9        | 15       |
| 1351751 | Parastagonospora          | 2        | 11       |
| 5073    | Penicillium               | 15       | 28       |
| 5296    | Puccinia                  | 7        | 13       |
| 48558   | Pyricularia               | 4        | 14       |
| 1322061 | Rhizoctonia               | 6        | 16       |
| 4930    | Saccharomyces             | 25       | 41       |
| 1890244 | Saitozyma                 | 1        | 15       |
| 5094    | Talaromyces               | 8        | 12       |
| 5543    | Trichoderma               | 36       | 61       |
| 1047167 | Zymoseptoria              | 5        | 12       |

## Groups and targets

Review `summary/collect.pass.csv` and `tree/groups.tsv`

## Collect proteins

### `all.pro.fa`

```shell
cd ~/data/Fungi

mkdir -p PROTEINS

cat summary/representative.lst | wc -l
#609

cat summary/representative.lst | grep -v -Fx -f ASSEMBLY/omit.lst | wc -l
#529

for STRAIN in $(cat summary/representative.lst | grep -v -Fx -f ASSEMBLY/omit.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz
done |
    pigz -p4 \
    > PROTEINS/all.pro.fa.gz

gzip -dcf PROTEINS/all.pro.fa.gz |
    perl -nl -e '
        BEGIN { our %seen; our $h; }

        if (/^>/) {
            $h = (split(" ", $_))[0];
            $seen{$h}++;
            $_ = $h;
        }
        print if $seen{$h} == 1;
    ' |
    pigz -p4 \
    > PROTEINS/all.uniq.fa.gz

# counting proteins
gzip -dcf PROTEINS/all.pro.fa.gz |
    grep "^>" |
    wc -l |
    numfmt --to=si
#5.4M

gzip -dcf PROTEINS/all.pro.fa.gz |
    grep "^>" |
    tsv-uniq |
    wc -l |
    numfmt --to=si
#5.4M

# annotations may be different
gzip -dcf PROTEINS/all.uniq.fa.gz |
    grep "^>" |
    wc -l |
    numfmt --to=si
#5.4M

```

### `all.replace.fa`

```shell
cd ~/data/Fungi

rm PROTEINS/all.strain.tsv PROTEINS/all.replace.fa.gz

for STRAIN in $(cat summary/representative.lst | grep -v -Fx -f ASSEMBLY/omit.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
        grep "^>" |
        cut -d" " -f 1 |
        sed "s/^>//" |
        STRAIN=${STRAIN} perl -nl -e '
            $n = $_;
            $s = $n;
            $s =~ s/\.\d+//;
            printf qq{%s\t%s_%s\t%s\n}, $n, $ENV{STRAIN}, $s, $ENV{STRAIN};
        ' \
    > PROTEINS/${STRAIN}.replace.tsv

    cut -f 2,3 PROTEINS/${STRAIN}.replace.tsv >> PROTEINS/all.strain.tsv

    faops replace -s \
        ASSEMBLY/${STRAIN}/*_protein.faa.gz \
        <(cut -f 1,2 PROTEINS/${STRAIN}.replace.tsv) \
        stdout |
        pigz -p4 \
        >> PROTEINS/all.replace.fa.gz

    rm PROTEINS/${STRAIN}.replace.tsv
done

gzip -dcf PROTEINS/all.replace.fa.gz |
    grep "^>" |
    wc -l |
    numfmt --to=si
#5.4M

(echo -e "#name\tstrain" && cat PROTEINS/all.strain.tsv)  \
    > temp &&
    mv temp PROTEINS/all.strain.tsv

faops size PROTEINS/all.replace.fa.gz > PROTEINS/all.replace.sizes

(echo -e "#name\tsize" && cat PROTEINS/all.replace.sizes) > PROTEINS/all.size.tsv

rm PROTEINS/all.replace.sizes

```

### `all.info.tsv`

```shell
cd ~/data/Fungi

for STRAIN in $(cat summary/representative.lst | grep -v -Fx -f ASSEMBLY/omit.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
        grep "^>" |
        sed "s/^>//" |
        perl -nl -e '/\[.+\[/ and s/\[/\(/; print' |
        perl -nl -e '/\].+\]/ and s/\]/\)/; print' |
        perl -nl -e 's/\s+\[.+?\]$//g; print' |
        perl -nl -e 's/MULTISPECIES: //g; print' |
        STRAIN=${STRAIN} perl -nl -e '
            /^(\w+)\.\d+\s+(.+)$/ or next;
            printf qq{%s_%s\t%s\n}, $ENV{STRAIN}, $1, $2;
        '
done \
    > PROTEINS/all.annotation.tsv

cat PROTEINS/all.annotation.tsv |
    wc -l |
    numfmt --to=si
#5.4M

(echo -e "#name\tannotation" && cat PROTEINS/all.annotation.tsv) \
    > temp &&
    mv temp PROTEINS/all.annotation.tsv

# check differences
cat PROTEINS/all.size.tsv |
    grep -F -f <(cut -f 1 PROTEINS/all.annotation.tsv) -v

tsv-join \
    PROTEINS/all.strain.tsv \
    --data-fields 1 \
    -f PROTEINS/all.size.tsv \
    --key-fields 1 \
    --append-fields 2 \
    > PROTEINS/all.strain_size.tsv

tsv-join \
    PROTEINS/all.strain_size.tsv \
    --data-fields 1 \
    -f PROTEINS/all.annotation.tsv \
    --key-fields 1 \
    --append-fields 2 \
    > PROTEINS/all.info.tsv

cat PROTEINS/all.info.tsv |
    wc -l |
    numfmt --to=si
#5.4M

```

## Phylogenetics with fungi61

### Find corresponding proteins by `hmmsearch`

* Download HMM models as described in [`HMM.md`](HMM.md)

* The `E_VALUE` was manually adjusted to 1e-20 to reach a balance between sensitivity and
  speciality.

```shell
E_VALUE=1e-20

cd ~/data/Fungi

# Find all genes
for marker in $(cat ~/data/HMM/fungi61/fungi61.lst); do
    >&2 echo "==> marker [${marker}]"

    mkdir -p PROTEINS/${marker}

    cat summary/representative.lst |
        grep -v -Fx -f ASSEMBLY/omit.lst |
        parallel --no-run-if-empty --linebuffer -k -j 8 "
            gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/HMM/fungi61/HMM/${marker}.HMM - |
                grep '>>' |
                perl -nl -e ' m{>>\s+(\S+)} and printf qq{%s\t%s\n}, \$1, {}; '
        " \
        > PROTEINS/${marker}/replace.tsv

    >&2 echo

done

```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Fungi

cat ~/data/HMM/fungi61/fungi61.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat PROTEINS/{}/replace.tsv |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
#523     543     1341

cat ~/data/HMM/fungi61/fungi61.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        cat PROTEINS/{}/replace.tsv |
            wc -l
    ' |
    paste - - |
    tsv-filter --invert --ge 2:400 --le 2:800 |
    cut -f 1 \
    > PROTEINS/fungi61.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process
cat ~/data/HMM/fungi61/fungi61.lst |
    grep -v -Fx -f PROTEINS/fungi61.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        cat PROTEINS/{}/replace.tsv \
            > PROTEINS/{}/{}.replace.tsv

        faops some PROTEINS/all.uniq.fa.gz <(
            cat PROTEINS/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > PROTEINS/{}/{}.pro.fa
    '

# Align each markers with muscle
cat ~/data/HMM/fungi61/fungi61.lst |
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

for marker in $(cat ~/data/HMM/fungi61/fungi61.lst); do
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
for marker in $(cat ~/data/HMM/fungi61/fungi61.lst); do
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
    > PROTEINS/fungi61.aln.fas

cat summary/representative.lst |
    grep -v -Fx -f ASSEMBLY/omit.lst |
    fasops concat PROTEINS/fungi61.aln.fas stdin -o PROTEINS/fungi61.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in PROTEINS/fungi61.aln.fa -out PROTEINS/fungi61.trim.fa -automated1

faops size PROTEINS/fungi61.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#111589
#14051

# To make it faster
FastTree -fastest -noml PROTEINS/fungi61.trim.fa > PROTEINS/fungi61.trim.newick

```

### Tweak the concat tree

```shell
cd ~/data/Fungi/tree

nw_reroot ../PROTEINS/fungi61.trim.newick Mit_dap_GCF_000760515_2 No_cera_GCF_000988165_1 |
    nw_order -c n - \
    > fungi61.reroot.newick

# rank::col
ARRAY=(
    'order::6'
    'family::5'
#    'genus::4'
#    'species::3'
)

rm fungi61.condensed.map
CUR_TREE=fungi61.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/genomes/bin/condense_tree.sh ${CUR_TREE} ../summary/strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick fungi61.${GROUP_NAME}.newick
    cat condense.map >> fungi61.condensed.map

    CUR_TREE=fungi61.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 fungi61.family.newick |
    rsvg-convert -o Fungi.fungi61.png

```

## InterProScan on all proteins of representative and typical strains

To be filled by other projects

```shell
mkdir -p ~/data/Fungi/STRAINS

```

## Divergence of Fungi

Ref.:

1. A genome-scale phylogeny of the kingdom Fungi. Cur Biology,
    2021. https://doi.org/10.1016/j.cub.2021.01.074

![1-s2.0-S0960982221001391-fx1_lrg.jpg](images%2F1-s2.0-S0960982221001391-fx1_lrg.jpg)
