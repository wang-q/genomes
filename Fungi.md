# Fungi

Download all genomes and analyze representative strains.

<!-- toc -->

## Taxon info

* [Fungi](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4751)

### List all ranks

```shell
nwr member Fungi |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat

```

| rank          | count |
|---------------|------:|
| kingdom       |     1 |
| no rank       |  4577 |
| species       | 64236 |
| subkingdom    |     1 |
| class         |    64 |
| order         |   240 |
| family        |   901 |
| genus         |  7410 |
| phylum        |    10 |
| subphylum     |    14 |
| strain        |  2237 |
| varietas      |  1053 |
| subspecies    |   161 |
| forma         |   207 |
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
        * GB3 - NOT 'contig'; genome_rep: 'Full'
    * '>= 2 genomes'
        * GB4 - assembly_level: 'Complete Genome', 'Chromosome'

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
#7410 genus.list

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
            AND genome_rep IN ('Full') -- fully representative
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
            AND genome_rep IN ('Full') -- fully representative
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
            AND genome_rep IN ('Full') -- fully representative
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
            AND assembly_level IN ('Complete Genome', 'Chromosome', 'Scaffold')
            AND genome_rep IN ('Full') -- fully representative
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
        HAVING count >= 2
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB4.tsv

wc -l RS*.tsv GB*.tsv
#   81 RS1.tsv
#  487 RS2.tsv
#    1 GB1.tsv
#    4 GB2.tsv
#   51 GB3.tsv
#   82 GB4.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e ${C}${N}.tsv ]; then
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#83
#492
#109
#750
#4473
#1082

```

* The names of some genera are abnormal

```shell
cd ~/data/Fungi/summary

cat RS*.tsv GB*.tsv |
    cut -f 1,2 |
    tsv-uniq |
    grep "\[" |
    nwr append stdin -r genus
#498019	[Candida] auris	Candida/Metschnikowiaceae
#1231522	[Candida] duobushaemulonis	Candida/Metschnikowiaceae
#45357	[Candida] haemuloni	Candida/Metschnikowiaceae
#418784	[Candida] pseudohaemulonii	Candida/Metschnikowiaceae
#561895	[Candida] subhashii	Spathaspora
#5477	[Candida] boidinii	Ogataea
#45354	[Candida] intermedia	Candida/Metschnikowiaceae

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

```

## Download all assemblies

### Create assembly.tsv

Five levels:

* 'RefSeq'
    * RS2 - genome_rep: 'Full'
* 'Genbank'
    * '>= 20 genomes'
        * GB1 - With strain ID; assembly_level: 'Complete Genome', 'Chromosome'
        * GB2 - assembly_level: 'Complete Genome', 'Chromosome'
        * GB3 - NOT 'contig'; genome_rep: 'Full'
    * '>= 2 genomes'
        * GB4 - assembly_level: 'Complete Genome', 'Chromosome'

If a refseq assembly is available, the corresponding genbank one is not downloaded

```shell
cd ~/data/Fungi/summary

cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,assembly_level,assembly_accession \
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
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND genome_rep IN ('Full') -- fully representative
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

# Avoid refseq in genbank
cat raw.tsv |
    cut -f 6 \
    > rs.acc.tsv

# GB1
SPECIES=$(
    cat GB1.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND tax_id NOT IN ($SPECIES) -- With strain ID
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 6 -e \
    >> raw.tsv

# GB2
SPECIES=$(
    cat GB2.tsv |
        tsv-join -f GB1.tsv -k 1 -e |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 6 -e \
    >> raw.tsv

# GB3
SPECIES=$(
    cat GB3.tsv |
        tsv-join -f GB1.tsv -k 1 -e |
        tsv-join -f GB2.tsv -k 1 -e |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND assembly_level IN ('Complete Genome', 'Chromosome', 'Scaffold')
        AND genome_rep IN ('Full') -- fully representative
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 6 -e \
    >> raw.tsv

# GB4
SPECIES=$(
    cat GB4.tsv |
        tsv-join -f GB1.tsv -k 1 -e |
        tsv-join -f GB2.tsv -k 1 -e |
        tsv-join -f GB3.tsv -k 1 -e |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND assembly_level IN ('Complete Genome', 'Chromosome', 'Scaffold')
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 6 -e \
    >> raw.tsv

datamash check < raw.tsv
#3602 lines, 6 fields

# Create abbr.
cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
    tsv-select -f 1-5 |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[5]}++; # abbr_name
        $seen{$F[5]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $F[5], $F[3], $F[1], $F[4];
        ' |
    keep-header -- sort -k3,3 -k1,1 \
    > Fungi.assembly.tsv

datamash check < Fungi.assembly.tsv
#3601 lines, 4 fields

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

### rsync and check

```shell
cd ~/data/Fungi

cat ~/Scripts/genomes/assembly/Fungi.assembly.tsv |
    tsv-filter -v --str-in-fld 2:http

nwr assembly ~/Scripts/genomes/assembly/Fungi.assembly.tsv \
    -o ASSEMBLY

# Remove dirs not in the list
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    tr "/" "\t" |
    cut -f 2 |
    tsv-join --exclude -k 1 -f ASSEMBLY/url.tsv -d 1 |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        echo Remove {}
        rm -fr ASSEMBLY/{}
    '

# Run
bash ASSEMBLY/rsync.sh

# Check md5
# rm ASSEMBLY/check.list
bash ASSEMBLY/check.sh

# Collect
bash ASSEMBLY/collect.sh

# Temporary files, possibly caused by an interrupted rsync process
find ASSEMBLY/ -type f -name ".*" > ASSEMBLY/temp.list

cat ASSEMBLY/temp.list |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        if [[ -f {} ]]; then
            echo Remove {}
            rm {}
        fi
    '

```

### Check N50 of assemblies

```shell
cd ~/data/Fungi

cat ASSEMBLY/collect.csv | # head |
    sed '1d' |
    tsv-select -d "," -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        1>&2 echo "==> {}"

        find ASSEMBLY/{} -type f -name "*_genomic.fna.gz" |
            grep -v "_from_" | # exclude CDS and rna
            xargs cat |
            faops n50 -C -S stdin |
            (echo -e "name\t{}" && cat) |
            datamash transpose
    ' |
    tsv-uniq | # keep the first header
    tee ASSEMBLY/n50.tsv

cat ASSEMBLY/n50.tsv |
    tsv-filter \
        -H --or \
        --le 4:1000 \
        --ge 2:100000 |
    tsv-filter -H --ge 3:1000000 |
    tr "\t" "," \
    > ASSEMBLY/n50.pass.csv

wc -l ASSEMBLY/n50* ASSEMBLY/collect.csv
#   3085 ASSEMBLY/n50.pass.csv
#   3601 ASSEMBLY/n50.tsv
#   3601 ASSEMBLY/collect.csv

# Strains without protein annotations
for STRAIN in $(cat ASSEMBLY/n50.pass.csv | cut -d, -f 1); do
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_protein.faa.gz" > /dev/null; then
        echo ${STRAIN}
    fi
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_cds_from_genomic.fna.gz" > /dev/null; then
        echo ${STRAIN}
    fi
done |
    tsv-uniq \
    > ASSEMBLY/omit.lst
wc -l ASSEMBLY/omit.lst
#1832 ASSEMBLY/omit.lst

tsv-join \
    ASSEMBLY/collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv \
    > summary/collect.pass.csv

wc -l summary/collect.pass.csv
#3085 summary/collect.pass.csv

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

# rsync -avP -e "ssh -T -c chacha20-poly1305@openssh.com -o Compression=no -x" \
#   wangq@202.119.37.251:data/Fungi/ASSEMBLY/ ~/data/Fungi/ASSEMBLY

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
# 3579

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
#520

cat summary/collect.pass.csv |
    tsv-select -H -d, -f Taxid |
    sed '1d' |
    tsv-uniq |
    nwr append stdin -r species |
    tsv-select -f 2 |
    tsv-uniq |
    wc -l
#519

cat summary/collect.pass.csv |
    tsv-filter -H -d, --or \
        --istr-eq "RefSeq_category:reference genome" --istr-eq "RefSeq_category:representative genome" |
    grep -v "symbiont " |
    tsv-select -H -d, -f name |
    sed '1d' \
    > summary/representative.lst

wc -l summary/representative.lst
#509 summary/representative.lst

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
|--------------------------|------:|
| Aspergillus flavus       |   144 |
| Aspergillus fumigatus    |    80 |
| Aspergillus niger        |    96 |
| Aspergillus oryzae       |    90 |
| Botryosphaeria dothidea  |   131 |
| Candida albicans         |    64 |
| Cryphonectria parasitica |    92 |
| Fusarium graminearum     |   114 |
| Komagataella phaffii     |   129 |
| Ophidiomyces ophidiicola |    63 |
| Parastagonospora nodorum |   171 |
| Penicillium chrysogenum  |    78 |
| Pyricularia oryzae       |   251 |
| Rhodotorula mucilaginosa |    66 |
| Saccharomyces cerevisiae |   111 |
| Saitozyma podzolica      |    56 |
| Venturia inaequalis      |    85 |

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
| 451869  | Botryosphaeriales | 4        | 131      |
| 5114    | Diaporthales      | 5        | 71       |
| 5042    | Eurotiales        | 97       | 675      |
| 5125    | Hypocreales       | 75       | 290      |
| 639021  | Magnaporthales    | 9        | 173      |
| 33183   | Onygenales        | 28       | 106      |
| 92860   | Pleosporales      | 44       | 271      |
| 4892    | Saccharomycetales | 240      | 727      |
| 231213  | Sporidiobolales   | 3        | 67       |
| 5234    | Tremellales       | 18       | 92       |

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
| 5598    | Alternaria                | 17       | 60       |
| 5052    | Aspergillus               | 66       | 454      |
| 45132   | Botryosphaeria            | 1        | 128      |
| 5475    | Candida                   | 51       | 137      |
| 2964429 | Candida/Metschnikowiaceae | 5        | 52       |
| 5115    | Cryphonectria             | 2        | 68       |
| 5506    | Fusarium                  | 40       | 224      |
| 460517  | Komagataella              | 3        | 127      |
| 27320   | Metschnikowia             | 6        | 53       |
| 1387562 | Ophidiomyces              | 1        | 63       |
| 1351751 | Parastagonospora          | 2        | 163      |
| 5073    | Penicillium               | 19       | 202      |
| 48558   | Pyricularia               | 8        | 172      |
| 5533    | Rhodotorula               | 3        | 67       |
| 4930    | Saccharomyces             | 121      | 203      |

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

# rank::col
ARRAY=(
    'order::6'
    'family::5'
#    'genus::4'
    'species::3'
)

rm mash.condensed.map
CUR_TREE=tree.nwk

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../summary/strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick mash.${GROUP_NAME}.newick
    cat condense.map >> mash.condensed.map

    CUR_TREE=mash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 mash.species.newick |
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
#105     2670

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
#87

find NR -name "redundant.lst" -empty | wc -l
#18

find NR -name "NR.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/NR.lst

wc -l summary/NR.lst
#504 summary/NR.lst

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

| #tax_id | genus          | #species | #strains |
|---------|----------------|----------|----------|
| 5598    | Alternaria     | 17       | 34       |
| 5052    | Aspergillus    | 56       | 67       |
| 5579    | Aureobasidium  | 6        | 16       |
| 5475    | Candida        | 18       | 21       |
| 5455    | Colletotrichum | 14       | 14       |
| 5506    | Fusarium       | 32       | 39       |
| 5073    | Penicillium    | 13       | 20       |
| 1322061 | Rhizoctonia    | 6        | 14       |
| 4930    | Saccharomyces  | 10       | 19       |
| 5543    | Trichoderma    | 10       | 10       |
| 1047167 | Zymoseptoria   | 5        | 10       |

## Groups and targets

Review `summary/collect.pass.csv` and `tree/groups.tsv`
