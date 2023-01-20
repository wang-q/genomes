# Bacteria

All genomes of *Bacteria* and *Archaea*, species by species

## Strain info

* [Bacteria](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2)
* [Archaea](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2157)

### `date`

Fri Jan 20 00:52:35 CST 2023

### List all ranks

```shell
nwr member Bacteria |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat

nwr member Archaea |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat

```

| rank             |  count |
|------------------|-------:|
| superkingdom     |      1 |
| phylum           |    173 |
| class            |    125 |
| order            |    293 |
| family           |    739 |
| no rank          |   6578 |
| species          | 110192 |
| genus            |   4959 |
| clade            |    133 |
| strain           |  41127 |
| varietas         |     24 |
| isolate          |    456 |
| subspecies       |    707 |
| subclass         |      6 |
| forma            |      4 |
| species group    |     99 |
| species subgroup |     31 |
| suborder         |      7 |
| biotype          |      7 |
| serotype         |    257 |
| serogroup        |    143 |
| subphylum        |      1 |
| subgenus         |      1 |
| tribe            |      2 |
| pathogroup       |      5 |
| subfamily        |      1 |

| rank          | count |
|---------------|------:|
| superkingdom  |     1 |
| phylum        |    42 |
| order         |    60 |
| no rank       |   327 |
| species       |  3414 |
| class         |    34 |
| family        |    82 |
| genus         |   257 |
| clade         |    44 |
| strain        |   353 |
| species group |     2 |
| isolate       |     6 |

### Species with assemblies

Four levels:

* '>= 100 genomes'
    * With strain ID
    * assembly_level: 'Complete Genome', 'Chromosome'
    * genome_rep: 'Full'
* '>= 2 genomes'
    * assembly_level: 'Complete Genome', 'Chromosome'

```shell
mkdir -p ~/data/Bacteria/summary
cd ~/data/Bacteria/summary

# should have a valid name of genus
# NCBI division Bacteria includes Bacteria and Archaea
nwr member Bacteria Archaea -r genus |
    grep -v -i "Candidatus " |
    grep -v -i "candidate " |
    grep -v " sp." |
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#4461 genus.list

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
    echo "
        SELECT
            species_id,
            species,
            COUNT(DISTINCT tax_id) AS count -- with strain ID
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
        GROUP BY species_id
        HAVING count >= 100
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > L1.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
        GROUP BY species_id
        HAVING count >= 100
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done  |
    tsv-sort -k2,2 \
    > L2.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND assembly_level IN ('Complete Genome', 'Chromosome', 'Scaffold')
            AND genome_rep IN ('Full') -- fully representative
        GROUP BY species_id
        HAVING count >= 100
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done  |
    tsv-sort -k2,2 \
    > L3.tsv

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
        GROUP BY species_id
        HAVING count >= 2
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done  |
    tsv-sort -k2,2 \
    > L4.tsv

wc -l L*.tsv
#    3 L1.tsv
#   45 L2.tsv
#  111 L3.tsv
# 1708 L4.tsv
# 1867 total

for L in L1 L2 L3 L4; do
    cat ${L}.tsv |
        tsv-summarize --sum 3
done
#817
#15391
#80067
#28548

cat L3.tsv |
    tsv-join -f L2.tsv -k 1 -e |
    tsv-summarize --sum 3
#11931

cat L4.tsv |
    tsv-join -f L2.tsv -k 1 -e |
    tsv-join -f L3.tsv -k 1 -e |
    tsv-summarize --sum 3
#10367

```

### Model organisms

```shell
cd ~/data/Bacteria/summary

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

cat reference.tsv |
    sed '1s/^/#/' |
    nwr append stdin -r phylum -r class |
    tsv-select -H -f 1,2,phylum,class |
    parallel --col-sep "\t" -j 1 '
        if [[ "{3}" != "Proteobacteria" ]]; then
            printf "%s\t%s\t%s\n" {1} {2} {3}
        else
            printf "%s\t%s\t%s\n" {1} {2} {4}
        fi
    ' |
    mlr --itsv --omd cat

```

| #tax_id | organism_name                                                    | phylum                |
|---------|------------------------------------------------------------------|-----------------------|
| 565050  | Caulobacter vibrioides NA1000                                    | Alphaproteobacteria   |
| 192222  | Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819      | Epsilonproteobacteria |
| 208964  | Pseudomonas aeruginosa PAO1                                      | Gammaproteobacteria   |
| 871585  | Acinetobacter pittii PHEA-2                                      | Gammaproteobacteria   |
| 511145  | Escherichia coli str. K-12 substr. MG1655                        | Gammaproteobacteria   |
| 386585  | Escherichia coli O157:H7 str. Sakai                              | Gammaproteobacteria   |
| 1125630 | Klebsiella pneumoniae subsp. pneumoniae HS11286                  | Gammaproteobacteria   |
| 99287   | Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 | Gammaproteobacteria   |
| 198214  | Shigella flexneri 2a str. 301                                    | Gammaproteobacteria   |
| 227377  | Coxiella burnetii RSA 493                                        | Gammaproteobacteria   |
| 272561  | Chlamydia trachomatis D/UW-3/CX                                  | Chlamydiae            |
| 93061   | Staphylococcus aureus subsp. aureus NCTC 8325                    | Bacillota             |
| 224308  | Bacillus subtilis subsp. subtilis str. 168                       | Bacillota             |
| 169963  | Listeria monocytogenes EGD-e                                     | Bacillota             |
| 83332   | Mycobacterium tuberculosis H37Rv                                 | Actinomycetota        |

## Download all assemblies

### Create assembly.tsv

Three levels:

* '>= 100 genomes'
    * assembly_level: 'Complete Genome', 'Chromosome'
    * assembly_level: NOT 'contig' AND genome_rep: 'Full'
* '>= 2 genomes'
    * assembly_level: 'Complete Genome', 'Chromosome'


```shell
cd ~/data/Bacteria/summary

cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,assembly_level \
    > raw.tsv

# L2
SPECIES=$(
    cat L2.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-filter --regex '2:^[A-Za-z ]+$' \
    >> raw.tsv

# L3
SPECIES=$(
    cat L3.tsv |
        tsv-join -f L2.tsv -k 1 -e |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome', 'Scaffold')
        AND genome_rep IN ('Full') -- fully representative
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-filter --regex '2:^[A-Za-z ]+$' \
    >> raw.tsv

# L4
SPECIES=$(
    cat L4.tsv |
        tsv-join -f L2.tsv -k 1 -e |
        tsv-join -f L3.tsv -k 1 -e |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-filter --regex '2:^[A-Za-z ]+$' \
    >> raw.tsv

datamash check < raw.tsv
#37552 lines, 5 fields

# Create abbr.
cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
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
    > Bacteria.assembly.tsv

datamash check < Bacteria.assembly.tsv
#37477 lines, 4 fields

# find potential duplicate strains or assemblies
cat Bacteria.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Bacteria.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Bacteria.assembly.tsv
# cp Bacteria.assembly.tsv ~/Scripts/genomes/assembly

# Comment out unneeded strains

# Cleaning
rm raw*.*sv

```

### rsync and check

```shell
cd ~/data/Bacteria

cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    tsv-filter -v --str-in-fld 2:http

nwr assembly ~/Scripts/genomes/assembly/Bacteria.assembly.tsv \
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
proxychains4 bash ASSEMBLY/rsync.sh

# md5
# rm ASSEMBLY/check.list
bash ASSEMBLY/check.sh

# collect
bash ASSEMBLY/collect.sh

```

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Bacteria/ \
    wangq@202.119.37.251:data/Bacteria

# rsync -avP wangq@202.119.37.251:data/Bacteria/ ~/data/Bacteria

# rsync -avP -e "ssh -T -c chacha20-poly1305@openssh.com -o Compression=no -x" \
#   wangq@202.119.37.251:data/Bacteria/ASSEMBLY/ ~/data/Bacteria/ASSEMBLY

```

## BioSample

```shell
cd ~/data/Bacteria

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

find biosample -name "SAM*.txt" | wc -l
# 1956

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
    > attributes.lst

cat attributes.lst |
    (echo -e "BioSample" && cat) |
    tr '\n' '\t' |
    sed 's/\t$/\n/' \
    > Pseudomonas.biosample.tsv

find biosample -name "SAM*.txt" |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        >&2 echo {/.}
        cat {} |
            perl -nl -MPath::Tiny -e '\''
                BEGIN {
                    our @keys = grep {/\S/} path(q{attributes.lst})->lines({chomp => 1});
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
    >> Pseudomonas.biosample.tsv

```
