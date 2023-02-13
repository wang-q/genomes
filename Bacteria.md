# Bacteria

All genomes of *Bacteria* and *Archaea*, species by species.

Download all genomes and analyze representative strains.

<!-- toc -->

- [Strain info](#strain-info)
  * [List all ranks](#list-all-ranks)
  * [Species with assemblies](#species-with-assemblies)
  * [Model organisms](#model-organisms)
- [Download all assemblies](#download-all-assemblies)
  * [Create assembly.tsv](#create-assemblytsv)
  * [rsync and check](#rsync-and-check)
  * [Rsync to hpcc](#rsync-to-hpcc)
- [BioSample](#biosample)
- [Count species and strains](#count-species-and-strains)
  * [Order](#order)
  * [Genus](#genus)
- [MinHash](#minhash)
  * [Compute MinHash](#compute-minhash)
  * [Raw phylo-tree of representative assemblies](#raw-phylo-tree-of-representative-assemblies)
  * [Tweak the mash tree](#tweak-the-mash-tree)
- [Non-redundant strains within species](#non-redundant-strains-within-species)
  * [Genus](#genus-1)
- [Collect proteins](#collect-proteins)
  * [`all.pro.fa`](#allprofa)
  * [`all.replace.fa`](#allreplacefa)
  * [`all.info.tsv`](#allinfotsv)
- [Phylogenetics with bac120](#phylogenetics-with-bac120)
  * [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
  * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
  * [Tweak the concat tree](#tweak-the-concat-tree)
- [Early divergence of Bacteria](#early-divergence-of-bacteria)
  * [Terrabacteria group](#terrabacteria-group)
  * [Proteobacteria](#proteobacteria)
- [InterProScan on all proteins of representative and typical strains](#interproscan-on-all-proteins-of-representative-and-typical-strains)
  * [`interproscan.sh`](#interproscansh)

<!-- tocstop -->

## Strain info

* [Bacteria](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2)
* [Archaea](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2157)

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
| order            |    295 |
| family           |    747 |
| no rank          |   6601 |
| species          | 110253 |
| genus            |   4972 |
| clade            |    133 |
| strain           |  41193 |
| varietas         |     24 |
| isolate          |    456 |
| subspecies       |    708 |
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
| no rank       |   328 |
| species       |  3416 |
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
    1. With strain ID
    2. assembly_level: 'Complete Genome', 'Chromosome'
    3. genome_rep: 'Full'
* '>= 2 genomes'
    4. assembly_level: 'Complete Genome', 'Chromosome'

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
#4470 genus.list

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
done |
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
done |
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
done |
    tsv-sort -k2,2 \
    > L4.tsv

wc -l L*.tsv
#    3 L1.tsv
#   43 L2.tsv
#  114 L3.tsv
# 1734 L4.tsv
# 1894 total

for L in L1 L2 L3 L4; do
    cat ${L}.tsv |
        tsv-summarize --sum 3
done
#817
#15299
#80300
#28732

cat L3.tsv |
    tsv-join -f L2.tsv -k 1 -e |
    tsv-summarize --sum 3
#13133

cat L4.tsv |
    tsv-join -f L2.tsv -k 1 -e |
    tsv-join -f L3.tsv -k 1 -e |
    tsv-summarize --sum 3
#10473

```

* Some species are divided into several separate species

```shell
cd ~/data/Bacteria/summary

nwr member Pseudomonas -r species |
    grep -v " sp." |
    grep -E "syringae|genomosp"
#251701  Pseudomonas syringae group genomosp. 3  species Bacteria
#251699  Pseudomonas syringae group genomosp. 7  species Bacteria
#317659  Pseudomonas syringae pv. coryli species Bacteria
#317     Pseudomonas syringae    species Bacteria

echo "
    SELECT
        *
    FROM ar
    WHERE 1=1
        AND species_id IN (251701,251699,317659)
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    GROUP BY species_id
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
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
        if [[ "{3}" == "Proteobacteria" || "{3}" == "Pseudomonadota" ]]; then
            printf "%s\t%s\t%s\n" {1} {2} {4}
        else
            printf "%s\t%s\t%s\n" {1} {2} {3}
        fi
    ' |
    mlr --itsv --omd cat

cp reference.tsv ~/Scripts/genomes/assembly/Bacteria.reference.tsv

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
    tsv-filter -H --regex '2:^[A-Z]' \
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
    tsv-filter -H --regex '2:^[A-Z]' \
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
    tsv-filter -H --regex '2:^[A-Z]' \
    >> raw.tsv

datamash check < raw.tsv
#38788 lines, 5 fields

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
    keep-header -- sort -k3,3 -k1,1 |
    tsv-filter -H --regex '1:^[A-Z]' \
    > Bacteria.assembly.tsv

datamash check < Bacteria.assembly.tsv
#38662 lines, 4 fields

# find potential duplicate strains or assemblies
cat Bacteria.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Bacteria.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

cat Bacteria.assembly.tsv |
    tsv-filter --or --str-in-fld 1:genomosp --str-in-fld 1:genomovar

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

### Rsync to hpcc

```bash
rsync -avP \
    ~/data/Bacteria/ \
    wangq@202.119.37.251:data/Bacteria

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Bacteria/ \
    wangq@58.213.64.36:data/Bacteria

# rsync -avP wangq@202.119.37.251:data/Bacteria/ ~/data/Bacteria

# rsync -avP -e "ssh -T -c chacha20-poly1305@openssh.com -o Compression=no -x" \
#   wangq@202.119.37.251:data/Bacteria/ASSEMBLY/ ~/data/Bacteria/ASSEMBLY

```

## BioSample

ENA's BioSample missed many strains, so NCBI's was used.

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

# Allowing samples not in the list
find biosample -name "SAM*.txt" | wc -l
# 39817

find biosample -name "SAM*.txt" |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat {} |
            perl -nl -e '\''
                print $1 if m{\s+\/([\w_ ]+)=};
            '\''
    ' |
    tsv-uniq --at-least 500 | # ignore rare attributes
    grep -v "^INSDC" |
    grep -v "^ENA" \
    > ASSEMBLY/attributes.lst

cat ASSEMBLY/attributes.lst |
    (echo -e "BioSample" && cat) |
    tr '\n' '\t' |
    sed 's/\t$/\n/' \
    > ASSEMBLY/biosample.tsv

find biosample -name "SAM*.txt" |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        >&2 echo {/.}
        cat {} |
            perl -nl -MPath::Tiny -e '\''
                BEGIN {
                    our @keys = grep {/\S/} path(q{ASSEMBLY/attributes.lst})->lines({chomp => 1});
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
    >> ASSEMBLY/biosample.tsv

```

## Count species and strains

```shell
cd ~/data/Bacteria

echo "
    SELECT
        COUNT(DISTINCT species_id)
    FROM ar
    WHERE 1=1
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
#19169

echo "
    SELECT
        COUNT(DISTINCT species_id)
    FROM ar
    WHERE 1=1
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
#6580

echo "
    SELECT
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '%symbiont %'
        AND refseq_category IN ('reference genome', 'representative genome')
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-filter --str-not-in-fld 1:"[" \
    > summary/assembly_accession.lst

wc -l summary/assembly_accession.lst
#5275 assembly_accession.lst

cat ASSEMBLY/collect.csv |
    tsv-select -H -d, -f Taxid |
    sed '1d' |
    tsv-uniq |
    nwr append stdin -r species |
    tsv-select -f 2 |
    tsv-uniq |
    wc -l
#1689

cat ASSEMBLY/collect.csv |
    tsv-filter -H -d, --or \
        --istr-eq "RefSeq_category:reference genome" --istr-eq "RefSeq_category:representative genome" |
    grep -v "symbiont " |
    tsv-select -H -d, -f name |
    sed '1d' \
    > summary/representative_assembly.lst

wc -l summary/representative_assembly.lst
#1432 summary/representative_assembly.lst

cat ASSEMBLY/collect.csv |
    sed -e '1d' |
    tr "," "\t" |
    tsv-select -f 1,3 |
    nwr append stdin -c 2 -r species -r genus -r family -r order \
    > summary/strains.taxon.tsv

cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    tsv-summarize -H -g 3 --count |
    tsv-filter -H --gt 2:200 |
    mlr --itsv --omd cat

```

| organism                      | count |
|-------------------------------|------:|
| Acinetobacter baumannii       |   531 |
| Bacillus cereus               |   413 |
| Bacillus subtilis             |   261 |
| Bacillus toyonensis           |   238 |
| Bacillus velezensis           |   278 |
| Bacteroides ovatus            |   201 |
| Bifidobacterium longum        |   224 |
| Bordetella pertussis          |   593 |
| Campylobacter coli            |   240 |
| Campylobacter jejuni          |   284 |
| Cronobacter sakazakii         |   208 |
| Enterobacter hormaechei       |   270 |
| Enterococcus faecium          |   313 |
| Erwinia amylovora             |   211 |
| Escherichia coli              |  2739 |
| Helicobacter pylori           |   395 |
| Klebsiella michiganensis      |   211 |
| Klebsiella pneumoniae         |  1611 |
| Klebsiella variicola          |   252 |
| Lactiplantibacillus plantarum |   211 |
| Limosilactobacillus reuteri   |   249 |
| Listeria monocytogenes        |   354 |
| Mycobacterium tuberculosis    |   546 |
| Mycobacteroides abscessus     |   946 |
| Pseudomonas aeruginosa        |   648 |
| Pseudomonas syringae          |   257 |
| Rhizobium leguminosarum       |   372 |
| Salmonella enterica           |  1377 |
| Shigella sonnei               |  1113 |
| Staphylococcus aureus         |  1124 |
| Staphylococcus haemolyticus   |   221 |
| Stenotrophomonas maltophilia  |   373 |
| Streptococcus pneumoniae      |   208 |
| Streptococcus pyogenes        |   263 |

### Order

```shell
cd ~/data/Bacteria

# Group by order
cat ASSEMBLY/collect.csv |
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
        n_species=$(cat ASSEMBLY/collect.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r order -r species |
            grep {} |
            tsv-select -f 1,3 |
            tsv-uniq |
            wc -l)

        n_strains=$(cat ASSEMBLY/collect.csv |
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
    tsv-filter --ge 4:100 |
    (echo -e '#tax_id\torder\t#species\t#strains' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | order               | #species | #strains |
|---------|---------------------|----------|----------|
| 135624  | Aeromonadales       | 27       | 277      |
| 135622  | Alteromonadales     | 56       | 135      |
| 1385    | Bacillales          | 503      | 4854     |
| 171549  | Bacteroidales       | 102      | 1127     |
| 85004   | Bifidobacteriales   | 89       | 797      |
| 80840   | Burkholderiales     | 233      | 2060     |
| 213849  | Campylobacterales   | 204      | 1090     |
| 51291   | Chlamydiales        | 123      | 262      |
| 84999   | Coriobacteriales    | 4        | 112      |
| 85007   | Corynebacteriales   | 350      | 2170     |
| 91347   | Enterobacterales    | 1128     | 10949    |
| 186802  | Eubacteriales       | 144      | 571      |
| 200644  | Flavobacteriales    | 79       | 500      |
| 356     | Hyphomicrobiales    | 260      | 1096     |
| 186826  | Lactobacillales     | 463      | 3980     |
| 118969  | Legionellales       | 32       | 177      |
| 1643688 | Leptospirales       | 104      | 201      |
| 85006   | Micrococcales       | 66       | 196      |
| 2887326 | Moraxellales        | 60       | 878      |
| 2085    | Mycoplasmatales     | 71       | 235      |
| 2790996 | Mycoplasmoidales    | 45       | 167      |
| 206351  | Neisseriales        | 54       | 350      |
| 135625  | Pasteurellales      | 91       | 633      |
| 85009   | Propionibacteriales | 90       | 163      |
| 72274   | Pseudomonadales     | 291      | 1727     |
| 204455  | Rhodobacterales     | 39       | 177      |
| 766     | Rickettsiales       | 97       | 219      |
| 136     | Spirochaetales      | 69       | 334      |
| 85011   | Streptomycetales    | 69       | 165      |
| 72273   | Thiotrichales       | 54       | 286      |
| 48461   | Verrucomicrobiales  | 2        | 175      |
| 135623  | Vibrionales         | 88       | 662      |
| 135614  | Xanthomonadales     | 141      | 926      |

### Genus

```shell
cd ~/data/Bacteria

# Group by genus
cat ASSEMBLY/collect.csv |
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
        n_species=$(cat ASSEMBLY/collect.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r genus -r species |
            grep {} |
            tsv-select -f 1,3 |
            tsv-uniq |
            wc -l)

        n_strains=$(cat ASSEMBLY/collect.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r genus |
            grep {} |
            wc -l)

        printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' |
    nwr append stdin --id |
    tsv-select -f 5,4,2,3 |
    tsv-sort -k2,2 |
    tsv-filter --ge 4:100 |
    (echo -e '#tax_id\tgenus\t#species\t#strains' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | genus               | #species | #strains |
|---------|---------------------|----------|----------|
| 222     | Achromobacter       | 10       | 125      |
| 469     | Acinetobacter       | 52       | 804      |
| 642     | Aeromonas           | 26       | 275      |
| 357     | Agrobacterium       | 16       | 160      |
| 239934  | Akkermansia         | 2        | 175      |
| 1386    | Bacillus            | 194      | 2068     |
| 816     | Bacteroides         | 42       | 714      |
| 1678    | Bifidobacterium     | 83       | 780      |
| 517     | Bordetella          | 23       | 807      |
| 64895   | Borreliella         | 15       | 120      |
| 234     | Brucella            | 95       | 278      |
| 32008   | Burkholderia        | 103      | 707      |
| 194     | Campylobacter       | 107      | 641      |
| 810     | Chlamydia           | 122      | 260      |
| 544     | Citrobacter         | 16       | 208      |
| 1870884 | Clostridioides      | 19       | 151      |
| 1485    | Clostridium         | 62       | 304      |
| 102106  | Collinsella         | 4        | 112      |
| 1716    | Corynebacterium     | 83       | 358      |
| 413496  | Cronobacter         | 8        | 210      |
| 1912216 | Cutibacterium       | 79       | 119      |
| 547     | Enterobacter        | 31       | 837      |
| 1350    | Enterococcus        | 32       | 675      |
| 551     | Erwinia             | 15       | 234      |
| 561     | Escherichia         | 191      | 2912     |
| 237     | Flavobacterium      | 11       | 119      |
| 262     | Francisella         | 50       | 202      |
| 724     | Haemophilus         | 16       | 217      |
| 209     | Helicobacter        | 82       | 423      |
| 570     | Klebsiella          | 72       | 2526     |
| 2759736 | Lacticaseibacillus  | 24       | 267      |
| 2767842 | Lactiplantibacillus | 14       | 228      |
| 1578    | Lactobacillus       | 64       | 537      |
| 1357    | Lactococcus         | 26       | 195      |
| 445     | Legionella          | 19       | 132      |
| 171     | Leptospira          | 104      | 201      |
| 2767887 | Ligilactobacillus   | 11       | 127      |
| 2742598 | Limosilactobacillus | 16       | 306      |
| 1637    | Listeria            | 103      | 388      |
| 75984   | Mannheimia          | 15       | 116      |
| 1763    | Mycobacterium       | 202      | 708      |
| 670516  | Mycobacteroides     | 18       | 957      |
| 2995234 | Mycoplasmoides      | 36       | 136      |
| 2767358 | Mycoplasmopsis      | 29       | 108      |
| 482     | Neisseria           | 35       | 310      |
| 53335   | Pantoea             | 17       | 158      |
| 375288  | Parabacteroides     | 9        | 145      |
| 745     | Pasteurella         | 17       | 145      |
| 909656  | Phocaeicola         | 6        | 153      |
| 2800373 | Priestia            | 11       | 136      |
| 583     | Proteus             | 7        | 132      |
| 286     | Pseudomonas         | 271      | 1672     |
| 48736   | Ralstonia           | 17       | 141      |
| 379     | Rhizobium           | 39       | 414      |
| 590     | Salmonella          | 479      | 1389     |
| 613     | Serratia            | 21       | 191      |
| 620     | Shigella            | 25       | 1243     |
| 1279    | Staphylococcus      | 98       | 1948     |
| 40323   | Stenotrophomonas    | 11       | 383      |
| 1301    | Streptococcus       | 183      | 1286     |
| 1883    | Streptomyces        | 68       | 163      |
| 104267  | Tenacibaculum       | 7        | 128      |
| 157     | Treponema           | 30       | 141      |
| 662     | Vibrio              | 76       | 627      |
| 338     | Xanthomonas         | 108      | 441      |
| 629     | Yersinia            | 65       | 368      |

## MinHash

### Compute MinHash

```shell
mkdir -p ~/data/Bacteria/mash
cd ~/data/Bacteria/mash

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
mkdir -p ~/data/Bacteria/tree
cd ~/data/Bacteria/tree

mash triangle -E -p 8 -l <(
    cat ../summary/representative_assembly.lst |
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

        group <- cutree(clusters, h=0.4) # k=5
        groups <- as.data.frame(group)
        groups$ids <- rownames(groups)
        rownames(groups) <- NULL
        groups <- groups[order(groups$group), ]
        write_tsv(groups, "groups.tsv")
    '

```

### Tweak the mash tree

```shell
cd ~/data/Bacteria/tree

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
    rsvg-convert -o Bacteria.mash.png

```

## Non-redundant strains within species

```shell
cd ~/data/Bacteria

mkdir NR

cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    tsv-summarize -H -g 3 --count |
    tsv-filter -H --ge 2:2 \
    > NR/species.tsv

tsv-summarize NR/species.tsv -H --count --sum count
#count   count_sum
#1689    38659

# each species
cat NR/species.tsv | sed '1d' | tsv-select -f 1 | #head -n 5 |
while read SPECIES; do
    1>&2 echo "==> ${SPECIES}"

    SPECIES_=$(
        echo "${SPECIES}" |
            tr " " "_"
    )

    mkdir -p NR/${SPECIES_}

    cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
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
        tsv-filter --ff-str-ne 1:2 --le 3:0.01 \
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
                our %rep = map { ($_, 1) } path(q{summary/representative_assembly.lst})->lines({chomp => 1});
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

find NR -name "NR.lst" | xargs cat | wc -l
#8934

find NR -name "redundant.lst" -size +0 | wc -l
#1125

find NR -name "redundant.lst" -empty | wc -l
#564

find NR -name "NR.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/NR.lst

```

### Genus

```shell
cd ~/data/Bacteria

cat summary/genus.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(
            cat ASSEMBLY/collect.csv |
                sed "1d" |
                grep -F -f <( cat summary/NR.lst summary/representative_assembly.lst | sort | uniq ) |
                tsv-select -d, -f 3 |
                nwr append stdin -r genus -r species |
                grep {} |
                tsv-select -f 1,3 |
                tsv-uniq |
                wc -l
        )

        n_strains=$(
            cat ASSEMBLY/collect.csv |
                sed "1d" |
                grep -F -f <( cat summary/NR.lst summary/representative_assembly.lst | sort | uniq ) |
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
    tsv-filter --ge 4:50 |
    (echo -e '#tax_id\tgenus\t#species\t#strains' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | genus               | #species | #strains |
|---------|---------------------|----------|----------|
| 469     | Acinetobacter       | 37       | 245      |
| 642     | Aeromonas           | 15       | 169      |
| 357     | Agrobacterium       | 14       | 79       |
| 1386    | Bacillus            | 97       | 389      |
| 816     | Bacteroides         | 30       | 221      |
| 1678    | Bifidobacterium     | 46       | 279      |
| 34098   | Blattabacterium     | 5        | 52       |
| 32199   | Buchnera            | 55       | 55       |
| 32008   | Burkholderia        | 33       | 124      |
| 194     | Campylobacter       | 54       | 107      |
| 544     | Citrobacter         | 15       | 91       |
| 1485    | Clostridium         | 48       | 134      |
| 102106  | Collinsella         | 3        | 55       |
| 1716    | Corynebacterium     | 51       | 96       |
| 547     | Enterobacter        | 18       | 255      |
| 1350    | Enterococcus        | 18       | 83       |
| 561     | Escherichia         | 21       | 197      |
| 724     | Haemophilus         | 12       | 79       |
| 209     | Helicobacter        | 55       | 289      |
| 570     | Klebsiella          | 18       | 130      |
| 1578    | Lactobacillus       | 36       | 112      |
| 1357    | Lactococcus         | 12       | 53       |
| 2767887 | Ligilactobacillus   | 9        | 50       |
| 2742598 | Limosilactobacillus | 7        | 96       |
| 909656  | Phocaeicola         | 5        | 56       |
| 286     | Pseudomonas         | 174      | 557      |
| 379     | Rhizobium           | 31       | 160      |
| 1827    | Rhodococcus         | 18       | 54       |
| 590     | Salmonella          | 37       | 73       |
| 613     | Serratia            | 14       | 69       |
| 1279    | Staphylococcus      | 41       | 100      |
| 40323   | Stenotrophomonas    | 9        | 159      |
| 1301    | Streptococcus       | 70       | 374      |
| 1883    | Streptomyces        | 57       | 88       |
| 662     | Vibrio              | 57       | 366      |
| 338     | Xanthomonas         | 55       | 83       |

## Collect proteins

### `all.pro.fa`

```shell script
cd ~/data/Bacteria

mkdir -p PROTEINS

for STRAIN in $(cat summary/representative_assembly.lst); do
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
#5.2M

gzip -dcf PROTEINS/all.pro.fa.gz |
    grep "^>" |
    tsv-uniq |
    wc -l |
    numfmt --to=si
#5.1M

# annotations may be different
gzip -dcf PROTEINS/all.uniq.fa.gz |
    grep "^>" |
    wc -l |
    numfmt --to=si
#5.1M

# ribonuclease
gzip -dcf PROTEINS/all.pro.fa.gz |
    grep "ribonuclease" |
    grep -v "deoxyribonuclease" |
    perl -nl -e 's/^>\w+\.\d+\s+//g; print' |
    perl -nl -e 's/\s+\[.+?\]$//g; print' |
    perl -nl -e 's/MULTISPECIES: //g; print' |
    sort |
    uniq -c |
    sort -nr

```

### `all.replace.fa`

```shell
cd ~/data/Bacteria

rm PROTEINS/all.strain.tsv PROTEINS/all.replace.fa.gz
for STRAIN in $(cat summary/representative_assembly.lst); do
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
#5.2M

(echo -e "#name\tstrain" && cat PROTEINS/all.strain.tsv)  \
    > temp &&
    mv temp PROTEINS/all.strain.tsv

faops size PROTEINS/all.replace.fa.gz > PROTEINS/all.replace.sizes

(echo -e "#name\tsize" && cat PROTEINS/all.replace.sizes) > PROTEINS/all.size.tsv

rm PROTEINS/all.replace.sizes

```

### `all.info.tsv`

```shell
cd ~/data/Bacteria

for STRAIN in $(cat summary/representative_assembly.lst); do
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
#5.2M

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
#5.2M

```

## Phylogenetics with bac120

### Find corresponding proteins by `hmmsearch`

* Download HMM models as described in [`HMM.md`](HMM.md)

* The `E_VALUE` was manually adjusted to 1e-20 to reach a balance between sensitivity and
  speciality.

```shell
E_VALUE=1e-20

cd ~/data/Bacteria

# Find all genes
for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    >&2 echo "==> marker [${marker}]"

    mkdir -p PROTEINS/${marker}

    cat summary/representative_assembly.lst |
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

```shell script
cd ~/data/Bacteria

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat PROTEINS/{}/replace.tsv |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
# 1385    1429    2168.25

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        cat PROTEINS/{}/replace.tsv |
            wc -l
    ' |
    paste - - |
    tsv-filter --invert --ge 2:1000 --le 2:2000 |
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
                tsv-uniq
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

fasops concat PROTEINS/bac120.aln.fas summary/representative_assembly.lst -o PROTEINS/bac120.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in PROTEINS/bac120.aln.fa -out PROTEINS/bac120.trim.fa -automated1

faops size PROTEINS/bac120.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#94432
#15688

# To make it faster
FastTree -fastest -noml PROTEINS/bac120.trim.fa > PROTEINS/bac120.trim.newick

```

### Tweak the concat tree

```shell script
cd ~/data/Bacteria/tree

cp ../PROTEINS/bac120.trim.newick .

# rank::col
ARRAY=(
    'order::6'
    'family::5'
#    'genus::4'
    'species::3'
)

rm mash.condensed.map
CUR_TREE=bac120.trim.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../summary/strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick bac120.${GROUP_NAME}.newick
    cat condense.map >> bac120.condensed.map

    CUR_TREE=bac120.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 bac120.species.newick |
    rsvg-convert -o Bacteria.bac120.png

```

## Early divergence of Bacteria

![molbiolevolmsn247f02_ht.jpeg](images%2Fmolbiolevolmsn247f02_ht.jpeg)

### Terrabacteria group

Bacillota == Firmicutes

Actinomycetota == Actinobacteria

```shell
mkdir -p ~/data/Bacteria/Terrabacteria
cd ~/data/Bacteria/Terrabacteria

FAMILY=$(
    nwr member "Terrabacteria group" -r family |
        sed '1d' |
        cut -f 1 |
        nwr append stdin -r order |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND family_id IN ($FAMILY)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    grep -v -i "symbiont " |
    tsv-filter --str-not-in-fld 1:"[" \
    > tmp.tsv

cat ../ASSEMBLY/collect.csv |
    grep -F -f <(cut -f 2 tmp.tsv) |
    grep -F -f <(cat ../summary/NR.lst ../summary/representative_assembly.lst | sort | uniq) |
    tsv-select -d, -f 1,3 |
    tr "," "\t" |
    nwr append stdin -c 2 -r species -r genus -r family -r order |
    sed 's/Pseudomonas paraeruginosa/Pseudomonas aeruginosa/g' \
    > strains.taxon.tsv

wc -l strains.taxon.tsv
#3106 strains.taxon.tsv

```

### Proteobacteria

Pseudomonadota == Proteobacteria

```shell
mkdir -p ~/data/Bacteria/Proteobacteria
cd ~/data/Bacteria/Proteobacteria

FAMILY=$(
    nwr member Proteobacteria -r family |
        sed '1d' |
        cut -f 1 |
        nwr append stdin -r order |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND family_id IN ($FAMILY)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    grep -v -i "symbiont " |
    tsv-filter --str-not-in-fld 1:"[" \
    > tmp.tsv

cat ../ASSEMBLY/collect.csv |
    grep -F -f <(cut -f 2 tmp.tsv) |
    grep -F -f <(cat ../summary/NR.lst ../summary/representative_assembly.lst | sort | uniq) |
    tsv-select -d, -f 1,3 |
    tr "," "\t" |
    nwr append stdin -c 2 -r species -r genus -r family -r order \
    > strains.taxon.tsv

wc -l strains.taxon.tsv
#4905 strains.taxon.tsv

```

## InterProScan on all proteins of representative and typical strains

```shell
mkdir -p ~/data/Bacteria/STRAINS
cd ~/data/Bacteria/STRAINS

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
    > typical.lst

# Target genus
GENUS=(
    Pseudomonas
    Halopseudomonas
    Stutzerimonas
    Azotobacter

    Acinetobacter
    Stenotrophomonas
    Burkholderia
    Bordetella
)
#GENUS=$(IFS=, ; echo "${GENUS[*]}")

for G in "${GENUS[@]}"; do
echo "
    SELECT
        organism_name,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND genus IN ('$GENUS')
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    grep -v -i "symbiont " |
    tsv-filter --str-not-in-fld 1:"[" \
    > tmp.tsv

cat ../ASSEMBLY/collect.csv |
    grep -F -f <(cut -f 2 tmp.tsv) |
    grep -F -f <(cat typical.lst ../summary/NR.lst ../summary/representative_assembly.lst | sort | uniq) |
    tsv-select -d, -f 1,3 |
    tr "," "\t" |
    tsv-uniq |
    nwr append stdin -c 2 -r species -r genus -r family -r order \
    > strains.taxon.tsv

wc -l strains.taxon.tsv

```

### `interproscan.sh`

```shell
cd ~/data/Bacteria

for S in $(cat STRAINS/strains.taxon.tsv | cut -f 1); do
    mkdir -p STRAINS/${S}
    faops split-about ASSEMBLY/${S}/*_protein.faa.gz 200000 STRAINS/${S}/
done

for S in $(cat STRAINS/strains.taxon.tsv | cut -f 1); do
    for f in $(find STRAINS/${S}/ -maxdepth 1 -type f -name "[0-9]*.fa" | sort); do
        >&2 echo "==> ${f}"
        if [ -e ${f}.tsv ]; then
            >&2 echo ${f}
            continue
        fi

        bsub -q mpi -n 24 -J "${f}" "
            interproscan.sh --cpu 24 -dp -f tsv,json -i ${f} --output-file-base ${f}
        "
    done
done

#for S in $(cat summary/typical.lst | head -n 1); do
#    for f in $(find STRAINS/${S}/ -maxdepth 1 -type f -name "[0-9]*.fa" | sort | head -n 1); do
#        >&2 echo "==> ${f}"
#        if [ -e ${f}.tsv ]; then
#            >&2 echo ${f}
#            continue
#        fi
#
#        interproscan.sh --cpu 12 -dp -f tsv,json -i ${f} --output-file-base ${f}
#    done
#done

find STRAINS -type f -name "*.json" | sort |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        if [ $(({#} % 10)) -eq "0" ]; then
            >&2 printf "."
        fi
        pigz -p 3 {}
    '

# same protein may have multiple families
for S in $(cat STRAINS/strains.taxon.tsv | cut -f 1); do
    for f in $(find STRAINS/${S} -maxdepth 1 -type f -name "[0-9]*.json.gz" | sort); do
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
            tsv-uniq
    done \
        > STRAINS/${S}/family.tsv
done

```
