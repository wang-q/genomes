# *Pseudomonas* HGT

<!-- toc -->

- [Strain info](#strain-info)
    * [List all ranks](#list-all-ranks)
    * [Species with assemblies](#species-with-assemblies)
- [Download all assemblies](#download-all-assemblies)
    * [Create assembly.tsv](#create-assemblytsv)
    * [rsync and check](#rsync-and-check)
    * [Rsync to hpcc](#rsync-to-hpcc)
- [BioSample](#biosample)
- [Count and group strains](#count-and-group-strains)
    * [Check N50 of assemblies](#check-n50-of-assemblies)
    * [Order](#order)
    * [Genus](#genus)
    * [Strains](#strains)
- [NCBI taxonomy](#ncbi-taxonomy)
- [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash)
    * [Tweak the mash tree](#tweak-the-mash-tree)
- [Collect proteins](#collect-proteins)
    * [`all.pro.fa`](#allprofa)
    * [`all.replace.fa`](#allreplacefa)
    * [`all.info.tsv`](#allinfotsv)

<!-- tocstop -->

## Strain info

* [Pseudomonas](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=286)
* [Acinetobacter](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=469)

According to a recent [paper](https://doi.org/10.1128/mSystems.00543-20), there are some order-level
changes in Gammaproteobacteria. We include both old and new orders.

* Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
* New ones: Moraxellales, Kangiellales, and Pseudomonadales

Another recent [paper](https://doi.org/10.1099/ijsem.0.005542) proposed that an outlier clade of P.
aeruginosa (containing PA7) split into a new species, Pseudomonas paraeruginosa.

Similarly, a recent [paper](https://doi.org/10.1016/j.syapm.2021.126289) suggested that some species
of Pseudomonas should be transferred to Halopseudomonas and Stutzerimonas.

### List all ranks

```shell
mkdir -p ~/data/Pseudomonas
cd ~/data/Pseudomonas

nwr member Pseudomonas |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

nwr member Acinetobacter |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

nwr member Pseudomonas Acinetobacter -r "species group" -r "species subgroup" |
    tsv-select -f 1-3 |
    keep-header -- tsv-sort -k3,3 -k2,2 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    mlr --itsv --omd cat

```

| rank             | count |
|------------------|-------|
| genus            | 1     |
| species          | 425   |
| strain           | 724   |
| subspecies       | 13    |
| no rank          | 121   |
| species group    | 5     |
| species subgroup | 4     |
| isolate          | 1     |

| rank             | count |
|------------------|-------|
| genus            | 1     |
| species          | 117   |
| species group    | 2     |
| species subgroup | 3     |
| strain           | 1110  |
| no rank          | 2     |
| subspecies       | 1     |
| isolate          | 2     |

| #tax_id | sci_name                                 | rank             |
|---------|------------------------------------------|------------------|
| 2839056 | A. Taxon 24                              | species group    |
| 909768  | A. calcoaceticus/baumannii complex       | species group    |
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

Check also the order Pseudomonadales.

* Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
* New ones: Moraxellales, Kangiellales, and Pseudomonadales

```shell
mkdir -p ~/data/Pseudomonas/summary
cd ~/data/Pseudomonas/summary

SPECIES=$(
    nwr member \
        Cellvibrionales Oceanospirillales Alteromonadales \
        Moraxellales Kangiellales Pseudomonadales \
        -r species |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        grep -v " sp." |
        grep -v -E "\bbacterium\b" |
        grep -v -E "\bsymbiont\b" |
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
    tsv-filter -H --ge CHR:5 |
    tsv-filter -H --invert --str-in-fld species:Pseudomonas --lt RS:30 |
    tsv-filter -H --invert --str-in-fld species:Acinetobacter --lt RS:30 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    mlr --itsv --omd cat

```

| species_id | species                         | RS   | CHR |
|------------|---------------------------------|------|-----|
| 287        | P. aeruginosa                   | 7458 | 640 |
| 470        | A. baumannii                    | 7191 | 528 |
| 33069      | P. viridiflava                  | 1541 | 9   |
| 317        | P. syringae                     | 612  | 50  |
| 48296      | A. pittii                       | 385  | 35  |
| 294        | P. fluorescens                  | 262  | 40  |
| 303        | P. putida                       | 223  | 67  |
| 480        | Moraxella catarrhalis           | 209  | 16  |
| 106654     | A. nosocomialis                 | 173  | 11  |
| 316        | Stutzerimonas stutzeri          | 141  | 28  |
| 29438      | P. savastanoi                   | 120  | 5   |
| 38313      | Shewanella algae                | 112  | 24  |
| 587753     | P. chlororaphis                 | 107  | 61  |
| 756892     | A. indicus                      | 102  | 21  |
| 47877      | P. amygdali                     | 90   | 10  |
| 380021     | P. protegens                    | 83   | 27  |
| 47880      | P. fulva                        | 75   | 5   |
| 40215      | A. junii                        | 71   | 11  |
| 108980     | A. ursingii                     | 68   | 6   |
| 106648     | A. bereziniae                   | 63   | 5   |
| 1530123    | A. seifertii                    | 60   | 25  |
| 40216      | A. radioresistens               | 60   | 5   |
| 29430      | A. haemolyticus                 | 57   | 14  |
| 62322      | Shewanella baltica              | 48   | 12  |
| 40214      | A. johnsonii                    | 47   | 19  |
| 28090      | A. lwoffii                      | 47   | 11  |
| 76759      | P. monteilii                    | 47   | 9   |
| 296        | P. fragi                        | 40   | 7   |
| 476        | Moraxella bovis                 | 39   | 37  |
| 202956     | A. towneri                      | 32   | 8   |
| 930166     | P. brassicacearum               | 31   | 9   |
| 78327      | P. mosselii                     | 30   | 6   |
| 332186     | Shewanella xiamenensis          | 30   | 5   |
| 28108      | Alteromonas macleodii           | 26   | 16  |
| 43657      | Pseudoalteromonas luteoviolacea | 25   | 5   |
| 34062      | Moraxella osloensis             | 23   | 10  |
| 43662      | Pseudoalteromonas piscicida     | 19   | 6   |
| 314275     | Alteromonas mediterranea        | 17   | 16  |
| 24         | Shewanella putrefaciens         | 15   | 9   |
| 2968969    | Stutzerimonas frequens          | 15   | 5   |
| 44935      | Halomonas venusta               | 10   | 7   |
| 386891     | Moraxella bovoculi              | 9    | 7   |
| 1697053    | Thiopseudomonas alkaliphila     | 7    | 7   |

## Download all assemblies

### Create assembly.tsv

```shell
cd ~/data/Pseudomonas/summary

cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    head -n 1 \
    > Pseudomonas.assembly.tsv

# Pseudomonas aeruginosa PAO1 is in the reference list
cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    grep -F -f <(
        cat ~/Scripts/genomes/assembly/Bacteria.reference.tsv |
            tsv-select -H -f assembly_accession
        ) \
    >> Pseudomonas.assembly.tsv

# Pseudomonadales and close relatives
# Species with 2 or more complete genomes were retained.
cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    sed '1d' |
    tsv-join -d 3 \
        -k 2 -f <(
            cat species.count.tsv |
                tsv-filter -H --ge CHR:2 |
                sed '1d'
        )\
    >> Pseudomonas.assembly.tsv

# Also includes representative strains of Gammaproteobacteria.
# families not in our orders
FAMILY=$(
    nwr member Gammaproteobacteria -r family |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        sed '1d' |
        cut -f 1 |
        nwr append stdin -r order |
        tsv-filter --str-ne 2:"Cellvibrionales" --str-ne 2:"Oceanospirillales" --str-ne 2:"Alteromonadales" |
        tsv-filter --str-ne 2:"Moraxellales" --str-ne 2:"Kangiellales" --str-ne 2:"Pseudomonadales" |
        tsv-filter --str-ne 2:"NA" |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        assembly_accession
    FROM ar
    WHERE 1=1
        AND family_id IN ($FAMILY)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND refseq_category IN ('representative genome')
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    grep -v -i "symbiont " |
    tsv-filter --str-not-in-fld 1:"[" \
    > tmp.tsv

cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    grep -F -f tmp.tsv \
    >> Pseudomonas.assembly.tsv

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Pseudomonas.assembly.tsv
# cp Pseudomonas.assembly.tsv ~/Scripts/genomes/assembly

# Comment out unneeded strains

# Cleaning
rm tmp*.*sv

```

### rsync and check

```shell
cd ~/data/Pseudomonas

cat ~/Scripts/genomes/assembly/Pseudomonas.assembly.tsv |
    tsv-filter -v --str-in-fld 2:http

nwr assembly ~/Scripts/genomes/assembly/Pseudomonas.assembly.tsv \
    -o ASSEMBLY

# Copy downloaded files
cat ASSEMBLY/url.tsv |
    tsv-join -f ASSEMBLY/check.list -k 1 -e |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -d ASSEMBLY/{1} ]; then
            if [ -d ~/data/Bacteria/ASSEMBLY/{1} ]; then
                echo >&2 "==> {1}"
                cp -r ~/data/Bacteria/ASSEMBLY/{1} ASSEMBLY/{1}
            fi
        fi
    '

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
    ~/data/Pseudomonas/ \
    wangq@202.119.37.251:data/Pseudomonas

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Pseudomonas/ \
    wangq@58.213.64.36:data/Pseudomonas

# rsync -avP wangq@202.119.37.251:data/Pseudomonas/ ~/data/Pseudomonas

# rsync -avP -e "ssh -T -c chacha20-poly1305@openssh.com -o Compression=no -x" \
#   wangq@202.119.37.251:data/Pseudomonas/ASSEMBLY/ ~/data/Pseudomonas/ASSEMBLY

```

## BioSample

```shell
mkdir -p ~/data/Pseudomonas/biosample
cd ~/data/Pseudomonas

ulimit -n `ulimit -Hn`

# Copy downloaded files
cat ASSEMBLY/collect.csv |
    tsv-select -H -d, -f BioSample |
    grep "^SAM" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -f biosample/{}.txt ]; then
            if [ -f ~/data/Bacteria/biosample/{}.txt ]; then
                echo >&2 "==> {}"
                cp ~/data/Bacteria/biosample/{}.txt biosample/
            fi
        fi
    '

# Download from NCBI
cat ASSEMBLY/collect.csv |
    tsv-select -H -d, -f BioSample |
    grep "^SAM" |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -s biosample/{}.txt ]; then
            >&2 echo {}
            curl -fsSL "https://www.ncbi.nlm.nih.gov/biosample/?term={}&report=full&format=text" -o biosample/{}.txt
        fi
    '

find biosample -name "SAM*.txt" | wc -l
# 2997

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
    > ASSEMBLY/attributes.lst

cat attributes.lst |
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

## Count and group strains

### Check N50 of assemblies

* Some strains were anomalously labeled and identified by the `mash` tree.
    * Pseudom_flu_GCF_900636635_1
    * Pseudom_chl_GCF_001023535_1
    * Pseudom_syr_GCF_004006335_1
    * Pseudom_puti_GCF_003228315_1 and Pseudom_puti_GCF_020172705_1

```shell
cd ~/data/Pseudomonas

for dir in $(find ASSEMBLY -maxdepth 1 -mindepth 1 -type d | sort); do
    1>&2 echo "==> ${dir}"
    name=$(basename ${dir})

    find ${dir} -type f -name "*_genomic.fna.gz" |
        grep -v "_from_" | # exclude CDS and rna
        xargs cat |
        faops n50 -C -S stdin |
        (echo -e "name\t${name}" && cat) |
        datamash transpose
done |
    tsv-uniq |
    tee ASSEMBLY/n50.tsv

cat ASSEMBLY/n50.tsv |
    tsv-filter \
        -H --or \
        --le 4:100 \
        --ge 2:100000 |
    tsv-filter -H --ge 3:1000000 |
    tr "\t" "," \
    > ASSEMBLY/n50.pass.csv

wc -l ASSEMBLY/n50*
#  2753 ASSEMBLY/n50.pass.csv
#  3001 ASSEMBLY/n50.tsv

tsv-join \
    ASSEMBLY/collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_900636635 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_001023535 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_004006335 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_003228315 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_020172705 \
    > ASSEMBLY/pass.csv

wc -l ASSEMBLY/*csv
#   3001 ASSEMBLY/collect.csv
#   2753 ASSEMBLY/n50.pass.csv
#   2748 ASSEMBLY/pass.csv

```

### Order

```shell
cd ~/data/Pseudomonas

# Group by order
cat ASSEMBLY/pass.csv |
    sed -e '1d' |
    tsv-select -d, -f 3 |
    tsv-uniq |
    nwr append stdin -r order |
    tsv-select -f 2 |
    tsv-uniq \
    > summary/order.lst

cat summary/order.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(cat ASSEMBLY/pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r order -r species |
            grep {} |
            tsv-select -f 1,3 |
            tsv-uniq |
            wc -l)

        n_strains=$(cat ASSEMBLY/pass.csv |
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
    (echo -e '#tax_id\torder\t#species\t#strains' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | order             | #species | #strains |
|---------|-------------------|----------|----------|
| 135624  | Aeromonadales     | 10       | 10       |
| 135622  | Alteromonadales   | 54       | 133      |
| 1385    | Bacillales        | 3        | 3        |
| 213849  | Campylobacterales | 1        | 1        |
| 135615  | Cardiobacteriales | 1        | 1        |
| 204458  | Caulobacterales   | 1        | 1        |
| 1706369 | Cellvibrionales   | 2        | 4        |
| 51291   | Chlamydiales      | 1        | 1        |
| 85007   | Corynebacteriales | 1        | 1        |
| 91347   | Enterobacterales  | 129      | 129      |
| 118969  | Legionellales     | 10       | 10       |
| 135618  | Methylococcales   | 2        | 2        |
| 2887326 | Moraxellales      | 59       | 865      |
| 135619  | Oceanospirillales | 14       | 28       |
| 1240482 | Orbales           | 1        | 1        |
| 135625  | Pasteurellales    | 21       | 21       |
| 72274   | Pseudomonadales   | 247      | 1461     |
| 72273   | Thiotrichales     | 12       | 12       |
| 135623  | Vibrionales       | 35       | 35       |
| 135614  | Xanthomonadales   | 28       | 28       |

### Genus

```shell
cd ~/data/Pseudomonas

# Group by order
cat ASSEMBLY/pass.csv |
    sed -e '1d' |
    tsv-select -d, -f 3 |
    tsv-uniq |
    nwr append stdin -r genus |
    tsv-select -f 2 |
    tsv-uniq \
    > summary/genus.lst

cat summary/genus.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(cat ASSEMBLY/pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r genus -r species |
            grep {} |
            tsv-select -f 1,3 |
            tsv-uniq |
            wc -l)

        n_strains=$(cat ASSEMBLY/pass.csv |
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
    tsv-filter --ge 4:10 |
    (echo -e '#tax_id\tgenus\t#species\t#strains' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | genus             | #species | #strains |
|---------|-------------------|----------|----------|
| 469     | Acinetobacter     | 51       | 791      |
| 226     | Alteromonas       | 16       | 38       |
| 544     | Citrobacter       | 10       | 10       |
| 2745    | Halomonas         | 7        | 19       |
| 570     | Klebsiella        | 10       | 10       |
| 475     | Moraxella         | 6        | 72       |
| 122277  | Pectobacterium    | 10       | 10       |
| 53246   | Pseudoalteromonas | 19       | 33       |
| 286     | Pseudomonas       | 227      | 1405     |
| 613     | Serratia          | 10       | 10       |
| 22      | Shewanella        | 16       | 58       |
| 2901164 | Stutzerimonas     | 9        | 37       |
| 662     | Vibrio            | 29       | 29       |
| 338     | Xanthomonas       | 15       | 15       |

### Strains

```shell
cd ~/data/Pseudomonas

# list strains
mkdir -p taxon

rm taxon/* strains.lst *.tmp
cat ASSEMBLY/pass.csv |
    sed -e '1d' |
    tr "," "\t" |
    tsv-select -f 1,2,3 |
    nwr append stdin -c 3 -r species -r "species group" -r genus -r family -r order |
    parallel --col-sep "\t" --no-run-if-empty --linebuffer -k -j 4 '
        echo {1} >> strains.lst

        echo {4} >> species.tmp
        echo {5} >> species_group.tmp
        echo {6} >> genus.tmp
        echo {7} >> family.tmp

        echo {8} >> order.tmp
        echo {1} >> taxon/{7}

        printf "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n" {1} {2} {3} {4} {5} {6} {7} {8}
    ' |
    sed 's/Pseudomonas paraeruginosa/Pseudomonas aeruginosa/g' \
    > strains.taxon.tsv

cat species.tmp | tsv-uniq > species.lst
cat species_group.tmp | tsv-uniq > species_group.lst
cat genus.tmp | tsv-uniq > genus.lst
cat family.tmp | tsv-uniq > family.lst
cat order.tmp | tsv-uniq > order.lst

# Omit strains without protein annotations
for STRAIN in $(cat strains.lst); do
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_protein.faa.gz" > /dev/null; then
        echo ${STRAIN}
    fi
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_cds_from_genomic.fna.gz" > /dev/null; then
        echo ${STRAIN}
    fi
done |
    tsv-uniq \
    > omit.lst
# All OK

rm *.tmp

```

## NCBI taxonomy

Done by `bp_taxonomy2tree.pl` from BioPerl.

```shell
mkdir -p ~/data/Pseudomonas/tree
cd ~/data/Pseudomonas/tree

bp_taxonomy2tree.pl -e \
    $(
        cat ../genus.lst |
            tr " " "_" |
            parallel echo '-s {}'
    ) \
    > ncbi.nwk

nw_display -s -b 'visibility:hidden' -w 600 -v 30 ncbi.nwk |
    rsvg-convert -o Pseudomonas.ncbi.png

```

## Raw phylogenetic tree by MinHash

```shell
mkdir -p ~/data/Pseudomonas/mash
cd ~/data/Pseudomonas/mash

for strain in $(cat ../strains.lst ); do
    2>&1 echo "==> ${strain}"

    if [[ -e ${strain}.msh ]]; then
        continue
    fi

    find ../ASSEMBLY/${strain} -name "*_genomic.fna.gz" |
        grep -v "_from_" |
        xargs cat |
        mash sketch -k 21 -s 100000 -p 8 - -I "${strain}" -o ${strain}
done

mash triangle -E -p 8 -l <(
    cat ../strains.lst | parallel echo "{}.msh"
    ) \
    > dist.tsv

# fill matrix with lower triangle
tsv-select -f 1-3 dist.tsv |
    (tsv-select -f 2,1,3 dist.tsv && cat) |
    (
        cut -f 1 dist.tsv |
            tsv-uniq |
            parallel -j 1 --keep-order 'echo -e "{}\t{}\t0"' &&
        cat
    ) \
    > dist_full.tsv

cat dist_full.tsv |
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
mkdir -p ~/data/Pseudomonas/tree
cd ~/data/Pseudomonas/tree

nw_reroot ../mash/tree.nwk Baci_subti_subtilis_168 Sta_aure_aureus_NCTC_8325 |
    nw_order -c n - \
    > mash.reroot.newick

# rank::col
ARRAY=(
    'order::8'
    'family::7'
    'genus::6'
    'species_group::5'
    'species::4'
)

rm mash.condensed.map
CUR_TREE=mash.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick mash.${GROUP_NAME}.newick
    cat condense.map >> mash.condensed.map

    CUR_TREE=mash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 20 mash.species.newick |
    rsvg-convert -o Pseudomonas.mash.png

```

## Collect proteins

### `all.pro.fa`

```shell script
cd ~/data/Pseudomonas

mkdir -p PROTEINS

find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    wc -l
# 3000

find ASSEMBLY -type f -name "*_protein.faa.gz" |
    wc -l
# 3000

cat strains.lst |
    wc -l
# 2747

for STRAIN in $(cat strains.lst); do
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
    wc -l
#12881967

gzip -dcf PROTEINS/all.pro.fa.gz |
    grep "^>" |
    tsv-uniq |
    wc -l
#4524852

# annotations may be different
gzip -dcf PROTEINS/all.uniq.fa.gz |
    grep "^>" |
    wc -l
#4436511

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
cd ~/data/Pseudomonas

rm PROTEINS/all.strain.tsv
for STRAIN in $(cat strains.lst); do
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

    faops replace -s ASSEMBLY/${STRAIN}/*_protein.faa.gz <(cut -f 1,2 PROTEINS/${STRAIN}.replace.tsv) stdout

    rm PROTEINS/${STRAIN}.replace.tsv
done |
    pigz -p4 \
    > PROTEINS/all.replace.fa.gz

gzip -dcf PROTEINS/all.replace.fa.gz |
    grep "^>" |
    wc -l
#12881967

(echo -e "#name\tstrain" && cat PROTEINS/all.strain.tsv)  \
    > temp &&
    mv temp PROTEINS/all.strain.tsv

faops size PROTEINS/all.replace.fa.gz > PROTEINS/all.replace.sizes

(echo -e "#name\tsize" && cat PROTEINS/all.replace.sizes) > PROTEINS/all.size.tsv

rm PROTEINS/all.replace.sizes

```

### `all.info.tsv`

```shell
cd ~/data/Pseudomonas

for STRAIN in $(cat strains.lst); do
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
    wc -l
#12881967

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
    wc -l
#12881968

```

## Phylogenetics with bac120

### Find corresponding proteins by `hmmsearch`

* Download HMM models as described in [`hmm/README.md`](../hmm/README.md)

* The `E_VALUE` was manually adjusted to 1e-20 to reach a balance between sensitivity and
  speciality.

```shell
E_VALUE=1e-20

cd ~/data/Pseudomonas

# Find all genes
for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    >&2 echo "==> marker [${marker}]"

    mkdir -p PROTEINS/${marker}

    for FAMILY in $(cat family.lst); do
        >&2 echo "==> FAMILY [${FAMILY}]"

        cat taxon/${FAMILY} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/HMM/bac120/HMM/${marker}.HMM - |
                    grep '>>' |
                    perl -nl -e ' m{>>\s+(\S+)} and printf qq{%s\t%s\n}, \$1, {}; '
            " \
            > PROTEINS/${marker}/${FAMILY}.replace.tsv
    done

    echo
done

```

### Align and concat marker genes to create species tree

```shell script
cd ~/data/Pseudomonas

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        for FAMILY in $(cat family.lst); do
            cat PROTEINS/{}/${FAMILY}.replace.tsv
        done |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
# 2741.75 2745    4019.75

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        for FAMILY in $(cat family.lst); do
            cat PROTEINS/{}/${FAMILY}.replace.tsv
        done |
            wc -l
    ' |
    paste - - |
    tsv-filter --invert --ge 2:2000 --le 2:3500 |
    cut -f 1 \
    > PROTEINS/bac120.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process
cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    grep -v -Fx -f PROTEINS/bac120.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        for FAMILY in $(cat family.lst); do
            cat PROTEINS/{}/${FAMILY}.replace.tsv
        done \
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

    # sequences in one line
    faops filter -l 0 PROTEINS/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > PROTEINS/bac120.aln.fas

fasops concat PROTEINS/bac120.aln.fas strains.lst -o PROTEINS/bac120.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in PROTEINS/bac120.aln.fa -out PROTEINS/bac120.trim.fa -automated1

faops size PROTEINS/bac120.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#49772
#25586

# To make it faster
FastTree -fastest -noml PROTEINS/bac120.trim.fa > PROTEINS/bac120.trim.newick

```

### Tweak the concat tree

```shell script
cd ~/data/Pseudomonas/tree

nw_reroot ../PROTEINS/bac120.trim.newick Baci_subti_subtilis_168 Sta_aure_aureus_NCTC_8325 |
    nw_order -c n - \
    > bac120.reroot.newick

rm bac120.condensed.map

# rank::col
ARRAY=(
#    'order::8'
#    'family::7'
    'genus::6'
    'species_group::5'
    'species::4'
)

rm bac120.condensed.map
CUR_TREE=bac120.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick bac120.${GROUP_NAME}.newick
    cat condense.map >> bac120.condensed.map

    CUR_TREE=bac120.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 20 bac120.species.newick |
    rsvg-convert -o Pseudomonas.bac120.png

```

## InterProScan on all proteins of typical strains

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

Pseudomonas fluorescens SBW25 was removed from refseq

```shell
cd ~/data/Pseudomonas

faops size ASSEMBLY/Pseudom_aeruginosa_PAO1/*_protein.faa.gz |
    wc -l
#5572

faops size ASSEMBLY/Pseudom_aeruginosa_PAO1/*_protein.faa.gz |
    tsv-summarize --sum 2
#1858983

mkdir -p STRAINS

#    Pseudom_fluo_SBW25_GCF_000009225_2 \
for S in \
    Pseudom_aeruginosa_PAO1 \
    Pseudom_putida_KT2440_GCF_000007565_2 \
    Pseudom_chl_aureofaciens_30_84_GCF_000281915_1 \
    Pseudom_entomophila_L48_GCF_000026105_1 \
    Pseudom_proteg_Pf_5_GCF_000012265_1 \
    Pseudom_savas_pv_phaseolicola_1448A_GCF_000012205_1 \
    Stu_stut_A1501_GCF_000013785_1 \
    Pseudom_syr_pv_syringae_B728a_GCF_000012245_1 \
    Pseudom_aeruginosa_UCBPP_PA14_GCF_000014625_1 \
    Pseudom_aeruginosa_PA7_GCF_000017205_1 \
    Pseudom_aeruginosa_PAK_GCF_000568855_2 \
    Pseudom_aeruginosa_LESB58_GCF_000026645_1 \
    ; do
    echo ${S}
done \
    > typical.lst

for S in $(cat typical.lst); do
    mkdir -p STRAINS/${S}
    faops split-about ASSEMBLY/${S}/*_protein.faa.gz 200000 STRAINS/${S}/
done

for S in $(cat typical.lst); do
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

for S in $(cat typical.lst); do
    for f in $(find STRAINS/${S}/ -maxdepth 1 -type f -name "[0-9]*.fa" | sort); do
        >&2 echo "==> ${f}"
        if [ -e ${f}.tsv ]; then
            >&2 echo ${f}
            continue
        fi

        interproscan.sh --cpu 12 -dp -f tsv,json -i ${f} --output-file-base ${f}
    done
done

find STRAINS -type f -name "*.json" | sort |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        if [ $(({#} % 10)) -eq "0" ]; then
            >&2 printf "."
        fi
        pigz -p 3 {}
    '

# same protein may have multiple families
for S in $(cat typical.lst); do
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

COUNT=
for S in $(cat typical.lst); do
    if [ ! -s STRAINS/${S}/family.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family.tsv |
        tsv-summarize -g 2,3 --count \
        > STRAINS/${S}/family-count.tsv

    COUNT=$((COUNT + 1))
done
echo $COUNT

# families in all strains
for S in $(cat typical.lst); do
    cat STRAINS/${S}/family-count.tsv
done |
    tsv-summarize -g 1,2 --count |
    tsv-filter -H --istr-not-in-fld 2:"probable" |
    tsv-filter -H --istr-not-in-fld 2:"putative" |
    tsv-filter -H --istr-not-in-fld 2:"Uncharacterised" |
    tsv-filter -H --istr-not-in-fld 2:" DUF" |
    tsv-filter --ge 3:$COUNT \
    > STRAINS/universal.tsv

# All other strains should have only 1 family member
cp STRAINS/universal.tsv STRAINS/family-1.tsv
for S in $(cat typical.lst | grep -v "_aeru_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/family-1.tsv |
        tsv-filter --eq 3:1 \
        > STRAINS/family-tmp.tsv

    mv STRAINS/family-tmp.tsv STRAINS/family-1.tsv
done

# All P_aeru strains should have multiple family members
cp STRAINS/family-1.tsv STRAINS/family-n.tsv
for S in $(cat typical.lst | grep "_aeru_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/family-n.tsv |
        tsv-filter --gt 3:1 \
        > STRAINS/family-tmp.tsv

    wc -l < STRAINS/family-tmp.tsv
    mv STRAINS/family-tmp.tsv STRAINS/family-n.tsv
done

wc -l STRAINS/Pseudom_aeru_PAO1/family.tsv STRAINS/universal.tsv STRAINS/family-1.tsv STRAINS/family-n.tsv
#  4084 STRAINS/Pseudom_aeru_PAO1/family.tsv
#  1567 STRAINS/universal.tsv
#   972 STRAINS/family-1.tsv
#    14 STRAINS/family-n.tsv

cat STRAINS/family-n.tsv |
    tsv-select -f 1,2 |
    (echo -e "#family\tcount" && cat) |
    mlr --itsv --omd cat

```
