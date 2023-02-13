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
- [Phylogenetics with bac120](#phylogenetics-with-bac120)
    * [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
    * [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
    * [Tweak the concat tree](#tweak-the-concat-tree)
- [InterProScan on all proteins of typical strains](#interproscan-on-all-proteins-of-typical-strains)
    * [Typical strains](#typical-strains)
    * [`interproscan.sh`](#interproscansh)
    * [P. aeruginosa](#p-aeruginosa)
    * [P. putida](#p-putida)
    * [P. syringae](#p-syringae)
    * [A. baumannii](#a-baumannii)
- [Protein families](#protein-families)
    * [IPR007416 - YggL 50S ribosome-binding protein](#ipr007416---yggl-50s-ribosome-binding-protein)

<!-- tocstop -->

## Strain info

* [Pseudomonas](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=286)
* [Acinetobacter](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=469)

According to a recent [paper](https://doi.org/10.1128/mSystems.00543-20), there are some order-level
changes in Gammaproteobacteria. We include both old and new orders.

* Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
* New ones: Moraxellales, Kangiellales, and Pseudomonadales

![msystems.00543-20-f0003.jpeg](..%2Fimages%2Fmsystems.00543-20-f0003.jpeg)

Another recent [paper](https://doi.org/10.1099/ijsem.0.005542) proposed that an outlier clade of P.
aeruginosa (containing PA7) split into a new species, Pseudomonas paraeruginosa.

Similarly, a recent [paper](https://doi.org/10.1016/j.syapm.2021.126289) suggested that some species
of Pseudomonas should be transferred to Halopseudomonas and Stutzerimonas.

### List all ranks

```shell
mkdir -p ~/data/Pseudomonas
cd ~/data/Pseudomonas

nwr member Pseudomonas Acinetobacter Stenotrophomonas Burkholderia Bordetella |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

nwr member Pseudomonas Acinetobacter Stenotrophomonas Burkholderia Bordetella -r "species group" -r "species subgroup" |
    tsv-select -f 1-3 |
    keep-header -- tsv-sort -k3,3 -k2,2 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    sed 's/Stenotrophomonas /S. /g' |
    sed 's/Burkholderia /Bu. /g' |
    sed 's/Bordetella /Bo. /g' |
    mlr --itsv --omd cat

```

| rank             | count |
|------------------|-------|
| genus            | 5     |
| species          | 655   |
| strain           | 2323  |
| subspecies       | 15    |
| no rank          | 141   |
| species group    | 10    |
| species subgroup | 7     |
| isolate          | 5     |

| #tax_id | sci_name                                 | rank             |
|---------|------------------------------------------|------------------|
| 2839056 | A. Taxon 24                              | species group    |
| 909768  | A. calcoaceticus/baumannii complex       | species group    |
| 87882   | Bu. cepacia complex                      | species group    |
| 136841  | P. aeruginosa group                      | species group    |
| 136842  | P. chlororaphis group                    | species group    |
| 136843  | P. fluorescens group                     | species group    |
| 136845  | P. putida group                          | species group    |
| 136849  | P. syringae group                        | species group    |
| 995085  | S. maltophilia group                     | species group    |
| 111527  | pseudomallei group                       | species group    |
| 2839060 | A. Taxon 24C                             | species subgroup |
| 2839057 | A. Taxon 24D                             | species subgroup |
| 2839061 | A. Taxon 24E                             | species subgroup |
| 627141  | P. nitroreducens/multiresinivorans group | species subgroup |
| 1232139 | P. oleovorans/pseudoalcaligenes group    | species subgroup |
| 251695  | P. syringae group genomosp. 1            | species subgroup |
| 251698  | P. syringae group genomosp. 2            | species subgroup |

### Species with assemblies

* Order Pseudomonadales
    * Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
    * New ones: Moraxellales, Kangiellales, and Pseudomonadales
    * Pseudomonas aeruginosa
    * Acinetobacter baumannii

* Order Xanthomonadales
    * Stenotrophomonas maltophilia

* Order Burkholderiales: It is classified as Betaproteobacteria, but is close Gammaproteobacteria
    * Burkholderia cepacia
    * Bordetella pertussis

* Order Legionellales
    * Coxiella burnetii

* Order Enterobacterales
    * Escherichia coli
    * Klebsiella pneumoniae
    * Salmonella enterica
    * Shigella flexneri

```shell
mkdir -p ~/data/Pseudomonas/summary
cd ~/data/Pseudomonas/summary

SPECIES=$(
    nwr member \
        Cellvibrionales Oceanospirillales Alteromonadales \
        Moraxellales Kangiellales Pseudomonadales \
        Xanthomonadales \
        Burkholderiales \
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
    tsv-filter -H --ge CHR:20 |
    tsv-filter -H --invert --str-in-fld species:Pseudomonas --lt RS:50 |
    tsv-filter -H --invert --str-in-fld species:Acinetobacter --lt RS:50 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    sed 's/Stenotrophomonas /S. /g' |
    sed 's/Burkholderia /Bu. /g' |
    sed 's/Bordetella /Bo. /g' |
    mlr --itsv --omd cat

```

| species_id | species                    | RS   | CHR |
|------------|----------------------------|------|-----|
| 287        | P. aeruginosa              | 7526 | 648 |
| 470        | A. baumannii               | 7313 | 531 |
| 28450      | Bu. pseudomallei           | 1782 | 154 |
| 520        | Bo. pertussis              | 878  | 593 |
| 40324      | S. maltophilia             | 731  | 89  |
| 317        | P. syringae                | 615  | 51  |
| 95486      | Bu. cenocepacia            | 522  | 71  |
| 87883      | Bu. multivorans            | 477  | 79  |
| 347        | Xanthomonas oryzae         | 413  | 140 |
| 48296      | A. pittii                  | 385  | 35  |
| 305        | Ralstonia solanacearum     | 305  | 123 |
| 294        | P. fluorescens             | 262  | 40  |
| 292        | Bu. cepacia                | 241  | 20  |
| 303        | P. putida                  | 226  | 70  |
| 346        | Xanthomonas citri          | 217  | 111 |
| 339        | Xanthomonas campestris     | 171  | 54  |
| 85698      | Achromobacter xylosoxidans | 147  | 23  |
| 56448      | Xanthomonas arboricola     | 146  | 20  |
| 2371       | Xylella fastidiosa         | 142  | 47  |
| 316        | Stutzerimonas stutzeri     | 140  | 27  |
| 38313      | Shewanella algae           | 112  | 24  |
| 587753     | P. chlororaphis            | 107  | 61  |
| 756892     | A. indicus                 | 102  | 21  |
| 13373      | Bu. mallei                 | 99   | 34  |
| 518        | Bo. bronchiseptica         | 97   | 25  |
| 519        | Bo. parapertussis          | 94   | 90  |
| 343        | Xanthomonas translucens    | 94   | 33  |
| 380021     | P. protegens               | 83   | 27  |
| 35814      | Bo. holmesii               | 80   | 66  |
| 57975      | Bu. thailandensis          | 73   | 25  |
| 337        | Bu. glumae                 | 67   | 51  |
| 1530123    | A. seifertii               | 60   | 25  |
| 164546     | Cupriavidus taiwanensis    | 43   | 38  |
| 476        | Moraxella bovis            | 39   | 37  |
| 1444770    | Xylella taiwanensis        | 32   | 22  |

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

# Also includes representative strains of Proteobacteria.
# families not in our orders
FAMILY=$(
    nwr member Proteobacteria -r family |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        sed '1d' |
        cut -f 1 |
        nwr append stdin -r order |
        tsv-filter --str-ne 2:"Cellvibrionales" --str-ne 2:"Oceanospirillales" --str-ne 2:"Alteromonadales" |
        tsv-filter --str-ne 2:"Moraxellales" --str-ne 2:"Kangiellales" --str-ne 2:"Pseudomonadales" |
        tsv-filter --str-ne 2:"Xanthomonadales" |
        tsv-filter --str-ne 2:"Burkholderiales" |
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
        AND organism_name NOT LIKE '%symbiont %'
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

# Check md5
# rm ASSEMBLY/check.list
bash ASSEMBLY/check.sh

# Run
proxychains4 bash ASSEMBLY/rsync.sh

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
# 6212

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
    >> summary/biosample.tsv

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
#  5789 ASSEMBLY/n50.pass.csv
#  6217 ASSEMBLY/n50.tsv

# Omit strains without protein annotations
for STRAIN in $(cat ASSEMBLY/n50.pass.csv | cut -d, -f 1); do
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_protein.faa.gz" > /dev/null; then
        echo ${STRAIN}
    fi
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_cds_from_genomic.fna.gz" > /dev/null; then
        echo ${STRAIN}
    fi
done |
    tsv-uniq
# All OK

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
#   6217 ASSEMBLY/collect.csv
#   5789 ASSEMBLY/n50.pass.csv
#   5784 ASSEMBLY/pass.csv

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
    tsv-uniq |
    sort \
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

| #tax_id | order               | #species | #strains |
|---------|---------------------|----------|----------|
| 225057  | Acidithiobacillales | 4        | 4        |
| 135624  | Aeromonadales       | 10       | 10       |
| 135622  | Alteromonadales     | 54       | 133      |
| 1385    | Bacillales          | 3        | 3        |
| 2024979 | Bacteriovoracales   | 2        | 2        |
| 213481  | Bdellovibrionales   | 2        | 2        |
| 1779134 | Bradymonadales      | 1        | 1        |
| 80840   | Burkholderiales     | 231      | 2007     |
| 213849  | Campylobacterales   | 40       | 40       |
| 135615  | Cardiobacteriales   | 1        | 1        |
| 204458  | Caulobacterales     | 7        | 7        |
| 1706369 | Cellvibrionales     | 2        | 4        |
| 51291   | Chlamydiales        | 1        | 1        |
| 85007   | Corynebacteriales   | 1        | 1        |
| 213118  | Desulfobacterales   | 2        | 2        |
| 213115  | Desulfovibrionales  | 4        | 4        |
| 69541   | Desulfuromonadales  | 5        | 5        |
| 91347   | Enterobacterales    | 131      | 131      |
| 356     | Hyphomicrobiales    | 64       | 64       |
| 118969  | Legionellales       | 10       | 10       |
| 135618  | Methylococcales     | 2        | 2        |
| 2887326 | Moraxellales        | 60       | 869      |
| 29      | Myxococcales        | 6        | 6        |
| 206351  | Neisseriales        | 21       | 21       |
| 32003   | Nitrosomonadales    | 2        | 2        |
| 135619  | Oceanospirillales   | 16       | 32       |
| 1240482 | Orbales             | 1        | 1        |
| 135625  | Pasteurellales      | 21       | 21       |
| 72274   | Pseudomonadales     | 254      | 1484     |
| 204455  | Rhodobacterales     | 21       | 21       |
| 206389  | Rhodocyclales       | 5        | 5        |
| 204441  | Rhodospirillales    | 19       | 19       |
| 766     | Rickettsiales       | 16       | 16       |
| 204457  | Sphingomonadales    | 12       | 12       |
| 72273   | Thiotrichales       | 12       | 12       |
| 135623  | Vibrionales         | 35       | 35       |
| 135614  | Xanthomonadales     | 141      | 792      |

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
    tsv-uniq |
    sort \
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
| 222     | Achromobacter     | 9        | 89       |
| 12916   | Acidovorax        | 5        | 28       |
| 469     | Acinetobacter     | 52       | 795      |
| 507     | Alcaligenes       | 2        | 23       |
| 226     | Alteromonas       | 16       | 38       |
| 517     | Bordetella        | 23       | 807      |
| 374     | Bradyrhizobium    | 13       | 13       |
| 32008   | Burkholderia      | 103      | 692      |
| 194     | Campylobacter     | 27       | 27       |
| 544     | Citrobacter       | 10       | 10       |
| 283     | Comamonas         | 5        | 19       |
| 106589  | Cupriavidus       | 14       | 70       |
| 80865   | Delftia           | 4        | 14       |
| 2745    | Halomonas         | 7        | 19       |
| 963     | Herbaspirillum    | 7        | 12       |
| 570     | Klebsiella        | 10       | 10       |
| 68      | Lysobacter        | 5        | 18       |
| 475     | Moraxella         | 6        | 72       |
| 846     | Oxalobacter       | 3        | 12       |
| 93217   | Pandoraea         | 5        | 16       |
| 1822464 | Paraburkholderia  | 10       | 24       |
| 122277  | Pectobacterium    | 10       | 10       |
| 44013   | Polynucleobacter  | 4        | 22       |
| 53246   | Pseudoalteromonas | 19       | 33       |
| 286     | Pseudomonas       | 234      | 1429     |
| 48736   | Ralstonia         | 17       | 141      |
| 613     | Serratia          | 10       | 10       |
| 22      | Shewanella        | 16       | 58       |
| 40323   | Stenotrophomonas  | 11       | 249      |
| 2901164 | Stutzerimonas     | 9        | 36       |
| 662     | Vibrio            | 29       | 29       |
| 338     | Xanthomonas       | 108      | 441      |
| 2370    | Xylella           | 13       | 69       |

### Strains

```shell
cd ~/data/Pseudomonas

# list strains
mkdir -p taxon

rm taxon/* summary/strains.lst *.tmp
cat ASSEMBLY/pass.csv |
    sed -e '1d' |
    tr "," "\t" |
    tsv-select -f 1,2,3 |
    nwr append stdin -c 3 -r species -r genus -r family -r order |
    parallel --col-sep "\t" --no-run-if-empty --linebuffer -k -j 4 '
        echo {1} >> summary/strains.lst

        echo {1} >> taxon/{7}

        printf "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n" {1} {2} {3} {4} {5} {6} {7}
    ' |
    sed 's/Pseudomonas paraeruginosa/Pseudomonas aeruginosa/g' \
    > strains.taxon.tsv

```

## NCBI taxonomy

Done by `bp_taxonomy2tree.pl` from BioPerl.

```shell
mkdir -p ~/data/Pseudomonas/tree
cd ~/data/Pseudomonas/tree

bp_taxonomy2tree.pl -e \
    $(
        cat ../summary/genus.lst |
            tr " " "_" |
            parallel echo '-s {}'
    ) \
    > ncbi.nwk

nw_display -s -b 'visibility:hidden' -w 1200 -v 20 ncbi.nwk |
    rsvg-convert -o Pseudomonas.ncbi.png

```

## Raw phylogenetic tree by MinHash

```shell
mkdir -p ~/data/Pseudomonas/tree
cd ~/data/Pseudomonas/tree

mash triangle -E -p 8 -l <(
    cat ../strains.taxon.tsv |
        cut -f 1 |
        (echo -e "Baci_subti_subtilis_168" && cat) |
        (echo -e "Sta_aure_aureus_NCTC_8325" && cat) |
        parallel --no-run-if-empty --linebuffer -k -j 1 '
            if [[ -e ../../Bacteria/mash/{}.msh ]]; then
                echo "../../Bacteria/mash/{}.msh"
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
mkdir -p ~/data/Pseudomonas/tree
cd ~/data/Pseudomonas/tree

nw_reroot ../mash/tree.nwk Baci_subti_subtilis_168 Sta_aure_aureus_NCTC_8325 |
    nw_order -c n - \
    > mash.reroot.newick

# rank::col
ARRAY=(
    'order::7'
    'family::6'
    'genus::5'
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
# 3034

find ASSEMBLY -type f -name "*_protein.faa.gz" |
    wc -l
# 3034

cat summary/strains.lst |
    wc -l
# 2781

for STRAIN in $(cat summary/strains.lst); do
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
#13058330

gzip -dcf PROTEINS/all.pro.fa.gz |
    grep "^>" |
    tsv-uniq |
    wc -l
#4567759

# annotations may be different
gzip -dcf PROTEINS/all.uniq.fa.gz |
    grep "^>" |
    wc -l
#4475736

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

rm PROTEINS/all.strain.tsv PROTEINS/all.replace.fa.gz
for STRAIN in $(cat summary/strains.lst); do
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
    wc -l
#13058330

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

for STRAIN in $(cat summary/strains.lst); do
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
#13058330

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
#13058331

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

    for ORDER in $(cat summary/order.lst); do
        >&2 echo "==> ORDER [${ORDER}]"

        cat taxon/${ORDER} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/HMM/bac120/HMM/${marker}.HMM - |
                    grep '>>' |
                    perl -nl -e ' m{>>\s+(\S+)} and printf qq{%s\t%s\n}, \$1, {}; '
            " \
            > PROTEINS/${marker}/${ORDER}.replace.tsv
    done

    echo
done

```

### Align and concat marker genes to create species tree

```shell script
cd ~/data/Pseudomonas

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        for ORDER in $(cat summary/order.lst); do
            cat PROTEINS/{}/${ORDER}.replace.tsv
        done |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
# 2776    2779    4068

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        for ORDER in $(cat summary/order.lst); do
            cat PROTEINS/{}/${ORDER}.replace.tsv
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

        for ORDER in $(cat summary/order.lst); do
            cat PROTEINS/{}/${ORDER}.replace.tsv
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

fasops concat PROTEINS/bac120.aln.fas summary/strains.lst -o PROTEINS/bac120.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in PROTEINS/bac120.aln.fa -out PROTEINS/bac120.trim.fa -automated1

faops size PROTEINS/bac120.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#50203
#25591

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
#    'order::7'
#    'family::6'
#    'genus::5'
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
nw_display -s -b 'visibility:hidden' -w 800 -v 20 bac120.species.newick |
    rsvg-convert -o Pseudomonas.bac120.png

```

## InterProScan on all proteins of typical strains

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

faops size ASSEMBLY/Pseudom_aeruginosa_PAO1/*_protein.faa.gz |
    wc -l
#5572

faops size ASSEMBLY/Pseudom_aeruginosa_PAO1/*_protein.faa.gz |
    tsv-summarize --sum 2
#1858983

cat summary/species.count.tsv |
    tsv-filter -H --or --ge RS:50 --ge CHR:10 |
    tsv-filter -H --or \
        --str-in-fld species:Pseudomonas \
        --str-in-fld species:Halopseudomonas \
        --str-in-fld species:Stutzerimonas \
        --str-in-fld species:Azotobacter \
        --str-in-fld species:Acinetobacter |
    mlr --itsv --omd cat

for ABBR in \
    Pseudom_aeruginosa_ \
    Pseudom_viridif_ \
    Pseudom_syringae_ \
    Pseudom_fluo_ \
    Pseudom_putida_ \
    Stu_stut_ \
    Pseudom_savas_ \
    Pseudom_chl_ \
    Pseudom_amyg_ \
    Pseudom_proteg_ \
    Pseudom_fulva_ \
    Pseudom_coro_ \
    Pseudom_synx_ \
    Pseudom_entomophila_ \
    Acin_bau_ \
    Acin_pit_ \
    Acin_nos_ \
    Acin_indicus_ \
    Acin_junii_ \
    Acin_urs_ \
    Acin_bere_ \
    Acin_sei_ \
    Acin_radior_ \
    Acin_haemolyticus_ \
    Acin_johnsonii_ \
    Acin_lwo_ \
    ; do
    cat ASSEMBLY/pass.csv |
        grep ${ABBR} |
        grep -E "Reference|Representative"

done |
    cut -d, -f 1

mkdir -p STRAINS

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
    Pseudom_syringae_pv_tagetis_GCF_022557255_1 \
    Pseudom_syringae_GCF_018394375_1 \
    Pseudom_fluo_GCF_900215245_1 \
    Pseudom_putida_KT2440_GCF_000007565_2 \
    Pseudom_putida_NBRC_14164_GCF_000412675_1 \
    Stu_stut_A1501_GCF_000013785_1 \
    Pseudom_savas_pv_phaseolicola_1448A_GCF_000012205_1 \
    Pseudom_savas_GCF_020917325_1 \
    Pseudom_chl_aureofaciens_30_84_GCF_000281915_1 \
    Pseudom_chl_GCF_014524625_1 \
    Pseudom_amyg_pv_tabaci_ATCC_11528_GCF_000145945_2 \
    Pseudom_proteg_Pf_5_GCF_000012265_1 \
    Pseudom_proteg_CHA0_GCF_900560965_1 \
    Pseudom_coro_pv_oryzae_1_6_GCF_000156995_2 \
    Pseudom_synx_GCF_003851555_1 \
    Pseudom_entomophila_L48_GCF_000026105_1 \
    Acin_bau_GCF_008632635_1 \
    Acin_bau_ATCC_17978_GCF_004794235_2 \
    Acin_bau_ATCC_19606_CIP_70_34_JCM_6841_GCF_019331655_1 \
    Acin_bau_GCF_016919505_2 \
    Acin_pit_PHEA_2 \
    Acin_nos_M2_GCF_005281455_1 \
    Acin_indicus_GCF_009914475_1 \
    Acin_junii_GCF_018336855_1 \
    Acin_bere_GCF_016576965_1 \
    Acin_sei_GCF_016064815_1 \
    Acin_radior_GCF_003258335_1 \
    Acin_haemolyticus_GCF_003323815_1 \
    Acin_lwo_GCF_019343495_1 \
    ; do
    echo ${S}
done \
    > summary/typical.lst

```

### `interproscan.sh`

```shell
cd ~/data/Pseudomonas

for S in $(cat summary/typical.lst); do
    mkdir -p STRAINS/${S}
    faops split-about ASSEMBLY/${S}/*_protein.faa.gz 200000 STRAINS/${S}/
done

for S in $(cat summary/typical.lst); do
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
for S in $(cat summary/typical.lst); do
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

# Pseudom
COUNT=
for S in $(cat summary/typical.lst | grep -E "^(Pseudom|Stu)"); do
    if [ ! -s STRAINS/${S}/family.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family.tsv |
        tsv-summarize -g 2,3 --count \
        > STRAINS/${S}/family-count.tsv

    COUNT=$((COUNT + 1))
done
echo $COUNT

# families in almost all strains
for S in $(cat summary/typical.lst | grep -E "^(Pseudom|Stu)"); do
    cat STRAINS/${S}/family-count.tsv
done |
    tsv-summarize -g 1,2 --count |
    tsv-filter -H --istr-not-in-fld 2:"probable" |
    tsv-filter -H --istr-not-in-fld 2:"putative" |
    tsv-filter -H --istr-not-in-fld 2:"Uncharacterised" |
    tsv-filter -H --istr-not-in-fld 2:" DUF" |
    tsv-filter --ge 3:$((COUNT - 2)) \
    > STRAINS/Pseudom-universal.tsv

# Acin
COUNT=
for S in $(cat summary/typical.lst | grep -E "^(Acin)"); do
    if [ ! -s STRAINS/${S}/family.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family.tsv |
        tsv-summarize -g 2,3 --count \
        > STRAINS/${S}/family-count.tsv

    COUNT=$((COUNT + 1))
done
echo $COUNT

# families in almost all strains
for S in $(cat summary/typical.lst | grep -E "^(Acin)"); do
    cat STRAINS/${S}/family-count.tsv
done |
    tsv-summarize -g 1,2 --count |
    tsv-filter -H --istr-not-in-fld 2:"probable" |
    tsv-filter -H --istr-not-in-fld 2:"putative" |
    tsv-filter -H --istr-not-in-fld 2:"Uncharacterised" |
    tsv-filter -H --istr-not-in-fld 2:" DUF" |
    tsv-filter --ge 3:$((COUNT - 1)) \
    > STRAINS/Acin-universal.tsv

```

### P. aeruginosa

```shell
PREFIX=Pseudom_aeruginosa

cd ~/data/Pseudomonas

# All other strains should have only 1 family member
cp STRAINS/Pseudom-universal.tsv STRAINS/${PREFIX}-1.tsv
for S in $(cat summary/typical.lst | grep -E "^(Pseudom|Stu)" | grep -v "${PREFIX}_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/${PREFIX}-1.tsv |
        tsv-filter --le 3:1 \
        > STRAINS/family-tmp.tsv

    mv STRAINS/family-tmp.tsv STRAINS/${PREFIX}-1.tsv
done

# All ${PREFIX} strains should have multiple family members
cp STRAINS/${PREFIX}-1.tsv STRAINS/${PREFIX}-n.tsv
for S in $(cat summary/typical.lst | grep -E "^(Pseudom|Stu)" | grep "${PREFIX}_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/${PREFIX}-n.tsv |
        tsv-filter --gt 3:1 \
        > STRAINS/family-tmp.tsv

    wc -l < STRAINS/family-tmp.tsv
    mv STRAINS/family-tmp.tsv STRAINS/${PREFIX}-n.tsv
done

wc -l STRAINS/Pseudom_aeruginosa_PAO1/family.tsv \
    STRAINS/Pseudom-universal.tsv STRAINS/${PREFIX}-1.tsv STRAINS/${PREFIX}-n.tsv
#  3971 STRAINS/Pseudom_aeruginosa_PAO1/family.tsv
#  1707 STRAINS/Pseudom-universal.tsv
#   912 STRAINS/Pseudom_aeruginosa-1.tsv
#     9 STRAINS/Pseudom_aeruginosa-n.tsv

cat STRAINS/${PREFIX}-n.tsv |
    tsv-select -f 1,2 |
    tsv-sort |
    (echo -e "#family\tcount" && cat) |
    mlr --itsv --omd cat

```

| #family   | count                                         |
|-----------|-----------------------------------------------|
| IPR001353 | Proteasome, subunit alpha/beta                |
| IPR001404 | Heat shock protein Hsp90 family               |
| IPR004361 | Glyoxalase I                                  |
| IPR005999 | Glycerol kinase                               |
| IPR006684 | Acyl-CoA thioester hydrolase YbgC/YbaW family |
| IPR007416 | YggL 50S ribosome-binding protein             |
| IPR011757 | Lytic transglycosylase MltB                   |
| IPR014311 | Guanine deaminase                             |
| IPR037532 | Peptidoglycan D,D-transpeptidase FtsI         |

### P. putida

```shell
PREFIX=Pseudom_putida

cd ~/data/Pseudomonas

# All other strains should have only 1 family member
cp STRAINS/Pseudom-universal.tsv STRAINS/${PREFIX}-1.tsv
for S in $(cat summary/typical.lst | grep -E "^(Pseudom|Stu)" | grep -v "${PREFIX}_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/${PREFIX}-1.tsv |
        tsv-filter --le 3:1 \
        > STRAINS/family-tmp.tsv

    mv STRAINS/family-tmp.tsv STRAINS/${PREFIX}-1.tsv
done

# All ${PREFIX} strains should have multiple family members
cp STRAINS/${PREFIX}-1.tsv STRAINS/${PREFIX}-n.tsv
for S in $(cat summary/typical.lst | grep -E "^(Pseudom|Stu)" | grep "${PREFIX}_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/${PREFIX}-n.tsv |
        tsv-filter --gt 3:1 \
        > STRAINS/family-tmp.tsv

    wc -l < STRAINS/family-tmp.tsv
    mv STRAINS/family-tmp.tsv STRAINS/${PREFIX}-n.tsv
done

wc -l STRAINS/Pseudom_putida_KT2440_GCF_000007565_2/family.tsv \
    STRAINS/Pseudom-universal.tsv STRAINS/${PREFIX}-1.tsv STRAINS/${PREFIX}-n.tsv
#  3731 STRAINS/Pseudom_putida_KT2440_GCF_000007565_2/family.tsv
#  1707 STRAINS/Pseudom-universal.tsv
#   860 STRAINS/Pseudom_putida-1.tsv
#     7 STRAINS/Pseudom_putida-n.tsv

cat STRAINS/${PREFIX}-n.tsv |
    tsv-select -f 1,2 |
    tsv-sort |
    (echo -e "#family\tcount" && cat) |
    mlr --itsv --omd cat

```

| #family   | count                                                      |
|-----------|------------------------------------------------------------|
| IPR001783 | Lumazine-binding protein                                   |
| IPR002033 | Sec-independent periplasmic protein translocase TatC       |
| IPR002446 | Lipocalin, bacterial                                       |
| IPR011342 | Shikimate dehydrogenase                                    |
| IPR012794 | Beta-ketoadipate transcriptional regulator, PcaR/PcaU/PobR |
| IPR018448 | Sec-independent protein translocase protein TatB           |
| IPR018550 | Lipid A 3-O-deacylase-related                              |

### P. protegens

```shell
PREFIX=Pseudom_proteg

cd ~/data/Pseudomonas

# All other strains should have only 1 family member
cp STRAINS/Pseudom-universal.tsv STRAINS/${PREFIX}-1.tsv
for S in $(cat summary/typical.lst | grep -E "^(Pseudom|Stu)" | grep -v "${PREFIX}_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/${PREFIX}-1.tsv |
        tsv-filter --le 3:1 \
        > STRAINS/family-tmp.tsv

    mv STRAINS/family-tmp.tsv STRAINS/${PREFIX}-1.tsv
done

# All ${PREFIX} strains should have multiple family members
cp STRAINS/${PREFIX}-1.tsv STRAINS/${PREFIX}-n.tsv
for S in $(cat summary/typical.lst | grep -E "^(Pseudom|Stu)" | grep "${PREFIX}_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/${PREFIX}-n.tsv |
        tsv-filter --gt 3:1 \
        > STRAINS/family-tmp.tsv

    wc -l < STRAINS/family-tmp.tsv
    mv STRAINS/family-tmp.tsv STRAINS/${PREFIX}-n.tsv
done

wc -l STRAINS/Pseudom_proteg_Pf_5_GCF_000012265_1/family.tsv \
    STRAINS/Pseudom-universal.tsv STRAINS/${PREFIX}-1.tsv STRAINS/${PREFIX}-n.tsv
#  4104 STRAINS/Pseudom_proteg_Pf_5_GCF_000012265_1/family.tsv
#  1707 STRAINS/Pseudom-universal.tsv
#   853 STRAINS/Pseudom_proteg-1.tsv
#     6 STRAINS/Pseudom_proteg-n.tsv

cat STRAINS/${PREFIX}-n.tsv |
    tsv-select -f 1,2 |
    tsv-sort |
    (echo -e "#family\tcount" && cat) |
    mlr --itsv --omd cat

```

| #family   | count                                    |
|-----------|------------------------------------------|
| IPR000926 | GTP cyclohydrolase II, RibA              |
| IPR004794 | Riboflavin biosynthesis protein RibD     |
| IPR009246 | Ethanolamine ammonia-lyase small subunit |
| IPR010099 | Epimerase family protein SDR39U1         |
| IPR010628 | Ethanolamine ammonia-lyase heavy chain   |
| IPR011842 | Coenzyme PQQ biosynthesis protein B      |

### P. syringae

### A. baumannii

```shell
PREFIX=Acin_bau

cd ~/data/Pseudomonas

# All other strains should have only 1 family member
cp STRAINS/Acin-universal.tsv STRAINS/${PREFIX}-1.tsv
for S in $(cat summary/typical.lst | grep -E "^(Acin)" | grep -v "${PREFIX}_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/${PREFIX}-1.tsv |
        tsv-filter --le 3:1 \
        > STRAINS/family-tmp.tsv

    mv STRAINS/family-tmp.tsv STRAINS/${PREFIX}-1.tsv
done

# All ${PREFIX} strains should have multiple family members
cp STRAINS/${PREFIX}-1.tsv STRAINS/${PREFIX}-n.tsv
for S in $(cat summary/typical.lst | grep -E "^(Acin)" | grep "${PREFIX}_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/${PREFIX}-n.tsv |
        tsv-filter --gt 3:1 \
        > STRAINS/family-tmp.tsv

    wc -l < STRAINS/family-tmp.tsv
    mv STRAINS/family-tmp.tsv STRAINS/${PREFIX}-n.tsv
done

wc -l STRAINS/Acin_bau_GCF_008632635_1/family.tsv \
    STRAINS/Acin-universal.tsv STRAINS/${PREFIX}-1.tsv STRAINS/${PREFIX}-n.tsv
#  2525 STRAINS/Acin_bau_GCF_008632635_1/family.tsv
#  1342 STRAINS/Acin-universal.tsv
#   890 STRAINS/Acin_bau-1.tsv
#     0 STRAINS/Acin_bau-n.tsv

cat STRAINS/${PREFIX}-n.tsv |
    tsv-select -f 1,2 |
    tsv-sort |
    (echo -e "#family\tcount" && cat) |
    mlr --itsv --omd cat

```

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
