# PDAN

<!-- toc -->

- [Strain info](#strain-info)
    * [Bacillota](#bacillota)
    * [Terrabacteria group](#terrabacteria-group)
    * [Pseudomonadota](#pseudomonadota)
    * [FCB group](#fcb-group)
    * [The rest](#the-rest)
    * [Extracting](#extracting)

<!-- tocstop -->

## Strain info

### Bacillota

Bacillota == Firmicutes

Mycoplasmatota == Tenericutes

Clostridia belongs to Bacillota

```shell
mkdir -p ~/data/Bacteria/Bacillota
cd ~/data/Bacteria/Bacillota

rm -f ASSEMBLY
rm -f Protein
ln -s ../ASSEMBLY ASSEMBLY
ln -s ../Protein Protein

mkdir -p summary

FAMILY=$(
    nwr member Bacillota Mycoplasmatota -r family |
        sed '1d' |
        cut -f 1
)
#FAMILY=$(IFS=, ; echo "${FAMILY[*]}")

```

### Terrabacteria group

Actinomycetota == Actinobacteria

Deinococcota == Deinococcus-Thermus

Cyanobacteriota == Cyanobacteria

Chloroflexi belongs to Chloroflexota

```shell
mkdir -p ~/data/Bacteria/Terrabacteria
cd ~/data/Bacteria/Terrabacteria

rm -f ASSEMBLY
rm -f Protein
ln -s ../ASSEMBLY ASSEMBLY
ln -s ../Protein Protein

mkdir -p summary

FAMILY=$(
    nwr member "Terrabacteria group" -r family |
        sed '1d' |
        nwr restrict Bacillota --exclude -f stdin -c 1 |
        nwr restrict Mycoplasmatota --exclude -f stdin -c 1 |
        cut -f 1
)

```

### Pseudomonadota

Pseudomonadota == Proteobacteria

```shell
mkdir -p ~/data/Bacteria/Pseudomonadota
cd ~/data/Bacteria/Pseudomonadota

rm -f ASSEMBLY
rm -f Protein
ln -s ../ASSEMBLY ASSEMBLY
ln -s ../Protein Protein

mkdir -p summary

FAMILY=$(
    nwr member Pseudomonadota -r family |
        sed '1d' |
        cut -f 1
)

```

### FCB group

```shell
mkdir -p ~/data/Bacteria/FCB
cd ~/data/Bacteria/FCB

rm -f ASSEMBLY
rm -f Protein
ln -s ../ASSEMBLY ASSEMBLY
ln -s ../Protein Protein

mkdir -p summary

FAMILY=$(
    nwr member "FCB group" -r family |
        sed '1d' |
        cut -f 1
)

```

### The rest

```shell
mkdir -p ~/data/Bacteria/TheRest
cd ~/data/Bacteria/TheRest

rm -f ASSEMBLY
rm -f Protein
ln -s ../ASSEMBLY ASSEMBLY
ln -s ../Protein Protein

mkdir -p summary

FAMILY=$(
    nwr member Bacteria -r family |
        sed '1d' |
        nwr restrict "Terrabacteria group" --exclude -f stdin -c 1 |
        nwr restrict Pseudomonadota --exclude -f stdin -c 1 |
        nwr restrict "FCB group" --exclude -f stdin -c 1 |
        cut -f 1
)

```

### Extracting

```shell
# cd ~/data/Bacteria/Bacillota

cat ../summary/collect.pass.tsv |
    head -n 1 \
    > summary/collect.pass.tsv

cat ../summary/collect.pass.tsv |
    sed '1d' | # 307775
    nwr restrict ${FAMILY[*]} -f stdin -c 3 | # restrict to these families 167899
    tsv-join -e -f ../ASSEMBLY/omit.lst -k 1 | # 167899
    sort \
    >> summary/collect.pass.tsv

cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    tsv-join -H -f summary/collect.pass.tsv -k 1 \
    > summary/assembly.tsv

# biosample.tsv
cp ../summary/attributes.lst summary/

cat ../summary/biosample.tsv |
    grep -Fw -f <(cat summary/collect.pass.tsv | tsv-select -H -f BioSample | sort | uniq) \
    > summary/biosample.tsv

```

```shell
cd ~/data/Bacteria

for GROUP in \
    Bacillota \
    Terrabacteria \
    Pseudomonadota \
    FCB \
    TheRest \
    ; do
    find ${GROUP}/summary -type f -name "collect.pass.tsv" | xargs wc -l
done
#90541 Bacillota/summary/collect.pass.tsv
#25210 Terrabacteria/summary/collect.pass.tsv
#167900 Pseudomonadota/summary/collect.pass.tsv
#10880 FCB/summary/collect.pass.tsv
#12845 TheRest/summary/collect.pass.tsv

```

## Collect proteins

```shell
cd ~/data/Bacteria/

ulimit -n `ulimit -Hn`

for GROUP in \
    Bacillota \
    Terrabacteria \
    Pseudomonadota \
    FCB \
    TheRest \
    ; do
    cd ~/data/Bacteria/${GROUP}

    nwr template summary/assembly.tsv \
        --pro \
        --parallel 8 \
        --clust-id 0.95 \
        --clust-cov 0.95

    # collect proteins and clustering
    # It may need to be run several times
    bash Protein/collect.sh Escherichia_albertii

done

fd -g strains.tsv ~/data/Bacteria/Protein --size +100k -l

cd ~/data/Bacteria/
rm -fr Protein/tmp/

for GROUP in \
    Bacillota \
    Terrabacteria \
    Pseudomonadota \
    FCB \
    TheRest \
    ; do
    cd ~/data/Bacteria/${GROUP}

    nwr template summary/assembly.tsv \
        --pro \
        --parallel 8 \
        --clust-id 0.95 \
        --clust-cov 0.95

    # info.tsv
    bash Protein/info.sh

done

cd ~/data/Bacteria/

find Protein -type f -name "pro.fa.gz" -size -4096c

nwr template ~/Scripts/genomes/assembly/Bacteria.assembly.tsv \
    --pro \
    --parallel 16 \
    --in ASSEMBLY/pass.lst \
    --not-in ASSEMBLY/omit.lst

# anything left
bash Protein/collect.sh
bash Protein/info.sh
# counts
bash Protein/count.sh

cat Protein/counts.tsv |
    tsv-summarize -H --count --sum 2-5 |
    sed 's/^count/species/' |
    datamash transpose |
    perl -nla -F"\t" -MNumber::Format -e '
        printf qq(%s\t%s\n), $F[0], Number::Format::format_number($F[1], 0,);
        ' |
    (echo -e "#item\tcount" && cat) |
    mlr --itsv --omd cat

```

| #item      | count       |
|------------|-------------|
| species    | 8,163       |
| strain_sum | 157,786     |
| total_sum  | 575,230,650 |
| dedup_sum  | 127,514,848 |
| rep_sum    | 59,532,317  |

## Rsync to hpcc

```bash
for GROUP in \
    Bacillota \
    Terrabacteria \
    Pseudomonadota \
    FCB \
    TheRest \
    ; do

    rsync -avP \
        ~/data/Bacteria/${GROUP}/ \
        wangq@202.119.37.251:data/Bacteria/${GROUP}

done

rsync -avP \
    ~/data/Bacteria/Protein/ \
    wangq@202.119.37.251:data/Bacteria/Protein

```
