# HMM related resources

<!-- toc -->

- [PFAM-A](#pfam-a)
- [TIGRFAM](#tigrfam)
- [40 single-copy genes](#40-single-copy-genes)
- [120 bacterial proteins `bac120`](#120-bacterial-proteins-bac120)
- [61 fungal marker genes](#61-fungal-marker-genes)

<!-- tocstop -->

## PFAM-A

```shell
mkdir -p ~/data/HMM/PFAM

cd ~/data/HMM/PFAM

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    proxychains4 wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
done

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    echo "==> ${basename}"
    gzip -dcf ${basename}.gz > ${basename}
done

```

## TIGRFAM

```shell
mkdir -p ~/data/HMM/TIGRFAM

cd ~/data/HMM/TIGRFAM

proxychains4 wget -N --content-disposition ftp://ftp.jcvi.org/data/TIGRFAMs/14.0_Release/TIGRFAMs_14.0_HMM.tar.gz

mkdir -p HMM
tar --directory HMM -xzvf TIGRFAMs_14.0_HMM.tar.gz TIGR02013.HMM
tar --directory HMM -xzvf TIGRFAMs_14.0_HMM.tar.gz TIGR00485.HMM

```

## 40 single-copy genes

Ref.:

1. Wu, D., Jospin, G. & Eisen, J. A. Systematic Identification of Gene Families for Use as “Markers”
   for Phylogenetic and Phylogeny-Driven Ecological Studies of Bacteria and Archaea and Their Major
   Subgroups. PLoS ONE 8, e77033 (2013).

2. Skennerton, C. T. et al. Phylogenomic analysis of Candidatus ‘Izimaplasma’ species: free-living
   representatives from a Tenericutes clade found in methane seeps. ISME J. 10, 2679–2692 (2016).

3. https://doi.org/10.6084/m9.figshare.722713.v1

4. `bacteria_and_archaea.tgz`: https://ndownloader.figshare.com/files/3093482

```shell
mkdir -p ~/data/HMM/scg40
cd ~/data/HMM/scg40

curl -LO https://ndownloader.figshare.com/files/3093482

tar xvfz 3093482

```

## 120 bacterial proteins `bac120`

Ref.:

1. Nature Microbiology volume 2, pages 1533–1542 (2017)
    * Save Supplementary Table 6 to `data/bac120.tsv`
    * https://doi.org/10.1038/s41564-017-0012-7
2. Nature Biotechnology volume 36, pages 996–1004 (2018)

```shell
mkdir -p ~/data/HMM/bac120
cd ~/data/HMM/bac120

cp ~/Scripts/genomes/data/bac120.tsv .

cat bac120.tsv |
    sed '1d' |
    tsv-select -f 1 \
    > bac120.lst

mkdir -p hmm

cat bac120.lst |
    grep '^TIGR' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        tar --directory hmm -xzvf ../TIGRFAM/TIGRFAMs_14.0_HMM.tar.gz {}.HMM
    '

cat bac120.lst |
    grep -v '^TIGR' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        curl -L http://pfam.xfam.org/family/{}/hmm > hmm/{}.HMM
    '

GZIP=-9 tar cvfz bac120.tar.gz \
    bac120.lst \
    bac120.tsv \
    hmm/

```

## 53 archaeal proteins `ar53`

Ref.:

1. Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life.
   Nat Microbiol 2, 1533–1542 (2017).
    * Supplementary Table 7 listed 122 marker genes
    * https://doi.org/10.1038/s41564-017-0012-7
2. A standardized archaeal taxonomy for the Genome Taxonomy Database.
   Nat Microbiol 6, 946–959 (2021).
    * https://doi.org/10.1038/s41564-021-00918-8

https://data.gtdb.ecogenomic.org/releases/release214/214.1/auxillary_files/ar53_msa_marker_info_r214.tsv

```shell
mkdir -p ~/data/HMM/ar53
cd ~/data/HMM/ar53

cp ~/Scripts/genomes/data/ar53.tsv .

cat ar53.tsv |
    sed '1d' |
    tsv-select -f 1 \
    > ar53.lst

mkdir -p hmm

cat ar53.lst |
    grep '^TIGR' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        tar --directory hmm -xzvf ../TIGRFAM/TIGRFAMs_14.0_HMM.tar.gz {}.HMM
    '

cat ar53.lst |
    grep -v '^TIGR' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        curl -L https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{}?annotation=hmm |
            gzip -dc \
            > hmm/{}.HMM
    '

GZIP=-9 tar cvfz ar53.tar.gz \
    ar53.lst \
    ar53.tsv \
    hmm/

```

## 61 fungal marker genes

Ref.: https://doi.org/10.1093/nar/gkac894

The link to the HMM file in their website is actually in the augustus prfl format.

We use hmmbuild to build hmm models from sequence alignments.

```shell
mkdir -p ~/data/HMM/fungi61
cd ~/data/HMM/fungi61

mkdir -p fatsa
mkdir -p hmm

curl -L https://ufcg.steineggerlab.com/ufcg/genes > genes.html

cat genes.html |
    pup 'table#genes tr td strong text{}' \
    > fungi61.lst

#https://ufcg.steineggerlab.workers.dev/msa/ACT1_aligned.fasta
cat fungi61.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        if [ -s fatsa/{}.fasta ]; then
            exit
        fi
        curl -L https://ufcg.steineggerlab.workers.dev/msa/{}_aligned.fasta > fatsa/{}.fasta
    '

cat fungi61.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        esl-reformat stockholm fatsa/{}.fasta > fatsa/{}.sto
        hmmbuild hmm/{}.HMM fatsa/{}.sto
    '

GZIP=-9 tar cvfz fungi61.tar.gz \
    fungi61.lst \
    genes.html \
    hmm/

```
