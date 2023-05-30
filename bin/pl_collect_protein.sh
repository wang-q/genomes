#!/usr/bin/env bash

#----------------------------#
# Colors in term
#----------------------------#
GREEN=
RED=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    NC='\033[0m' # No Color
fi

log_warn () {
    echo >&2 -e "${RED}==> $@ <==${NC}"
}

log_info () {
    echo >&2 -e "${GREEN}==> $@${NC}"
}

log_debug () {
    echo >&2 -e "==> $@"
}

#----------------------------#
# USAGE
#----------------------------#
USAGE="
This is a script that belongs to the pipeline for collecting proteins.

Expecting dirs/files:
    ASSEMBLY/
    summary/strains.lst

Creating files:
    PROTEINS/all.all.fa.gz
    PROTEINS/all.uniq.fa.gz

$ cd ~/data/Fungi
$ bash ~/Scripts/genome/bin/pl_collect_protein.sh

"

# Check dirs/files
if [[ ! -d "ASSEMBLY/" ]]; then
    echo >&2 "ASSEMBLY/ doesn't exist"
    echo >&2 "$USAGE"
    exit 1
fi
if [[ ! -s "summary/strains.lst" ]]; then
    echo >&2 "summary/strains.lst doesn't exist or is empty"
    echo >&2 "$USAGE"
    exit 1
fi

#----------------------------#
# all.pro.fa
#----------------------------#
log_info "all.pro.fa"

mkdir -p PROTEINS

log_debug "PROTEINS/all.pro.fa.gz"
cat summary/strains.lst | grep -v -Fw -f ASSEMBLY/omit.lst |
while read STRAIN; do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz
done |
    pigz -p4 \
    > PROTEINS/all.pro.fa.gz

log_debug "PROTEINS/all.uniq.fa.gz"
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

#----------------------------#
# all.replace.fa
#----------------------------#
log_info "all.replace.fa"

rm -f PROTEINS/all.strain.tsv PROTEINS/all.replace.fa.gz

log_debug "PROTEINS/all.replace.fa.gz"
cat summary/strains.lst | grep -v -Fw -f ASSEMBLY/omit.lst |
while read STRAIN; do
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

log_debug "PROTEINS/PROTEINS/all.size.tsv"
(echo -e "#name\tstrain" && cat PROTEINS/all.strain.tsv)  \
    > temp &&
    mv temp PROTEINS/all.strain.tsv

faops size PROTEINS/all.replace.fa.gz > PROTEINS/all.replace.sizes

(echo -e "#name\tsize" && cat PROTEINS/all.replace.sizes) > PROTEINS/all.size.tsv

rm PROTEINS/all.replace.sizes

#----------------------------#
# `all.info.tsv`
#----------------------------#
log_info "all.info.tsv"

log_debug "PROTEINS/PROTEINS/all.annotation.tsv"
cat summary/strains.lst | grep -v -Fw -f ASSEMBLY/omit.lst |
while read STRAIN; do
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

(echo -e "#name\tannotation" && cat PROTEINS/all.annotation.tsv) \
    > temp &&
    mv temp PROTEINS/all.annotation.tsv

# check differences
#cat PROTEINS/all.size.tsv |
#    grep -v -F -f <(cut -f 1 PROTEINS/all.annotation.tsv)

log_debug "PROTEINS/PROTEINS/all.info.tsv"
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

#----------------------------#
# Counts
#----------------------------#
log_info "Counts"

printf "#item\tcount\n" \
    > PROTEINS/counts.tsv

printf "Proteins\t%s\n" $(
    gzip -dcf PROTEINS/all.pro.fa.gz |
        grep "^>" |
        wc -l |
        numfmt --to=si
    ) \
    >> PROTEINS/counts.tsv

printf "Unique headers and annotations\t%s\n" $(
    gzip -dcf PROTEINS/all.pro.fa.gz |
        grep "^>" |
        tsv-uniq |
        wc -l |
        numfmt --to=si
    ) \
    >> PROTEINS/counts.tsv

printf "Unique proteins\t%s\n" $(
    gzip -dcf PROTEINS/all.uniq.fa.gz |
        grep "^>" |
        wc -l |
        numfmt --to=si
    ) \
    >> PROTEINS/counts.tsv

printf "Unique proteins\t%s\n" $(
    gzip -dcf PROTEINS/all.uniq.fa.gz |
        grep "^>" |
        wc -l |
        numfmt --to=si
    ) \
    >> PROTEINS/counts.tsv

printf "all.replace.fa\t%s\n" $(
    gzip -dcf PROTEINS/all.replace.fa.gz |
        grep "^>" |
        wc -l |
        numfmt --to=si
    ) \
    >> PROTEINS/counts.tsv

printf "all.replace.fa\t%s\n" $(
    gzip -dcf PROTEINS/all.replace.fa.gz |
        grep "^>" |
        wc -l |
        numfmt --to=si
    ) \
    >> PROTEINS/counts.tsv

printf "all.annotation.tsv\t%s\n" $(
    cat PROTEINS/all.annotation.tsv |
        wc -l |
        numfmt --to=si
    ) \
    >> PROTEINS/counts.tsv

printf "all.info.tsv\t%s\n" $(
    cat PROTEINS/all.annotation.tsv |
        wc -l |
        numfmt --to=si
    ) \
    >> PROTEINS/counts.tsv
