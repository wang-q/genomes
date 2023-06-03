#!/usr/bin/env bash

#----------------------------#
# USAGE
#----------------------------#
USAGE="
Usage: $0 TSV_FILE MIN_COUNT

Default values:
    MIN_COUNT   1

TSV headers:
* strain
* tax_id
* species
* genus
* family
* order

$ bash ~/Scripts/genomes/bin/taxon_count.sh summary/strains.taxon.tsv 2

"

if [ "$#" -lt 1 ]; then
    echo >&2 "$USAGE"
    exit 1
fi

# Set default parameters
TSV_FILE=$1
MIN_COUNT=${2:-1}

#----------------------------#
# Run
#----------------------------#
cat ${TSV_FILE} |
    tsv-select -f 5,4,3 |
    tsv-sort -k1,1 -k2,2 -k3,3 |
    tsv-summarize -g 1,2,3 --count |
    tsv-filter --ge "4:${MIN_COUNT}" |
    perl -nla -F'\t' -e '
            BEGIN { our $family = q(); our $genus = q(); }

            # record the current family
            if ($F[0] eq $family) {
                printf qq(\t);
            } else {
                $family = $F[0];
                printf qq($family\t);
            }
            # record the current genus
            if ($F[1] eq $genus) {
                printf qq(\t);
            } else {
                $genus = $F[1];
                printf qq($genus\t);
            }

            print join qq(\t), ($F[2], $F[3]);
        ' |
    (echo -e '#family\tgenus\tspecies\tcount' && cat) |
    mlr --itsv --omd cat |
    sed 's/-\s*|$/-:|/'
