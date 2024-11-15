
```shell
nwr member Stenotrophomonas Serratia Burkholderia Bordetella |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    rgr md stdin --num

nwr member Stenotrophomonas Burkholderia Bordetella -r "species group" -r "species subgroup" |
    tsv-select -f 1-3 |
    keep-header -- tsv-sort -k3,3 -k2,2 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    sed 's/Stenotrophomonas /St. /g' |
    sed 's/Serratia /Se. /g' |
    sed 's/Burkholderia /Bu. /g' |
    sed 's/Bordetella /Bo. /g' |
    rgr md stdin

```

* Order Xanthomonadales
    * Stenotrophomonas maltophilia

* Order Enterobacterales
    * Serratia marcescens

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
