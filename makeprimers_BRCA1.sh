#!/bin/sh
primerway.py -o BRCA1_primers -n BRCA1 -p NP_009231 \
            -R GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
            -G GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz \
            -V 00-common_all.vcf.gz -c /opt/primer-way/primerway.cfg