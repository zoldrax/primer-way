#!/bin/sh
python /opt/primer-way/primerway.py -o /data/BRCA1_primers -n BRCA1 -p NP_009231 \
            -R /data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
            -G /data/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz \
            -V /data/00-common_all.vcf.gz -c /opt/primer-way/primerway.cfg