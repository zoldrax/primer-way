# PrimerWay
PrimerWay is a pipeline, which would facilitate the process of large-scale primer design for the analysis of entire coding sequence of the gene in an automatic mode. It obtains the sequence of interest with flanking regions using protein_id or genome coordinates from the reference genome. It uses Primer3 to design set of primer pair candidates. The choice is optimized to avoid frequent variants in the primer sequence extracted from the common subset of the dbSNP database. ThermonucleotideBLAST is used for the check of pair specificity. As result, the optimal set of overlapped pairs is generated.

## Installation
PrimerWay is a simple Python script which requires installed [Primer3-py][1], [ThermonucleotideBLAST][2], and [HTSlib][3].
#### Get PrimerWay
```
wget http://github.com/zoldrax/primer-way/raw/master/primerway.py
chmod a+x primerway.py
```
#### Get and Install [Primer3-py][1]
```
pip install primer3-py
```
#### Get and Install [ThermonucleotideBLAST][2]
```
wget http://public.lanl.gov/jgans/tntblast/tntblast-2.04.tar.gz
tar zvfx tntblast-2.04.tar.gz
cd tntblast-2.04
./configure --enable-mpi
make
sudo make install
```
#### Get and Install [HTSlib][3]
```
wget http://github.com/samtools/htslib/archive/1.6.tar.gz
tar zvfx 1.6.tar.gz
cd htslib-1.6
make
sudo make install
```


## Getting started
Design primer pairs for entire coding sequence of the BRCA1 gene:
```
./primerway.py -o BRCA1_primers -n BRCA1 -p NP_009231 \
-R GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
-G GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz \
-V common_all_20170403_GRCh38p7.vcf.gz
```

## Acknowledgements
The work has been supported by RSF grant [#16-45-02011][9].

[1]: http://github.com/libnano/primer3-py
[2]: http://public.lanl.gov/jgans/tntblast
[3]: http://www.htslib.org
[9]: http://rscf.ru/en/enprjcard?rid=16-45-02011
