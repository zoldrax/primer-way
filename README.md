# PrimerWay
PrimerWay is a pipeline, which would facilitate the process of large-scale primer design for the analysis of entire coding sequence of the gene in an automatic mode. It obtains the sequence of interest with flanking regions using protein_id or genome coordinates from the reference genome. It uses Primer3 to design set of primer pair candidates. The choice is optimized to avoid frequent variants in the primer sequence extracted from the common subset of the dbSNP database. ThermonucleotideBLAST is used for the check of pair specificity. As result, the optimal set of overlapped pairs is generated.

Please cite us: [Ekaterina Sh Kuligina, Anna P. Sokolenko, Ilya V. Bizin. PCR-based validation of novel cancer-predisposing mutations identified upon exome sequencing: An emphasis on primer design [abstract]. In: Proceedings of the American Association for Cancer Research Annual Meeting 2018; 2018 Apr 14-18; Chicago, IL. Philadelphia (PA): AACR; Cancer Res 2018;78(13 Suppl):Abstract nr 2283.][9]

## Installation
PrimerWay is a simple Python script which requires installed [Primer3-py][1], [ThermonucleotideBLAST][2], [pysam][3], and [PyVCF][4].

You can easily install PrimerWay with [docker][8]:
```
mkdir data
docker run -it --rm -v $PWD/data:/data zoldrax/primer-way download_GRCh38.sh
docker run -it --rm -v $PWD/data:/data zoldrax/primer-way protein_id.sh BRCA1 NP_009231
```
The download is needed only once, and about 5 Gb of free disk space is needed.
"BRCA1" is the prefix you may choose.

or do subsequent steps without docker instead:

#### Get PrimerWay
```
wget http://github.com/zoldrax/primer-way/raw/master/primerway.py
chmod +x primerway.py
```

#### Get and Install [Primer3-py][1] [pysam][3] and [PyVCF][4]
```
pip install primer3-py pysam pyvcf
```

#### Get and Install [ThermonucleotideBLAST][2]
```
wget http://public.lanl.gov/jgans/tntblast/tntblast-2.04.tar.gz
tar zvfx tntblast-2.04.tar.gz
cd tntblast-2.04
./configure --enable-openmp
make CXXFLAGS='-std=c++03 -O1'
sudo make install
cd ..
```

#### Get Reference, Annotation, and Variants from [The Genome Reference Consortium][4] and [dbSNP][5]
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
gzip -d GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/00-common_all.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/00-common_all.vcf.gz.tbi
```

## Getting started
Design primer pairs for entire coding sequence of the BRCA1 gene:
```
./primerway.py -o BRCA1_primers -n BRCA1 -p NP_009231 \
-R GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
-G GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz \
-V 00-common_all.vcf.gz
```

## Acknowledgements
The work has been supported by RSF grant [#16-45-02011][7].

[1]: http://github.com/libnano/primer3-py
[2]: http://public.lanl.gov/jgans/tntblast
[3]: http://github.com/pysam-developers/pysam
[4]: http://github.com/jamescasbon/PyVCF
[5]: http://www.ncbi.nlm.nih.gov/grc/human
[6]: http://www.ncbi.nlm.nih.gov/SNP
[7]: http://rscf.ru/en/enprjcard?rid=16-45-02011
[8]: http://hub.docker.com/r/zoldrax/primer-way
[9]: https://cancerres.aacrjournals.org/content/78/13_Supplement/2283
