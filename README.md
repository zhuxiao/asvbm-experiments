# SV_STAT experiment
In our experiments, we used the high accuracy of the PacBio sequencing platform in CCS mode, combined with NGMLR's alignment capability for long-read data and 30X the coverage depth. This approach provides a well-balanced strategy for detecting variants. The identification results of HG002 CCS data were benchmarked using ASVCLR(v1.4.1), SVDSS (v2.0.0), Debreak (v2.0.3), Sniffles2 (v2.0.2), pbsv (v2.9.0), cuteSV (v2.0.3), and SVIM(v2.0.0), respectively. The benchmark dataset was the high-confidence HG002 dataset created by the Genome in a Bottle Consortium (GIAB). More specific experimental information was shown as follows.
## Prerequisites

### Tools

We used  [SV_STAT](https://github.com/zhuxiao/sv_stat) to benchmark variant calling results.

```sh
$ wget -c https://github.com/zhuxiao/sv_stat/releases/download/1.0.1/sv_stat_1.0.1.tar.xz
$ tar -xf sv_stat_1.0.1.tar.xz
$ cd sv_stat_1.0.1/
$ ./autogen.sh
```

And the binary file `sv_stat` will be output into the folder `bin` in this package directory.

We used the following detection methods to variant calling. In addition to the SV detection method included in this experiment, we also introduced a new SV detection tool ASVCLR and benchmarked its SV identification results.

```sh
# Get ASVCLR 
$ wget -c https://github.com/zhuxiao/asvclr/releases/download/1.4.1/asvclr_1.4.1.tar.xz
$ tar -xf asvclr_1.4.1.tar.xz
$ cd asvclr_1.4.1/
$ ./auto_gen.sh
# Or get from github
$ git clone https://github.com/zhuxiao/asvclr.git
$ tar -xf asvclr_1.4.1.tar.xz
$ cd asvclr_1.4.1/
$ ./auto_gen.sh
```

And the binary file `asvclr` will be output into the folder `bin` in this package directory.

```sh
# Get SVDSS v2.0.0-alpha.1
$ wget -c https://github.com/Parsoa/SVDSS/archive/refs/tags/v2.0.0-alpha.1.tar.gz
$ tar -zxvf SVDSS-2.0.0-alpha.1.tar.gz
$ cd SVDSS
$ mkdir build ; cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make
# Get DeBreak cuteSV pbsv sniffles svim and samtools
$ conda install debreak=2.0.3 cuteSV=2.0.3 pbsv=2.9.0 sniffles2=2.0.2 svim=2.0.0 samtools  
# We need ngmlr v0.2.7 to align fasta or fastq with reference
$ wget https://github.com/philres/ngmlr/releases/download/v0.2.7/ngmlr-0.2.7-linux-x86_64.tar.gz
$ tar xvzf ngmlr-0.2.7-linux-x86_64.tar.gz
$ cd ngmlr-0.2.7/
$ mkdir build ; cd build
$ cmake ..
$ make
# We also need sratoolkit to download PacBio CCS data.
# centOS
$ wget -c https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.10/sratoolkit.3.0.10-centos_linux64.tar.gz
$ tar -zxvf sratoolkit.3.0.10-centos_linux64.tar.gz
$ cd sratoolkit.3.0.10-centos/bin/
$ ln -s prefetch /usr/local/bin/prefetch
$ ln -s fastq-dump /usr/local/bin/fastq-dump
$ ln -s fasterq-dump /usr/local/bin/fasterq-dump
# ubuntu
$ wget -c https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.10/sratoolkit.3.0.10-ubuntu64.tar.gz
$ tar -zxvf sratoolkit.3.0.10-ubuntu64.tar.gz
$ cd sratoolkit.3.0.10-ubuntu64
$ ln -s prefetch /usr/local/bin/prefetch
$ ln -s fastq-dump /usr/local/bin/fastq-dump
$ ln -s fasterq-dump /usr/local/bin/fasterq-dump
```

And the binary file `prefetch`、 `fastq-dump`  and `fasterq-dump` will be output into the folder `bin` in this package directory.

### Data

In our experimental benchmarking, we used hg19.

#### Download reference

```sh
# Reference
$ wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
$ gunzip hs37d5.fa.gz
# Extract chromosomes from 1 to 22 and X and Y
$ samtools faidx hs37d5.fa 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y > hs37d5.chroms.fa
```

## HG002

Download [HG002 PacBio CCS](https://www.ncbi.nlm.nih.gov/sra/SRX5327410) data.  After aligning with ngmlr, convert them into bam files using samtools, sort and create index. For convenience, we provide a shell script and a list of accessions to help you obtain the fastq file (see `script` folder). Significantly, you need to ensure that the file and the script are in the same folder.

```sh
$ ./prefetch_fastq.sh SRR_Acc_List.txt SRR885_whole.fastq
$ ngmlr -t 12 --rg-id na24385_pb_ccs -r hs37d5.chroms.fa -q SRR885_whole.fastq -o HG002_pacbio_ccs.sam
$ samtools view -bSh -@ 12 HG002_pacbio_ccs.sam > HG002_pacbio_ccs.bam
$ samtools sort -@ 12 -o HG002_pacbio_ccs_sorted.bam HG002_pacbio_ccs.bam
$ samtools index -@ 12 HG002_pacbio_ccs_sorted.bam HG002_pacbio_ccs_sorted.bai
# remove fastq to save storage space
$ rm -rf HG002_pacbio_ccs.sam HG002_pacbio_ccs.bam
```

## Variant calling

```sh
# ASVCLR
$ asvclr all -m 20 -p hg002_hs37d5_ngmlr_asvclr -o output_debug_hs37d5_q0_20240408 hs37d5.fa ngmlr_HG002_hs37d5_pacbio_ccs_sorted.bam
$ mv hg002_hs37d5_ngmlr_asvclr_variants.vcf output_ASVCLR.vcf
```

More  detailed usage of ASVCLR can be obtained from Github ([ASVCLR](https://github.com/zhuxiao/asvclr)).

You can get variant detection results in folder `4_results` and variant detection results are reported in VCF file format in this file folder: `genome_variants.vcf`.

```sh
#SVDSS
$ SVDSS index --threads 12 --reference hs37d5.fa --index hs37d5.fmd
$ SVDSS smooth --threads 12 --reference hs37d5.fa --bam ngmlr_HG002_hs37d5_pacbio_ccs_ash_sorted.bam > ngmlr_HG002_hs37d5_pacbio_ccs_ash_sorted_smoothed.bam
$ samtools index ngmlr_HG002_hs37d5_pacbio_ccs_ash_sorted_smoothed.bam
$ SVDSS search --threads 12 --index hs37d5.fmd --bam ngmlr_HG002_hs37d5_pacbio_ccs_ash_sorted_smoothed.bam > specifics.txt
$ SVDSS call --threads 12 --reference hs37d5.fa --bam ngmlr_HG002_hs37d5_pacbio_ccs_ash_sorted_smoothed.bam --sfs specifics.txt > output_SVDSS.vcf
# DeBreak
$ debreak --thread 12 --min_size 20 --bam ngmlr_HG002_hs37d5_pacbio_ccs_ash_sorted.bam --outpath output_debreak_hs37d5_ngmlr --rescue_large_ins --poa --ref hs37d5.fa 
$ cd output_debreak_hs37d5_ngmlr && mv debreak.vcf output_DebBeak.vcf
# cuteSV
$ cuteSV -t 12 -l 20 --genotype HG002_pacbio_ccs_sorted.bam hs37d5.chroms.fa output_cuteSV.vcf $PWD
# pbsv
$ pbsv discover -s HG002_30X_CCS HG002_pacbio_ccs_sorted.bam ref.out.svsig.gz
$ pbsv call -m 20 --ccs hs37d5.chroms.fa ref.out.svsig.gz output_pbsv.vcf
# Sniffles2
$ sniffles --minsvlen 20 -i HG002_pacbio_ccs_sorted.bam -v output_Sniffles.vcf
# svim
$ svim alignment --min_sv_size 20 --sample HG002_CCS_30X HG002_pacbio_ccs_sorted.bam hs37d5.chroms.fa
$ mv variants.vcf output_SVIM.vcf
```

You can get the following seven results. Identification results are saved within the `Identification_result` directory:

* **ASVCLR** : `output_ASVCLR.vcf`
* **SVDSS** : `output_SVDSS.vcf`
* **DeBreak** : `output_DeBreak.vcf`
* **cuteSV** : `output_cuteSV.vcf`
* **pbsv** : `output_pbsv.vcf`
* **Sniffles2** : `output_sniffles2.vcf`
* **SVIM** : `output_SVIM.vcf`

## GIAB analysis

In this experiment, we used HG002_SVs_Tier1_v0.6.vcf as the gold benchmark set required for the benchmarking of user-called sets, which can be downloaded as follows:

```sh
# Get GIAB VCF Tier 1 
$ wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
$ gunzip HG002_SVs_Tier1_v0.6.vcf.gz

# Run sv_stat against the Tier1 callset and SV_STAT can benchmark multiple user-called sets simultaneously.
$ sv_stat -m 50000 -T "ASVCLR;SVDSS;DeBreak;Sniffles2;pbsv;cuteSV;SVIM" -C "1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;X;Y" genome_variants.vcf output_SVDSS.vcf output_DeBreak.vcf output_sniffles2.vcf output_pbsv.vcf output_cuteSV.vcf output_SVIM.vcf hs37d5.chroms.fa -o Tier1_eval
```

We used -T option to specify the name of the detection method and -C option to specify the set of chromosomes to be benchmarked.  In this experiment, we benchmarked the autosomes and sex chromosomes X/Y in the SV Identification results of different detection methods.

## Benchmarking results

### Performance benchmarking

Typically, benchmarking results are saved within the `Tier1_eval` directory, with each tool's results saved in the subfolder named after the respective tool.  Additionally, a local HTML file (`sv_stat_reports.html`) is generated to store the benchmarking results for each user-called set. Review the comprehensive benchmarking results more conveniently through the `sv_stat_reports.pdf` file.

The benchmarking results are shown in the table:

|   Tool   |  SVs   | TP_benchmark | TP_user |  FP   |  FN   |  recall  | precision | F1  score |  Identity  |
| :------: | :----: | :----------: | :-----: | :---: | :---: | :------: | :-------: | :-------: | :-------: |
|  ASVCLR  | 54423  |    45986     |  45807  | 8616  | 28026 | 0.621332 | 0.841685  | 0.714914  | 0.973949  |
|  SVDSS   | 45787  |    34689     |  37221  | 8566  | 39328 | 0.468627 | 0.812916  | 0.594524  | 0.977113  |
| DeBreak  | 49868  |    43644     |  41248  | 7565  | 30368 | 0.589688 | 0.845021  | 0.694634  | 0.936123  |
| Sniffles2| 54545  |    44973     |  43106  | 10168 | 29039 | 0.607645 | 0.809138  | 0.694063  | 0.973332  |
|   pbsv   | 52807  |    44492     |  42927  | 9253  | 29520 | 0.601146 | 0.822672  | 0.694676  | 0.983104  |
|  cuteSV  | 44937  |    39438     |  36952  | 6416  | 34574 | 0.532860 | 0.852057  | 0.655674  | 0.975145  |
|   SVIM   | 116615 |    48022     |  47230  | 30995 | 25990 | 0.648841 | 0.603771  | 0.625459  | 0.979929  |

The figure below displays the benchmarking results of different detection methods, including two categories of basic metrics, where Identity represents the sequence identity calculated for matched SVs containing sequences.  Detailed statistics can be found in the corresponding text files within the respective folders.

<div style="text-align: center;">
    <img src="img/evaluation_result.png" alt="Performance comparison between different tools" style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\result_classification.png" alt="Benchmark results between different tools" style="display: inline-block;" width="400"/>
</div>



### Statistical results of deviation of overlapping variants

Moreover, for regions with overlapping variant, the quantities of  region size ratio and center distance were statistically analyzed to  provide a more intuitive presentation of benchmarking information. The statistical results are as follows:

(1) Deviation statistics of center distance

The statistical results of the deviation of center distance are as follows：

|   Tool   | -200~-151 | -150~-101 | -100~-51 | -50~-1 | 0~50  | 51~100 | 101~150 | 151~200 |
| :------: | :-------: | :-------: | :------: | :----: | :---: | :----: | :-----: | :-----: |
|  ASVCLR  |    241    |    358    |   756    |  6408  | 31364 |  4160  |  1924   |  1313   |
|  SVDSS   |    281    |    373    |   902    |  5463  | 27417 |  2188  |  1520   |  1122   |
|  DeBreak |    278    |    451    |   849    |  3435  | 31649 |  3741  |  1546   |   963   |
| Sniffles2|    276    |    355    |   646    |  3154  | 32890 |  4559  |  2006   |  1394   |
|   pbsv   |    288    |    438    |   812    |  5027  | 32552 |  2314  |  1476   |  1071   |
|  cuteSV  |    182    |    276    |   604    |  4750  | 27134 |  3410  |  1541   |  1107   |
|   SVIM   |    553    |    686    |   1159   |  5988  | 35797 |  3855  |  2291   |  1737   |

(2) Deviation statistics of the region size ratio

The statistical results of the deviation of the region size ratio are as follows：

|   Tool   | 0.0~0.5 | 0.5~0.7 | 0.7~1.2 | 1.2~2.0 | 2.0~5.0 | 5.0~10.0 | 10.0~50.0 | 50.0~100.0 | >100.0 |
| :------: | :-----: | :-----: | :-----: | :-----: | :-----: | :------: | :-------: | :--------: | :----: |
|  ASVCLR  |  2885   |   773   |  46869  |   1287   |   1222   |    306    |    248    |    31     |  23  |
|   SVDSS  |  2874   |   737   |  39357  |   1186   |   1330   |    283    |    154    |    6     |  6  |
|  DeBreak |  2559   |   657   |  40170  |   2723   |   1841   |    456    |    292    |    21     |  40  |
| Sniffles2 |  3257   |   752   |  46170  |   1244   |   1236   |    342    |    231    |    31     |  43  |
|   pbsv   |  3500   |   840   |  43660  |   1724   |   1489   |   401    |   321    |    40    |  62  |
| cuteSV |  2916   |   615   |  39079  |   1039   |   934   |   213    |    175    |    16     |  24  |
|   SVIM   |  7184   |   1303   |  51061  |   1730   |   1752   |   607    |   524    |    85     |  173  |


### Benchmarking results for metrics of different SV size regions

Additionally, basic metrics for different structural variant (SV) size ranges were computed, primarily categorized into the following seven intervals. The results are shown as follows:

(1) Benchmarking results for metrics of different SV size regions with different methods

SVs are categorized into seven size regions and metrics are computed for comprehensive benchmarking of different detection methods within each region. The benchmarking results are as follows:

<div style="text-align: center;">
    <img src="img\different_range.png" alt="Evaluation results of different SV size regions" style="display: inline-block; margin: 0 auto; text-align: center;" width="800"/>
</div>


(2) Statistics of the count of different SV lengths in the user-called set (ASVCLR):

|    region    | TP_bench | TP_user |  FP  |  FN   |  recall  | precision | F1 score |  Identity  |
| :----------: | :------: | :-----: | :--: | :---: | :------: | :-------: | :------: | :-------: |
|   1-100bp    |  38704   |  38369  | 7279 | 25440 | 0.603392 | 0.840635  | 0.702525 | 0.972428  |
|  101-250bp   |   2050   |  1889   | 968  | 1908  | 0.517938 | 0.661183  | 0.580860 | 0.992270  |
|  251-500bp   |   1707   |  1669   | 1481 | 1265  | 0.572342 | 0.524953  | 0.547624 | 0.995352  |
|  501-1000bp  |   364    |   358   | 870  |  748  | 0.327338 | 0.291531  | 0.308399 | 0.996147  |
| 1001-2500bp  |   387    |   462   | 312  |  448  | 0.463473 | 0.596899  | 0.521792 | 0.998879  |
| 2501-5000bp  |   187    |   288   | 112  |  216  | 0.464020 | 0.720000  | 0.564339 | 0.999939  |
| 5001-10000bp |   138    |   138   |  116  |  198  | 0.410714 | 0.543307 | 0.467797 | 1.000000  |
|   >10001bp   |    45    |   44    |  42  |  207  | 0.174603 | 0.488889  | 0.257310 | 1.000000  |

The benchmarking results of ASVCLR in different SV size regions are shown as follows with figures:

<div style="text-align: center;">
    <img src="img\evaluation_result_ASVCLR.png" alt="Performance comparison between different tools" style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\result_classification_ASVCLR.png" alt="Benchmark results between different tools" style="display: inline-block;" width="400"/>
</div>

(3) Statistics of the count of different SV lengths in the user-called set (SVDSS):

|    region    | TP_bench | TP_user |  FP  |  FN   |  recall  | precision | F1 score |  Identity  |
| :----------: | :------: | :-----: | :--: | :---: | :------: | :-------: | :------: | :-------: |
|   1-100bp    |  29479   |  31916  | 8150 | 34665 | 0.459575 | 0.796586  | 0.582873 | 0.973558  |
|  101-250bp   |   2086   |  2126   | 428  | 1872  | 0.527034 | 0.832420  | 0.645426 | 1.000000  |
|  251-500bp   |   1728   |  1869   | 159 | 1244  | 0.581427 | 0.921598  | 0.713018 | 1.000000  |
|  501-1000bp  |   406    |   433   | 75  |  706  | 0.365108 | 0.852362  | 0.511231 | 1.000000  |
| 1001-2500bp  |   395    |   424   | 35  |  440  | 0.473054 | 0.923747  | 0.625690 | 1.000000 |
| 2501-5000bp  |   158    |   164   | 8  |  245  | 0.392060 | 0.953488  | 0.555646 | 1.000000 |
| 5001-10000bp |   0    |   0   |  0  |  0  | 0 | 0 | 0 | 0  |
|   >10001bp   |     0    |   0   |  0  |  0  | 0 | 0 | 0 | 0  |

The benchmarking results of SVDSS in different SV size regions are shown as follows with figures:

<div style="text-align: center;">
    <img src="img\evaluation_result_SVDSS.png" alt="Performance comparison between different tools" style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\result_classification_SVDSS.png" alt="Benchmark results between different tools" style="display: inline-block;" width="400"/>
</div>

(4) Statistics of the count of different SV lengths in the user-called set (DeBreak):

|    region    | TP_bench | TP_user |  FP  |  FN   |  recall  | precision | F1 score |  Identity  |
| :----------: | :------: | :-----: | :--: | :---: | :------: | :-------: | :------: | :-------: |
|   1-100bp    |  31287   |  29770  | 5756 | 32857 | 0.487762 | 0.837978  | 0.616612 | 0.930387  |
|  101-250bp   |   2014   |  1843   | 3658  | 1944  | 0.508843 | 0.335030  | 0.404036 | 0.999851  |
|  251-500bp   |   1743   |  1699   | 2963 | 1229  | 0.586474 | 0.364436  | 0.449532 | 0.999830  |
|  501-1000bp  |   421    |   410   | 1165  |  691  | 0.378597 | 0.260317  | 0.308509 | 0.998064  |
| 1001-2500bp  |   408    |   398   | 510  |  427  | 0.488623 | 0.438326  | 0.462110 | 1.000000 |
| 2501-5000bp  |   198    |   196   | 119  |  205  | 0.491315 | 0.622222  | 0.549074 | 1.000000 |
| 5001-10000bp |   137    |   138   |  33  |  199  | 0.407738 | 0.807018 | 0.541758 | 1.000000  |
|   >10001bp   |     58    |   59   |  96  |  194  | 0.230159 | 0.380645 | 0.286864 | 1.000000  |

The benchmarking results of DeBreak in different SV size regions are shown as follows with figures:

<div style="text-align: center;">
    <img src="img\evaluation_result_DeBreak.png" alt="Performance comparison between different tools" style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\result_classification_DeBreak.png" alt="Benchmark results between different tools" style="display: inline-block;" width="400"/>
</div>

(5) Statistics of the count of different SV lengths in the user-called set (Sniffles2):

|    region    | TP_bench | TP_user |  FP  |  FN   |  recall  | precision | F1 score |  Identity  |
| :----------: | :------: | :-----: | :--: | :---: | :------: | :-------: | :------: | :-------: |
|   1-100bp    |  38459   |  36520  | 9296 | 25685 | 0.599573 | 0.797101  | 0.684369 | 0.914315  |
|  101-250bp   |   2096   |  1926   | 592  | 1862  | 0.529560 | 0.764893  | 0.625835 | 1.000000  |
|  251-500bp   |   1737   |  1699   | 1072 | 1235  | 0.584455 | 0.613136  | 0.598452 | 1.000000  |
|  501-1000bp  |   416    |   413   | 524  |  696  | 0.374101 | 0.440768  | 0.404707 | 1.000000  |
| 1001-2500bp  |   404    |   417   | 173  |  431  | 0.483832 | 0.706780  | 0.574432 | 1.000000  |
| 2501-5000bp  |   205    |   216   | 68  |  198  | 0.508685 | 0.760563  | 0.609632 | 1.000000  |
| 5001-10000bp |   142    |   144   |  66  |  194  | 0.422619 | 0.685714  | 0.522940 | 1.000000  |
|   >10001bp   |    52    |   55    |  93  |  200  | 0.206349 | 0.371622  | 0.265355 | 1.000000  |

The benchmarking results of Sniffles2 in different SV size regions are shown as follows with figures:

<div style="text-align: center;">
    <img src="img\evaluation_result_Sniffles.png" alt="Performance comparison between different tools" style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\result_classification_Sniffles.png" alt="Benchmark results between different tools" style="display: inline-block;" width="400"/>
</div>


(6) Statistics of the count of different SV lengths in the user-called set (pbsv):

|    region    | TP_bench | TP_user |  FP  |  FN   |  recall  | precision | F1 score |  Identity  |
| :----------: | :------: | :-----: | :--: | :---: | :------: | :-------: | :------: | :-------: |
|   1-100bp    |  37343   |  35360  | 7614 | 26801 | 0.582174 | 0.822823  | 0.681890 | 0.962946  |
|  101-250bp   |   1970   |  1792   | 1855 | 1988  | 0.497726 | 0.491363  | 0.494524 | 1.000000  |
|  251-500bp   |   1656   |  1618   | 1609 | 1316  | 0.557201 | 0.501395  | 0.527827 | 1.000000  |
|  501-1000bp  |   355    |   351   | 738  |  757  | 0.319245 | 0.322314  | 0.320772 | 1.000000  |
| 1001-2500bp  |   377    |   373   | 210  |  458  | 0.451497 | 0.639794  | 0.529401 | 1.000000  |
| 2501-5000bp  |   207    |   204   | 83  |  196  | 0.513648 | 0.710801  | 0.596352 | 1.000000  |
| 5001-10000bp |   139    |   140   |  67  |  197  | 0.413690 | 0.676328  | 0.513368 | 1.000000  |
|   >10001bp   |    57    |   59    | 107  |  195  | 0.226190 | 0.355422  | 0.276449 | 1.000000  |

The benchmarking results of pbsv in different SV size regions are shown as follows with figures:

<div style="text-align: center;">
    <img src="img\evaluation_result_pbsv.png" alt="Performance comparison between different tools" style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\result_classification_pbsv.png" alt="Benchmark results between different tools" style="display: inline-block;" width="400"/>
</div>

(7) Statistics of the count of different SV lengths in the user-called set (cuteSV):

|    region    | TP_bench | TP_user |  FP  |  FN   |  recall  | precision | F1 score |  Identity  |
| :----------: | :------: | :-----: | :--: | :---: | :------: | :-------: | :------: | :-------: |
|   1-100bp    |  33438   |  31461  | 5502 | 30706 | 0.521296 | 0.851148  | 0.646584 | 0.912347  |
|  101-250bp   |   1865   |  1683   | 549  | 2093  | 0.471198 | 0.754032  | 0.579970 | 1.000000  |
|  251-500bp   |   1618   |  1576   | 892  | 1354  | 0.544415 | 0.638574  | 0.587747 | 1.000000  |
|  501-1000bp  |   349    |   341   | 472  |  763  | 0.313849 | 0.419434  | 0.359040 | 1.000000  |
| 1001-2500bp  |   349    |   342   | 129  |  486  | 0.417964 | 0.726115  | 0.530540 | 1.000000  |
| 2501-5000bp  |   169    |   169   | 38  |  234  | 0.419355 | 0.816425  | 0.554098 | 1.000000  |
| 5001-10000bp |   113    |   113   |  18  |  223  | 0.336310 | 0.862595  | 0.483940 | 1.000000  |
|   >10001bp   |    39    |   40    |  43  |  213  | 0.154762 | 0.481925  | 0.235287 | 1.000000  |

The benchmarking results of cuteSV in different SV size regions are shown as follow with figures:

<div style="text-align: center;">
    <img src="img\evaluation_result_cuteSV.png" alt="Performance comparison between different tools" style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\result_classification_cuteSV.png" alt="Benchmark results between different tools" style="display: inline-block;" width="400"/>
</div>



(8) Statistics of the count of different SV lengths in the user-called set (SVIM):

|    region    | TP_bench | TP_user |  FP   |  FN   |  recall  | precision | F1 score |  Identity  |
| :----------: | :------: | :-----: | :---: | :---: | :------: | :-------: | :------: | :-------: |
|   1-100bp    |  41223   |  39159  | 26170 | 22921 | 0.642663 | 0.599412  | 0.620285 | 0.952275  |
|  101-250bp   |   2138   |  1954   | 2506  | 1820  | 0.540172 | 0.438117  | 0.483821 | 1.000000  |
|  251-500bp   |   1776   |  1742   | 2507  | 1196  | 0.597577 | 0.409979  | 0.486313 | 1.000000  |
|  501-1000bp  |   442    |   438   | 1312  |  670  | 0.397482 | 0.250286  | 0.307160 | 1.000000  |
| 1001-2500bp  |   426    |   436   |  677  |  409  | 0.510180 | 0.391734  | 0.443179 | 1.000000  |
| 2501-5000bp  |   213    |   215   |  281  |  190  | 0.528536 | 0.433468  | 0.476304 | 1.000000  |
| 5001-10000bp |   143    |   146   |  218  |  193  | 0.425595 | 0.401099  | 0.412984 | 1.000000  |
|   >10001bp   |    62    |   71    |  393  |  190  | 0.246032 | 0.153017  | 0.188684 | 1.000000  |

The benchmarking results of SVIM in different SV size regions are shown as follows with figures:

<div style="text-align: center;">
    <img src="img\evaluation_result_SVIM.png" alt="Performance comparison between different tools" style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\result_classification_SVIM.png" alt="Benchmark results between different tools" style="display: inline-block;" width="400"/>
</div>

### Quantitative distribution statistics of variants

Additionally, the distribution of the number of variants in the reference set and  user-called set was statistically analyzed, as shown below:

(1) Statistics of the count of different SV lengths in the benchmark set:
The SV reference region size statistics for benchmark set: Total SVs number：74012
<div align="center">
<img src="img\ref_reg_size_benchmark.png" alt= "ref_reg_size_benchmark" style="display: inline-block; margin: 0 auto;" width="450"> 
</div>

(2) Statistics of the count of different SV lengths in the user-called set (ASVCLR):

The SV reference region size statistics before filtering for the user-called set (ASVCLR): Total SVs number：50674
The SV reference region size statistics after filtering for the user-called set (ASVCLR): Total SVs number：50674                                                                                          
The result statistics before filtering are shown in the left figure, and the result statistics after filtering are shown in the right figure:

<div style="text-align: center;">
    <img src="img\ref_reg_size_ASVCLR_before.png" alt="Performance comparison between different tools"  style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\ref_reg_size_ASVCLR_after.png" alt="Benchmark results between different tools" width="400"/>
</div> 

(3) Statistics of the count of different SV lengths in the user-called set (SVDSS):

The SV reference region size statistics before filtering for the user-called set (ASVCLR): Total SVs number：50674
The SV reference region size statistics after filtering for the user-called set (ASVCLR): Total SVs number：50674                                                                                          
The result statistics before filtering are shown in the left figure, and the result statistics after filtering are shown in the right figure:

<div style="text-align: center;">
    <img src="img\ref_reg_size_SVDSS_before.png" alt="Performance comparison between different tools"  style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\ref_reg_size_SVDSS_after.png" alt="Benchmark results between different tools" width="400"/>
</div>

(4) Statistics of the count of different SV lengths in the user-called set (DeBreak):

The SV reference region size statistics before filtering for the user-called set (ASVCLR): Total SVs number：50674
The SV reference region size statistics after filtering for the user-called set (ASVCLR): Total SVs number：50674                                                                                          
The result statistics before filtering are shown in the left figure, and the result statistics after filtering are shown in the right figure:

<div style="text-align: center;">
    <img src="img\ref_reg_size_DeBreak_before.png" alt="Performance comparison between different tools"  style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\ref_reg_size_DeBreak_after.png" alt="Benchmark results between different tools" width="400"/>
</div>

(5) Statistics of the count of different SV lengths in the user-called set (Sniffles2):

The SV reference region size statistics before filtering for the user-called set (Sniffles2): Total SVs number：54545
The SV reference region size statistics after filtering for the user-called set (Sniffles2): Total SVs number：54458
The result statistics before filtering are shown in the left figure, and the result statistics after filtering are shown in the right figure:

<div style="text-align: center;">
    <img src="img\ref_reg_size_Sniffles_before.png" alt="Performance comparison between different tools" style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\ref_reg_size_Sniffles_after.png" alt="Benchmark results between different tools" width="400"/>
</div>


(6) Statistics of the count of different SV lengths in the user-called set (pbsv):

The SV reference region size statistics before filtering for the user-called set (pbsv): Total SVs number：52807
The SV reference region size statistics after filtering for the user-called set (pbsv): Total SVs number：52741
The result statistics before filtering are shown in the left figure, and the result statistics after filtering are shown in the right figure:

<div style="text-align: center;">
    <img src="img\ref_reg_size_pbsv_before.png" alt="Performance comparison between different tools" style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\ref_reg_size_pbsv_after.png" alt="Benchmark results between different tools" width="400"/>

</div>

(7) Statistics of the count of different SV lengths in the user-called set (cuteSV):

The SV reference region size statistics before filtering for the user-called set (cuteSV): Total SVs number：44937
The SV reference region size statistics after filtering for the user-called set (cuteSV): Total SVs number：44928
The result statistics before filtering are shown in the left figure, and the result statistics after filtering are shown in the right figure:

<div style="text-align: center;">
    <img src="img\ref_reg_size_cuteSV_before.png" alt="Performance comparison between different tools" style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\ref_reg_size_cuteSV_after.png" alt="Benchmark results between different tools" width="400"/>
</div>




(8) Statistics of the count of different SV lengths in the user-called set (SVIM):

The SV reference region size statistics before filtering for the user-called set (SVIM): Total SVs number：116615
The SV reference region size statistics after filtering for the user-called set (SVIM): Total SVs number：116427
The result statistics before filtering are shown in the left figure, and the result statistics after filtering are shown in the right figure:

<div style="text-align: center;">
    <img src="img\ref_reg_size_SVIM_before.png" alt="Performance comparison between different tools" style="display: inline-block; margin-right: 20px;" width="400"/>
    <img src="img\ref_reg_size_SVIM_after.png" alt="Benchmark results between different tools" width="400"/>
</div>

 More detailed result information can be found in the `sv_stat_reports.html` after completing the assessment following the above steps.

## Contact

If you have problems or some suggestions, please contact: zhuxiao_hit@yeah.net without hesitation. 

---- Enjoy !!! -----
