# somaticContamination

IonServer plugin for the assessment of contamination in somatic samples. This plugin creates plots of allele frequencies across samples so that contamination can be visually assessed.

We wrote this plugin because we have:
* unpaired somatic samples
* potentially very low allele frequency mutations

The plugin 1st gets a combined list of variants for all samples and then, using pysam, looks at samples bam to look at raw support for that variant. Using seaborn it then plots a heatmap of allele frequencies computed from raw coverage in each bam.

## Requirements

You *possibly* need to install the following packages on your IonServer

libpng-dev
libjpeg8-dev
libfreetype6-dev
build-dep
python-dev
python-numpy

## Installation

Besides the code in this repo you also need exac vcf file which is downloaded by `deploy.sh` from amazon s3. This script also makes vcfanno executable.

1. Clone the git repo into `/results/plugins` on your IonServer
2. `cd /results/plugins`
3. `bash deploy.sh`

