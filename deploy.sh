#!/usr/bin/env bash

#make exac dir

mkdir -p /results/sdgs-resources/hg19

#get exac

cd /results/sdgs-resources/hg19
wget https://s3.eu-west-2.amazonaws.com/sheffieldgenomics/reference_files/hg19/ExAC.r0.3.1.sites.vep.vcf.gz
wget https://s3.eu-west-2.amazonaws.com/sheffieldgenomics/reference_files/hg19/ExAC.r0.3.1.sites.vep.vcf.gz.tbi


#make vcfanno executale

chmod +x /results/plugins/somaticContamination/bin/vcfanno

#make? & activate venv and install requirements

