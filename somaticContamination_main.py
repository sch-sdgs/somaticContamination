import pysam
import argparse
import vcf
import os
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
import progressbar
import seaborn as sns
import json

def is_non_zero_file(fpath):
    """
    checks if a file exists and  is not empty

    :param fpath: full path to the file
    :return: True or False
    """
    if os.path.isfile(fpath) and os.path.getsize(fpath) > 0:
        return True
    else:
        return False

def get_coverage(list_of_variants,bams):
    """

    gets coverage for a variant in a list of bam files

    :param list_of_variants: list of dictionaries containing CHROM POS ALT & REF as keys
    :param listofbams: plain list of the full path of the bam files to interrogate coverage
    :return: final result file - list of dicts now containing information from bams
    """

    bam_number = len(bams)
    print bam_number
    #with progressbar.ProgressBar(maxval=bam_number) as bar:
        #for i in range(bam_number):
    for bam in bams:
            #bam=bams[i]
        if is_non_zero_file(bam):
            samfile = pysam.AlignmentFile(bam, "rb")
            for j in range(len(list_of_variants)):
                record = {}
                record["BAM"] = bam
                record["REF_COUNT"] = 0
                record["ALT_COUNT"] = 0
                variant = list_of_variants[j]
                for pileupcolumn in samfile.pileup(variant["CHROM"], variant["POS"]-1,  variant["POS"]):
                    if pileupcolumn.pos == variant["POS"]-1:
                        for pileupread in pileupcolumn.pileups:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                if pileupread.alignment.query_sequence[pileupread.query_position] == variant["REF"]:
                                    record["REF_COUNT"]+=1
                                if pileupread.alignment.query_sequence[pileupread.query_position] == variant["ALT"]:
                                    record["ALT_COUNT"]+=1
                if record["REF_COUNT"]+record["ALT_COUNT"] != 0:
                    record["FREQ"]=record["ALT_COUNT"]/float(record["REF_COUNT"]+record["ALT_COUNT"])
                else:
                    record["FREQ"]=0 #todo change this to NA and deal with that in the data frame

                list_of_variants[j]["BAMS"].append(record)
        else:
            print "WARNING: " + bam.rstrip() + " is empty or not found"
     #       bar.update(i)
    return list_of_variants



def make_variant_list(listofvcfs):
    """
    makes a list of variants for analysis from a vcf file

    :param listofvcfs:
    :return: list of dicts containing key information about each variant
    """
    result=[]
    print "Loading VCF files...."
    for vcf_file in listofvcfs:
        for i in vcf.VCFReader(filename=vcf_file):
            if "AF_exac" in i.INFO:
                af = i.INFO["AF_exac"]
            else:
                af=[]
                for j in i.ALT:
                    af.append(0)
            if not i.is_indel:
                for count,alt in enumerate(i.ALT):
                        record = {}
                        record["CHROM"] = i.CHROM
                        record["POS"] = i.POS
                        record["REF"] = i.REF
                        record["ALT"] = alt
                        record["BAMS"] = list()
                        if af[count] >= 0.01:
                            record["EXAC"] = 1
                        else:
                            record["EXAC"] = 0
                        id = i.ID
                        if id is None:
                            id = "."
                        record["ID"] = id
                        if record not in result:
                            result.append(record)

    return result


def output_as_matrix(result,type):
    """
    reformats the result into a list (rows) of lists (colums)

    :param result: the result of the analysis
    :param type: type of result to output - frequency or count data
    :return: list of lists
    """
    output = []

    header = ["variant","id","exac"]
    for variant in result:
        variantid = variant["CHROM"] + ":" + str(variant["POS"]) + "." + variant["REF"] + "/" + str(variant["ALT"])
        row = [variantid,variant["ID"],str(variant["EXAC"])]
        for b in variant["BAMS"]:
            bam = os.path.basename(b["BAM"])
            bam = bam.rstrip("\.bam")
            if bam not in header:
                header.append(bam)

            if type == 'freq':
                freq = str(b["FREQ"])
                row.append(freq)
            if type == 'count':
                count = str(b["ALT_COUNT"])+"|"+str(b["REF_COUNT"]+b["ALT_COUNT"])
                row.append(count)
        output.append(row)
    output.append(header)
    return output


def write_output_to_file(output,output_dir):
    """
    writes the result of the analysis to a file

    :param output: the analysis in matric format (list of lists)
    :param file: output file
    """
    for_writing='\n'.join('\t'.join(inner) for inner in output)
    #print for_writing
    target = open(output_dir + "/contamination.tsv", 'w')
    target.write(for_writing)
    target.close()


def heatmap(output,output_dir):
    """
    draws 2 heatmaps of the data - one for presumed germline variants and one for presumed somatic

    :param output: the data in a matrix format (list of lists)
    :param output_dir: the output directory to write the images to
    """

    length = len(output)
    headers = output.pop(length-1)
    df = pd.DataFrame(output,columns=headers,dtype=float)

    df.set_index("variant",inplace=True,drop=True)
    df.drop('id',axis=1,inplace=True)

    somatic = df.loc[(df.exac == 0.0)]
    germline = df.loc[(df.exac == 1.0)]
    somatic.drop('exac', axis=1, inplace=True)
    germline.drop('exac', axis=1, inplace=True)

    a4_dims = (8.27, 11)

    fig, ax = plt.subplots(figsize=a4_dims)
    ax.set_title("Somatic or Rare/De Novo Variant (ExAC < 0.01)")
    plot = sns.heatmap(somatic, annot=False, vmin=0, vmax=0.1, square=True, ax=ax,cmap="YlGnBu",linewidths=.5)
    plt.subplots_adjust(bottom=0.2,left=0.3)
    fig = plot.get_figure()
    fig.savefig(output_dir+"somatic.png")
    fig, ax = plt.subplots(figsize=a4_dims)
    ax.set_title("ExAC >= 0.01")
    plot = sns.heatmap(germline, annot=False, vmin=0, vmax=1, square=True, ax=ax,cmap="YlGnBu",linewidths=.5)
    plt.subplots_adjust(bottom=0.2,left=0.3)
    fig = plot.get_figure()
    fig.savefig(output_dir+"germline.png")


def main():
    parser = argparse.ArgumentParser(description='Calculates cross sample comtamination - somatic')
    parser.add_argument('--listofbams', metavar='listofbams', type=str, help='list of bam files to check')
    parser.add_argument('--listofvcfs', metavar='listofvcfs', type=str, help='list of vcf files')
    parser.add_argument('--output_type', metavar='output_type', type=str, help='type freq or count')
    parser.add_argument('--output_dir', metavar='output_dir', type=str, help='output dir')
    args = parser.parse_args()

    out = open(args.output_dir+"/log.txt","wb")
    out.write("hello")
    out.close()


    listofvcfs = args.listofvcfs.split(",")
    listofbams = args.listofbams.split(",")

    variant_list = make_variant_list(listofvcfs)
    print variant_list
    result = get_coverage(variant_list,listofbams)
    print result
    output = output_as_matrix(result,args.output_type)
    write_output_to_file(output, args.output_dir)
    print len(variant_list)
    if len(variant_list) < 100:
        heatmap(output,args.output_dir)


if __name__ == '__main__':
    main()
