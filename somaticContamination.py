#!/usr/bin/env/python
# somaticContamination plugin
import os
import glob
import json
import traceback
import subprocess
import requests
from ion.utils import blockprocessing
from ion.plugin import *
from django.conf import settings
from django import template
from django.template.loader import render_to_string
from django.conf import global_settings
from subprocess import *
# global_settings.LOGGING_CONFIG=None


class somaticContamination(IonPlugin):
    """
    somaticContamination
    """

    version = '1.0'
    allow_autorun = True
    major_block = True

    def createreport(self, reportName, reportTemplate, reportData, data):
        print "Creating Report\n"
        if not settings.configured:
            plugin_dir = data['runinfo']['plugin']['path']
            settings.configure(DEBUG=False, TEMPLATE_DEBUG=False,
                               INSTALLED_APPS=('django.contrib.humanize',),
                               TEMPLATE_DIRS=(os.path.join(plugin_dir, 'templates'),))
        print reportData
        with open(reportName, 'w') as report:
            report.write(render_to_string(reportTemplate, { 'data': reportData }))

        report.close()

    def load_startpluginjson(self):
        with open('startplugin.json','r') as file:
            return json.load(file)

    def launch(self, data=None):
        data = self.load_startpluginjson()

        print json.dumps(data, indent=4)
        info = {}
        for sample_name in data["plan"]["barcodedSamples"]:
            for barcode in data["plan"]["barcodedSamples"][sample_name]["barcodeSampleInfo"]:
                info[barcode] = {}
                info[barcode]["sample"] = sample_name
                info[barcode]["bed"] = data["plan"]["barcodedSamples"][sample_name]["barcodeSampleInfo"][barcode][
                    "targetRegionBedFile"]

        bam_dir = data["runinfo"]["alignment_dir"]
        bams = glob.glob(bam_dir + "/*.bam")
        for bam in bams:
            barcode = os.path.basename(bam).replace("_rawlib.bam", "")
            info[barcode]["bam"] = bam

        # find vcf files

        main_plugin_dir = os.path.dirname(data['runinfo']['results_dir'])
        variant_caller_dirs = glob.glob(main_plugin_dir + "/variant*")
        vcfs = []
        for barcode in info:
            vcf = variant_caller_dirs[0]+"/"+barcode+"/TSVC_variants.vcf"
            if os.path.isfile(vcf):
                vcfs.append(vcf)


        bam_list = ",".join(bams)

        print bam_list

        print data['runinfo']['results_dir']+"/"

        #run annotation
        annotated_vcfs = []
        for vcf in vcfs:
            print "annotaing " + vcf
            barcode = os.path.basename(os.path.dirname(vcf))
            out_vcf = data['runinfo']['results_dir']+"/"+barcode+".annotated.vcf"

            with open(out_vcf,"wb") as out:
                out.write(subprocess.check_output([
                    '%s/bin/vcfanno' % os.environ['DIRNAME'], '%s/conf/conf.toml' % os.environ['DIRNAME'], vcf
                ]))
                out.close()
            annotated_vcfs.append(out_vcf)

        vcf_list = ",".join(annotated_vcfs)

        #plot contamination

        print "getting coverage & plotting contamination"

        #print '%s/somaticContamination_venv/bin/python' % os.environ['DIRNAME'],'%s/somaticContamination_main.py' % os.environ['DIRNAME'], '--listofbams', bam_list, '--listofvcfs', vcf_list, '--output_dir', data['runinfo']['results_dir']+"/", '--output_type', 'freq'

        plugin = Popen([
            '%s/somaticContamination_venv/bin/python' % os.environ['DIRNAME'],'%s/somaticContamination_main.py' % os.environ['DIRNAME'], '--listofbams', bam_list, '--listofvcfs', vcf_list, '--output_dir', data['runinfo']['results_dir']+"/", '--output_type', 'freq'
        ], stdout=PIPE, shell=False)
        plugin.communicate()

        block_file = 'somaticContamination_block.html'

        result = {"somaticImage": data['runinfo']['results_dir'].replace("/results/analysis","") + "/somatic.png",
                  "germlineImage": data['runinfo']['results_dir'].replace("/results/analysis","") + "/germline.png"}

        self.createreport(block_file, 'report_block.html', result, data)

        #sys.exit(plugin.poll())

if __name__ == '__main__':
    PluginCLI(somaticContamination)