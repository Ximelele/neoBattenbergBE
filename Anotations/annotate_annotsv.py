from AnnotSV import start_analysis_annotsv
from helper_functions import safe_run
import os


def annotate_with_annotcv(vcf_input: str):
    annot_output = vcf_input + "_annotSV.tsv"
    annotsv_command = f"/app/AnnotSV/bin/AnnotSV -SVinputFile {vcf_input} -outputFile {annot_output} -outputDir ."
    if not os.path.exists(annot_output):
        safe_run(annotsv_command)
    start_analysis_annotsv(annot_output)
