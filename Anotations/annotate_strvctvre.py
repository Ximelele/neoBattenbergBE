from StrVCTVRE_annot import start_analysis
from helper_functions import safe_run
import os


def annotate_with_strvctvre(vcf_input: str):
    annotated_vcf_output = vcf_input + "_StrVCTVRE.vcf"

    strvctvre_command = f"python /app/StrVCTVRE/StrVCTVRE.py -i {vcf_input} -o {annotated_vcf_output}"
    if not os.path.exists(annotated_vcf_output):
        safe_run(strvctvre_command)
    print(f"VCF file annotated using StrVCTVRE. Annotated VCF file: {annotated_vcf_output}")
    start_analysis(annotated_vcf_output)
