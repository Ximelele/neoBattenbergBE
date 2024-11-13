import subprocess
import os
from pathvalidate import sanitize_filepath
from helper_functions import safe_run


def process_bam_file(bam_file, reference_genome="hg38.fa"):
    # Sanitize file paths
    bam_file = sanitize_filepath(bam_file)
    reference_genome = sanitize_filepath(reference_genome)

    # Generate output prefix based on the bam file name
    output_prefix_parts = bam_file.split("/")[-1].split(".")
    output_prefix_parts.pop()  # Remove the file extension
    output_prefix = '_'.join(output_prefix_parts)

    # Define the output files
    bcf_output = f"{output_prefix}.bcf"
    vcf_output = f"{output_prefix}.vcf"

    # Check if the VCF file already exists
    if not os.path.exists(vcf_output):

        try:
            # Step 1: Run delly call to create BCF file
            delly_command = f"delly call -o {bcf_output} -g {reference_genome} {bam_file}"
            safe_run(delly_command)
            print(f"Delly call completed. BCF file generated: {bcf_output}")

            # Step 2: Convert BCF to VCF using bcftools
            bcftools_command = f"bcftools view {bcf_output} -o {vcf_output}"
            safe_run(bcftools_command)
            print(f"BCF file converted to VCF: {vcf_output}")

            # Step 3: Delete the BCF file and its index file if it exists
            os.remove(bcf_output)
            if os.path.exists(bcf_output + ".csi"):
                os.remove(bcf_output + ".csi")
            print(f"BCF file and index file deleted: {bcf_output}")

        except subprocess.CalledProcessError as e:
            print("An error occurred during the execution of a command:", e)
        except FileNotFoundError as e:
            print("An error occurred while trying to delete the BCF file:", e)

    return vcf_output


