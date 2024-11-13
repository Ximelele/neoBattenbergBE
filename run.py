import StrVCTVRE_annot
import AnnotSV
import threading

from annotate_annotsv import annotate_with_annotcv
from annotate_strvctvre import annotate_with_strvctvre
from annotated_vcf import process_bam_file
from helper_functions import safe_run


def run_battenberg():
    r_command = f"Rscript batt.R"
    safe_run(r_command)


def main():
    # Set up command-line argument parsing
    print("Running delly")
    vcf_file_name = process_bam_file("/data/Resources/Whitespring.HG38/result/mapping/bqsr/Lynch.1827.06.N.bam",
                                     "/data/Results/MartinD/hg38.fa")

    threads = []
    print("Running strvctvre")
    str_thread = threading.Thread(target=annotate_with_strvctvre, args=(vcf_file_name,))
    threads.append(str_thread)
    print("Running annotsv")
    ann_thread = threading.Thread(target=annotate_with_annotcv, args=(vcf_file_name,))
    print("Running batt")
    batt_thread = threading.Thread(target=run_battenberg, args=())
    threads.append(batt_thread)
    threads.append(ann_thread)

    print("Starting threads")
    for thread in threads:
        thread.start()

    for thread in threads:
        thread.join()
    print("Finnished")


if __name__ == "__main__":
    main()
