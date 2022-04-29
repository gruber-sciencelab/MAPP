"""
##############################################################################
#
#   Parser for the snakemake config template.
#   Generates the real config.yml file for the MAPP pipeline.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 04-07-2020
#   LICENSE: Apache_2.0
#
##############################################################################
"""

# imports
import time
import logging
import logging.handlers
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import yaml
import glob
import pandas as pd


def parse_arguments():
    """Parser of the command-line arguments."""
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-v",
        "--verbosity",
        dest="verbosity",
        choices=("DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"),
        default="ERROR",
        help="Verbosity/Log level. Defaults to ERROR",
    )
    parser.add_argument(
        "-l", "--logfile", dest="logfile", help="Store log to this file."
    )
    parser.add_argument(
        "--config-template",
        dest="template",
        required=True,
        help="Meta-config file for the snakemake pipeline.",
    )
    parser.add_argument(
        "--pipeline-configfile",
        dest="pipeline_configfile",
        required=True,
        help="Path for the output YAML config file of the MAPP pipeline.",
    )
    return parser


##############################################################################


def generate_windows(size, up, down):
    """Generate encoding for the sliding windows."""
    windows = []
    s = -1 * up
    e = s + size
    while e <= down:
        windows.append((s, e))
        s = s + int(float(size) / 2.0)
        e = e + int(float(size) / 2.0)
    # convert to dirextory names
    windows_directories = []
    for win in windows:
        if win[0] < 0 and win[1] <= 0:  # upstream ss
            windows_directories.append(
                "u" + str(-1 * win[0]) + "to" + str(-1 * win[1]) + ".d0to0"
            )
        elif win[0] >= 0 and win[1] > 0:  # downstream ss
            windows_directories.append("u0to0.d" + str(win[0]) + "to" + str(win[1]))
        else:  # window goes through ss
            windows_directories.append(
                "u" + str(-1 * win[0]) + "to0.d0to" + str(win[1])
            )
    return windows_directories


def check_config_paths(config_template):
    """Check if all the paths have been provided correctly."""

    if not os.path.isdir(config_template["MAPP_directory"]):
        print("### ERROR ###")
        print("Please provide a correct path for the \"MAPP_directory\" parameter")
        exit()
    
    if not os.path.exists(config_template["genomic_annotation"]):
        print("### ERROR ###")
        print("Please provide a correct path for the \"genomic_annotation\" parameter")
        exit()
    
    if not os.path.exists(config_template["genomic_sequence"]):
        print("### ERROR ###")
        print("Please provide a correct path for the \"genomic_sequence\" parameter")
        exit()
    
    if not os.path.exists(config_template["analysis_design_table"]):
        print("### ERROR ###")
        print("Please provide a correct path for the \"analysis_design_table\" parameter")
        exit()
    
    if config_template["matrix_type"] == "kmers":

        if config_template["PWM_directory"] != "":
            print("### ERROR ###")
            print("Invalid value for \"PWM_directory\" parameter: please provide \"\" for a k-mer based analysis")
            exit()
        
        if config_template["seqlogo_directory"] != "":
            print("### ERROR ###")
            print("Invalid value for \"seqlogo_directory\" parameter: please provide \"\" for a k-mer based analysis")
            exit()

    if config_template["matrix_type"] == "pwms":

        if not os.path.isdir(config_template["PWM_directory"]):
            print("### ERROR ###")
            print("Invalid value for \"PWM_directory\" parameter: please provide a correct path for a PWM based analysis")
            exit()
        
        if not os.path.isdir(config_template["seqlogo_directory"]):
            print("### ERROR ###")
            print("Invalid value for \"seqlogo_directory\" parameter: please provide a correct path for a PWM based analysis")
            exit()

    if not os.path.exists(config_template["PAS_atlas"]):
        print("### ERROR ###")
        print("Please provide a correct path for the \"PAS_atlas\" parameter")
        exit()

    if config_template["genomic_index"] == "":
        home_index = config_template["genomic_sequence"].split(os.sep)[-1]
        suffix = "__index_sjdbOverhang" + str(config_template["sjdbOverhang"])
        home_index = os.path.join(
            os.path.expanduser("~"),
            home_index + suffix
        )
        config_template["genomic_index"] = home_index
        print("### WARNING ###")
        print("No value provided for the \"genomic_index\" parameter")
        print("Setting a default value of: " + home_index)


def main():
    """Main body of the script."""

    # read the YAML config template
    with open(options.template, "r") as f:
        template = yaml.safe_load(f)

    # check paths in the config template
    check_config_paths(template)

    # parse the samples design table
    design_table = pd.read_csv(template["analysis_design_table"], sep="\t")
    sample_IDs = design_table["sample"]
    conditions = design_table["condition"]

    # default directory name for the per-modules output
    default_output_dir_name = "output"

    # set the pipeline configfile path
    if os.path.isabs(options.pipeline_configfile):
        pipeline_configfile = options.pipeline_configfile
    else:
        pipeline_configfile = os.path.join(
            template["MAPP_directory"], options.pipeline_configfile
        )

    # dynamicaly created entries:
    ##########################################################################
    QEI_sample_bam_dict = []
    for s in sample_IDs:
        bam_path = os.path.join(
            template["MAPP_directory"],
            "modules",
            "PREPROCESSING",
            default_output_dir_name,
            "alignments",
            s,
            s + ".Aligned.toTranscriptome.out.bam",
        )
        QEI_sample_bam_dict.append('  "' + s + '": "' + bam_path + '"')
    QEI_sample_bam_dict = os.linesep.join(QEI_sample_bam_dict)

    PAQ_sample_bam_dict = []
    for s in sample_IDs:
        bam_path = os.path.join(
            template["MAPP_directory"],
            "modules",
            "PREPROCESSING",
            default_output_dir_name,
            "alignments",
            s,
            s + ".Aligned.out.sorted.bam",
        )
        PAQ_sample_bam_dict.append('  "' + s + '": "' + bam_path + '"')
    PAQ_sample_bam_dict = os.linesep.join(PAQ_sample_bam_dict)

    PAQ_sample_bai_dict = []
    for s in sample_IDs:
        bai_path = os.path.join(
            template["MAPP_directory"],
            "modules",
            "PREPROCESSING",
            default_output_dir_name,
            "alignments",
            s,
            s + ".Aligned.out.sorted.bam.bai",
        )
        PAQ_sample_bai_dict.append('  "' + s + '": "' + bai_path + '"')

    PAQ_sample_bai_dict = os.linesep.join(PAQ_sample_bai_dict)

    PAQ_sample_condition_dict = []
    for s, condition in zip(sample_IDs, conditions,):
        PAQ_sample_condition_dict.append('  "' + s + '": "' + condition + '"')
    PAQ_sample_condition_dict = os.linesep.join(PAQ_sample_condition_dict)

    if template["library_type"] == "stranded":
        unstranded_string = "no"
    else:
        unstranded_string = "yes"

    if template["quality_check"]:
        updated_design_table = os.path.join(
            template["MAPP_directory"],
            "modules",
            "PREPROCESSING",
            default_output_dir_name,
            "design_table_quality_filtered.tsv",
        )
    else:
        updated_design_table = template["analysis_design_table"]

    MAE_windows_3ss = []
    if template["sitecount_matrices_directories"]["3ss"] == "":
        for w in generate_windows(
            template["window_size"],
            template["3ss_region_up"],
            template["3ss_region_down"],
        ):
            matrix_path = os.path.join(
                template["MAPP_directory"],
                "modules",
                "CREATE_SITECOUNT_MATRICES",
                default_output_dir_name,
                "3ss",
                w,
                "matrix.tsv",
            )
            MAE_windows_3ss.append('    "' + w + '": "' + matrix_path + '"')
    else:
        for window_path in glob.glob(
            os.path.join(template["sitecount_matrices_directories"]["3ss"], "*")
        ):
            w = window_path.split(os.sep)[-1]
            MAE_windows_3ss.append(
                '    "' + w + '": "' + window_path + os.sep + 'matrix.tsv"'
            )
            MAE_windows_3ss = sorted(MAE_windows_3ss)
    MAE_windows_3ss = os.linesep.join(MAE_windows_3ss)

    MAE_windows_5ss = []
    if template["sitecount_matrices_directories"]["5ss"] == "":
        for w in generate_windows(
            template["window_size"],
            template["5ss_region_up"],
            template["5ss_region_down"],
        ):
            matrix_path = os.path.join(
                template["MAPP_directory"],
                "modules",
                "CREATE_SITECOUNT_MATRICES",
                default_output_dir_name,
                "5ss",
                w,
                "matrix.tsv",
            )
            MAE_windows_5ss.append('    "' + w + '": "' + matrix_path + '"')
    else:
        for window_path in glob.glob(
            os.path.join(template["sitecount_matrices_directories"]["5ss"], "*")
        ):
            w = window_path.split(os.sep)[-1]
            MAE_windows_5ss.append(
                '    "' + w + '": "' + window_path + os.sep + 'matrix.tsv"'
            )
            MAE_windows_5ss = sorted(MAE_windows_5ss)
    MAE_windows_5ss = os.linesep.join(MAE_windows_5ss)

    KPC_windows_pas = []
    if template["sitecount_matrices_directories"]["pas"] == "":
        for w in generate_windows(
            template["window_size"],
            template["pas_region_up"],
            template["pas_region_down"],
        ):
            matrix_path = os.path.join(
                template["MAPP_directory"],
                "modules",
                "CREATE_SITECOUNT_MATRICES",
                default_output_dir_name,
                "pas",
                w,
                "matrix.tsv",
            )
            KPC_windows_pas.append('    "' + w + '": "' + matrix_path + '"')
    else:
        for window_path in glob.glob(
            os.path.join(template["sitecount_matrices_directories"]["pas"], "*")
        ):
            w = window_path.split(os.sep)[-1]
            KPC_windows_pas.append(
                '    "' + w + '": "' + window_path + os.sep + 'matrix.tsv"'
            )
            KPC_windows_pas = sorted(KPC_windows_pas)
    KPC_windows_pas = os.linesep.join(KPC_windows_pas)

    if template["tandem_pas_coordinates"] == "":
        tandem_pas_coordinates = os.path.join(
            template["MAPP_directory"],
            "modules",
            "PREPARE_TANDEM_PAS",
            default_output_dir_name,
            "tandem_pas.terminal_exons.clean.bed",
        )
    else:
        tandem_pas_coordinates = template["tandem_pas_coordinates"]

    if template["tandem_pas_representative_sites_coordinates"] == "":
        tandem_pas_representative_sites_coordinates = os.path.join(
            template["MAPP_directory"],
            "modules",
            "PREPARE_TANDEM_PAS",
            default_output_dir_name,
            "tandem_pas.terminal_exons.representative_sites.bed",
        )
    else:
        tandem_pas_representative_sites_coordinates = template[
            "tandem_pas_representative_sites_coordinates"
        ]

    if template["3ss_coordinates"] == "":
        _3ss_coordinates = os.path.join(
            template["MAPP_directory"],
            "modules",
            "EXTRACT_AS_EXONS",
            default_output_dir_name,
            "3ss.bed",
        )
    else:
        _3ss_coordinates = template["3ss_coordinates"]

    if template["5ss_coordinates"] == "":
        _5ss_coordinates = os.path.join(
            template["MAPP_directory"],
            "modules",
            "EXTRACT_AS_EXONS",
            default_output_dir_name,
            "5ss.bed",
        )
    else:
        _5ss_coordinates = template["5ss_coordinates"]

    if template["skipped_exon_events"] == "":
        skipped_exon_events = os.path.join(
            template["MAPP_directory"],
            "modules",
            "EXTRACT_AS_EXONS",
            default_output_dir_name,
            "events_SE_selected.tsv",
        )
    else:
        skipped_exon_events = template["skipped_exon_events"]

    ##########################################################################

    config_string = f"""---

MAPP_analysis_name: "{template["analysis_name"]}"
MAPP_pipeline_directory: "{template["MAPP_directory"]}"
MAPP_pipeline_configfile: "{pipeline_configfile}"
MAPP_seqlogos_directory: "{template["seqlogo_directory"]}"

### module: PREPROCESSING
PQA_scripts_dir: "{template["MAPP_directory"]}/modules/PREPROCESSING/scripts"
PQA_outdir: "{template["MAPP_directory"]}/modules/PREPROCESSING/{default_output_dir_name}"
PQA_adapters_sequences: "{template["MAPP_directory"]}/modules/PREPROCESSING/resources/adapters.txt"
PQA_genomic_sequence: "{template["genomic_sequence"]}"
PQA_genomic_annotation: "{template["genomic_annotation"]}"
PQA_index: "{template["genomic_index"]}"
PQA_design_file: "{template["analysis_design_table"]}"
PQA_sjdbOverhang: {template["sjdbOverhang"]}
PQA_transcript_biotypes: "{template["transcript_biotypes"]}"
PQA_min_median_TIN_score: {template["min_median_TIN_score"]}
PQA_RNASeQC_min_mapping_rate: {template["RNASeQC_min_mapping_rate"]}
PQA_RNASeQC_min_unique_rate_of_mapped: {template["RNASeQC_min_unique_rate_of_mapped"]}
PQA_RNASeQC_min_high_quality_rate: {template["RNASeQC_min_high_quality_rate"]}
PQA_RNASeQC_max_intergenic_rate: {template["RNASeQC_max_intergenic_rate"]}
PQA_RNASeQC_max_rRNA_rate: {template["RNASeQC_max_rRNA_rate"]}
PQA_storage_efficient: {template["storage_efficient"]}

### module: EXTRACT_AS_EXONS
ASE_scripts_dir: "{template["MAPP_directory"]}/modules/EXTRACT_AS_EXONS/scripts"
ASE_outdir: "{template["MAPP_directory"]}/modules/EXTRACT_AS_EXONS/{default_output_dir_name}"
ASE_genomic_annotation: "{template["genomic_annotation"]}"
ASE_genomic_sequence: "{template["genomic_sequence"]}"
ASE_3ss_region_up: {template["3ss_region_up"]}
ASE_3ss_region_down: {template["3ss_region_down"]}
ASE_5ss_region_up: {template["5ss_region_up"]}
ASE_5ss_region_down: {template["5ss_region_down"]}
ASE_transcript_biotypes: "{template["transcript_biotypes"]}"
ASE_clustering_max_distance: {template["exon_clustering_max_coverage_distance"]}

### module: PREPARE_TANDEM_PAS
TPA_scripts_dir: "{template["MAPP_directory"]}/modules/PREPARE_TANDEM_PAS/scripts"
TPA_outdir: "{template["MAPP_directory"]}/modules/PREPARE_TANDEM_PAS/{default_output_dir_name}"
TPA_pas_atlas: "{template["PAS_atlas"]}"
TPA_genomic_annotation: "{template["genomic_annotation"]}"
TPA_library_type: "{template["library_type"]}"
TPA_min_protocols: {template["polyA_sites_min_no_protocol"]}
TPA_biotype_key: "transcript_biotype"
TPA_biotype_values:
    - "{template["transcript_biotypes"]}"
TPA_three_prime_end_offset: {template["three_prime_end_offset"]}
TPA_transcript_locus_offset: {template["transcript_locus_offset"]}
TPA_downstream_extend: {template["tandem_polyA_exon_extension"]}

### module: CREATE_SITECOUNT_MATRICES
CSM_scripts_dir: "{template["MAPP_directory"]}/modules/CREATE_SITECOUNT_MATRICES/scripts"
CSM_outdir: "{template["MAPP_directory"]}/modules/CREATE_SITECOUNT_MATRICES/{default_output_dir_name}"
CSM_genomic_sequence: "{template["genomic_sequence"]}"
CSM_window_size: {template["window_size"]}
CSM_regions_ranges:
  "pas": "{template["pas_region_up"]}-{template["pas_region_down"]}"
  "3ss": "{template["3ss_region_up"]}-{template["3ss_region_down"]}"
  "5ss": "{template["5ss_region_up"]}-{template["5ss_region_down"]}"
CSM_regions_files:
  "pas": "{tandem_pas_representative_sites_coordinates}"
  "3ss": "{_3ss_coordinates}"
  "5ss": "{_5ss_coordinates}"
CSM_window_step: {int(template["window_size"]/2)}
CSM_matrix_type: "{template["matrix_type"]}"
CSM_kmer_min: {template["k_min"]}
CSM_kmer_max: {template["k_max"]}
CSM_pwm_directory: "{template["PWM_directory"]}"
CSM_MotEvo_bg_binding_prior: {template["MotEvo_bg_binding_prior"]}
CSM_MotEvo_min_binding_posterior: {template["MotEvo_min_binding_posterior"]}
CSM_MotEvo_Markov_chain_order: {template["MotEvo_Markov_chain_order"]}

### module: QUANTIFICATION
QEI_scripts_dir: "{template["MAPP_directory"]}/modules/QUANTIFICATION/scripts"
QEI_outdir: "{template["MAPP_directory"]}/modules/QUANTIFICATION/{default_output_dir_name}"
QEI_genomic_annotation: "{template["genomic_annotation"]}"
QEI_genomic_sequence: "{template["genomic_sequence"]}"
QEI_alignment_files:
{QEI_sample_bam_dict}
QEI_design_table: "{updated_design_table}"
QEI_exon_set: "{skipped_exon_events}"

### module: PAQR
PAQ_scripts_dir: "{template["MAPP_directory"]}/modules/PAQR/scripts"
PAQ_outdir: "{template["MAPP_directory"]}/modules/PAQR/{default_output_dir_name}"
PAQ_tandem_pas: "{tandem_pas_coordinates}"
PAQ_mapped_samples:
{PAQ_sample_bam_dict}
PAQ_mapped_samples_indices:
{PAQ_sample_bai_dict}
PAQ_mapped_samples_conditions:
{PAQ_sample_condition_dict}
PAQ_design_file: "{updated_design_table}"
PAQ_coverage_unstranded: "{unstranded_string}"
PAQ_coverage_downstream_extension: {template["coverage_extension"]}
PAQ_min_distance_start_to_proximal: {template["cvg_start2prox_minDist"]}
PAQ_read_length: {template["readLength"]}
PAQ_min_length_mean_coverage: {template["relUse_minLength_meanCvg"]}
PAQ_min_mean_exon_coverage: {template["relUse_minMeanCvg_perSample"]}
PAQ_distal_downstream_extension: {template["relUse_distal_ds"]}
PAQ_max_mean_coverage: {template["relUse_distal_ds_maxCvg"]}
PAQ_cluster_distance: {template["relUse_min_cluster_distance"]}
PAQ_upstream_cluster_extension: {template["relUse_us_reg_for_best_breakPoint"]}
PAQ_coverage_mse_ratio_limit: {template["relUse_mse_ratio_threshold"]}
PAQ_expression_pseudocount: {template["pas_expression_pseudocount"]}

### module: MAEI
MAE_scripts_dir: "{template["MAPP_directory"]}/modules/MAEI/scripts"
MAE_outdir: "{template["MAPP_directory"]}/modules/MAEI/{default_output_dir_name}"
MAE_inclusion_table: "{template["MAPP_directory"]}/modules/QUANTIFICATION/{default_output_dir_name}/inclusion_table.tsv"
MAE_design_file: "{updated_design_table}"
MAE_sitecount_matrices:
  "3ss":
{MAE_windows_3ss}
  "5ss":
{MAE_windows_5ss}
MAE_min_motif_fraction: {template["min_motif_fraction"]}
MAE_min_expression: {template["min_transcript_expression"]}
MAE_sorting_strategy: "{template["sorting_strategy"]}"

### module: KAPAC
KPC_scripts_dir: "{template["MAPP_directory"]}/modules/KAPAC/scripts"
KPC_outdir: "{template["MAPP_directory"]}/modules/KAPAC/{default_output_dir_name}"
KPC_tandem_polyA_sites: "{template["MAPP_directory"]}/modules/PAQR/{default_output_dir_name}/tandem_pas_expression_normalized_pseudocount.tsv"
KPC_design_file: "{updated_design_table}"
KPC_sitecount_matrices:
{KPC_windows_pas}
KPC_window_size: {template["window_size"]}
KPC_min_motif_fraction: {template["min_motif_fraction"]}
KPC_row_center_expression: True
KPC_consider_only_excess_counts: False
KPC_nr_random_runs: 0
KPC_upstream_window: {template["unique_region_upstream_pas"]}
KPC_downstream_window: {template["unique_region_downstream_pas"]}
KPC_sorting_strategy: "{template["sorting_strategy"]}"

### module: REPORT_RESULTS
RES_scripts_dir: "{template["MAPP_directory"]}/modules/REPORT_RESULTS/scripts"
RES_outdir: "{template["MAPP_directory"]}/modules/REPORT_RESULTS/{default_output_dir_name}"
RES_splicing_results_table_3ss: "{template["MAPP_directory"]}/modules/MAEI/{default_output_dir_name}/3ss/collapsed_results.tsv"
RES_splicing_results_table_5ss: "{template["MAPP_directory"]}/modules/MAEI/{default_output_dir_name}/5ss/collapsed_results.tsv"
RES_polyadenylation_results_table: "{template["MAPP_directory"]}/modules/KAPAC/{default_output_dir_name}/model_results/collapsed_results.tsv"
RES_design_file: "{updated_design_table}"
RES_max_pval: {template["max_pval"]}
RES_sorting_strategy: "{template["sorting_strategy"]}"

...
"""

    with open(options.pipeline_configfile, "w") as outfile:
        outfile.write(config_string)


##############################################################################

if __name__ == "__main__":

    try:
        # parse the command-line arguments
        options = parse_arguments().parse_args()

        # set up logging during the execution
        formatter = logging.Formatter(
            fmt="[%(asctime)s] %(levelname)s - %(message)s",
            datefmt="%d-%b-%Y %H:%M:%S",
        )
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger = logging.getLogger("logger")
        logger.setLevel(logging.getLevelName(options.verbosity))
        logger.addHandler(console_handler)
        if options.logfile is not None:
            logfile_handler = logging.handlers.RotatingFileHandler(
                options.logfile, maxBytes=50000, backupCount=2
            )
            logfile_handler.setFormatter(formatter)
            logger.addHandler(logfile_handler)

        # execute the body of the script
        start_time = time.time()
        logger.info("Starting script")
        main()
        seconds = time.time() - start_time

        # log the execution time
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        logger.info(
            "Successfully finished in {hours}h:{minutes}m:{seconds}s",
            hours=int(hours),
            minutes=int(minutes),
            seconds=int(seconds) if seconds > 1.0 else 1,
        )
    # log the exception in case it happens
    except Exception as e:
        logger.exception(str(e))
        raise e
