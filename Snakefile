##############################################################################
#
#   Snakemake pipeline: MAPP
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 30-06-2020
#   LICENSE: Apache_2.0
#
##############################################################################

# imports
import sys
import os

# local rules
localrules: all

# create logging directory for the top-level rules of MAPP summary
os.makedirs(
    os.path.join(
        config["MAPP_pipeline_directory"],
        "logs",
    ),
    exist_ok = True
)

##############################################################################
### Target rule with final output of the pipeline
##############################################################################

rule all:
    """
    Gathering all output
    """
    input:
        TARGZ_summary = expand(
            os.path.join(
                "{MAPP_pipeline_directory}",
                "summary.tar.gz"
            ),
            MAPP_pipeline_directory = config["MAPP_pipeline_directory"]
        )

##############################################################################
### Include all MAPP modules:
##############################################################################

include: "modules/PREPROCESSING/Snakefile"
include: "modules/CREATE_SITECOUNT_MATRICES/Snakefile"
include: "modules/PREPARE_TANDEM_PAS/Snakefile"
include: "modules/PAQR/Snakefile"
include: "modules/KAPAC/Snakefile"
include: "modules/EXTRACT_AS_EXONS/Snakefile"
include: "modules/QUANTIFICATION/Snakefile"
include: "modules/MAEI/Snakefile"
include: "modules/REPORT_RESULTS/Snakefile"

##############################################################################
### Collect files for summary:
##############################################################################

rule MAPP_collect_summary:
    """
    Copy relevant text files and Motif Activity Heatmaps
    """
    input:
        DIR_Zscores_heatmaps = expand(
            os.path.join(
                "{RES_output_dir}",
                "Zscores_heatmaps"
            ),
            RES_output_dir = config["RES_outdir"]
        )

    output:
        DIR_pipeline_summary = directory(
            os.path.join(
                "{MAPP_pipeline_directory}",
                "summary"
            )
        )

    params:
        YML_pipeline_config = config["MAPP_pipeline_configfile"],
        BED_3ss_list = config["CSM_regions_files"]["3ss"],
        BED_5ss_list = config["CSM_regions_files"]["5ss"],
        BED_pas_list = config["CSM_regions_files"]["pas"],
        DIR_PAQ_module_outdir = config["PAQ_outdir"],
        DIR_QEI_module_outdir = config["QEI_outdir"],
        DIR_RES_module_outdir = config["RES_outdir"],
        STR_pwm_seqlogos_copy_command = lambda wildcards: \
            "" if config["CSM_matrix_type"] == "kmers" \
            else \
                "cp -r " + config["MAPP_seqlogos_directory"] + " " + os.path.join(
                    wildcards.MAPP_pipeline_directory,
                    "summary",
                    "sequence-logos"
                ),
        LOG_cluster_log = os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_collect_summary.log"
        )

    threads: 1

    log:
        LOG_local_stdout = os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_collect_summary.stdout.log"
        ),
        LOG_local_stderr = os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_collect_summary.stderr.log"
        )

    benchmark:
        os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_collect_summary.benchmark.log"
        )

    conda:
        "env/bash.yml"

    singularity:
        "docker://bash:4.4.18"

    shell:
        """
        (mkdir {output.DIR_pipeline_summary}
        mkdir {output.DIR_pipeline_summary}/results
        cp \
            {params.YML_pipeline_config} \
            {output.DIR_pipeline_summary}/pipeline-config.yml
        cp \
            {params.BED_3ss_list} \
            {output.DIR_pipeline_summary}/results/coordinates-3ss.bed
        cp \
            {params.BED_5ss_list} \
            {output.DIR_pipeline_summary}/results/coordinates-5ss.bed
        cp \
            {params.BED_pas_list} \
            {output.DIR_pipeline_summary}/results/coordinates-pas.bed
        cp \
            {params.DIR_PAQ_module_outdir}/filtered_pas_expression.tsv \
            {output.DIR_pipeline_summary}/results/pas-expression-table.tsv
        cp \
            {params.DIR_PAQ_module_outdir}/relative_pas_positions.tsv \
            {output.DIR_pipeline_summary}/results/pas-relative-positions.tsv
        cp \
            {params.DIR_QEI_module_outdir}/filtered_merged_quantification_genes.tsv \
            {output.DIR_pipeline_summary}/results/quantification-genes.tsv
        cp \
            {params.DIR_QEI_module_outdir}/filtered_merged_quantification_transcripts.tsv \
            {output.DIR_pipeline_summary}/results/quantification-transcripts.tsv
        cp \
            {params.DIR_QEI_module_outdir}/inclusion_table.tsv \
            {output.DIR_pipeline_summary}/results/exon-inclusion-table.tsv
        cp \
            {params.DIR_RES_module_outdir}/results_updated_3ss.tsv \
            {output.DIR_pipeline_summary}/results/results-3ss.tsv
        cp \
            {params.DIR_RES_module_outdir}/results_updated_5ss.tsv \
            {output.DIR_pipeline_summary}/results/results-5ss.tsv
        cp \
            {params.DIR_RES_module_outdir}/results_updated_pas.tsv \
            {output.DIR_pipeline_summary}/results/results-pas.tsv
        cp \
            {params.DIR_RES_module_outdir}/motifs_3ss.txt \
            {output.DIR_pipeline_summary}/results/motifs-list-3ss.txt
        cp \
            {params.DIR_RES_module_outdir}/motifs_5ss.txt \
            {output.DIR_pipeline_summary}/results/motifs-list-5ss.txt
        cp \
            {params.DIR_RES_module_outdir}/motifs_pas.txt \
            {output.DIR_pipeline_summary}/results/motifs-list-pas.txt
        cp -r \
            {params.DIR_RES_module_outdir}/Zscores_heatmaps \
            {output.DIR_pipeline_summary}/results/heatmaps
        {params.STR_pwm_seqlogos_copy_command})
        1> {log.LOG_local_stdout} 2> {log.LOG_local_stderr}
        """

##############################################################################
### Prepare final summary table for the report in TSV format
##############################################################################

rule MAPP_prepare_summary_table:
    """
    Prepare the main table with results for the HTML report
    in text format first
    """
    input:
        DIR_pipeline_summary = os.path.join(
            "{MAPP_pipeline_directory}",
            "summary"
        ),
        SCRIPT_ = os.path.join(
            "{MAPP_pipeline_directory}",
            "scripts",
            "prepare-final-summary-table.py"
        )

    output:
        TSV_summary_table = os.path.join(
            "{MAPP_pipeline_directory}",
            "summary",
            "main-table.tsv"
        )

    params:
        STR_ranking_score = "Avg",
        LOG_cluster_log = os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_prepare_summary_table.log"
        )

    threads: 1

    log:
        LOG_local_stdout = os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_prepare_summary_table.stdout.log"
        ),
        LOG_local_stderr = os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_prepare_summary_table.stderr.log"
        )

    benchmark:
        os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_prepare_summary_table.benchmark.log"
        )

    conda:
        "env/python.yml"

    singularity:
        "docker://zavolab/mapp_base_python:1.1.1"

    shell:
        """
        python {input.SCRIPT_} \
        --ranking-score {params.STR_ranking_score} \
        --summary-directory {input.DIR_pipeline_summary} \
        1> {log.LOG_local_stdout} 2> {log.LOG_local_stderr}
        """

##############################################################################
### Generate MAPP report:
##############################################################################

rule MAPP_generate_report:
    """
    Create a simple HTML report with a summary of MAPP run
    """
    input:
        TSV_summary_table = os.path.join(
            "{MAPP_pipeline_directory}",
            "summary",
            "main-table.tsv"
        ),
        SCRIPT_ = os.path.join(
            "{MAPP_pipeline_directory}",
            "scripts",
            "generate-final-report.py"
        )

    output:
        HTML_report = os.path.join(
            "{MAPP_pipeline_directory}",
            "summary",
            "report.html"
        )

    params:
        DIR_pipeline_summary = os.path.join(
            "{MAPP_pipeline_directory}",
            "summary"
        ),
        DIR_resources = os.path.join(
            "{MAPP_pipeline_directory}",
            "resources"
        ),
        LOG_cluster_log = os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_generate_report.log"
        )

    threads: 1

    log:
        LOG_local_stdout = os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_generate_report.stdout.log"
        ),
        LOG_local_stderr = os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_generate_report.stderr.log"
        )

    benchmark:
        os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_generate_report.benchmark.log"
        )

    conda:
        "env/python.yml"

    singularity:
        "docker://zavolab/mapp_base_python:1.1.1"

    shell:
        """
        python {input.SCRIPT_} \
        --summary-directory {params.DIR_pipeline_summary} \
        --resources-directory {params.DIR_resources} \
        1> {log.LOG_local_stdout} 2> {log.LOG_local_stderr}
        """

##############################################################################
### Summary Tarball:
##############################################################################

rule MAPP_summary_tarball:
    """
    Archive and compress the whole summary directory
    """
    input:
        HTML_report = os.path.join(
            "{MAPP_pipeline_directory}",
            "summary",
            "report.html"
        )

    output:
        TARGZ_summary = os.path.join(
            "{MAPP_pipeline_directory}",
            "summary.tar.gz"
        )

    params:
        DIR_pipeline_directory = "{MAPP_pipeline_directory}",
        LOG_cluster_log = os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_summary_tarball.log"
        )

    threads: 1

    log:
        LOG_local_stdout = os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_summary_tarball.stdout.log"
        ),
        LOG_local_stderr = os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_summary_tarball.stderr.log"
        )

    benchmark:
        os.path.join(
            "{MAPP_pipeline_directory}",
            "logs",
            "MAPP_summary_tarball.benchmark.log"
        )

    conda:
        "env/tar.yml"

    singularity:
        "docker://ubuntu:20.04"

    shell:
        """
        tar \
        --exclude=.snakemake_timestamp \
        -czvf \
        {output.TARGZ_summary} \
        -C {params.DIR_pipeline_directory} \
        summary \
        1> {log.LOG_local_stdout} 2> {log.LOG_local_stderr}
        """
