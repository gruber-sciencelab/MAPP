"""
##############################################################################
#
#   Generate final MAPP report in html with jinja2 templating mechanism.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 03-06-2021
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
import shutil
import yaml
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
        "--summary-directory",
        dest="summary_directory",
        required=True,
        help="Path to the directory with MAPP summary.",
    )
    parser.add_argument(
        "--resources-directory",
        dest="resources_directory",
        required=True,
        help="Path to the internal directory with resources for the report.",
    )
    return parser


##############################################################################


def main():
    """Main body of the script."""

    # read in MAPP config file
    with open(os.path.join(options.summary_directory, "pipeline-config.yml")) as f:
        configfile = yaml.safe_load(f)

    # load report template
    with open(
        os.path.join(options.resources_directory, "template.html")
    ) as template_html_file:
        template_html = template_html_file.read().splitlines()

    # load the main summary table
    summary = pd.read_csv(
        os.path.join(options.summary_directory, "main-table.tsv"), sep="\t"
    )
    summary = summary.round(3).fillna("NA")
    print(summary)

    # set distinct columns width parameters for the table
    # depending on if we run on kmers or pwms
    if configfile["CSM_matrix_type"] == "kmers":
        text_col_width = "10"
        heatmap_col_width = "50"
        heatmap_img_scaling = "30"
    else:
        assert configfile["CSM_matrix_type"] == "pwms"
        text_col_width = "10"
        heatmap_col_width = "40"
        seqlogo_col_width = "10"
        heatmap_img_scaling = "30"
        seqlogo_img_scaling = "50"

    # render the HTML report template MANUALLY (for now)
    rendered_html = []
    for _ in range(11):
        rendered_html.append(template_html[_])
    rendered_html[8] = "<h3>" + configfile["MAPP_analysis_name"] + "</h3>"
    rendered_html.append("<table>")
    rendered_html.append("  <colgroup>")
    if configfile["CSM_matrix_type"] == "pwms":
        rendered_html.append('    <col style="width:' + seqlogo_col_width + '%">')
    rendered_html.append('    <col style="width:' + text_col_width + '%">')
    rendered_html.append('    <col style="width:' + text_col_width + '%">')
    rendered_html.append('    <col style="width:' + text_col_width + '%">')
    rendered_html.append('    <col style="width:' + text_col_width + '%">')
    rendered_html.append('    <col style="width:' + text_col_width + '%">')
    rendered_html.append('    <col style="width:' + heatmap_col_width + '%">')
    rendered_html.append("  </colgroup>")
    rendered_html.append("  <tr>")
    for h in summary.columns.values:
        rendered_html.append("    <th>" + h + "</th>")
    rendered_html.append("  </tr>")
    for i, row in summary.iterrows():
        rendered_html.append("  <tr>")
        if configfile["CSM_matrix_type"] == "pwms":
            rendered_html.append(
                '    <td style="text-align:center"><a href='
                + row["Sequence Logos"]
                + "><img src="
                + row["Sequence Logos"]
                + ' width="'
                + seqlogo_img_scaling
                + '%"></a></td>'
            )
        rendered_html.append(
            '    <td style="text-align:center">' + row["motif"] + "</td>"
        )
        rendered_html.append(
            '    <td style="text-align:center">' + str(row["Ranking score"]) + "</td>"
        )
        rendered_html.append(
            '    <td style="text-align:center">' + str(row["3'ss Zscore"]) + "</td>"
        )
        rendered_html.append(
            '    <td style="text-align:center">' + str(row["5'ss Zscore"]) + "</td>"
        )
        rendered_html.append(
            '    <td style="text-align:center">' + str(row["pas Zscore"]) + "</td>"
        )
        rendered_html.append(
            '    <td style="text-align:center"><a href='
            + row["Activity Map"]
            + "><img src="
            + row["Activity Map"]
            + ' width="'
            + heatmap_img_scaling
            + '%"></a></td>'
        )
        rendered_html.append("  </tr>")
    rendered_html.append("</table>")
    for _ in range(12, len(template_html)):
        rendered_html.append(template_html[_])

    with open(os.path.join(options.summary_directory, "report.html"), "w") as report:
        for line in rendered_html:
            report.write(line + "\n")

    # include font-awesome in the report dir
    shutil.copytree(
        os.path.join(options.resources_directory, "font-awesome-4.7.0"),
        os.path.join(options.summary_directory, "font-awesome-4.7.0"),
    )


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
