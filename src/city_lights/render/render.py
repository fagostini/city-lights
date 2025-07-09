import argparse
import json
import logging
import pathlib
import re
import warnings

import altair as alt
import polars
import scanpy

from city_lights import utils
from city_lights.setup import setup_logging

warnings.filterwarnings("ignore", category=UserWarning, module="anndata")
logging.getLogger("anndata").setLevel(logging.WARNING)


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description="FASTA, FASTQ or plain text sequence generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Generate random nucleotide sequences in FASTA, FASTQ or plain text format",
    )
    parser.add_argument(
        "--input-path",
        type=pathlib.Path,
        help="Path to the cell2stats output folder.",
    )
    parser.add_argument(
        "--stats-json",
        type=pathlib.Path,
        help="Path to the cell2stats 'RunStats.json' file.",
    )
    parser.add_argument(
        "--panel-json",
        type=pathlib.Path,
        help="Path to the cell2stats 'Panel.json' file.",
    )
    parser.add_argument(
        "--raw-parquet",
        type=pathlib.Path,
        help="Path to the cell2stats 'RawCellStats.parquet' file.",
    )
    parser.add_argument(
        "--wells",
        type=str,
        help="""Comma-separated list of wells to include in the analysis.
        By default, all wells will be used.""",
    )
    parser.add_argument(
        "--batches",
        type=str,
        help="""Comma-separated list of batches to include in the analysis.
        By default, all batches will be used.""",
    )
    parser.add_argument(
        "--plot-list",
        type=str,
        help="""Comma-separated list of plots to generate. Accepted values include 'All',
        'BatchWell', 'Well', 'Count', 'Correlation', 'UMAP'. Default: All.""",
        default="all",
    )
    parser.add_argument(
        "--format", type=str, choices=["html", "png", "pdf"], default="png"
    )
    parser.add_argument(
        "--output-path",
        type=pathlib.Path,
        help="Path to the output folder where the generated plots will be saved.",
        default="ngi_plots",
    )
    parser.set_defaults(parse=validate_args, run=main)

    return parser


def validate_args(args) -> argparse.Namespace:
    """
    Validate the command line arguments.
    """
    if not args.input_path and not (
        args.stats_json or args.panel_json or args.raw_parquet
    ):
        logging.error(
            "At least one of --input-path, --stats-json, --panel-json, or --raw-parquet must be provided."
        )
        raise argparse.ArgumentTypeError
    if args.input_path and not args.input_path.is_dir():
        logging.error(f"Input path '{args.input_path}' is not a directory.")
        raise argparse.ArgumentTypeError
    if not args.stats_json:
        args.stats_json = args.input_path.joinpath("RunStats.json")
    if not args.stats_json.is_file():
        logging.error(f"Stats JSON file '{args.stats_json}' does not exist.")
        raise argparse.ArgumentTypeError
    if not args.panel_json:
        args.panel_json = args.input_path.joinpath("Panel.json")
    if not args.panel_json.is_file():
        logging.error(f"Panel JSON file '{args.panel_json}' does not exist.")
        raise argparse.ArgumentTypeError
    if not args.raw_parquet:
        args.raw_parquet = args.input_path.joinpath("RawCellStats.parquet")
    if not args.raw_parquet.is_file():
        logging.error(f"Raw Parquet file '{args.raw_parquet}' does not exist.")
        raise argparse.ArgumentTypeError

    if args.wells:
        args.wells = [x.strip().upper() for x in args.wells.split(",")]
        invalid_wells = [x for x in args.wells if re.match("^[A-Z][1-2]$")]
        if invalid_wells:
            logging.error(f"Well value(s) {invalid_wells} are not valid well names!")
            raise argparse.ArgumentError

    if args.batches:
        args.batches = [x.strip().upper() for x in args.batches.split(",")]
        invalid_batches = [x for x in args.batches if re.match("^B0[1-9]$")]
        if invalid_batches:
            logging.error(
                f"Batch value(s) {invalid_batches} are not valid batch names!"
            )
            raise argparse.ArgumentError

    args.plot_list = [x.strip().lower() for x in args.plot_list.split(",")]
    accepted_values = ["all", "batchwell", "well", "count", "correlation", "umap"]
    raise_error = False
    for x in args.plot_list:
        if x not in accepted_values:
            logging.error(f"Plot value {x} cannot be found among the accepted values!")
            raise_error = True
    if raise_error:
        raise argparse.ArgumentError

    if not args.output_path.is_dir():
        args.output_path.mkdir(parents=True, exist_ok=True)
    else:
        logging.warning(
            f"Output path '{args.output_path}' already exists. Files may be overwritten."
        )

    return args


def plot_per_batch_and_well(wells, batch2well_values, label, run_name):
    return (
        alt.Chart(
            polars.from_dict(batch2well_values)
            .with_columns(polars.Series("Well", wells))
            .unpivot(index="Well", variable_name="Batch", value_name=label),
            title=run_name,
        )
        .mark_bar()
        .encode(
            x=alt.X("Well:N"),
            y=alt.Y(f"{label}:Q"),
            xOffset="Batch:N",
            color="Batch:N",
        )
        .properties(height=480, width=1200)
    )


def plot_per_well(wells, values, label, run_name):
    return (
        alt.Chart(
            polars.DataFrame({"Well": wells, f"{label}": values}),
            title=run_name,
        )
        .mark_bar()
        .encode(
            x=alt.X("Well:N"),
            y=alt.Y(f"{label}:Q"),
        )
        .properties(height=480, width=1200)
    )


def plot_correlation(data: polars.DataFrame, name: str, label: str):
    return (
        alt.Chart(
            data,
            title=f"{name} well correlation - {label}",
        )
        .mark_rect()
        .encode(
            alt.X("X:O").title(None),
            alt.Y("Label:O").title(None),
            alt.Color("Z:Q")
            .scale(scheme="lightgreyred", type="linear", zero=True)
            .title("\u03c1"),
        )
        .configure_scale(bandPaddingInner=0.1)
        .properties(height=480, width=480)
    )


def plot_umap(anndata, anntype, outname: str = None):
    # create the plot object
    if outname:
        scanpy.pl.umap(anndata, color=anntype, show=False).figure.savefig(outname)
    else:
        scanpy.pl.umap(anndata, color=anntype)


def main(args: argparse.Namespace) -> None:
    setup_logging(args)
    logging.debug("Run parameters:")
    logging.debug(f"   Input Directory: '{args.input_path}'")
    logging.debug(f"   Stats JSON File: '{args.stats_json}'")
    logging.debug(f"   Panel JSON File: '{args.panel_json}'")
    logging.debug(f"   Raw Parquet File: '{args.raw_parquet}'")
    logging.debug(f"   Wells: {args.wells}")
    logging.debug(f"   Batches: {args.batches}")
    logging.debug(f"   Plots List: {args.plot_list}")
    logging.debug(f"   Plots format: '{args.format}'")
    logging.debug(f"   Plots Output Path: '{args.output_path}'")
    logging.info("Running render module...")
    run_stats_json = args.stats_json
    panel_json = args.panel_json
    raw_cell_stats_parquet = args.raw_parquet

    with open(run_stats_json) as f:
        run_stats = json.load(f)

    with open(panel_json) as f:
        panel = json.load(f)

    wells = (
        args.wells
        if args.wells
        else sorted(
            [well_data["WellLocation"] for well_data in run_stats["CytoStats"]["Wells"]]
        )
    )
    batches = (
        args.batches
        if args.batches
        else sorted(
            [
                batch_data["BatchName"]
                for batch_data in run_stats["DemuxStats"]["Batches"]
            ]
        )
    )

    if "all" in args.plot_list or "batchwell" in args.plot_list:
        for demux_stat in ["PercentAssignedReads", "PercentMismatch"]:
            values = utils.extract_demux_stat(demux_stat, wells, run_stats)
            plot = plot_per_batch_and_well(
                wells, values, demux_stat, run_stats["RunName"]
            )
            plot.save(f"{args.output_path.joinpath(demux_stat)}.{args.format}")

    if "all" in args.plot_list or "well" in args.plot_list:
        for segmentation_metric in [
            "PercentConfluency",
            "CellCount",
            "MedianCellDiameter",
        ]:
            values = utils.extract_segmentation_metric(
                segmentation_metric, wells, run_stats
            )
            plot = plot_per_well(
                wells, values, segmentation_metric, run_stats["RunName"]
            )
            plot.save(f"{args.output_path.joinpath(segmentation_metric)}.{args.format}")

    if "all" in args.plot_list or "count" in args.plot_list:
        for cyto_stat in [
            "AssignedCountsPerMM2",
        ]:
            values = utils.extract_cyto_stat(cyto_stat, run_stats)
            plot = plot_per_batch_and_well(
                wells, values, cyto_stat, run_stats["RunName"]
            )
            plot.save(f"{args.output_path.joinpath(cyto_stat)}.{args.format}")

    if (
        "all" in args.plot_list
        or "correlation" in args.plot_list
        or "umap" in args.plot_list
    ):
        # Read and filter the raw cell stats data
        df = polars.read_parquet(raw_cell_stats_parquet)
        df = utils.filter_cells(df, batches, wells).to_pandas()

        if "all" in args.plot_list or "correlation" in args.plot_list:
            for data_type in ["RNA", "Protein"]:
                corr_distances = utils.generate_distances_2d(df, panel, data_type)
                plot = plot_correlation(corr_distances, run_stats["RunName"], data_type)
                plot.save(
                    f"{args.output_path.joinpath(data_type)}_correlation.{args.format}"
                )

        if "all" in args.plot_list or "umap" in args.plot_list:
            anndata = utils.calculate_umap(df, panel)
            for umap_type in ["WellLabel", "phase"]:
                plot_umap(
                    anndata,
                    umap_type,
                    outname=f"{args.output_path.joinpath(umap_type)}_umap.{args.format}",
                )
