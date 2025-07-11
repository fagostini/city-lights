"""Inspect module main script."""

import argparse
import json
import logging
import pathlib
import re
import warnings

import polars

from city_lights import utils
from city_lights.setup import setup_logging

warnings.filterwarnings("ignore", category=UserWarning, module="anndata")
logging.getLogger("anndata").setLevel(logging.WARNING)


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Initialise module subparser."""
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description="Module to extract useful metrics from AVITI Teton output folder",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Extract useful metrics from AVITI Teton output folder",
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
        "--stats-list",
        type=str,
        help="""Comma-separated list of metrics to produce. Accepted values include 'All',
        'BatchWell', 'Well', 'Count', 'Correlation'. Default: All.""",
        default="all",
    )
    parser.add_argument(
        "--format", type=str, choices=["csv", "parquet", "stdout"], default="csv"
    )
    parser.add_argument(
        "--output-path",
        type=pathlib.Path,
        help="Path to the output folder where the metrics will be saved.",
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

    args.stats_list = [x.strip().lower() for x in args.stats_list.split(",")]
    accepted_values = ["all", "batchwell", "well", "count", "correlation", "umap"]
    raise_error = False
    for x in args.stats_list:
        if x not in accepted_values:
            logging.error(
                f"Metric value {x} cannot be found among the accepted values!"
            )
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


def save_dataframe(df, output, format):
    """Function to save the polars dataframe in different formats."""
    if format == "csv":
        df.write_csv(output)
    elif format == "parquet":
        df.write_parquet(output)
    else:
        print(df.to_pandas())


def main(args: argparse.Namespace) -> None:
    """Main function."""
    if args.format == "stdout":
        args.quiet = True
    setup_logging(args)
    logging.debug("Run parameters:")
    logging.debug(f"   Input Directory: '{args.input_path}'")
    logging.debug(f"   Stats JSON File: '{args.stats_json}'")
    logging.debug(f"   Panel JSON File: '{args.panel_json}'")
    logging.debug(f"   Raw Parquet File: '{args.raw_parquet}'")
    logging.debug(f"   Wells: {args.wells}")
    logging.debug(f"   Batches: {args.batches}")
    logging.debug(f"   Metrics List: {args.stats_list}")
    logging.debug(f"   Outputs Format: '{args.format}'")
    logging.debug(f"   Outputs Path: '{args.output_path}'")
    logging.info(f"Running {__name__.split('.')[-1]} module...")
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

    if "all" in args.stats_list or "batchwell" in args.stats_list:
        df_list = [
            polars.from_dict(
                utils.extract_demux_stat(demux_stat, wells, run_stats)
            ).with_columns(
                [
                    polars.Series("Well", wells),
                    polars.lit(demux_stat).alias("DemuxStats"),
                ]
            )
            for demux_stat in ["PercentAssignedReads", "PercentMismatch"]
        ]
        save_dataframe(
            polars.concat(df_list).unpivot(
                index=["Well", "DemuxStats"], variable_name="Batch", value_name="Value"
            ),
            f"{args.output_path.joinpath('DemuxStats')}.{args.format}",
            args.format,
        )

    if "all" in args.stats_list or "well" in args.stats_list:
        df = (
            polars.DataFrame(
                {
                    segmentation_metric: utils.extract_segmentation_metric(
                        segmentation_metric, wells, run_stats
                    )
                    for segmentation_metric in [
                        "PercentConfluency",
                        "CellCount",
                        "MedianCellDiameter",
                    ]
                }
            )
            .with_columns(polars.Series("Well", wells))
            .select(["Well", "PercentConfluency", "CellCount", "MedianCellDiameter"])
        )
        save_dataframe(
            df,
            f"{args.output_path.joinpath('SegmentationMetrics')}.{args.format}",
            args.format,
        )

    if "all" in args.stats_list or "count" in args.stats_list:
        df_list = [
            polars.from_dict(
                utils.extract_cyto_stat(cyto_stat, run_stats)
            ).with_columns(
                [polars.Series("Well", wells), polars.lit(cyto_stat).alias("CytoStats")]
            )
            for cyto_stat in [
                "AssignedCountsPerMM2",
            ]
        ]
    save_dataframe(
        polars.concat(df_list).unpivot(
            index=["Well", "CytoStats"], variable_name="Batch", value_name="Value"
        ),
        f"{args.output_path.joinpath('CytoStats')}.{args.format}",
        args.format,
    )

    if "all" in args.stats_list or "correlation" in args.stats_list:
        # Read and filter the raw cell stats data
        df = polars.read_parquet(raw_cell_stats_parquet)
        df = utils.filter_cells(df, batches, wells)

        if args.format in ["parquet"]:
            save_dataframe(
                df,
                f"{args.output_path.joinpath('FilteredCellsTable')}.{args.format}",
                args.format,
            )

        df = df.to_pandas()

        for data_type in ["RNA", "Protein"]:
            corr_distances = utils.generate_distances_2d(df, panel, data_type)
            corr_distances = corr_distances.join(
                corr_distances.select(["Y", "Label"]).unique(),
                left_on="X",
                right_on="Y",
                suffix="_X",
                how="left",
            ).rename({"Label": "Label_Y"})
            save_dataframe(
                corr_distances,
                f"{args.output_path.joinpath(data_type)}DistCorr.{args.format}",
                args.format,
            )
