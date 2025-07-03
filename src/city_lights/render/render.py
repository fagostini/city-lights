import argparse
import json
import logging
import pathlib
import warnings
from collections import defaultdict

import altair as alt
import numpy
import polars
import scanpy

from city_lights.setup import setup_logging
from cytoprofiling import assign_cell_phase, cytoprofiling_to_anndata

warnings.filterwarnings("ignore", category=UserWarning, module="anndata")


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
        "--output-path",
        type=pathlib.Path,
        help="Path to the output folder where the generated sequences will be saved.",
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

    if not args.output_path.parent.is_dir():
        args.output_path.parent.mkdir(parents=True, exist_ok=True)
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


def extract_demux_stat(label, wells, run_stats) -> dict:
    values = {}
    for batch_data in run_stats["DemuxStats"]["Batches"]:
        batch_values = []
        for well in wells:
            well_value = float("nan")
            well_data = [
                well_data
                for well_data in batch_data["Wells"]
                if well_data["WellLocation"] == well
            ]
            if len(well_data) == 1:
                well_value = well_data[0][label]
                if abs(float(well_value) + 999) < 0.0001:
                    well_value = float("nan")
            batch_values.append(well_value)
        values[batch_data["BatchName"]] = batch_values
    return values


def extract_segmentation_metric(label, wells, run_stats) -> list:
    values = []
    for well in wells:
        well_value = float("nan")
        well_data = [
            well_data
            for well_data in run_stats["CytoStats"]["Wells"]
            if well_data["WellLocation"] == well
        ]
        if len(well_data) == 1:
            well_value = (
                float("nan")
                if abs(well_data[0][label] + 999) < 0.0001
                else well_data[0][label]
            )
        values.append(well_value)
    return values


def extract_cyto_stat(label, run_stats) -> dict:
    values = {}
    for well_data in run_stats["CytoStats"]["Wells"]:
        for batch_data in well_data["Batches"]:
            batch_name = batch_data["BatchName"]
            if batch_name.startswith("CP"):
                continue
            batch_value = batch_data.get(label, -999)
            if batch_value + 999 < 0.0001:
                batch_value = float("nan")
            values.setdefault(batch_data["BatchName"], []).append(batch_value)
    return values


def filter_cells(
    df: polars.DataFrame,
    batch_names: list[str] = None,
    well_names: list[str] = None,
    stats=None,
) -> polars.DataFrame:
    """Perform default filtering of cells (filter by area, assigned counts, and assigned rate)

    Args:
      df : Input dataframe with cells
      batch_names : List of batches to use for normalization
      well_names : List of wells to include in normalization
      stats : Dictionary to store statistics about filtering under the key 'filter_cells'

    Returns:
      Filtered data frame
    """

    def get_wells(df: polars.DataFrame) -> list[str]:
        """Get unique well names from the DataFrame."""
        return df.get_column("Well").unique().sort().to_list()

    def get_barcoding_batches(df: polars.DataFrame) -> list[str]:
        """Get barcoding batch names from the DataFrame columns."""
        return [
            val.split(".")[1]
            for val in df.columns
            if "Nuclear" not in val and "Unassigned" in val
        ]

    def filter_cells_by_area_percentile(
        df: polars.DataFrame,
        well_names: list[str] = None,
        min_quantile: float = 0.03,
        max_quantile: float = 0.97,
    ) -> polars.DataFrame:
        """Filter cells by area percentile for each well."""
        subset = df.filter(polars.col("Well").is_in(well_names)) if well_names else df
        return (
            subset.with_columns(
                polars.col("Area")
                .quantile(min_quantile)
                .over("Well")
                .alias("MinThreshold"),
                polars.col("Area")
                .quantile(max_quantile)
                .over("Well")
                .alias("MaxThreshold"),
            )
            .filter(
                (polars.col("Area") >= polars.col("MinThreshold"))
                & (polars.col("Area") <= polars.col("MaxThreshold"))
            )
            .drop(["MinThreshold", "MaxThreshold"])
        )

    def filter_cells_by_assigned_counts_percentile(
        df: polars.DataFrame,
        batch_names: list[str] = None,
        well_names: list[str] = None,
        min_quantile: float = 0.03,
        max_quantile: float = 0.97,
    ) -> polars.DataFrame:
        """Filter cells by assigned counts percentile for each well and batch."""
        subset = (
            df.lazy()
            .unpivot(
                index=["Cell", "Well", "WellLabel", "Tile", "X", "Y", "Area"],
                variable_name="Target",
                value_name="Count",
            )
            .with_columns(polars.col("Target").str.extract(r"\.(\w+)$").alias("Batch"))
        )

        if batch_name:
            subset = subset.filter(polars.col("Batch").is_in(batch_names))

        if well_names:
            subset = subset.filter(polars.col("Well").is_in(well_names))

        selected_cells = (
            subset.filter(
                [
                    polars.col("Batch").str.contains("Unassigned").not_()
                    & polars.col("Batch").str.starts_with("Nuclear").not_()
                ]
            )
            .filter(
                [
                    polars.col("Target").str.contains("_Decoy").not_()
                    & polars.col("Target").str.starts_with("NSB").not_()
                    & polars.col("Target").str.starts_with("Unassigned").not_()
                    & polars.col("Target").str.contains("_Nuclear.").not_()
                ]
            )
            .with_columns(
                polars.col("Count")
                .sum()
                .over(["Cell", "Well", "Batch"])
                .alias("TotalCount")
            )
            .with_columns(
                polars.col("TotalCount")
                .filter(polars.col("TotalCount") > 0)
                .quantile(min_quantile)
                .over(["Batch", "Well"])
                .alias("MinThreshold"),
                polars.col("TotalCount")
                .filter(polars.col("TotalCount") > 0)
                .quantile(max_quantile)
                .over(["Batch", "Well"])
                .alias("MaxThreshold"),
            )
            .filter(
                (polars.col("TotalCount") >= polars.col("MinThreshold"))
                & (polars.col("TotalCount") <= polars.col("MaxThreshold"))
            )
            .drop(["MinThreshold", "MaxThreshold"])
            .group_by(["Cell", "Well"])
            .agg(polars.col("Batch").unique().len().alias("BatchCount"))
            .filter(polars.col("BatchCount") == polars.col("BatchCount").max())
            .sort(["Well", "Cell"])
            .collect()
        )

        return df.join(
            selected_cells,
            on=["Cell", "Well"],
            how="semi",
        )

    def filter_cells_by_assigned_rate(
        df: polars.DataFrame,
        batch_names: list[str] = None,
        well_names: list[str] = None,
        min_rate: float = 0.5,
    ) -> polars.DataFrame:
        """Filter cells by assigned rate for each well and batch."""
        subset = (
            df.lazy()
            .unpivot(
                index=["Cell", "Well", "WellLabel", "Tile", "X", "Y", "Area"],
                variable_name="Target",
                value_name="Count",
            )
            .with_columns(polars.col("Target").str.extract(r"\.(\w+)$").alias("Batch"))
        )

        if batch_name:
            subset = subset.filter(polars.col("Batch").is_in(batch_names))

        if well_names:
            subset = subset.filter(polars.col("Well").is_in(well_names))

        unassigned = subset.filter(
            [
                polars.col("Target").str.contains("Unassigned")
                & polars.col("Target").str.contains("Unassigned_").not_()
            ]
        ).select(["Cell", "Well", "Batch", "Count"])

        selected_cells = (
            subset.filter(
                [
                    polars.col("Batch").str.contains("Unassigned").not_()
                    & polars.col("Batch").str.starts_with("Nuclear").not_()
                ]
            )
            .filter(
                [
                    polars.col("Target").str.contains("_Decoy").not_()
                    & polars.col("Target").str.starts_with("NSB").not_()
                    & polars.col("Target").str.starts_with("Unassigned").not_()
                    & polars.col("Target").str.contains("_Nuclear.").not_()
                ]
            )
            .with_columns(
                polars.col("Count")
                .sum()
                .over(["Cell", "Well", "Batch"])
                .alias("TotalCount")
            )
            .join(unassigned, on=["Cell", "Well", "Batch"], how="inner")
            .filter(
                (polars.col("TotalCount") > 0)
                & (
                    polars.col("TotalCount")
                    / (polars.col("TotalCount") + polars.col("Count"))
                    > min_rate
                )
            )
            .select(["Cell", "Well"])
            .unique()
            .collect()
        )

        return df.join(
            selected_cells,
            on=["Cell", "Well"],
            how="semi",
        )

    if well_names is None:
        well_names = get_wells(df)
        print(f"Using wells: {well_names}")
    if batch_names is None:
        batch_names = get_barcoding_batches(df)
        print(f"Using batches: {batch_names}")
    if stats is None:
        stats = df.group_by("Well").agg(polars.len().alias("Total")).sort("Well")

    result = filter_cells_by_area_percentile(df, well_names)
    stats = stats.join(
        result.group_by("Well").agg(polars.len().alias("Area")), on="Well", how="left"
    )

    result = filter_cells_by_assigned_counts_percentile(result, batch_names, well_names)
    stats = stats.join(
        result.group_by("Well").agg(polars.len().alias("AssignedCounts")),
        on="Well",
        how="left",
    )

    result = filter_cells_by_assigned_rate(result, batch_names, well_names, 0.5)
    stats = stats.join(
        result.group_by("Well").agg(polars.len().alias("Passing")),
        on="Well",
        how="left",
    )

    stats = (
        stats.with_columns(
            (polars.col("AssignedCounts") - polars.col("Passing")).alias("AssignedRate")
        )
        .with_columns(
            (polars.col("Area") - polars.col("AssignedCounts")).alias("AssignedCounts")
        )
        .with_columns((polars.col("Total") - polars.col("Area")).alias("Area"))
    )

    stats = stats.select(["Area", "AssignedCounts", "AssignedRate", "Passing"])

    return result


def generate_distances_2d(df, panel, data_type):
    # Convert dataframe to anndata
    adata = cytoprofiling_to_anndata(df, panel)

    # filter data columns to only include simple counts for data_type
    adata = adata[
        :,
        (~adata.var["is_unassigned"])
        & (~adata.var["is_nuclear"])
        & numpy.isin(
            adata.var["measurement_type"],
            [
                data_type,
            ],
        ),
    ]

    wells_x = adata.obs["Well"].unique()

    # average anndata object adata by well observation
    grouped_adata = scanpy.get.aggregate(adata, func="mean", by="Well")

    well2label_x = {}
    for well in wells_x:
        well2label_x[well] = df["WellLabel"][df["Well"] == well].unique()[0]

    distances_2d = []
    distances = defaultdict(list)
    for well_idx in range(len(wells_x)):
        current_distances = []
        for well_idy in range(len(wells_x)):
            well_x = wells_x[well_idx]
            well_y = wells_x[well_idy]

            r2 = (
                numpy.corrcoef(
                    numpy.log2(
                        grouped_adata[grouped_adata.obs["Well"] == well_x]
                        .layers["mean"]
                        .flatten()
                    ),
                    numpy.log2(
                        grouped_adata[grouped_adata.obs["Well"] == well_y]
                        .layers["mean"]
                        .flatten()
                    ),
                )[0][1]
                ** 2
            )
            current_distances.append(r2)
            if well_idy == well_idx:
                r2 = 1.0
            distances.setdefault("X", []).append(well_x)
            distances.setdefault("Y", []).append(well_y)
            distances.setdefault("Z", []).append(r2)
        distances_2d.append(current_distances)
    return (
        polars.DataFrame(distances)
        .with_columns(
            polars.Series(
                "Label",
                [
                    df["WellLabel"][df["Well"] == well].unique()[0]
                    for well in distances["Y"]
                ],
            ),
        )
        .with_columns(
            polars.concat_str(
                [polars.col("Y"), polars.col("Label")], separator=" "
            ).alias("Label")
        )
    )


def plot_correlation(data: polars.DataFrame, name: str, label: str):
    return (
        alt.Chart(
            data,
            title=f"{name} well correlation - {label}",
        )
        .mark_rect(grid=False)
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


def calculate_and_plot_umap(df: polars.DataFrame, panel: dict):
    # Convert dataframe to anndata
    adata = cytoprofiling_to_anndata(df, panel)

    # filter data columns to only include simple counts for protein and RNA
    adata = adata[
        :,
        (~adata.var["is_unassigned"])
        & (~adata.var["is_nuclear"])
        & numpy.isin(adata.var["measurement_type"], ["RNA", "Protein"]),
    ]

    # convert column names to gene names and remove any resulting duplicates
    adata.var_names = adata.var["gene"]
    adata = adata[:, ~adata.var_names.duplicated()].copy()

    # do processing of data to prepare for UMAP and cell cycle determination
    n_comps = 10
    scanpy.pp.normalize_total(adata, target_sum=1e4)
    scanpy.pp.log1p(adata)
    scanpy.tl.pca(adata, n_comps=n_comps)
    scanpy.pp.neighbors(adata, n_pcs=n_comps)

    # assign cell phase
    assign_cell_phase(adata)

    # calculate UMAP
    scanpy.tl.umap(adata)

    # plot UMAP with well label
    scanpy.pl.umap(adata, color="WellLabel")

    # plot UMAP with calculated cell phase
    scanpy.pl.umap(adata, color="phase")


def main(args: argparse.Namespace) -> None:
    setup_logging(args)
    logging.info("Running render module...")
    run_stats_json = args.stats_json
    panel_json = args.panel_json
    raw_cell_stats_parquet = args.raw_parquet

    with open(run_stats_json) as f:
        run_stats = json.load(f)

    with open(panel_json) as f:
        panel = json.load(f)

    wells = sorted(
        [well_data["WellLocation"] for well_data in run_stats["CytoStats"]["Wells"]]
    )
    batches = sorted(
        [batch_data["BatchName"] for batch_data in run_stats["DemuxStats"]["Batches"]]
    )

    for demux_stat in ["PercentAssignedReads", "PercentMismatch"]:
        values = extract_demux_stat(demux_stat, wells, run_stats)
        plot_per_batch_and_well(wells, values, demux_stat, run_stats["RunName"])

    for segmentation_metric in [
        "PercentConfluency",
        "CellCount",
        "MedianCellDiameter",
    ]:
        values = extract_segmentation_metric(segmentation_metric, wells, run_stats)
        plot_per_well(wells, values, segmentation_metric, run_stats["RunName"])

    for cyto_stat in [
        "AssignedCountsPerMM2",
    ]:
        values = extract_cyto_stat(cyto_stat, run_stats)
        plot_per_batch_and_well(wells, values, cyto_stat, run_stats["RunName"])

    # Read and filter the raw cell stats data
    df = polars.read_parquet(raw_cell_stats_parquet)
    df = filter_cells(df).to_pandas()

    for data_type in ["RNA", "Protein"]:
        corr_distances = generate_distances_2d(df, panel, data_type)
        plot_correlation(corr_distances, run_stats["RunName"], data_type)

    calculate_and_plot_umap(df, panel)
