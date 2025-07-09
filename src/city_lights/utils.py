from collections import defaultdict

import numpy
import polars
import scanpy

from cytoprofiling import assign_cell_phase, cytoprofiling_to_anndata


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

        if batch_names:
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

        if batch_names:
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


def calculate_umap(df: polars.DataFrame, panel: dict):
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

    return adata
