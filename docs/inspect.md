# Inspect

This module performs the extraction of useful metrics that can be either shown on the standard output or saved them to file, in CSV or Parquet format.

## Usage

The minimal command is `cl inspect --input-path <PATH> --output-path <PATH>`, where `input-path` is the path to the AVITI data folder and output-path is the directory where the outputs will be placed.

## Options

- `--stats-json`: Path to the cell2stats 'RunStats.json' file. Required only if different from `<INPUT_PATH>/RunStats.json`.
- `--panel-json`: Path to the cell2stats 'Panel.json' file. Required only if different from `<INPUT_PATH>/Panel.json`.
- `--raw-parquet`: Path to the cell2stats 'RawCellStats.parquet' file. Required only if different from `<INPUT_PATH>/RawCellStats.parquet`.
- `--wells`: Comma-separated list of wells to include in the analysis. By default, all wells will be used.
- `--batches`: Comma-separated list of batches to include in the analysis. By default, all batches will be used.
- `--stats-list`: Comma-separated list of metrics to produce. Accepted values include 'All', 'BatchWell', 'Well', 'Count', 'Correlation'. Default: All.
- `--format`: It can be one of 'csv', 'parquet' or 'stdout'. Default is 'csv'.
