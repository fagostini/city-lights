# Welcome to City-lights

City-lights is a Python package for inspecting and visualising AVITI Teton 24 data.

To check the package information and available sub-commands, you can run `cl --help`.

## Available tools

- [Inspect](inspect.md) - Extract useful metrics that can be either shown on the terminal or saved to file.
- [Render](render.md) - Extract useful metrics and use them to generate static or interactive plots.

## General options

In addition to the `--version` option, which prints the current version of the city-lights package, the following options can be provided before any sub-command (_e.g._ `cl <OPTION> <SUBCOMMAND>`), and they will be propagated to the sub-command execution:

- `--verbose`: Increase verbosity of the output by setting the logging level to `DEBUG`. This will print additional information about the execution of the sub-command.
- `--quiet`: Decrease verbosity of the output by setting the logging level to `WARNING`. This will suppress most of the output, only printing warnings and errors.

> _**Note:** By default, the logging level is set to `INFO`, which will print general information about the execution of the sub-command, but not the detailed debug information._

> _**Important:** In some cases (e.g. when the output will be printed to the standard output), the `--quiet` option will be applied automatically to avoid cluttering the output with unnecessary information. In such cases, the logging level will be set to `WARNING` regardless of the `--verbose` or `--quiet` options._
