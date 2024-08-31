# Makefile

This Makefile is designed to simplify and automate routine actions for developing, maintaining, and deploying Python packages.
It provides a standardized set of commands for managing your project's environment, dependencies, formatting, testing, building, and documentation.

## Prerequisites

You must have the following packages installed on your system:

-   [GNU Make](https://www.gnu.org/software/make/), and
-   [Conda](https://docs.conda.io/en/latest/).

All other dependencies will be installed within a conda environment.

## Environment variables

The Makefile for this project utilizes a set of environmental variables to streamline the build process and maintain consistency across different development environments.
These parameters, defined at the beginning of the Makefile, play crucial roles in various tasks throughout the development lifecycle.
By centralizing these configuration parameters, we've created a flexible and maintainable build system that can easily adapt to different project requirements while maintaining consistency in our development and build processes.

!!! example
    ```make
    SHELL := /usr/bin/env bash
    PYTHON_VERSION := 3.12
    PYTHON_VERSION_CONDENSED := 312
    PACKAGE_NAME := simlify
    PACKAGE_PATH := $(PACKAGE_NAME)/
    TESTS_PATH := tests/
    CONDA_NAME := $(PACKAGE_NAME)-dev
    CONDA := conda run -n $(CONDA_NAME)
    CONDA_LOCK_OPTIONS := -p linux-64 -p osx-64 --channel conda-forge
    ```

At the core of our configuration is the specification of the Python environment.
We currently target Python 3.12, as indicated by the `PYTHON_VERSION` parameter.
This version number is also represented in a condensed format (`PYTHON_VERSION_CONDENSED`) for compatibility with certain tools and naming conventions.

For testing and quality assurance, we've designated a specific directory (`TESTS_PATH`) where all test files reside.
This separation of concerns helps maintain a clean and organized project structure.

Our development environment leverages Conda for package management and environment isolation.
The `CONDA_NAME` parameter defines the name of our Conda environment, while the `CONDA` parameter provides a convenient shorthand for running commands within this environment.

To ensure reproducibility across different platforms, we use conda-lock with specific options defined in `CONDA_LOCK_OPTIONS`. These options target multiple platforms (`linux-64`, `osx-64`) and specify conda-forge as our primary channel for package sourcing.

## Environment Management

Managing our project's environment is crucial for reproducibility and consistency across different development setups.

!!! tip
    To automatically build our entire development environment using:

    ```bash
    make environment
    ```

    This command runs a series of other commands (that follow) to set up the Conda environment, install all dependencies, and configure any necessary settings.

### Conda Packages

Conda is used as the initial package and environment manager for installing non-Python or conda-only packages.
Using PyPI is often tremendously easier than conda, so you should only use pip once you finish installing conda.

The following commands help set up and manage your Conda environment.

-   **`make conda-create`**: Creates a new, minimal Conda environment for your project.
    Only Python and the `conda-lock` packages are installed.
-   **`make conda-setup`**: Installs basic conda packages needed to ensure the remaining environment build pipeline works correctly; namely `poetry`, `pre-commit`, and `conda-poetry-liaison`.
-   **`make conda-dependencies`**: Install all project-specific conda packages that are not available on PyPI.
    By default, this installs nothing.
    If you wanted to add a package, you would add a line like `$(CONDA) conda install -y conda-forge::gcc`.
-   **`make nodejs-dependencies`**: Installs Node.js dependencies.
    This is useful if your project includes JavaScript components or uses Node.js-based tools.

### Pre-commit Hooks

To set up pre-commit hooks for your project, run:

```bash
make pre-commit-install
```

This installs pre-commit hooks defined in your .pre-commit-config.yaml file, which can help maintain code quality by running checks before each commit.

### PyPI Packages

For packages not available through Conda, we use poetry within your Conda environment:

1.  Activate your Conda environment:
    ```bash
    conda activate simlify-dev
    ```
2.  Install the package using poetry:
    ```bash
    poetry add <package_name>
    ```

### Updating Dependencies

To update all dependency lock files:

```bash
make locks
```

This command updates both `conda-lock.yml` and `poetry.lock`, ensuring all dependencies are at their latest compatible versions.

## Code Formatting and Linting

Maintaining consistent code style and quality is important for project maintainability.

1.  `make formatting`: Applies automatic code formatting.
2.  `make validate`: Validates code formatting without making changes.

## Testing

To run your project's test suite:

```bash
make test
```

## Building

To build your project (e.g., creating a distributable package):

```bash
make build
```

## Cleaning

To clean up temporary files and build artifacts:

```bash
make cleanup
```

This helps maintain a clean project directory by removing generated files, caches, and other temporary data.

## Documentation

To serve your project's documentation locally:

```bash
make serve
```
