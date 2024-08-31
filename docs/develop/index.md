# Developer guide

Welcome to the **simlify** codebase!
This document will help you navigate the code, set up your development environment, and understand the build process.

## Directory structure

-   **simlify/**: The main source code of the project.
-   **tests/**: Unit and integration tests for the project.
-   **docs/**: Documentation files.
-   **.github/**: GitHub-specific files for workflows and actions.

## Important files

-   **Makefile**: Contains commands to build, test, and clean the project.
-   **pyproject.toml**: Project configuration file for building and packaging.
-   **mkdocs.yml**: Configuration file for generating documentation with MkDocs.

## DevOps

We use a [Makefile](./makefile.md) to automate our DevOps process.
This fill will automate the many stages of development such as running tests, building documentation, installing dependencies, etc.
