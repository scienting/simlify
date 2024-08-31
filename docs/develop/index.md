# Developer guide

Welcome to the **simlify** codebase!
This document will help you navigate the code, set up your development environment, and understand the build process.

## Contributing

We welcome contributions from the community!
If you'd like to contribute to simlify, please follow these steps:

1.  Fork the repository.
2.  Create a new branch (`git checkout -b feature-branch`).
3.  Make your changes.
4.  Commit your changes (`git commit -m 'Add new feature'`).
5.  Push to the branch (`git push origin feature-branch`).
6.  Open a pull request.

Please make sure to follow the Code of Conduct guidelines.

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
