<h1 align="center">simlify</h1>
<h4 align="center">Simplify your molecular simulation workflow.</h4>
<p align="center">
    <a href="https://gitlab.com/oasci/software/simlify/-/pipelines">
        <img src="https://gitlab.com/oasci/software/simlify/badges/main/pipeline.svg" alt="Build Status ">
    </a>
    <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/simlify">
    <a href="https://codecov.io/gl/oasci:software/simlify">
        <img src="https://codecov.io/gl/oasci:software/simlify/graph/badge.svg?token=KVGB7NU117" alt="codecov">
    </a>
    <a href="https://github.com/oasci/simlify/releases">
        <img src="https://img.shields.io/github/v/release/oasci/simlify" alt="GitHub release (latest by date)">
    </a>
    <a href="https://github.com/oasci/simlify/blob/main/LICENSE" target="_blank">
        <img src="https://img.shields.io/github/license/oasci/simlify" alt="License">
    </a>
    <a href="https://github.com/oasci/simlify/" target="_blank">
        <img src="https://img.shields.io/github/repo-size/oasci/simlify" alt="GitHub repo size">
    </a>
    <a href="https://github.com/psf/black" target="_blank">
        <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Black style">
    </a>
    <a href="https://github.com/PyCQA/pylint" target="_blank">
        <img src="https://img.shields.io/badge/linting-pylint-yellowgreen" alt="Black style">
    </a>
    <a href="https://github.com/astral-sh/ruff" target="_blank">
        <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json" alt="Black style">
    </a>
</p>
<h4 align="center"><a href="https://simlify.oasci.org">Documentation</a></h4>

Add information about simlify here.

## Deploying

We use [bump-my-version](https://github.com/callowayproject/bump-my-version) to release a new version.
This will create a git tag that is used by [poetry-dynamic-version](https://github.com/mtkennerly/poetry-dynamic-versioning) to generate version strings and update `CHANGELOG.md`.

For example, to bump the `minor` version you would run the following command.

```bash
poetry run bump-my-version bump minor
```

After releasing a new version, you need to push and include all tags.

```bash
git push --follow-tags
```

## License

This project is released under the Apache-2.0 License as specified in `LICENSE.md`.
