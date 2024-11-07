depmap_omics_upload
===

This repo contains a Python module for various DepMap post-workflow processing steps.

# Installation

1. Install the required system dependencies:

    - [pyenv](https://github.com/pyenv/pyenv)
    - [Poetry](https://python-poetry.org/)

2. Install the required Python version (3.9.18):

   ```bash
   pyenv install "$(cat .python-version)"
   ```

3. Confirm that `python` maps to the correct version:

   ```
   python --version
   ```

4. Set the Poetry interpreter and install the Python dependencies:

   ```bash
   poetry env use "$(pyenv which python)"
   poetry install
   ```
