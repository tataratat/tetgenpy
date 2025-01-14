# Contributing
tetgenpy welcomes and appreciates discussions, issues and pull requests!

## Quick start
Once the repo is forked, one possible starting point would be creating a new python environments, for example, using [conda](https://docs.conda.io/en/latest/miniconda.html) with `python=3.9`
```bash
conda create -n tetgenpyenv python=3.9
conda activate tetgenpyenv
git clone git@github.com:<path-to-your-fork>
cd tetgenpy  # or <forkname>
git submodule update --init --recursive
git checkout -b new-feature0
pip install -e .
```

## Automatic formatting / style check
To check the format and style of your code use the following commands at tetgenpy root:
```bash
pip install pre-commit
precommit run -a
```
