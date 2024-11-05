# Docking

## Description
Collection of tools for docking and some extras.
The tool requires CCDC python API. Parts work with Schrödinger python API or are written in MOE SVL, though not having those will not make the package unusable.

## Project setup
First create the python package by running the command:
```
python -m build
```

Then move to `dist` subfolder and run:
```
pip install docking-0.0.2-py3-none-any.whl
```

Package requirements:
- bipython
- numpy
- pandas
- ccdc
- MOE for using the included SVL scripts
- optionally schrodinger
