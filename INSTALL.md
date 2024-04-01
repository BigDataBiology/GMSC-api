# Install

## Install dependencies

```bash
conda create -n gmsc-api python=3.11
conda activate gmsc-api
conda install -c default -c conda-forge flask numpy pandas jug requests python-xz
```

## Download data

```bash
jug execute download-data.py
```

## Create indices

```bash
jug execute make-indices.py
```

