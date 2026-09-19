# Sophiark Deployment Package

This directory is a deployment copy of Sophiark. The source workspace is not modified.

## Included

* Supported Streamlit entry point: `app.py`
* Application packages, UI assets, local browser libraries, and configuration
* Identity, OmniPath, evidence, structural-context, and interpretation lookups
* Human and mouse baseline annotation seeds required by the frozen engines
* Local ENSP-symbol and enrichment caches required for bounded/offline startup
* Tests and acceptance scripts for pre-release verification

The runtime database package is distributed separately through Zenodo.

**Sophiark Database v1.0**
DOI: https://doi.org/10.5281/zenodo.22849093

Large build-only sources are intentionally excluded, including:

* raw BindingDB
* AlphaFold archives and structures
* BioGRID archives
* raw STRING build inputs
* staging outputs
* previous analysis outputs
* virtual environments

## Local Installation

Clone the Sophiark repository and download the corresponding database release from Zenodo.

After extracting the database package, the project should contain the required `data/` directory.

Run the local verification checks:

```powershell
python scripts/data_health.py
python scripts/verify_cftr_lung_engine_flow.py
python -m pytest -q
streamlit run app.py
```

## Docker

The Docker image contains the application but does not contain the runtime database package.

Mount the local `data/` directory into the container as a read-only volume. Docker Compose uses persistent named volumes for human and mouse outputs.

`docker-entrypoint.sh` initializes a new empty output volume with the two frozen-engine baseline annotation tables. It never overwrites an existing output volume.

The entrypoint also seeds the local symbol/enrichment cache so normal analysis does not fan out into thousands of MyGene lookups.

```bash
docker compose build
docker compose up -d
docker compose ps
```

Open:

```text
http://SERVER:8501
```

Put an HTTPS reverse proxy and authentication in front of the service before exposing it publicly.

## Production Notes

* Recommended initial host: 4 vCPU, 8 GB RAM, 30+ GB SSD.
* Keep the runtime data snapshot and `config/engine-freeze.json` together.
* Run the data-health, engine-flow, and regression checks after every release.
* Do not enable unrestricted concurrent heavy analyses without per-job output isolation and a concurrency limit.
* Review redistribution and licensing terms for every bundled scientific dataset before public distribution.
* Third-party datasets remain subject to their respective licenses and terms of use.

## Release Reproducibility

A Sophiark release consists of the application source code and its corresponding database snapshot.

For reproducible deployment, use the database release associated with the same Sophiark version rather than mixing database snapshots from different releases.
