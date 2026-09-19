# Sophiark deployment package

This directory is a deployment copy. The source workspace is not modified.

## Included

- Supported Streamlit entry point: `app.py`
- Application packages, UI assets, local browser libraries and configuration
- Human and mouse runtime databases
- Identity, OmniPath, evidence, structural-context and interpretation lookups
- Human and mouse baseline annotation seeds required by the frozen engines
- Local ENSP-symbol and enrichment caches needed for bounded/offline startup
- Tests and acceptance scripts for pre-release verification

Large build-only sources are intentionally excluded: raw BindingDB, AlphaFold
archives/structures, BioGRID archives, raw STRING build inputs, staging outputs,
previous analysis outputs and virtual environments.

## Local verification

```powershell
python scripts/data_health.py
python scripts/verify_cftr_lung_engine_flow.py
python -m pytest -q
streamlit run app.py
```

## Docker

The Docker image contains the application but not the runtime datasets. Compose
mounts the local `data/` directory read-only and uses persistent named volumes
for human and mouse outputs. `docker-entrypoint.sh` initializes a new empty
output volume with the two frozen-engine baseline annotation tables; it never
overwrites an existing volume. It also seeds the local symbol/enrichment cache
so normal analysis does not fan out into thousands of MyGene lookups.

```bash
docker compose build
docker compose up -d
docker compose ps
```

Open `http://SERVER:8501`. Put an HTTPS reverse proxy and authentication in
front of the service before exposing it publicly.

## Production notes

- Recommended initial host: 4 vCPU, 8 GB RAM, 30+ GB SSD.
- Keep the runtime data snapshot and `config/engine-freeze.json` together.
- Run the data-health, engine-flow and regression checks after every release.
- Do not enable unrestricted concurrent heavy analyses without per-job output
  isolation and a concurrency limit.
- Review redistribution/licensing terms for every bundled scientific dataset
  before public distribution.
