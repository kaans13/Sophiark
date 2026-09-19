# Sophiark

**Deterministic protein-network perturbation simulation for computational biology and bioinformatics.**

Sophiark is a research-oriented Streamlit application for simulating how a **protein-protein interaction (PPI) network changes after target-focused perturbation**.

Rather than treating a target as an isolated gene, Sophiark examines how perturbing its network connections changes topology, centrality, alternative connectivity, local efficiency, and related network properties within a tissue-filtered biological context.

> **Sophiark does not use machine learning for its core network simulation.**
>
> It is a deterministic, evidence-bounded network analysis tool. Its outputs describe **network topology and computationally defined perturbation responses**, not gene expression, biological causality, clinical outcomes, or treatment efficacy.

[![Dataset DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22849093.svg)](https://doi.org/10.5281/zenodo.22849093)

## What Sophiark Does

Sophiark connects biological evidence layers with network-topology analysis to investigate questions such as:

* How does a target-focused perturbation redistribute network connectivity?
* Which proteins become more structurally important after perturbation?
* Where does local network efficiency decrease?
* Which alternative routes remain available?
* How do these responses differ between tissues?
* How do human and mouse network contexts differ?
* Which biological processes or pathways are associated with the affected network regions?

The goal is not to claim that a computationally observed topology change is automatically a biological effect. Instead, Sophiark provides a reproducible computational framework for **generating and investigating network-level hypotheses**.

## Core Approach

The main workflow is:

```text
Biological data
      ↓
PPI network
      ↓
Tissue-specific filtering / weighting
      ↓
Baseline network metrics
      ↓
Target-focused edge perturbation
      ↓
Perturbed network
      ↓
Topology comparison
      ↓
Redistribution / compensation / structural analysis
      ↓
Optional biological interpretation
```

The core calculation is separated from optional interpretation layers. Enrichment, regulatory evidence, and research-context information can help interpret results, but they do not redefine the underlying PPI calculation.

## Capabilities

### Network simulation

* Human (*Homo sapiens*) and mouse (*Mus musculus*) network paths
* Tissue-specific PPI subnetworks
* Single- and multi-target perturbation
* Target-focused edge attenuation
* Baseline vs. perturbed network comparison
* PageRank
* Betweenness centrality
* Network Importance / Hinterland analysis
* Community analysis
* Local and global efficiency analysis
* Alternative connectivity and redistribution analysis

### Biological context

* Tissue-expression filtering
* GO enrichment
* KEGG enrichment
* TRRUST evidence
* OmniPath directional evidence
* Structural and binding-related context
* Human–mouse analysis support
* Research-context and export workflows

### Reproducibility

Sophiark is designed around explicit data snapshots, deterministic calculations, documented parameters, and reproducibility metadata.

## Example Research Question

A typical Sophiark analysis can be framed as:

> **What happens to the topology of a tissue-specific PPI network when the connections surrounding a selected target are computationally attenuated?**

The resulting changes can then be examined for:

* loss of local connectivity,
* redistribution toward alternative routes,
* changes in network importance,
* bottleneck changes,
* efficiency changes,
* affected communities,
* and associated biological annotations.

These observations can be used to develop hypotheses for further biological investigation.

## Data

The Sophiark v1.0 database package is distributed separately from the source code.

**Dataset:** Sophiark Database v1.0
**Version:** 1.0.0
**DOI:** [10.5281/zenodo.22849093](https://doi.org/10.5281/zenodo.22849093)

The database package contains the biological data snapshots required by the corresponding Sophiark release.

> The dataset contains third-party biological resources. Their respective sources and applicable terms should be considered independently of the Sophiark software license.

## Required Snapshots

| Path                                          | Purpose                                           |
| --------------------------------------------- | ------------------------------------------------- |
| `data/raw/hinterland_core.db`                 | Human PPI interactions and tissue-expression data |
| `data/raw/mouse_hinterland_core.db`           | Mouse PPI interactions and tissue-expression data |
| `data/processed/ensp_with_symbols.csv`        | Human protein-ID to gene-symbol mapping           |
| `data/raw/10090.protein.info.v12.0.txt.gz`    | Mouse STRING-ID to symbol mapping                 |
| `data/regulatory/omnipath/omnipath_9606.tsv`  | Human OmniPath directional-evidence snapshot      |
| `data/regulatory/omnipath/omnipath_10090.tsv` | Mouse OmniPath directional-evidence snapshot      |
| `data/processed/structural_binding.db`        | Structural and binding-related lookup data        |

## Installation

### Requirements

* Python 3.10+
* Local Sophiark database snapshots
* Windows, Linux, or macOS environment capable of running Streamlit

### Clone

```powershell
git clone https://github.com/kaans13/Sophiark.git
cd Sophiark
```

### Create environment

```powershell
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### Add the database

Download the **Sophiark Database v1.0** package from Zenodo and extract its contents into the project directory so that the required `data/` paths are available.

### Launch

```powershell
streamlit run app.py
```

The interface opens in English by default.

Use:

**System → Language**

to switch to Turkish.

## Documentation

* [Researcher Guide](docs/USER_GUIDE.md) — workflow, components, interpretation, limitations, glossary, and FAQ.
* [Parameter Reference](docs/PARAMETERS.md) — source-backed defaults, ranges, options, and practical effects.
* [Product Architecture](docs/product_architecture.md) — application and data-layer architecture.

## Interpretation Boundaries

Sophiark deliberately distinguishes computational topology from biological claims.

* STRING confidence is **not** causality or biological-effect strength.
* Edge attenuation is **not** equivalent to a gene knockout, drug dose, or biophysical inhibition.
* Centrality changes are **not** measurements of gene or protein expression.
* Network redistribution is **not proof of biological compensation**.
* Enrichment does **not** establish pathway activation or inhibition.
* Missing directional evidence is **not** negative biological evidence.
* Computationally identified network changes should be treated as **hypothesis-generating observations**, not experimental validation.

## Scientific Scope

Sophiark is intended for:

* computational biology research,
* bioinformatics,
* systems biology,
* protein-network analysis,
* perturbation hypothesis generation,
* comparative network analysis,
* and research-oriented exploration of potential network responses.

It is **not intended to provide clinical recommendations, therapeutic decisions, or claims of biological causality**.

## License

The Sophiark software is distributed under the **Sophiark Non-Commercial Research License v1.0**.

Academic, research, educational, and personal use is permitted under the license terms. Commercial use requires separate written permission from the copyright holder.

See [`LICENSE`](LICENSE) for the complete terms.

Third-party datasets included in the Sophiark database package remain subject to their respective licenses, terms of use, and attribution requirements.

## Citation

If you use Sophiark in research, please cite the software and the corresponding dataset release.

**Software**

> Temizer, Metin Kaan. *Sophiark — deterministic protein-network perturbation simulation.* 2026.
> https://github.com/kaans13/Sophiark

**Dataset**

> Temizer, Metin Kaan. *Sophiark Database v1.0.* Zenodo, 2026.
> https://doi.org/10.5281/zenodo.22849093

## Project Status

Sophiark is an actively developing research project.

The current public release represents a reproducible research-oriented deployment snapshot. Future releases may extend network models, biological evidence integration, perturbation methods, and comparative analyses.

---

**Repository:** https://github.com/kaans13/Sophiark
**Dataset:** https://doi.org/10.5281/zenodo.22849093
