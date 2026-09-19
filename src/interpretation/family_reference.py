"""Local, curated reference readers used only by the interpretation layer."""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Mapping

import pandas as pd


_ROOT = Path(__file__).resolve().parents[2]
HGNC_GROUPS_PATH = _ROOT / "data" / "interpretation_reference" / "hgnc" / "hgnc_gene_groups.tsv"
COMPLEX_PORTAL_PATH = _ROOT / "data" / "interpretation_reference" / "complex_portal" / "complex_portal_9606.tsv"

# Groups which describe a domain or a generic technical category are not surfaced
# as a family claim. The underlying HGNC rows remain available for audit.
_LOW_SPECIFICITY_TERMS = (
    "protein coding host genes",
    "small nucleolar rna",
)
_LOW_SPECIFICITY_NAMES = {
    # This broad cross-class label would turn unrelated METTL symbols into a
    # purported family. Specific HGNC methyltransferase groups remain valid.
    "dna/rna methyltransferases",
    # A shared adaptor-domain annotation spans AP-1 through AP-5. It is not
    # evidence that AP2 and AP3 form one complex or one paralog family.
    "clathrin/coatomer adaptor, adaptin-like, n-terminal domain containing",
}


@lru_cache(maxsize=1)
def hgnc_groups() -> Mapping[str, Mapping[str, object]]:
    """Return stable group records keyed by HGNC group id.

    The source is intentionally read from a versioned local download. No web
    call is performed during analysis or rendering.
    """
    if not HGNC_GROUPS_PATH.exists():
        return {}
    frame = pd.read_csv(HGNC_GROUPS_PATH, sep="\t", dtype=str, keep_default_na=False)
    required = {"Group ID", "Group name", "Approved symbol"}
    if not required.issubset(frame.columns):
        return {}
    result: dict[str, Mapping[str, object]] = {}
    for group_id, part in frame.groupby("Group ID", sort=False):
        name = str(part["Group name"].iloc[0]).strip()
        if (not name or name.casefold() in _LOW_SPECIFICITY_NAMES
                or any(term in name.casefold() for term in _LOW_SPECIFICITY_TERMS)):
            continue
        result[str(group_id)] = {
            "id": str(group_id),
            "name": name,
            "members": frozenset(
                symbol.strip().upper()
                for symbol in part["Approved symbol"].tolist()
                if str(symbol).strip()
            ),
        }
    return result


@lru_cache(maxsize=1)
def hgnc_groups_by_symbol() -> Mapping[str, tuple[str, ...]]:
    """Index curated group memberships by approved human symbol."""
    index: dict[str, list[str]] = defaultdict(list)
    for group_id, group in hgnc_groups().items():
        for symbol in group["members"]:
            index[str(symbol)].append(group_id)
    return {symbol: tuple(ids) for symbol, ids in index.items()}


def complex_portal_available() -> bool:
    """Signal data availability without guessing UniProt-to-symbol mappings."""
    return COMPLEX_PORTAL_PATH.exists() and COMPLEX_PORTAL_PATH.stat().st_size > 0
