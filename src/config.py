# ---------------------------------------------------------------------------
# config.py
# Genel ayarlar, sabitler, logging ve TF/Kompleks lookup + sembol cache.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# PROJE KÖKÜ — tek yerden yönetilen path (migrate.py tarafından eklendi)
# ---------------------------------------------------------------------------
from pathlib import Path as _Path
BASE_DIR = _Path(__file__).resolve().parents[1]   # src/../  →  genom/

import logging
import pickle
import sys
import warnings
from pathlib import Path
from typing import Optional

import requests

from .runtime_logging import ensure_utf8_standard_streams

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# LOGGING
# ---------------------------------------------------------------------------
ensure_utf8_standard_streams()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger("HinterlandMapping")

# ---------------------------------------------------------------------------
# HASAR TİPİ SABİTLERİ
# ---------------------------------------------------------------------------
HASAR_BIRINCIL = "Birincil_Hasar_Hedef"
HASAR_IKINCIL  = "İkincil_Hasar_Dökülen"
HASAR_UCUNCUL  = "Üçüncül_Hasar_Stresli"

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
DB_PATH  = str(BASE_DIR / "data" / "raw" / "hinterland_core.db")
OUT_CSV  = str(BASE_DIR / "outputs" / "reports" / "genom_analiz_sonuc_v5.csv")
OUT_DIR  = BASE_DIR / "outputs"
OMNIPATH_SNAPSHOT_PATH = BASE_DIR / "data" / "regulatory" / "omnipath" / "omnipath_9606.tsv"
ENABLE_DIRECTED_SIGNALING = True

HIGH_CONF_THRESHOLD = 500  # StringDB combined_score eşik değeri (300-999 arası önerilir)
PAGERANK_DAMPING    = 0.85
PAGERANK_MAX_ITER   = 200
LOUVAIN_SEED        = 42

NOISE_MAX_DEGREE        = 4
NOISE_MAX_SCORE         = 700
NOISE_MAX_PR_PERCENTILE = 5

THRESHOLD_EMPEROR   = 99
THRESHOLD_STRATEGIC = 90
THRESHOLD_BRIDGE    = 50

EFFICIENCY_SAMPLE_SIZE = 3000

SWEEP_MIN  = 400
SWEEP_MAX  = 900
SWEEP_STEP = 100

RESCUE_THRESHOLDS   = [450, 400, 350, 300, 250, 200]
RESCUE_MIN_PARTNERS = 1

FORCED_DISTANCE_SCALE = 1.5
BLOCK_WEIGHT_FRACTION = 0.001

WBI_EXT_RATIO_MIN  = 0.30
WBI_BC_PERCENTILE  = 60.0

# ---------------------------------------------------------------------------
# DOZ-YANITLA İLGİLİ SABİTLER
# ---------------------------------------------------------------------------
# DOZLAR = "kalan aktivite fraksiyonu" (survival fraction, 0–1 arası).
#   0.90 → %10 inhibisyon, %90 aktivite kaldı.
#   0.01 → %99 inhibisyon, %1 aktivite kaldı.
#   IC50 → sinyal kaybının %50 olduğu survival_fraction noktası.
#   IC50 büyükse (ör. 0.80) → az inhibisyon yeterli → sistem hassas.
#   IC50 küçükse (ör. 0.05) → çok inhibisyon gerekiyor → sistem dirençli.
DOZLAR                = [0.90, 0.75, 0.50, 0.25, 0.10, 0.01]
DOSE_RESPONSE_DAMPING = 0.50
HILL_COEFFICIENT      = 2
PARALOG_BOOST_FACTOR  = 0.10
BLAST_SCORE_MIN       = 200
NULL_MODEL_ITER       = 200
NULL_MODEL_Z_CUTOFF   = 2.0

HUB_PENALTY_DAMPING = 0.5

# ===========================================================================
# BİYOLOJİK BONUS SİSTEMİ - TF-Hedef + Kompleks (v5.1)
# ===========================================================================
_BASE    = BASE_DIR / "data" / "processed"
_CACHE_PATH = BASE_DIR / "outputs" / "cache" / "ensp_symbol_cache.pkl"

# ── TF & Kompleks lookup tablolarını yükle ──
try:
    with open(BASE_DIR / "data" / "processed" / "tf_targets.pkl", "rb") as f:
        TF_TARGETS = pickle.load(f)
    with open(BASE_DIR / "data" / "processed" / "complex_members.pkl", "rb") as f:
        COMPLEX_MEMBERS = pickle.load(f)
    BONUS_ENABLED = True
except FileNotFoundError:
    TF_TARGETS = {}
    COMPLEX_MEMBERS = {}
    BONUS_ENABLED = False

# ── Kalıcı ENSP→Sembol cache (dosya tabanlı) ──
def _load_symbol_cache() -> dict:
    if _CACHE_PATH.exists():
        try:
            with open(_CACHE_PATH, "rb") as f:
                return pickle.load(f)
        except Exception:
            pass
    return {}

def _save_symbol_cache(cache: dict) -> None:
    try:
        with open(_CACHE_PATH, "wb") as f:
            pickle.dump(cache, f)
    except Exception:
        pass

_ENSPSYMBOL_CACHE = _load_symbol_cache()

def _ensp_to_symbol(ensp: str) -> Optional[str]:
    """ENSP ID'sinden gen sembolüne çevir (kalıcı dosya cache'li)"""
    if ensp in _ENSPSYMBOL_CACHE:
        return _ENSPSYMBOL_CACHE[ensp]

    try:
        clean = str(ensp).replace("9606.", "")
        r = requests.post(
            "https://mygene.info/v3/query",
            json={"q": clean, "scopes": "ensembl.protein", "fields": "symbol", "species": "human"},
            timeout=5
        )
        if r.status_code == 200:
            hits = r.json()
            if hits and "symbol" in hits[0]:
                symbol = hits[0]["symbol"]
                _ENSPSYMBOL_CACHE[ensp] = symbol
                # Her 100 yeni sembolde bir dosyaya kaydet
                if len(_ENSPSYMBOL_CACHE) % 100 == 0:
                    _save_symbol_cache(_ENSPSYMBOL_CACHE)
                return symbol
    except Exception:
        pass

    _ENSPSYMBOL_CACHE[ensp] = None
    _save_symbol_cache(_ENSPSYMBOL_CACHE)
    return None


# ---------------------------------------------------------------------------
# GO:CC HİYERARŞİ HARİTASI
# ---------------------------------------------------------------------------
GO_CC_HIERARCHY: list[tuple[str, str]] = [
    ("nucleus",                "Nucleus"),
    ("nucleol",                "Nucleus"),
    ("nuclear",                "Nucleus"),
    ("chromatin",              "Nucleus"),
    ("chromosome",             "Nucleus"),
    ("nuclear envelope",       "Nucleus"),
    ("nuclear pore",           "Nucleus"),
    ("mitochondri",            "Mitochondrion"),
    ("endoplasmic reticulum",  "Endoplasmic_Reticulum"),
    ("er lumen",               "Endoplasmic_Reticulum"),
    ("er membrane",            "Endoplasmic_Reticulum"),
    ("sarcoplasmic reticulum", "Endoplasmic_Reticulum"),
    ("golgi",                  "Golgi_Apparatus"),
    ("trans-golgi",            "Golgi_Apparatus"),
    ("cis-golgi",              "Golgi_Apparatus"),
    # PLAZMA MEMBRANI - sadece spesifik terimler
    ("plasma membrane",        "Cell_Membrane"),
    ("cell surface",           "Cell_Membrane"),
    ("cell periphery",         "Cell_Membrane"),
    ("basolateral plasma membrane", "Cell_Membrane"),
    ("apical plasma membrane", "Cell_Membrane"),
    ("integral component of plasma membrane", "Cell_Membrane"),
    ("intrinsic component of plasma membrane", "Cell_Membrane"),
    # GENEL "membrane" KALDIRILDI - çok geniş kapsamlı
    ("extracellular",          "Extracellular"),
    ("secreted",               "Extracellular"),
    ("collagen",               "Extracellular"),
    ("extracellular matrix",   "Extracellular"),
    ("extracellular space",    "Extracellular"),
    ("lysosom",                "Lysosome_Endosome"),
    ("endosom",                "Lysosome_Endosome"),
    ("vacuol",                 "Lysosome_Endosome"),
    ("phagosom",               "Lysosome_Endosome"),
    ("autophagosom",           "Lysosome_Endosome"),
    ("cytoplasm",              "Cytoplasm"),
    ("cytosol",                "Cytoplasm"),
    ("cytoskeleton",           "Cytoplasm"),
    ("ribosom",                "Cytoplasm"),
    ("proteasom",              "Cytoplasm"),
    ("peroxisom",              "Peroxisome"),
]

MYGENE_URL     = "https://mygene.info/v3/query"
MYGENE_FIELDS = "go.CC,go.MF"
MYGENE_SPECIES = "human"
API_CHUNK_SIZE = 500
API_RETRY_MAX  = 3
API_RETRY_WAIT = 5
API_TIMEOUT    = 30

DB_FIX_MAP = {"ENSP00000426309": "ENSP00000003084"}

_DARK_BG  = "#0F0F1A"
_PANEL_BG = "#13132A"


# UI/kanıt katmanı uyumluluğu: bu tablolar simülasyon matematiğini
# değiştirmez; yalnızca mevcut TRRUST ve essentiality açıklamalarını besler.
try:
    with open(BASE_DIR / "data" / "processed" / "target_to_regulators.pkl", "rb") as f:
        TARGET_TO_REGULATORS = pickle.load(f)
    REGULATORS_ENABLED = True
except FileNotFoundError:
    TARGET_TO_REGULATORS = {}
    REGULATORS_ENABLED = False


def format_regulators(symbol: str, max_show: int | None = None) -> str:
    """Tüm doğrulanmış düzenleyici TF'leri, istenirse açık bir sınırla biçimlendir."""
    if not symbol:
        return "—"
    regs = TARGET_TO_REGULATORS.get(symbol, [])
    if not regs:
        return "—"
    visible = regs if max_show is None else regs[:max_show]
    parts = [f"{tf} ({yon})" for tf, yon in visible]
    extra = f" +{len(regs) - max_show} daha" if max_show is not None and len(regs) > max_show else ""
    return ", ".join(parts) + extra


try:
    with open(BASE_DIR / "data" / "processed" / "biogrid_essentiality.pkl", "rb") as f:
        GENE_ESSENTIALITY = pickle.load(f)
    ESSENTIALITY_ENABLED = True
except FileNotFoundError:
    GENE_ESSENTIALITY = {}
    ESSENTIALITY_ENABLED = False


def format_essentiality(symbol: str, min_test: int = 1) -> str:
    if not symbol:
        return "—"
    veri = GENE_ESSENTIALITY.get(symbol)
    if not veri or veri["test_sayisi"] < min_test:
        return "—"
    return f"{veri['hit_sayisi']}/{veri['test_sayisi']} ekran (%{veri['hit_orani']*100:.0f})"


def get_essentiality_ratio(symbol: str) -> Optional[float]:
    veri = GENE_ESSENTIALITY.get(symbol)
    return veri["hit_orani"] if veri else None



