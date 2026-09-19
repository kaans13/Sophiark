# ---------------------------------------------------------------------------
# mouse_config.py
# Fare (Mus musculus) için ayrı konfigürasyon.
# İnsan config.py'den tamamen bağımsızdır, ayrı veritabanı ve ayrı dosyalar kullanır.
# ---------------------------------------------------------------------------

from pathlib import Path as _Path
from ..runtime_logging import ensure_utf8_standard_streams

ensure_utf8_standard_streams()
BASE_DIR = _Path(__file__).resolve().parents[2]  # src/mouse/../.. → proje kökü
_STRING_INFO_PATH = BASE_DIR / "data" / "raw" / "10090.protein.info.v12.0.txt.gz"

def _load_mouse_symbol_map() -> dict:
    """STRING protein.info dosyasından protein_id -> sembol sözlüğü yükler."""
    import pandas as pd  # ← LOKAL IMPORT GARANTİ

    try:
        _STRING_INFO_PATH = BASE_DIR / "data" / "raw" / "10090.protein.info.v12.0.txt.gz"
        with open(_STRING_INFO_PATH, "rb") as f:
            magic = f.read(2)
        compression = "gzip" if magic == b"\x1f\x8b" else None

        df = pd.read_csv(_STRING_INFO_PATH, sep="\t", compression=compression)
        df.columns = df.columns.str.replace("#", "", regex=False)
        df["string_protein_id"] = df["string_protein_id"].str.replace("10090.", "", regex=False)
        return dict(zip(df["string_protein_id"], df["preferred_name"]))
    except Exception as e:
        print(f"STRING info okunamadı: {e}")
        return {}

MOUSE_SYMBOL_MAP = _load_mouse_symbol_map()
print(f"MOUSE_SYMBOL_MAP yüklendi: {len(MOUSE_SYMBOL_MAP)} protein")
import logging
import pickle
import sys
import warnings
from pathlib import Path
from typing import Optional
import pandas as pd
import requests

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# LOGGING
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger("HinterlandMapping_Mouse")

# ---------------------------------------------------------------------------
# HASAR TİPİ SABİTLERİ
# ---------------------------------------------------------------------------
HASAR_BIRINCIL = "Birincil_Hasar_Hedef"
HASAR_IKINCIL  = "İkincil_Hasar_Dökülen"
HASAR_UCUNCUL  = "Üçüncül_Hasar_Stresli"

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
DB_PATH  = str(BASE_DIR / "data" / "raw" / "mouse_hinterland_core.db")
OUT_DIR  = BASE_DIR / "outputs_mouse"
# Fare OmniPath snapshot'Ä± ayrÄ± kurulmalÄ±dÄ±r; insan evidence'Ä± sessizce kullanÄ±lmaz.
OMNIPATH_SNAPSHOT_PATH = BASE_DIR / "data" / "regulatory" / "omnipath" / "omnipath_10090.tsv"
ENABLE_DIRECTED_SIGNALING = True
OUT_CSV  = str(OUT_DIR / "reports" / "mouse_genom_analiz_sonuc_v5.csv")

HIGH_CONF_THRESHOLD = 500
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
# DOZ-YANIT SABİTLERİ
# ---------------------------------------------------------------------------
DOZLAR                = [0.90, 0.75, 0.50, 0.25, 0.10, 0.01]
DOSE_RESPONSE_DAMPING = 0.50
HILL_COEFFICIENT      = 2
PARALOG_BOOST_FACTOR  = 0.10
BLAST_SCORE_MIN       = 200
NULL_MODEL_ITER       = 200
NULL_MODEL_Z_CUTOFF   = 2.0

HUB_PENALTY_DAMPING = 0.5

# ---------------------------------------------------------------------------
# BİYOLOJİK BONUS SİSTEMİ — Fare TRRUST + CORUM
# ---------------------------------------------------------------------------
try:
    with open(BASE_DIR / "data" / "processed" / "mouse_tf_targets.pkl", "rb") as f:
        TF_TARGETS = pickle.load(f)
    with open(BASE_DIR / "data" / "processed" / "mouse_complex_members.pkl", "rb") as f:
        COMPLEX_MEMBERS = pickle.load(f)
    BONUS_ENABLED = True
except FileNotFoundError:
    TF_TARGETS = {}
    COMPLEX_MEMBERS = {}
    BONUS_ENABLED = False

# Kalıcı ENSP→Sembol cache (fare için ayrı)
_CACHE_PATH = BASE_DIR / "outputs_mouse" / "cache" / "ensmusp_symbol_cache.pkl"

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

def _ensmusp_to_symbol(ensmusp: str) -> Optional[str]:
    """ENSMUSP ID'sinden fare gen sembolüne çevir (yerel STRING haritasından, API'siz)."""
    clean = str(ensmusp).replace("10090.", "")
    return MOUSE_SYMBOL_MAP.get(clean)

# ---------------------------------------------------------------------------
# GO:CC HİYERARŞİ HARİTASI (insanla aynı)
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
    ("plasma membrane",        "Cell_Membrane"),
    ("cell surface",           "Cell_Membrane"),
    ("cell periphery",         "Cell_Membrane"),
    ("basolateral plasma membrane", "Cell_Membrane"),
    ("apical plasma membrane", "Cell_Membrane"),
    ("integral component of plasma membrane", "Cell_Membrane"),
    ("intrinsic component of plasma membrane", "Cell_Membrane"),
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
MYGENE_SPECIES = "mouse"
API_CHUNK_SIZE = 500
API_RETRY_MAX  = 3
API_RETRY_WAIT = 5
API_TIMEOUT    = 30

# Fare için şimdilik boş; gerekirse doldurulur
DB_FIX_MAP = {}

_DARK_BG  = "#0F0F1A"
_PANEL_BG = "#13132A"

# ---------------------------------------------------------------------------
# TRRUST DÜZENLEYİCİ LOOKUP — fare hedef gen → onu düzenleyen TF'ler
# ---------------------------------------------------------------------------
try:
    with open(BASE_DIR / "data" / "processed" / "mouse_target_to_regulators.pkl", "rb") as f:
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

# ---------------------------------------------------------------------------
# BİOGRID ORCS ESSENTIALITY LOOKUP — farede şimdilik boş
# ---------------------------------------------------------------------------
def _load_mouse_ortholog_essentiality() -> dict:
    """İnsan CRISPR essentiality kanıtını fare ortoloğuna aktar.

    Bir fare geni birden çok insan ortoloğuna eşleşirse ekran sayıları ve hit
    sayıları toplanır. Bu alan fare deney verisi değil, translasyonel kanıttır.
    """
    try:
        with open(BASE_DIR / "data" / "processed" / "biogrid_essentiality.pkl", "rb") as f:
            human_scores = pickle.load(f)
        with open(BASE_DIR / "data" / "processed" / "mouse_ortholog_map.pkl", "rb") as f:
            human_to_mouse = pickle.load(f)
    except (FileNotFoundError, pickle.UnpicklingError):
        return {}

    aggregated: dict[str, dict[str, int]] = {}
    for human_symbol, mouse_symbol in human_to_mouse.items():
        record = human_scores.get(human_symbol)
        if not record or not mouse_symbol:
            continue
        item = aggregated.setdefault(mouse_symbol, {"test_sayisi": 0, "hit_sayisi": 0})
        item["test_sayisi"] += int(record.get("test_sayisi", 0))
        item["hit_sayisi"] += int(record.get("hit_sayisi", 0))

    for item in aggregated.values():
        tests = item["test_sayisi"]
        item["hit_orani"] = round(item["hit_sayisi"] / tests, 4) if tests else 0.0
    return aggregated


GENE_ESSENTIALITY = _load_mouse_ortholog_essentiality()
ESSENTIALITY_ENABLED = bool(GENE_ESSENTIALITY)

def format_essentiality(symbol: str, min_test: int = 1) -> str:
    record = GENE_ESSENTIALITY.get(symbol)
    if not record or int(record.get("test_sayisi", 0)) < min_test:
        return "—"
    return (
        f"Translasyonel ORCS: %{float(record.get('hit_orani', 0.0)) * 100:.0f} "
        f"({int(record.get('hit_sayisi', 0))}/{int(record.get('test_sayisi', 0))})"
    )


try:
    with open(BASE_DIR / "data" / "processed" / "mouse_mygene_cache.pkl", "rb") as f:
        MOUSE_MYGENE_CACHE = pickle.load(f)
except (FileNotFoundError, pickle.UnpicklingError):
    MOUSE_MYGENE_CACHE = {}
