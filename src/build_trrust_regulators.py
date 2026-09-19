"""
build_trrust_regulators.py
============================
TRRUST ham verisinden (trrust_rawdata.human.tsv) iki adet pkl üretir:

1. tf_targets.pkl          -> {TF_symbol: [target_symbol, ...]}
   (config.py'deki mevcut TF_TARGETS ile AYNI FORMAT — biology_logic.py'deki
    apply_biological_bonus fonksiyonu bunu zaten kullanıyor. Bu dosyayı
    üzerine yazarsan bonus sistemi otomatik olarak TRRUST'ı kullanmaya başlar.)

2. target_to_regulators.pkl -> {target_symbol: [(TF_symbol, yön), ...]}
   (YENİ — stresli/hedef gen tablolarında "Düzenleyici TF'ler" sütunu için.
    Yön: "Aktivasyon" / "Baskılama" / "Bilinmiyor")

KULLANIM:
    python -m src.build_trrust_regulators
    (SOPHIARK_ROOT / SOPHIARK_DATA_DIR ile path override edilebilir.)

Unknown kayıtlar ÇIKARILMAZ — "etkileşim var, yön belirsiz" anlamına gelir
ve hem tf_targets.pkl (bonus sistemi topolojik ilişkiyi zaten yön-bağımsız
kullanıyor) hem de target_to_regulators.pkl'de (yön "Bilinmiyor" olarak
etiketlenerek) korunur.
"""

import pickle
from collections import defaultdict
from pathlib import Path
import sys

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.product.paths import ProjectPaths

# ---------------------------------------------------------------------------
# YOLLAR — merkezi, taşınabilir resolver
# ---------------------------------------------------------------------------
_PATHS = ProjectPaths.discover()
RAW_TSV_PATH   = _PATHS.raw / "trrust_rawdata.human.tsv"
PROCESSED_DIR  = _PATHS.processed

TF_TARGETS_OUT       = PROCESSED_DIR / "tf_targets.pkl"
REGULATORS_OUT       = PROCESSED_DIR / "target_to_regulators.pkl"

YON_MAP = {
    "Activation": "Aktivasyon",
    "Repression": "Baskılama",
    "Unknown":    "Bilinmiyor",
}


def main():
    if not RAW_TSV_PATH.exists():
        raise FileNotFoundError(f"TRRUST ham dosyası bulunamadı: {RAW_TSV_PATH}")

    print(f"[1/4] TRRUST okunuyor: {RAW_TSV_PATH}")
    df = pd.read_csv(
        RAW_TSV_PATH, sep="\t", header=None,
        names=["TF", "Target", "Yon", "PMID"],
    )
    print(f"      Toplam kayıt: {len(df):,}")
    print(f"      Yön dağılımı:\n{df['Yon'].value_counts().to_string()}")

    print("[2/4] tf_to_targets ve target_to_regulators oluşturuluyor…")
    tf_to_targets = defaultdict(set)          # set -> tekrarları otomatik eler
    target_to_regulators = defaultdict(list)

    for tf, target, yon in zip(df["TF"], df["Target"], df["Yon"]):
        tf_to_targets[tf].add(target)
        target_to_regulators[target].append((tf, YON_MAP.get(yon, yon)))

    # set -> list (pickle/JSON uyumluluğu ve mevcut TF_TARGETS formatıyla birebir eşleşme için)
    tf_to_targets = {tf: sorted(targets) for tf, targets in tf_to_targets.items()}
    target_to_regulators = dict(target_to_regulators)

    print(f"      {len(tf_to_targets):,} benzersiz TF")
    print(f"      {len(target_to_regulators):,} benzersiz hedef gen")

    print("[3/4] Kaydediliyor…")
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    with open(TF_TARGETS_OUT, "wb") as f:
        pickle.dump(tf_to_targets, f)
    print(f"      → {TF_TARGETS_OUT}")

    with open(REGULATORS_OUT, "wb") as f:
        pickle.dump(target_to_regulators, f)
    print(f"      → {REGULATORS_OUT}")

    print("[4/4] Doğrulama örneği:")
    ornek_gen = next(iter(target_to_regulators))
    print(f"      '{ornek_gen}' geninin düzenleyicileri: {target_to_regulators[ornek_gen][:5]}")

    print("\n✔ Tamamlandı.")


if __name__ == "__main__":
    main()
