import gzip
import zipfile
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[1]
RAW = BASE_DIR / "data" / "raw"

files = [
    "10090.protein.links.v12.0.txt.gz",
    "10090.protein.links.v12.0.txt",
    "10090.protein.info.v12.0.txt.gz",
    "10090.protein.info.v12.0.txt",
    "Mus_musculus_expr_simple.tsv",
]

def sniff_file(path: Path):
    if not path.exists():
        print(f"❌ YOK: {path.name}")
        return

    size = path.stat().st_size
    print(f"\n📄 {path.name}  ({size:,} byte)")

    # İlk 8 byte'ı oku
    with open(path, "rb") as f:
        magic = f.read(8)

    if magic.startswith(b"\x1f\x8b"):
        print("   ✅ Gerçek gzip dosyası")
    elif magic.startswith(b"PK"):
        print("   ⚠️ ZIP dosyası (gzip değil)")
    else:
        print("   ❌ Düz metin / sıkıştırılmamış (gzip DEĞİL)")

    # İlk 3 satırı göster (gzip ise açarak)
    try:
        if magic.startswith(b"\x1f\x8b"):
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
                lines = [next(f) for _ in range(3)]
        else:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                lines = [next(f) for _ in range(3)]
        print("   İlk satırlar:")
        for line in lines:
            print(f"     {line.rstrip()[:120]}")
    except Exception as e:
        print(f"   [okuma hatası] {e}")

def main():
    print("=== SOPHIARK RAW DOSYA KONTROLÜ ===\n")
    for fname in files:
        sniff_file(RAW / fname)

if __name__ == "__main__":
    main()