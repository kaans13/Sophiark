import sqlite3
from pathlib import Path

import pandas as pd

BASE_DIR = Path(__file__).resolve().parents[1]
DB_PATH = BASE_DIR / "data" / "raw" / "mouse_hinterland_core.db"

STRING_PPI_FILE = BASE_DIR / "data" / "raw" / "10090.protein.links.v12.0.txt"
STRING_INFO_FILE = BASE_DIR / "data" / "raw" / "10090.protein.info.v12.0.txt.gz"
BGEE_EXPR_FILE = BASE_DIR / "data" / "raw" / "Mus_musculus_expr_simple.tsv"


def is_gzip(path: Path) -> bool:
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def read_csv_smart(path: Path, sep: str):
    if not path.exists():
        raise FileNotFoundError(f"Dosya bulunamadı: {path}")
    compression = "gzip" if is_gzip(path) else None
    return pd.read_csv(path, sep=sep, compression=compression)


def create_mouse_ppi_table(con):
    print("[1/3] STRING fare PPI dosyası okunuyor...")
    df = read_csv_smart(STRING_PPI_FILE, sep=" ")
    df["protein1"] = df["protein1"].str.replace("10090.", "", regex=False)
    df["protein2"] = df["protein2"].str.replace("10090.", "", regex=False)
    df = df[df["combined_score"] >= 500]
    df.to_sql("mouse_interactions", con, if_exists="replace", index=False)
    print(f"      → {len(df)} yüksek güvenli kenar yüklendi.")


def create_mouse_tissue_table(con):
    print("[2/3] STRING protein info okunuyor...")
    info_df = read_csv_smart(STRING_INFO_FILE, sep="\t")
    info_df.columns = info_df.columns.str.replace("#", "", regex=False)
    info_df["string_protein_id"] = info_df["string_protein_id"].str.replace("10090.", "", regex=False)
    symbol_to_protein = dict(zip(info_df["preferred_name"], info_df["string_protein_id"]))
    print(f"      {len(symbol_to_protein)} sembol → protein eşleşmesi hazır.")

    print("Bgee doku ifadesi okunuyor...")
    df = read_csv_smart(BGEE_EXPR_FILE, sep="\t")
    print(f"      Toplam ifade kaydı: {len(df)}")

    df = df[["Gene name", "Anatomical entity name", "Expression score"]]
    df.columns = ["symbol", "tissue", "expression_level"]
    df = df[df["expression_level"] > 0]

    df["protein_id"] = df["symbol"].map(symbol_to_protein)
    df = df.dropna(subset=["protein_id"]).copy()
    df = df[["protein_id", "tissue", "expression_level"]]

    df.to_sql("mouse_tissue_expression", con, if_exists="replace", index=False)
    print(f"      → {len(df)} kayıt, {df['protein_id'].nunique()} protein ile yüklendi.")


def main():
    print("=== Sophiark: Ayrı fare veritabanı oluşturuluyor ===\n")
    if DB_PATH.exists():
        DB_PATH.unlink()

    con = sqlite3.connect(str(DB_PATH))
    try:
        create_mouse_ppi_table(con)
        create_mouse_tissue_table(con)

        print("[3/3] İndeksler oluşturuluyor...")
        con.execute("CREATE INDEX idx_mouse_int_p1 ON mouse_interactions(protein1)")
        con.execute("CREATE INDEX idx_mouse_int_p2 ON mouse_interactions(protein2)")
        con.execute("CREATE INDEX idx_mouse_tissue ON mouse_tissue_expression(tissue, protein_id)")
        con.commit()
        print("\n✔ Veritabanı hazır:", DB_PATH)
    finally:
        con.close()


if __name__ == "__main__":
    main()