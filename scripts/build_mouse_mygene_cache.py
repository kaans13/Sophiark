import sqlite3
import pickle
import time
import pandas as pd
import requests
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[1]
DB_PATH = BASE_DIR / "data" / "raw" / "mouse_hinterland_core.db"
CACHE_PATH = BASE_DIR / "data" / "processed" / "mouse_mygene_cache.pkl"
MYGENE_URL = "https://mygene.info/v3/query"
CHUNK_SIZE = 500

def fetch_batch(ids):
    clean = [i.replace("10090.", "") for i in ids]
    payload = {
        "q": ",".join(clean),
        "scopes": "ensembl.protein",
        "fields": "symbol,name,go.BP,go.CC,go.MF",
        "species": "mouse",
        "size": str(len(clean))
    }
    r = requests.post(MYGENE_URL, data=payload, timeout=60)
    r.raise_for_status()
    data = r.json()
    return data.get("hits", []) if isinstance(data, dict) else data

def main():
    if not DB_PATH.exists():
        print(f"Veritabanı bulunamadı: {DB_PATH}")
        return

    con = sqlite3.connect(str(DB_PATH))
    # Benzersiz protein ID'lerini al
    rows = con.execute("""
        SELECT DISTINCT p FROM (
            SELECT protein1 AS p FROM mouse_interactions
            UNION
            SELECT protein2 AS p FROM mouse_interactions
        )
    """).fetchall()
    con.close()

    all_ids = [r[0] for r in rows]
    print(f"Toplam benzersiz protein: {len(all_ids)}")

    cache = {}
    for i in range(0, len(all_ids), CHUNK_SIZE):
        chunk = all_ids[i:i+CHUNK_SIZE]
        print(f"İşleniyor: {i+1}-{min(i+CHUNK_SIZE, len(all_ids))} / {len(all_ids)}")
        try:
            hits = fetch_batch(chunk)
        except Exception as e:
            print(f"  Hata: {e}, 5 sn bekleniyor...")
            time.sleep(5)
            continue

        for hit in hits:
            if not isinstance(hit, dict) or "notfound" in hit and hit["notfound"]:
                continue
            ensp = hit.get("query")
            if not ensp:
                continue
            go_bp = hit.get("go", {}).get("BP", [])
            go_cc = hit.get("go", {}).get("CC", [])
            go_mf = hit.get("go", {}).get("MF", [])

            def _term_list(go_data):
                if isinstance(go_data, dict):
                    go_data = [go_data]
                return [t.get("term", "") for t in go_data if isinstance(t, dict)]

            cache[ensp] = {
                "symbol": hit.get("symbol", "—"),
                "name": hit.get("name", "—"),
                "go_bp": _term_list(go_bp),
                "go_cc": _term_list(go_cc),
                "go_mf": _term_list(go_mf),
            }
        time.sleep(0.3)

    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(CACHE_PATH, "wb") as f:
        pickle.dump(cache, f)
    print(f"✓ {len(cache)} protein kaydı {CACHE_PATH} dosyasına yazıldı.")

if __name__ == "__main__":
    main()