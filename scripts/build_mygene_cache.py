"""MyGene API'den tüm ENSP'ler için veri çekip yerel pickle dosyası oluşturur."""
import pickle
import time
import pandas as pd
import requests
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[1]
MYGENE_URL = "https://mygene.info/v3/query"
CHUNK_SIZE = 500  # Her seferde 500 gen sorgula

def fetch_batch(ensp_list):
    clean_ids = [eid.replace("9606.", "") for eid in ensp_list]
    payload = {
        "q": ",".join(clean_ids),
        "scopes": "ensembl.protein",
        "fields": "symbol,name,go.CC,go.MF,go.BP",
        "species": "human",
        "size": str(len(clean_ids))
    }
    r = requests.post(MYGENE_URL, data=payload, timeout=60)
    r.raise_for_status()
    return r.json()

def main():
    ensp_df = pd.read_csv(BASE_DIR / "data" / "processed" / "ensp_with_symbols.csv")
    all_ensp = ensp_df["gene"].unique().tolist()
    total = len(all_ensp)
    
    cache = {}
    for i in range(0, total, CHUNK_SIZE):
        chunk = all_ensp[i:i+CHUNK_SIZE]
        print(f"Çekiliyor: {i+1}-{min(i+CHUNK_SIZE, total)} / {total}")
        
        try:
            hits = fetch_batch(chunk)
            if isinstance(hits, dict):
                hits = hits.get("hits", [])
            for hit in hits:
                if not isinstance(hit, dict) or hit.get("notfound"):
                    continue
                ensp = hit.get("query")
                cache[ensp] = {
                    "symbol": hit.get("symbol", "—"),
                    "name": hit.get("name", "—"),
                    "go_cc": hit.get("go", {}).get("CC", []),
                    "go_mf": hit.get("go", {}).get("MF", []),
                    "go_bp": hit.get("go", {}).get("BP", []),
                }
        except Exception as e:
            print(f"Hata: {e}, 10 saniye bekleniyor...")
            time.sleep(10)
        
        time.sleep(0.3)
    
    cache_path = BASE_DIR / "data" / "processed" / "mygene_cache.pkl"
    with open(cache_path, "wb") as f:
        pickle.dump(cache, f)
    print(f"✓ {len(cache)} gen için veri kaydedildi: {cache_path}")

if __name__ == "__main__":
    main()