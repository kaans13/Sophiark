# Sophiark ürün mimarisi ve veri güncelleme sözleşmesi

Bilimsel motorlar bu mimarinin iç bileşeni değil, dondurulmuş hesaplama
referanslarıdır. `config/engine-freeze.json` Classic, Directed ve Evidence
hesaplama kaynaklarının SHA-256 değerlerini tutar. Test paketi bu dosyalarda
onaysız drift olduğunda durur. Veri, sunum veya orkestrasyon ihtiyacı motor
matematiği değiştirilerek çözülmez.

```text
data/raw → validate → normalize/map → data/processed + data/derived
         → fingerprint → capability/preflight → frozen engines
         → context/synthesis → presentation → versioned export
```

## Sınırlar

- `src/product/paths.py`: cwd ve kullanıcı dizininden bağımsız yollar.
- `config/datasets.json`: required/optional veri kataloğu, roller ve şemalar.
- `src/product/datasets.py`: registry, şema kontrolü ve health durumu.
- `src/product/capabilities.py`: UI ve orchestration için ortak availability.
- `src/product/entities.py`: canonical ENSP normalizasyonu ve mapping coverage.
- `src/product/build.py`: deterministik, atomik STRING fixture/build sözleşmesi.
- `src/product/fingerprints.py`: dataset/config duyarlı cache kimliği.
- `src/product/provenance.py`: engine, data, preprocessing ve export sürümlerini
  ayrı taşıyan run manifest.
- `src/product/preflight.py`: Streamlit'ten bağımsız Classic/Directed/Evidence
  veri önkoşulu kontrolü.

`python scripts/data_health.py` hızlı availability ve schema raporu üretir.
`--fingerprint` tüm dataset dosyalarını içerik bazında SHA-256 ile doğrular;
büyük SQLite snapshot'larında bu seçenek doğal olarak daha yavaştır.

HPA, OmniPath, CORUM ve TRRUST shadow-build/update akışları ile authoritative
source envanteri için [`data_ecosystem_inventory.md`](data_ecosystem_inventory.md)
belgesine bakın. `python scripts/data_build.py --dataset all` yalnız staging
artifact üretir; production promotion yapmaz ve internetten veri indirmez.

## Dataset güncelleme

Yeni ham dosya önce `data/raw` altında immutable snapshot olarak saklanır.
Source-specific builder şemayı, organizmayı, identifier biçimini, duplicate
anahtarları ve sayı aralıklarını doğrular. İşlenmiş çıktı `data/processed`
altına atomik olarak yazılır; build report source/valid/mapped/unmapped/dropped
sayılarını, coverage ve fingerprint'i kaydeder. Derived indexler
`data/derived`, runtime cache'leri `outputs/cache` altında tutulur.

STRING kimlik dönüşümü yalnız gerçek human STRING biçimi için
`9606.ENSP… → ENSP…` uygular. Zaten canonical ENSP değişmez; yabancı taxon ve
bilinmeyen biçim açık hatadır. `combined_score` 0–1000 kaynak ölçeğinde kalır.

## Veri rolleri

- STRING core: PPI ve Classic tissue-context omurgası (required).
- ENSP mapping: canonical internal identity (required).
- OmniPath: Directed yön politikası (optional, engine-specific).
- STRING full channels/physical: Evidence provenance (optional).
- CORUM-derived membership: hesaplama sonrası complex context (optional).

Optional kaynak yokluğu başka motor için sessiz fallback üretmez. Örneğin
OmniPath yokken Classic preflight READY, Directed DATA_UNAVAILABLE döner.
