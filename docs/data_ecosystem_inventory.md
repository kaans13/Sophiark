# Sophiark aktif veri ekosistemi

Bu envanter runtime importları, loader/query yolları ve mevcut production
artifact'ları üzerinden hazırlanmıştır. Bilinmeyen release bilgileri tahmin
edilmemiştir.

| Kaynak | Aktif rol | Authoritative runtime temsili | Kimlik / normalizasyon | Durum |
|---|---|---|---|---|
| STRING v12 | PPI ve evidence channels | `data/raw/hinterland_core.db`; `data/processed/string_evidence_v12_9606.sqlite` | STRING prefix → ENSP | Standard builder mevcut |
| HPA | Human tissue context | `hinterland_core.db:tissue_expression` | ENSG → ENSP; tissue title-case; runtime `MAX(expression_level)` | Shadow builder mevcut; production tablo legacy discrete ve TPM kayıtlarını birlikte içeriyor |
| OmniPath | Directed yön/sign context | raw TSV + `omnipath_signaling_9606.sqlite` | symbol → unambiguous ENSP; yön/sign metadata aynen korunur | Shadow builder exact-compatible |
| CORUM | Classic legacy complex modifier; Unified/Evidence context | Classic: `complex_members.pkl`; context: `coreComplexes.txt` | Classic: all-organism primary symbol uppercase; context: Human symbol ve audit için symbol → ENSP | Classic legacy semantiği reproducible; clean Human-only builder ayrı tutulur |
| TRRUST | TF/target modifier ve regulatory context | raw `trrust_human.tsv`; iki processed pickle | symbol; ayrıca unambiguous symbol → ENSP audit tablosu | Shadow builder exact-compatible |
| MyGene-derived cache | GO/name/local annotation | `mygene_cache.pkl`, mouse karşılığı | canonical ENSP keyed cache | `DEFERRED_SOURCE_BUILDER`: cache üretim scripti var; external-response provenance/version lifecycle ayrı çalışma ister |
| HGNC Gene Groups | family/context mapping | interpretation reference TSV + metadata JSON | ENSG/symbol → ENSP, ambiguity rejected | `DEFERRED_SOURCE_BUILDER`: mevcut immutable snapshot/provider zaten schema kontrollü |
| Complex Portal | research complex context | interpretation reference TSV + metadata JSON | canonical provider mapping | `DEFERRED_SOURCE_BUILDER`: CORUM semantics ile birleştirilmemeli |
| SIGNOR | regulatory evidence | raw `signor.tsv`, processed `signor_relations.tsv` | source-specific symbols/relations | `DEFERRED_SOURCE_BUILDER`: iki aktif temsilin authority/migration semantics'i ayrıca doğrulanmalı |
| BioGRID ORCS | essentiality annotation | `biogrid_essentiality.pkl` | symbol keyed annotation | `DEFERRED_SOURCE_BUILDER`: archival raw script ve provenance temizliği gerekir |
| AlphaFold / BindingDB | structural/binding context | `structural_binding.db` | ENSP → UniProt | `DEFERRED_SOURCE_BUILDER`: çok büyük multi-source pipeline; bu turda güvenli shadow rebuild yapılmadı |
| Ensembl orthology | human/mouse comparison | orthology TSV + metadata JSON | ENSP/ENSMUSP pairs | Mevcut dedicated fetch/build akışı; bu turun human priority kapsamı dışında |
| GO/KEGG Enrichr libraries | optional enrichment | runtime external library names/results | gene symbols; species-specific library | Local immutable raw dataset değildir; automatic download builder eklenmedi |

## HPA legacy davranışı

Production `tissue_expression` tablosu aynı protein/tissue için hem continuous
HPA TPM hem de eski `High/Medium/Low/Not detected → 3/2/1/0.1` kayıtları
içerebilir. Runtime bunların `MAX` değerini kullanır. Raw RNA HPA shadow build
bu legacy discrete kaynağı yeniden üretmez. Bu davranış
`LEGACY_SCIENTIFIC_DATA_BEHAVIOR` olarak korunur; shadow HPA production'a
promote edilmez.

## CORUM authority sınırı

Unified/Evidence context doğrudan `coreComplexes.txt` içindeki Human kayıtları
okur. Classic pickle aynı dosyanın organism filtresiz, primary-symbol uppercase
legacy dönüşümüdür ve membership-set semantiği tam yeniden üretilebilmiştir.
Clean Human-only builder Classic'e ayrı scientific migration kararı olmadan
geçirilmez. Ayrıntı: `docs/legacy_data_provenance.md`.

## Güncelleme komutu

Builderlar yalnız `data/derived/staging` altına yazar ve promotion yapmaz:

```powershell
python scripts/data_build.py --dataset hpa
python scripts/data_build.py --dataset omnipath
python scripts/data_build.py --dataset corum
python scripts/data_build.py --dataset trrust
python scripts/data_build.py --dataset all
```

Alternatif raw kaynak `--source`, izole staging kökü `--output-root` ile
verilebilir. İnternetten otomatik indirme yapılmaz.
