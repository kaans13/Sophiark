# Pre-UI/UX teknik hardening

Bu geçiş bilimsel algoritmaları değiştirmez. Reuse yalnız hedeften bağımsız,
immutable product kaynaklarında yapılır; motorlar detached çalışma grafikleri
almaya devam eder.

## Ölçülen performans

Tek Windows/Python 3.11 sürecinde, Lung dokusu ve `bc_sample_sources=4` ile:

| Senaryo | Önce | Sonra | Sonuç |
|---|---:|---:|---:|
| Cold CFTR | 33.95 sn | 33.42 sn | %1.5; ölçüm gürültüsü düzeyinde |
| Aynı doku ikinci hedef, TGFBR1 | 34.75 sn | 15.78 sn | %54.6 daha hızlı |
| Aynı doku üçüncü hedef, ANO2 | 35.17 sn | 15.75 sn | %55.2 daha hızlı |
| Unified cache miss, hazırlanmış kaynaklarla | 32.30 sn | 30.01 sn | %7.1 daha hızlı |
| Unified cache hit | 0.0011 sn | 0.0012 sn | eşdeğer |
| `.sophiark` load | 0.52–0.55 sn | 0.50–0.52 sn | eşdeğer |

En büyük maliyetler doku grafı hazırlığı (~14 sn), Classic (~7–8 sn) ve
Directed source/overlay/calculation (~8–10 sn) idi. Aynı process içinde doku
grafı hazırlığı ikinci hedefte kaldırıldı. Dataset fingerprint provenance
maliyeti ~2.6 saniyeden yaklaşık 0.003–0.005 saniyeye, OmniPath tekrar okuması
~1.3–1.7 saniyeden cache lookup düzeyine indi. Üç hedef sonrası gözlenen RSS
yaklaşık 1226 MB yerine 935 MB oldu; target-specific Directed graph LRU tek
öğeyle sınırlandı.

Parallel motor execution eklenmedi. Classic/Directed algoritmalarını,
yakınsama kriterlerini, BC yaklaşımını veya graph matematiğini değiştirecek
optimizasyonlar ayrı bilimsel doğrulama gerektirir.

## Presentation contract

`src.ui.metric_presentation` canonical alan adı → Türkçe label, açıklama,
birim, precision, semantic group, sıra, missing ve boolean gösterimini tutar.
`render_dataframe` yalnız gösterilen satırlar üzerinde bu contract'ı uygular;
export DataFrame'i, `.sophiark` içeriği ve raw bilimsel değerler değişmez.

Kolon düzeni identity → response → rank → comparison → evidence → context →
engine-specific şeklindedir. Yüzde olarak saklanan yanıtlar ile 0–1 fraction
alanlar metadata üzerinden ayrılır; büyüklükten tahmin edilmez. `False`
"Hayır", kaynağı olmayan boolean ise "Kullanılamıyor" gösterilir.

## Bağımlılık ve proje hijyeni

- Core, UI, analysis, benchmark ve test bağımlılıkları `pyproject.toml`
  gruplarında beyan edilmiştir. `requirements.txt` aynı uyumlu aralıkları
  sunar.
- Python 3.11 clean venv içinde `.[test]` kurulumu, `pip check`, dataset
  health, Directed preflight ve targeted testler doğrulanmıştır.
- `venv/`, `__pycache__/`, `outputs/`, cache, log ve backup kalıpları ignore
  edilir. Kullanıcının yaklaşık 0.93 GB mevcut `venv` dizini silinmemiştir.
- Raw/processed bilimsel veri, legacy HPA/CORUM, `.sophiark` geçmişi ve
  belirsiz root artefact'ları korunmuştur.
- Aktif runtime/config kodunda `C:\Users\` veya `OneDrive\` sabiti yoktur.
- Product facade lazy-load olur; focused product/Unified/presentation
  importları Streamlit'i veya motorları başlatmaz.

Kalan legacy presentation alanları özellikle eski cross-species/Research
Explorer tablolarında Türkçe display anahtarlarını ara DataFrame kolonları
olarak kullanabilir. Bunlar mevcut export tüketicilerini kırmamak için bu
geçişte zorla canonical schema'ya çevrilmemiştir ve UI/UX fazında adapter ile
ele alınmalıdır.
