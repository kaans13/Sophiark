# Engine call provenance audit

## Sonuç

**PASS.** Kararlı referans `genom-29-agustos` ile güncel uygulamanın çağrı
zincirleri karşılaştırıldı. UI yeniden düzenlenmiş olsa da ana bilimsel
giriş noktaları korunuyor; sonuçlar tek bir tamamlanmış çalıştırma
bağlamından downstream yüzeylere aktarılıyor.

## Mod bazında gerçek çağrı zinciri

| UI modu | Bilimsel çağrı | Ekrana bağlanan kaynak | Durum |
|---|---|---|---|
| Human Classic | `analysis_service.execute_simulation` → `src.main.run_infection_simulation` | Classic rapor + Classic signed response | PASS |
| Human Directed | Classic baz çağrısı → `directed_analysis.run_optional_directed_engine` | Directed presentation + Directed full response | PASS |
| Human Evidence | `analysis_service.execute_simulation` → `src.evidence.main.run_infection_simulation` → `run_evidence_simulation` | Evidence rapor + Evidence signed response | PASS |
| Human Compare | Aynı girdilerle Classic baz çağrısı → Directed çağrısı → karşılaştırma | Classic ana yüzey + ayrı Classic/Directed karşılaştırması | PASS |
| Mouse Classic | `analysis_service.execute_simulation` → `src.mouse.main.run_infection_simulation` | Mouse Classic rapor + Mouse signed response | PASS |
| Unified | `src.unified.run_unified_analysis` → Classic + Directed | Bundle içindeki ayrı, karıştırılmamış sonuçlar | PASS |

Directed ve Compare modlarında Classic'in önce çalışması yanlış
yönlendirme değildir: Directed yeniden-dağılımı ve aynı girdiyle
karşılaştırma için kayıtlı baz sonuçtur. Evidence ise ayrı adaptör ve
motor üzerinden çalışır.

## Kararlı sürüm karşılaştırması

- 12 korumalı bilimsel kaynağın 11'i kararlı referansla byte-identical.
- Tek fark `src/biology_logic.py` içindeki doz-yanıt fonksiyonunda kararlı
  sürümde hesaplamadan önce dönen erken `return` satırının kaldırılmış
  olmasıdır. Bu formül değişikliği değil, ulaşılamayan mevcut hesaplamayı
  çalıştıran kontrol-akışı onarımıdır.
- Kararlı CFTR · Human · Lung · Classic CSV'si ile güncel kabul
  çalıştırması aynı 21 gen kaydını verdi. Kontrol edilen temel sayısal
  alanlarda en büyük mutlak fark `1.95e-11` (yalnızca kayan nokta toleransı).
- Güncel engine-freeze manifesti: **12/12 PASS**.

## Canlı CFTR · Human · Lung kanıtı

- Classic: 21 rapor satırı, 20 yeniden-dağılım adayı; ekran kaynağı `classic`.
- Directed: 21 rapor satırı, 20 Directed adayı; yedi aşamalı Directed
  pipeline tamamlandı ve ekran kaynağı `directed` oldu.
- Evidence: 101 rapor satırı, 100 Evidence adayı; zenginleştirme,
  yorumlama ve Research snapshot kaynaklarının tamamı `evidence` oldu.
- Compare: aynı hedef/doku/parametrelerle Classic ve Directed çalıştı;
  ayrı karşılaştırma tablosu oluştu.

## Yapılan sunum/izlenebilirlik düzeltmeleri

- Evidence çağrısını yanlış `CLASSIC` diye etiketleyen log düzeltildi.
- Evidence adaylarının zenginleştirme fingerprint'i artık `evidence`
  kimliğiyle saklanıyor; Biyoloji yüzeyi doğru zenginleştirmeyi gösteriyor.
- Streamlit sıcak kaynak yenilemesinden kalan eski Research service nesnesi
  güvenli biçimde yeniden kuruluyor; simülasyon sonucu değiştirilmiyor.

## Regresyon

- Tam test matrisi: **502 passed**.
- Aktif simülasyon aileleri kabul matrisi: **16/16 PASS**.
- Bilimsel motor freeze: **12/12 PASS**.

