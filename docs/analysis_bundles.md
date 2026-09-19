# Portable `.sophiark` analysis bundles

Bir `.sophiark` dosyası tek bir tamamlanmış Unified Research Report'un kalıcı,
paylaşılabilir bilimsel kaydıdır. ZIP-compatible container içindeki şeffaf v1
yapı şunları içerir:

- `manifest.json`: analysis UUID, target/tissue/status, engine/config sürümleri,
  dataset sürüm ve fingerprint'leri, capability ve legacy HPA/CORUM kimliği;
- `report.json`: agreement, provenance, availability, limitation ve hata durumu;
- `tables/*.json`: Unified candidates, Classic raw, Directed raw ve CORUM
  complex context tabloları; ham sayısal değerler ve dtype şemasıyla;
- `integrity.json`: her internal dosyanın SHA-256 değeri.

Public format pickle/joblib içermez. Loader arşiv üyelerini doğrudan okur;
dosya çıkartmaz, code çalıştırmaz ve Classic/Directed/Evidence motorlarını
çağırmaz. Path traversal, duplicate archive isimleri, aşırı açılmış boyut,
eksik/hash'i bozuk içerik ve desteklenmeyen schema güvenli biçimde reddedilir.

## History ve cache farkı

Bundle durable result'tır; scientific computation cache değildir. Yeni analiz
isteklerinde görünmez cache hit olarak kullanılmaz. Varsayılan history dizini
`outputs/history` olup `SOPHIARK_HISTORY_DIR` ile taşınabilir. SQLite
`history.sqlite3` yalnız hızlı listeleme indexidir ve scientific tabloları
kopyalamaz. Silinir veya bozulursa `.sophiark` dosyaları taranarak
`HistoryService.rebuild_history_index()` ile yeniden oluşturulur.

Unified UI başarılı veya açıkça partial bir sonucu, mevcut `UnifiedResearchReport`
nesnesinden otomatik kaydeder; motorları save için yeniden çalıştırmaz. Bundle
önce geçici dosyaya yazılır, doğrulanır, atomik olarak finalize edilir ve ancak
sonra indexlenir. Indexleme başarısız olsa da bundle source-of-truth olarak kalır.

## Açma, import ve eski veriler

External bundle kopyalanmadan açılabilir veya managed history'ye import
edilebilir. Aynı `analysis_id` + aynı content `ALREADY_IMPORTED`, aynı ID + farklı
content `ID_CONFLICT` üretir; hiçbir dosya sessizce overwrite edilmez.

Bundle'daki dataset fingerprint'i current registry'den farklıysa kayıt bozuk
sayılmaz. `SavedUnifiedResult.is_historical` işaretlenir ve saved result aynen,
yeniden hesaplanmadan ortak Unified renderer ile gösterilir. Gelecekteki schema
geçişleri `migrate_to_current_schema` sınırına eklenecektir; current
`bundle_schema_version` değeri `1`'dir.

Eski `session_state` history satırları tam scientific tabloları taşımadığından
`LEGACY_HISTORY_NOT_BUNDLE_COMPATIBLE` kabul edilir ve silinmez.
