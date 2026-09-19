# Deterministic Biological Interpretation V2 — İlerleme

- Tamamlanan faz: Interpretation V2 gerçek UI/runtime debug düzeltmesi.
- Değiştirilen dosyalar: `src/interpretation/ontology_normalizer.py`, `src/interpretation/biological_facts.py`, `src/interpretation/theme_detector_v2.py`, `src/interpretation/sentence_bank.py`, `src/interpretation/candidate_interpreter.py`, `src/interpretation/system_interpreter.py`, `src/interpretation/biological_interpreter.py`, `src/ui/biological_interpretation.py`, `tests/test_interpretation_production_shape.py`.
- Mimari karar: Modern ham forward raporunda eksik olan GO/MyGene fact'leri yalnızca mevcut yerel `mygene_cache.pkl` içinden read-only tamamlanır. Geniş “membran taşınması” teması baskın ekseni gölgeleyemez; ana/ikincil eksen, pozitif ΔPageRank ağırlığı ve anotasyon desteğiyle seçilir. Aile tekrarları bağımsız kanıt sayılmaması için özetin içinde uyarılır.
- Test: `compileall` başarılı; 56/56 test başarılı. Production-shaped test gerçek ham CFTR CSV'siyle 382 yerel-MyGene annotation fallback'i, tema ve üst seviye Türkçe özeti doğruladı. Geçici read-only Streamlit harness'inde gerçek UI DOM'u; ana/ikincil eksen, destek adayları, aile uyarısı, sistem yanıtı, sınırlama ve doğrulama önerilerini görünür olarak doğruladı. Testler sırasında ham CSV hash'i değişmedi.
- Scientific engine: Değiştirilmedi; yeni veri, LLM veya harici runtime isteği kullanılmadı.
- Sıradaki faz: Tamamlandı; gerçek simülasyon çalıştırılmadan mevcut CFTR çıktısı yorumlandı ve UI render yolu doğrulandı.

- Tamamlanan faz: Signed ΔPageRank / Redistribution Losses görünümü.
- Değiştirilen dosyalar: `src/biology_logic.py`, `src/mouse/biology_logic.py`, `src/services/analysis_service.py`, `src/interpretation/biological_interpreter.py`, `src/interpretation/system_interpreter.py`, `src/interpretation/theme_detector_v2.py`, `src/ui/analysis_results.py`, `src/ui/biological_interpretation.py`, `src/ui/pages/04_ayarlar.py`, `app.py`, `tests/test_signed_redistribution_losses.py`.
- Mimari karar: Mevcut `classify_redistribution` çıktısının pozitif geriye-uyumluluk filtresinden önceki iki-yönlü kopyası yalnızca DataFrame `attrs` üzerinden UI/yorumlama katmanına taşınır. Disk CSV, pozitif ana tablo, seçimi ve export davranışı değişmez; negatif görünüm ayrı eşik/Top-N ile read-only seçilir.
- Test: `compileall` başarılı; 11/11 hedefli test geçti: signed seçim/sıralama, pozitif tablo değişmezliği, negative safe-language, canonical bileşen kontrolü, ayrı export ve UI handoff.
- Scientific engine: Formül, ağırlık, PageRank çağrısı ve pozitif seçim değişmedi; yalnız hesaplanmış signed sonucun UI handoff metadatası eklendi.
- Sıradaki faz: Import ve signed-loss unit testleri.
