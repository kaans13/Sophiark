# Runtime execution acceptance

## Sonuç

**PASS.** Aktif simülasyon aileleri gerçek üretim verileriyle tamamlandı; motor koruması ve final regresyon geçti.

## Kanıt özeti

- `Unified Research Report (Classic + Directed)` — COMPLETE / PASS — unified_panel/{CFTR,TGFBR1,ANO2}
- `Classic/Directed standalone parity panel` — COMPLETE / PASS — ../runtime_acceptance_closure/bounded_parity.json
- `Multi-target perturbation` — COMPLETE / PASS — secondary_smoke/secondary_smoke.json
- `Null/FDR acceptance smoke` — COMPLETE / PASS — secondary_smoke/secondary_smoke.json
- `Dose response` — COMPLETE / PASS — ../runtime_acceptance_closure/dose_response/validation.json
- `Threshold Sweep` — COMPLETE / PASS — threshold_smoke/validation.json
- `Target Stress Search` — COMPLETE / PASS — secondary_smoke/secondary_smoke.json
- `Compensation BETA` — COMPLETE / PASS — secondary_smoke/secondary_smoke.json
- `Severity sensitivity` — COMPLETE / PASS — tests/test_forward_smoke.py
- `Multi-target network interaction` — COMPLETE / PASS — tests/test_forward_smoke.py
- `Propagation trace` — COMPLETE / PASS — tests/test_propagation_trace.py
- `Tissue Differential` — COMPLETE / PASS — tissue_smoke/validation.json
- `Research Explorer` — COMPLETE / PASS — research_validation.json
- `Human–Mouse comparison` — COMPLETE / PASS — ../../outputs_mouse/audits/mouse_cross_species_benchmark.json
- `Evidence BETA` — COMPLETE / PASS — evidence_validation.json
- `.sophiark bundle and history` — COMPLETE / PASS — bundle_history_check.json

## Güvenlik

- Engine freeze: 12/12 exact.
- Data health: 7 capability READY; başlangıç/bitiş fingerprint kayıtları saklandı.
- Final regression: 502 passed in 31.73s.

## Düzeltmeler

- Gen ve ENSP kimlikleri ayrı sunuluyor; mevcut bilimsel sıralama ve formüller değiştirilmedi.
- Farklı doku grafı hazırlıkları igraph'ın süreç-geneli RNG durumunu paylaşamayacak şekilde sıralandı.
- Türler arası tablo İnsan ENSP, Fare protein kimliği ve ilişki durumunu ayrı sunuyor.
- Evidence log, zenginleştirme ve downstream sonuç kimliği gerçek Evidence
  çağrısıyla eşleştirildi; hesaplama formülleri değiştirilmedi.
- Sıcak kaynak yenilemesinden kalan eski Research service nesnesi yeniden
  kurularak Araştırma yüzeyindeki sonuç kaybı engellendi.
