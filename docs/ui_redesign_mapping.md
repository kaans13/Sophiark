# Sophiark UI redesign: old → new mapping

This mapping is the presentation-preservation contract for the redesigned root
workspace. Scientific calculations, ranking semantics, thresholds, and source
data are intentionally outside its scope.

| Previous surface | New destination | Preservation note |
|---|---|---|
| Permanent sidebar: workspace, organism, tissue, engine, target/mode, language, execution actions | `Analizi değiştir · Ayarlar` popover | Same widget keys and callbacks; context changes still clear incompatible active results. |
| Sidebar advanced simulation parameters, local/resource status, diagnostics, stress/compensation triggers, reset and threshold sweep | Popover → `Gelişmiş ayarlar` and existing secondary disclosures | Controls and defaults retained; the panel no longer consumes analysis width. |
| Result identity, target, tissue, organism, attenuation and engine status | Header + `Genel Bakış` | Current context remains visible while controls are closed. |
| Analysis summary, main finding, top candidates, network verdict, limitations | `Genel Bakış` | Four primary readouts are shown first; extended metrics remain in network/full data. |
| Target cards and candidate cards | `Genel Bakış` → closed `Hedef ve öne çıkan aday kartları` | Existing cards retained without an extra bordered wrapper. |
| Candidate rankings, positive/negative redistribution, PageRank, BC, gate, efficiency, response matrices | `Adaylar ve Ağ` | Existing source tables and network views retained. |
| Candidate annotations, pathway membership, network neighborhood and technical fields | `Adaylar ve Ağ` → candidate lens | One selected candidate is inspected at a time; `Adayın tüm kaynak alanları` preserves every field. |
| 2D/3D network visualization and display controls | `Adaylar ve Ağ` → `Sistem yanıtı` | Visualization limits affect display only, not calculations. |
| Tissue expression, localization, GO, KEGG, CORUM, functional and deep interpretation | `Biyolojik Bağlam` | Evidence/statistics stay beside the relevant interpretation; no separate Evidence tab exists. |
| Null/FDR support, validation and selection-mode explanations | `Biyolojik Bağlam` plus full source fields in `Tüm Sonuçlar` | Conditional unrun/unavailable states remain explicit. |
| Research Explorer, families, functions, relationships, local/external reference context and disease context | `Araştırma` | Snapshot binding, provider restrictions and complete research export retained. |
| Classic/Directed/Evidence comparisons and engine-specific full tables | `Tüm Sonuçlar` | Full schemas and exports retained; no candidate-set blending. |
| Main report, complete signed response, raw rows, provenance and technical fields | `Tüm Sonuçlar` | Positional pagination reaches all rows; source order and precision are unchanged. |
| Active target/context and advanced parameter values | `Tüm Sonuçlar` → `Analiz bağlamı ve kullanılan parametreler` | Read-only audit table; does not alter calculations. |
| Dose response, threshold sweep, propagation trace, tissue differential, severity sensitivity, target stress and predicted compensation | `İleri Analizler` | Existing implementations and conditional availability retained. |
| Human–Mouse comparison and Unified/bundle history | `İleri Analizler` and the on-demand workspace selector | Existing bundle validation and transient-state invalidation retained. |
| Per-table CSV/XLSX actions and complete research/evidence exports | Their relevant result section and `Tüm Sonuçlar` | No export action is removed. |
| Help text, interpretation guides and technical explanations | Contextual captions and closed `Nasıl yorumlanır?`/detail disclosures | Educational content remains available without dominating the first view. |

Conditional outputs keep their original guards. `app_state.clear_active_workspace`,
the completed-analysis context match, tool-context checks, research snapshot IDs,
and bundle validation continue to prevent stale results from being presented as
belonging to a changed configuration.
