"""Streamlit presentation for the read-only biological interpretation object."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from src.interpretation import BiologicalInterpretation
from src.interpretation.biological_rules import fmt_percent
from src.ui.components import render_scientific_note, render_section_header
from src.ui.localization import t
from src.ui.tables import render_dataframe


_STATISTICAL_VALIDATION_MESSAGES = {
    "null_fdr seçilmiş ancak gösterilecek istatistik alanı yok.": "statistical_validation_fields_unavailable",
}


def _statistical_validation_message(message: object) -> str:
    key = _STATISTICAL_VALIDATION_MESSAGES.get(str(message))
    return t(key) if key else str(message)


def _display_location(value: object) -> object:
    return " ".join(str(value).replace("_", " ").split()) if value is not None else value


def _candidate_table(rows: list[dict[str, object]]) -> pd.DataFrame:
    return pd.DataFrame([{
        "Gen": row["gene"],
        "ENSP": row.get("entity_id", "—"),
        "Ağ önemi değişimi (%)": fmt_percent(row.get("delta_pagerank_pct"), signed=True),
        "Başlangıç ağ önemi": row.get("pagerank_baseline") if row.get("pagerank_baseline") is not None else "veri yok",
        "Pertürbasyon sonrası ağ önemi": row.get("pagerank_perturbed") if row.get("pagerank_perturbed") is not None else "veri yok",
        "Ağdaki Önemi": row.get("hinterland_baseline") if row.get("hinterland_baseline") is not None else "veri yok",
        "Geçiş Merkeziliği": row.get("bc_baseline") if row.get("bc_baseline") is not None else "veri yok",
        "Hücresel konum": _display_location(row.get("localization")),
        "Yorum kapsamı": row.get("interpretation_confidence_label", "veri yok"),
    } for row in rows])


def _group_table(rows: list[dict[str, object]]) -> pd.DataFrame:
    return pd.DataFrame([{
        "Gen": row["gene"],
        "ENSP": row.get("entity_id", "—"),
        "Ağ önemi değişimi (%)": fmt_percent(row.get("delta_pagerank_pct"), signed=True),
        "Ağdaki Önemi": row.get("hinterland_baseline") if row.get("hinterland_baseline") is not None else "veri yok",
    } for row in rows])


def render_biological_interpretation(interpretation: BiologicalInterpretation) -> None:
    """Render all explanation sections without recalculating scientific outputs."""
    render_section_header(
        "Biyolojik yorum",
        "Mevcut ağ sonucu ve anotasyonları okuyan, sıralama/formül değiştirmeyen açıklama katmanı.",
        label="BİYOLOJİK YORUM · READ-ONLY",
    )
    render_scientific_note("Bu bölüm yeni bir bilimsel skor üretmez. Ağ motorunun oluşturduğu mevcut sonuçları ve mevcut anotasyon/evidence alanlarını açıklanabilir biçimde bir araya getirir.")

    summary = interpretation.analysis_summary
    st.markdown("### Ana yeniden dağılım")
    st.write(summary["primary"])
    if summary.get("secondary"):
        st.markdown("### İkincil biyolojik eksen")
        st.write(summary["secondary"])
    balance = interpretation.redistribution_balance
    st.markdown("### Yeniden dağılım dengesi")
    balance_cols = st.columns(4)
    balance_cols[0].metric("Pozitif aday", balance.get("positive_candidates", 0))
    balance_cols[1].metric("Negatif aday", balance.get("negative_candidates", 0))
    strongest_gain = balance.get("strongest_gain")
    strongest_loss = balance.get("strongest_loss")
    balance_cols[2].metric("En güçlü artış", f"{strongest_gain.get('gene')} {fmt_percent(strongest_gain.get('delta_pagerank_pct'), signed=True)}" if strongest_gain else "veri yok")
    balance_cols[3].metric("En güçlü kayıp", f"{strongest_loss.get('gene')} {fmt_percent(strongest_loss.get('delta_pagerank_pct'), signed=True)}" if strongest_loss else "veri yok")
    if balance.get("dominant_gain_themes"):
        st.caption("Baskın artış temaları: " + ", ".join(balance["dominant_gain_themes"]))
    if balance.get("dominant_loss_themes"):
        st.caption("Baskın kayıp temaları: " + ", ".join(balance["dominant_loss_themes"]))
    if balance.get("shared_themes"):
        render_scientific_note("Aynı anotasyon teması hem artış hem kayıp tarafında görülüyor: " + ", ".join(balance["shared_themes"]) + ". Bu, tema içi yeniden dağılımın heterojen olabileceğini gösterir; yolak durumu çıkarımı değildir.")
    if balance.get("consistency_warning"):
        st.warning(str(balance["consistency_warning"]))
    if summary.get("loss_summary"):
        st.markdown("### Negatif yeniden dağılım")
        st.write(summary["loss_summary"])
        if interpretation.loss_candidates:
            render_dataframe(
                _candidate_table(interpretation.loss_candidates)[["Gen", "Ağ önemi değişimi (%)", "Ağdaki Önemi", "Geçiş Merkeziliği", "Yorum kapsamı"]],
                label="biyolojik_negatif_ozet", key="bio_loss_summary", width="stretch", hide_index=True,
            )
    if summary.get("supporting_candidates"):
        supporting = set(summary["supporting_candidates"])
        rows = [row for row in interpretation.top_candidates if str(row["gene"]).upper() in supporting]
        st.markdown("### Bu yorumu destekleyen başlıca adaylar")
        render_dataframe(
            _candidate_table(rows)[["Gen", "ENSP", "Ağ önemi değişimi (%)", "Ağdaki Önemi", "Geçiş Merkeziliği", "Yorum kapsamı"]],
            label="biyolojik_destekleyen_adaylar", key="bio_supporting_candidates", width="stretch", hide_index=True,
        )
    if summary.get("family_warning"):
        render_scientific_note(str(summary["family_warning"]))
    if interpretation.canonical_component_check:
        with st.expander("Canonical bileşen kontrolü", expanded=False):
            render_dataframe(pd.DataFrame(interpretation.canonical_component_check).rename(columns={
                "gene": "Gen", "side": "Signed sonuç tarafı", "delta_pagerank_pct": "Ağ önemi değişimi (%)", "message": "Yorum",
            }), label="canonical_bilesen_kontrolu", key="bio_canonical_components", width="stretch", hide_index=True)
    if summary.get("tissue_note"):
        st.caption(summary["tissue_note"])
    st.markdown("### Sistem düzeyi ağ yanıtı")
    st.write(interpretation.system_response)
    st.markdown("### Sonraki doğrulama")
    for item in summary.get("validation", []):
        st.write(f"- {item}")
    render_scientific_note(interpretation.limitations)

    if st.session_state.get("ui_developer_mode", False):
        with st.expander("Interpretation debug", expanded=False):
            st.json(interpretation.debug, expanded=False)

    sections = st.tabs([
        "Analiz özeti", "Sistem yanıtı", "Öne çıkan adaylar", "Negatif adaylar", "Başlangıç / perturbasyon",
        "İşlevsel temalar", "Gen aileleri", "Yönlü kanıt", "Zenginleştirme", "Sınırlar",
    ])
    with sections[0]:
        st.write(interpretation.overview)
    with sections[1]:
        st.write(interpretation.system_response)
    with sections[2]:
        if interpretation.top_candidates:
            render_dataframe(
                _candidate_table(interpretation.top_candidates),
                label="biyolojik_one_cikan_adaylar", key="bio_top_candidates", width="stretch", hide_index=True, height=500,
            )
            with st.expander("Aday açıklamaları, anotasyon ve doğrulama önerileri", expanded=False):
                for row in interpretation.top_candidates:
                    st.markdown(f"**{row['gene']}** — {row['interpretation']}")
                    suggestions = row.get("validation_suggestions", [])
                    if suggestions:
                        st.caption(" ".join(str(item) for item in suggestions))
                render_dataframe(pd.DataFrame([{
                    "Gen": row["gene"], "İşlevsel anotasyon": row["annotation"], "Köprü etiketi (başlangıç)": row["bridge_label"],
                    "İlaçlanabilirlik skoru (yardımcı bağlam)": row["drug_score"] if row["drug_score"] is not None else "veri yok",
                    "Dış kanıt bağlamı": row["external_evidence"],
                } for row in interpretation.top_candidates]),
                    label="biyolojik_aday_anotasyonlari", key="bio_top_annotations", width="stretch", hide_index=True,
                )
        else:
            st.info("Ağ önemi değişimi alanı olmadığı için aday yorumu üretilemedi.")
    with sections[3]:
        if interpretation.loss_candidates:
            render_dataframe(
                _candidate_table(interpretation.loss_candidates),
                label="biyolojik_negatif_adaylar", key="bio_loss_candidates", width="stretch", hide_index=True, height=500,
            )
            with st.expander("Negatif aday açıklamaları", expanded=False):
                for row in interpretation.loss_candidates:
                    st.markdown(f"**{row['gene']}** — {row['interpretation']}")
        else:
            st.info("Mevcut signed sonuçta seçili eşik altında negatif aday yok.")
    with sections[4]:
        left, right = st.columns(2)
        with left:
            st.markdown("#### Perturbasyona duyarlı, başlangıçta baskın olmayan adaylar")
            st.caption("Aday kümesi içindeki yüzdelik gruplamadır; bilimsel sıralamayı değiştirmez.")
            sensitive = interpretation.baseline_vs_perturbation["perturbation_sensitive_baseline_nondominant"]
            render_dataframe(
                _group_table(sensitive), label="perturbasyona_duyarli_adaylar", key="bio_sensitive_candidates", width="stretch", hide_index=True,
            ) if sensitive else st.info("Bu gruba giren aday yok.")
        with right:
            st.markdown("#### Başlangıçta merkezi, perturbasyona zayıf yanıt veren adaylar")
            st.caption("Yüksek başlangıç merkeziliği ile düşük/orta mutlak ağ önemi değişimi birlikteliği aranır.")
            central = interpretation.baseline_vs_perturbation["baseline_central_weakly_responsive"]
            render_dataframe(
                _group_table(central), label="baslangicta_merkezi_adaylar", key="bio_central_candidates", width="stretch", hide_index=True,
            ) if central else st.info("Bu gruba giren aday yok.")
    with sections[5]:
        if interpretation.functional_themes:
            for theme in interpretation.functional_themes:
                st.markdown(f"**{theme['theme']}** · {theme['count']} aday")
                st.caption(", ".join(theme["genes"]))
            render_scientific_note("Temalar, yalnızca mevcut adayların GO/protein anotasyonlarında tekrarlanan ifadelerden çıkarılır; yolak durumu çıkarımı değildir.")
        else:
            st.info("Mevcut aday anotasyonlarında tekrar eden, tanımlı işlevsel tema bulunamadı.")
    with sections[6]:
        if interpretation.family_clusters:
            clusters = pd.DataFrame([{
                "Küme": item["name"], "Tür": item["cluster_type"], "Üyeler": ", ".join(item["members"]),
                "Güven": item["confidence_label"], "Kanıt": item["evidence"],
                "Topoloji çoğalma riski": item["topology_amplification_risk"],
            } for item in interpretation.family_clusters])
            render_dataframe(clusters, label="biyolojik_aile_kumeleri", key="bio_family_clusters", width="stretch", hide_index=True)
            for item in interpretation.family_clusters:
                if item.get("warning"):
                    st.caption(f"{item['name']}: {item['warning']}")
            render_scientific_note("Aile, paralog, protein kompleksi ve üst işlevsel sistem aynı kavram değildir. Bu kümeler yalnızca mevcut adayların kaynak-destekli bağlamını gösterir; sıralamaya geri beslenmez.")
        elif interpretation.family_patterns:
            for pattern in interpretation.family_patterns:
                st.markdown(f"**{pattern['family']} ailesi** · {pattern['count']} aday")
                st.caption(", ".join(pattern["genes"]))
            render_scientific_note("Tekrarlayan aile üyeleri ağ düzeyi bir ilişki örüntüsünü gösterir; aktivasyon veya telafi kanıtı değildir.")
        else:
            st.info("Öne çıkan adaylarda en az iki üyeli tanımlı bir gen ailesi örüntüsü yok.")
    with sections[7]:
        if interpretation.directed_evidence:
            render_dataframe(pd.DataFrame(interpretation.directed_evidence).rename(columns={
                "source": "Kaynak", "target": "Hedef", "effect": "Yönlü etki", "provenance": "Kaynak veri", "references": "Referanslar",
            }), label="yonlu_kanit", key="bio_directed_evidence", width="stretch", hide_index=True)
            render_scientific_note("Bu bölüm yalnızca TRRUST/OmniPath gibi gerçekten yönlü evidence kaydı bulunan ilişkileri gösterir. STRING kenarlarından yön türetilmez.")
        else:
            st.info("Ana hedef veya öne çıkan adaylarla eşleşen yönlü TRRUST/OmniPath kaydı yok. Bu, negatif biyolojik kanıt değildir.")
    with sections[8]:
        enrichment = interpretation.enrichment
        st.write(enrichment["message"])
        if enrichment.get("available") and enrichment.get("rows"):
            source = pd.DataFrame(enrichment["rows"])
            display = source.drop(columns=["adjusted_p"], errors="ignore").rename(columns={"term": "Terim", "genes": "Temsilci genler"})
            render_dataframe(
                display, label="biyolojik_zenginlestirme", key="bio_enrichment",
                width="stretch", hide_index=True, export_df=source,
            )
    with sections[9]:
        render_scientific_note(interpretation.limitations)
