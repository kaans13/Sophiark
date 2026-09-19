"""Researcher-facing, read-only views over completed scientific outputs.

No engine imports, inferred rankings or biological conclusions belong here.
"""
from __future__ import annotations

from html import escape
import re
import pandas as pd
import streamlit as st

from src.ui.metric_presentation import format_metric, metric_label, is_missing
from src.ui.localization import locale
from src.ui.tables import render_dataframe


def matching_enrichment(enrichment: pd.DataFrame | None, symbol: str) -> pd.DataFrame:
    """Exact recorded membership only: substring proximity is not evidence."""
    if not isinstance(enrichment, pd.DataFrame) or "Genes" not in enrichment or not symbol:
        return pd.DataFrame()
    def matches(value):
        tokens = value if isinstance(value, (list, tuple, set)) else re.split(r"[;,|\s]+", str(value))
        return symbol.casefold() in {str(token).strip().casefold() for token in tokens}
    return enrichment.loc[enrichment["Genes"].map(matches)].copy(deep=True)


def search_results(frame: pd.DataFrame, query: str) -> pd.DataFrame:
    """Literal search across recorded cells, preserving full schema and row order."""
    if not query.strip():
        return frame.copy(deep=True)
    mask = pd.Series(False, index=frame.index)
    for column in frame.columns:
        mask |= frame[column].astype(str).str.contains(query.strip(), case=False, regex=False, na=False)
    return frame.loc[mask].copy(deep=True)


def align_candidate_evidence(frame: pd.DataFrame, signed: pd.DataFrame | None) -> pd.DataFrame:
    """Fill missing candidate fields from exact identity matches in signed output.

    Recorded primary values win on conflicts. Neither input nor row order changes.
    """
    if not isinstance(signed, pd.DataFrame) or signed.empty or "gene" not in signed or "gene" not in frame:
        return frame.copy(deep=True)
    detail = signed.drop_duplicates("gene", keep="first").set_index("gene")
    result = frame.copy(deep=True)
    for column in detail:
        values = result["gene"].map(detail[column])
        if column not in result:
            result[column] = values
        else:
            result[column] = result[column].where(result[column].notna(), values)
    return result


def result_destinations():
    """Task-oriented reading route; switching destinations does not run Python."""
    if locale() == "en":
        st.markdown(
            '''<section class="sk-result-route" aria-label="Result reading route">
<header><span>RESULT READING ROUTE</span><small>Three tasks provide different levels of examination for the same computed result.</small></header>
<div class="sk-result-route-grid">
  <div class="primary"><b>01</b><span><strong>Understand the result</strong><small>First review what changed and the response balance.</small></span></div>
  <div class="primary"><b>02</b><span><strong>Explore the network</strong><small>Inspect affected genes and their biological context.</small></span></div>
  <div><b>03</b><span><strong>Validate evidence</strong><small>Continue to sources, recorded data, and advanced tools.</small></span></div>
</div></section>''', unsafe_allow_html=True,
        )
        summary_stage, explore_stage, verify_stage = st.tabs([
            "1 · Understand the result", "2 · Explore the network", "3 · Validate evidence",
        ])
        with summary_stage:
            st.markdown('<div class="sk-tab-intent"><b>Starting point</b><span>Review the response magnitude, directional balance, and key findings.</span></div>', unsafe_allow_html=True)
        with explore_stage:
            st.markdown('<div class="sk-tab-intent"><b>Exploration</b><span>Review candidates and the network response first, then inspect their biological context.</span></div>', unsafe_allow_html=True)
            network_tab, biology_tab = st.tabs(["Network Response", "Biological Context"])
        with verify_stage:
            st.markdown('<div class="sk-tab-intent"><b>Validation</b><span>Move from interpretation to recorded sources; exports and advanced tools are available here.</span></div>', unsafe_allow_html=True)
            research_tab, data_tab, tools_tab = st.tabs(["Evidence & sources", "Raw data & download", "Advanced tools"])
        return summary_stage, network_tab, biology_tab, research_tab, data_tab, tools_tab
    st.markdown(
        '''<section class="sk-result-route" aria-label="Sonuçları okuma rotası">
<header><span>SONUÇLARI OKUMA ROTASI</span><small>Üç ana görev, aynı hesaplanmış sonucun farklı inceleme derinlikleridir.</small></header>
<div class="sk-result-route-grid">
  <div class="primary"><b>01</b><span><strong>Sonucu anla</strong><small>Önce ne değiştiğini ve yanıt dengesini görün.</small></span></div>
  <div class="primary"><b>02</b><span><strong>Ağı keşfet</strong><small>Etkilenen genleri ve biyolojik bağlamı inceleyin.</small></span></div>
  <div><b>03</b><span><strong>Kanıtı doğrula</strong><small>Kaynaklara, ham verilere ve ileri araçlara geçin.</small></span></div>
</div></section>''',
        unsafe_allow_html=True,
    )
    summary_stage, explore_stage, verify_stage = st.tabs([
        "1 · Sonucu Anla", "2 · Ağı Keşfet", "3 · Kanıtı Doğrula",
    ])
    with summary_stage:
        st.markdown(
            '<div class="sk-tab-intent"><b>Başlangıç noktası</b><span>Yanıtın büyüklüğünü, yön dengesini ve en önemli bulguları okuyun.</span></div>',
            unsafe_allow_html=True,
        )
    with explore_stage:
        st.markdown(
            '<div class="sk-tab-intent"><b>Keşif</b><span>Önce aday ve ağ yanıtını, ardından bu kayıtların biyolojik bağlamını inceleyin.</span></div>',
            unsafe_allow_html=True,
        )
        network_tab, biology_tab = st.tabs(["Ağ Yanıtı", "Biyolojik Bağlam"])
    with verify_stage:
        st.markdown(
            '<div class="sk-tab-intent"><b>Doğrulama</b><span>Yorumdan ham kaynağa ilerleyin; dışa aktarma ve ileri araçlar burada bulunur.</span></div>',
            unsafe_allow_html=True,
        )
        research_tab, data_tab, tools_tab = st.tabs([
            "Kanıt & Kaynaklar", "Ham Veri & İndirme", "İleri Araçlar",
        ])
    return summary_stage, network_tab, biology_tab, research_tab, data_tab, tools_tab


def _value(row: pd.Series, field: str) -> str:
    value = row.get(field)
    if is_missing(value):
        return "Bu sonuçta mevcut değil"
    if isinstance(value, (list, tuple, set)):
        return "; ".join(map(str, value))
    if field == "Lokalizasyon":
        return str(value).replace("_", " ")
    if field == "Hasar_Tipi":
        return {"Birincil_Hasar_Hedef": "Pertürbasyon hedefi", "Üçüncül_Hasar_Stresli": "Yeniden-dağılım adayı"}.get(str(value), str(value))
    return format_metric(field, value)


def evidence_html(row: pd.Series, fields: tuple[str, ...]) -> str:
    entries = []
    for field in fields:
        if field in row and not is_missing(row[field]):
            label = {"Hasar_Tipi": "Sonuç sınıfı", "Selection_Reason": "Seçilme nedeni", "Response_Direction": "Yanıt yönü", "Düzenleyici_TFler": "Kayıtlı düzenleyiciler", "Tissue": "Doku"}.get(field, metric_label(field))
            entries.append(f'<dt>{escape(label)}</dt><dd>{escape(_value(row, field))}</dd>')
    return '<dl class="sk-evidence-list">' + ''.join(entries) + '</dl>' if entries else '<p>Bu sonuçta kayıtlı kanıt yok.</p>'


def _candidate_rows_html(frame: pd.DataFrame, *, identity: str, symbol_field: str, limit: int = 10) -> str:
    delta_field = "Delta_PageRank_Pct" if "Delta_PageRank_Pct" in frame else "Directed_Redistribution_Pct"
    deltas = pd.to_numeric(frame.get(delta_field, pd.Series(index=frame.index, dtype=float)), errors="coerce")
    finite = deltas[deltas.notna()]
    scale = max(float(finite.abs().max()), 1e-12) if not finite.empty else 1.0
    rows = []
    for position, (_, row) in enumerate(frame.head(limit).iterrows(), start=1):
        delta = pd.to_numeric(pd.Series([row.get(delta_field)]), errors="coerce").iloc[0]
        tone = "gain" if pd.notna(delta) and delta > 0 else "loss" if pd.notna(delta) and delta < 0 else "neutral"
        delta_label = "—" if pd.isna(delta) else f"{float(delta):+.2f}%"
        bar = 0 if pd.isna(delta) else max(4.0, min(100.0, abs(float(delta)) / scale * 100.0))
        importance = pd.to_numeric(pd.Series([row.get("Hinterland_Skoru")]), errors="coerce").iloc[0]
        location = str(row.get("Lokalizasyon", "—")).replace("_", " ").replace("|", " · ")
        gateway = "Compartment Bottleneck" if str(row.get("Gümrük_Kapisi", "")).casefold() in {"true", "1", "✓"} else "Standart düğüm"
        rows.append(
            f'<div class="sk-candidate-board-row sk-object-{tone}">'
            f'<span class="sk-board-rank">{position:02d}</span>'
            f'<span class="sk-board-gene"><b>{escape(str(row.get(symbol_field, row.get(identity, "—"))))}</b><small>{escape(str(row.get(identity, "—")))}</small></span>'
            f'<span class="sk-board-response"><b>{escape(delta_label)}</b><i><em style="width:{bar:.1f}%"></em></i></span>'
            f'<span class="sk-board-baseline"><small>AĞDAKİ ÖNEMİ</small><b>{"—" if pd.isna(importance) else f"{float(importance):.1f}"}</b></span>'
            f'<span class="sk-board-location">{escape(location)}</span>'
            f'<span class="sk-board-gateway">{escape(gateway)}</span>'
            '</div>'
        )
    return (
        '<div class="sk-candidate-board"><div class="sk-candidate-board-head">'
        '<span>#</span><span>ADAY</span><span>AĞ YANITI</span><span>BAŞLANGIÇ</span><span>HÜCRESEL KONUM</span><span>ROL</span>'
        '</div>' + "".join(rows) + '</div>'
    )


@st.fragment
def render_candidate_lens(frame: pd.DataFrame, enrichment: pd.DataFrame | None = None, *, key="candidate_lens"):
    """Progressive candidate browser over recorded values; no engine reruns."""
    if frame is None or frame.empty:
        st.info("Bu sonuçta incelenecek aday yok.")
        return
    identity = next((field for field in ("gene", "entity_id", "entity", "Protein ID", "Aday") if field in frame), None)
    if identity is None:
        st.info("Bu sonuçta aday kimliği bulunmuyor; tüm değerlere tam tablodan ulaşabilirsiniz.")
        return
    candidate_rows = frame
    if "Hasar_Tipi" in frame:
        mask = frame["Hasar_Tipi"].astype(str).eq("Üçüncül_Hasar_Stresli")
        if mask.any():
            candidate_rows = frame.loc[mask]
    symbol_field = next((field for field in ("Symbol", "symbol", "gene_symbol") if field in candidate_rows), identity)

    st.markdown("### Aday gezgini")
    st.caption("Önce yanıt büyüklüğünü ve başlangıç bağlamını tarayın; yalnız gerekirse kayıt ayrıntısını açın.")
    st.markdown(
        _candidate_rows_html(candidate_rows, identity=identity, symbol_field=symbol_field),
        unsafe_allow_html=True,
    )
    if len(candidate_rows) > 10:
        st.caption(f"İlk 10 / {len(candidate_rows)} aday · Tam sonuç aşağıdaki ağ yanıtında ve Kanıtı Doğrula › Ham Veri & İndirme bölümünde korunur.")

    with st.expander("Seçili adayın kayıtlı ayrıntıları", expanded=False):
        ids = frame[identity].astype(str).tolist()
        symbols = frame.get("Symbol", frame.get("symbol", frame.get("gene_symbol", frame[identity]))).astype(str).tolist()
        selected = st.selectbox(
            "İncelenecek kayıt", range(len(frame)),
            format_func=lambda index: f"{symbols[index]} · {ids[index]}", key=f"{key}_selected",
        )
        row = frame.iloc[selected]
        metric_cols = st.columns(3)
        metric_cols[0].metric("Gen", symbols[selected])
        metric_cols[1].metric("Ağ yanıtı", _value(row, "Delta_PageRank_Pct") if "Delta_PageRank_Pct" in row else _value(row, "Directed_Redistribution_Pct"))
        metric_cols[2].metric("Ağdaki Önemi", _value(row, "Hinterland_Skoru"))

        response_tab, biology_tab, evidence_tab = st.tabs(["Hesaplanan yanıt", "Biyolojik bağlam", "Kanıt ve kaynak"])
        with response_tab:
            st.markdown(evidence_html(row, ("Hasar_Tipi", "Delta_PageRank_Pct", "Directed_Redistribution_Pct", "PageRank_Baseline", "PageRank_Perturbed", "Selection_Reason", "Response_Direction", "Hinterland_Skoru", "BC_Skoru", "Directed_BC", "Topluluk_ID", "Gümrük_Kapisi", "Lokalizasyon")), unsafe_allow_html=True)
            st.caption("Göreli ağ önemi değişimi; ekspresyon, aktivasyon veya biyolojik işlev ölçümü değildir.")
        with biology_tab:
            st.markdown(evidence_html(row, ("Protein_Adi", "MyGene Adı", "GO Biyolojik Süreç", "GO Moleküler İşlev", "GO Hücresel Bileşen", "GO_CC_Terimleri", "GO_MF_Terimleri", "Düzenleyici_TFler", "Essentiality", "Hedef Gen Etkisi", "Tissue")), unsafe_allow_html=True)
            pathways = matching_enrichment(enrichment, symbols[selected])
            if not pathways.empty:
                show_pathways = st.checkbox(
                    "Bu adayı içeren zenginleştirme kayıtlarını göster",
                    key=f"{key}_show_pathways",
                )
                if show_pathways:
                    render_dataframe(pathways, label="aday_yolak_kaniti", key=f"{key}_pathways", row_limit=100, paginate=True)
            else:
                st.caption("Bu aday için mevcut zenginleştirme tablosunda eşleşen üyelik yok.")
        with evidence_tab:
            st.markdown(evidence_html(row, ("empirical_p", "q_value", "significant_redistribution", "Evidence_Provenance_Status", "Non_Text_Evidence", "Yapisal_Kanit", "Relation", "Support")), unsafe_allow_html=True)
            st.caption("Eksik kayıt, ilişkinin yokluğuna kanıt değildir. Anotasyon bir pertürbasyon sonucu değildir.")
            show_all_fields = st.checkbox(
                "Tüm kaynak alanlarını göster",
                key=f"{key}_show_all_fields",
            )
            if show_all_fields:
                render_dataframe(
                    pd.DataFrame({"Kaynak alan": row.index, "Kaydedilen değer": [str(v) for v in row.values]}),
                    label="aday_tum_alanlar", key=f"{key}_all_fields", row_limit=200, paginate=True,
                )


def render_reading_guide():
    if locale() == "en":
        st.markdown('''<div class="sk-reading-guide">
<div><b>COMPUTED</b><p>Network response after target edges are attenuated.</p></div>
<div><b>CONTEXT / EVIDENCE</b><p>Tissue, annotations, enrichment, and available statistical support.</p></div>
<div><b>HYPOTHESIS</b><p>Compensation or bypass candidates require follow-up; they do not validate a mechanism.</p></div>
</div>''', unsafe_allow_html=True)
        return
    st.markdown('''<div class="sk-reading-guide">
<div><b>HESAPLANAN</b><p>Hedef kenarları zayıflatıldıktan sonraki ağ yanıtı.</p></div>
<div><b>BAĞLAM / KANIT</b><p>Doku, anotasyon, zenginleştirme ve varsa istatistiksel destek.</p></div>
<div><b>HİPOTEZ</b><p>Telafi veya bypass adayları ileri inceleme gerektirir; mekanizma doğrulaması değildir.</p></div>
</div>''', unsafe_allow_html=True)
@st.fragment
def render_research_download(bundle, *, config=None):
    """Prepare the complete existing export only on explicit demand."""
    from src.research.ui import ResearchExportFormat, research_explorer_download_payload
    selected = st.radio("Tüm araştırma çıktıları için format",
        tuple(item.value for item in ResearchExportFormat), horizontal=True,
        key="desktop_research_export_format")
    identity = (bundle.snapshot_id, selected)
    cached = st.session_state.get("desktop_research_export")
    if st.button("Tam araştırma dosyasını hazırla", key="desktop_prepare_research_export"):
        with st.spinner("Tüm araştırma tabloları dışa aktarılıyor…"):
            try:
                payload = research_explorer_download_payload(bundle, ResearchExportFormat(selected), config=config)
                cached = (identity, payload)
                st.session_state["desktop_research_export"] = cached
            except RuntimeError:
                st.warning("Bu format hazırlanamadı; CSV paketi seçerek yeniden deneyin.")
    if cached is not None and cached[0] == identity:
        filename, mime, data = cached[1]
        st.download_button("Tüm araştırma çıktılarını indir", data=data,
            file_name=filename, mime=mime, key="desktop_download_research", on_click="ignore")
    else:
        st.caption("Tam dışa aktarım yalnızca istendiğinde hazırlanır; tüm araştırma tablolarını içerir.")
