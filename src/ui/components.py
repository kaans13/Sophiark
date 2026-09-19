"""Reusable product-facing UI components, independent from scientific logic."""

from __future__ import annotations

import base64
from html import escape
from functools import lru_cache
from pathlib import Path

import streamlit as st

from src.ui.design_tokens import BG_SURFACE, BORDER, render_badge
from src.ui.guides import render_quick_start
from src.ui.localization import locale, t


@lru_cache(maxsize=1)
def _brand_logo_data_uri() -> str:
    logo = Path(__file__).resolve().parents[2] / "assets" / "logo.png"
    if not logo.exists():
        return ""
    return "data:image/png;base64," + base64.b64encode(logo.read_bytes()).decode("ascii")


def render_command_identity(*, context: str = "", ready: bool = False) -> None:
    """Render the non-interactive portion of the compact application command bar."""
    context_html = (
        f'<div class="sk-command-context">{escape(context)}</div>' if context else ""
    )
    status_html = (
        f'<div class="sk-command-status"><span></span>{"Ready" if locale() == "en" else "Hazır"}</div>' if ready else ""
    )
    logo_uri = _brand_logo_data_uri()
    logo_html = (
        f'<span class="sk-brand-mark"><img src="{logo_uri}" alt="Sophiark logo"></span>'
        if logo_uri else '<span class="sk-brand-mark">S</span>'
    )
    st.markdown(
        '<div class="sk-command-identity">'
        f'<div class="sk-command-brand">{logo_html}'
        '<span><b>Sophiark</b><small>NETWORK BIOLOGY WORKSPACE</small></span></div>'
        f'{context_html}{status_html}</div>',
        unsafe_allow_html=True,
    )


def render_analysis_hero(*, target: str, species: str, tissue: str, engine: str) -> None:
    """Present the completed analysis context without changing any result source."""
    if locale() == "en":
        st.markdown(
            '<section class="sk-analysis-hero"><div class="sk-hero-copy">'
            '<div class="sk-hero-kicker"><span></span>ANALYSIS WORKSPACE</div>'
            f'<h1>{escape(target)}</h1><p>Network redistribution, biological context, and traceable source results.</p>'
            '</div><div class="sk-hero-context">'
            f'<div><small>ORGANISM</small><strong>{escape(species)}</strong></div>'
            f'<div><small>TISSUE</small><strong>{escape(tissue)}</strong></div>'
            f'<div><small>ENGINE</small><strong>{escape(engine)}</strong></div>'
            '</div></section>', unsafe_allow_html=True,
        )
        return
    st.markdown(
        '<section class="sk-analysis-hero">'
        '<div class="sk-hero-copy">'
        '<div class="sk-hero-kicker"><span></span>ANALİZ ÇALIŞMA ALANI</div>'
        f'<h1>{escape(target)}</h1>'
        '<p>Ağ yeniden dağılımı, biyolojik bağlam ve doğrulanabilir kaynak sonuçları.</p>'
        '</div>'
        '<div class="sk-hero-context">'
        f'<div><small>ORGANİZMA</small><strong>{escape(species)}</strong></div>'
        f'<div><small>DOKU</small><strong>{escape(tissue)}</strong></div>'
        f'<div><small>MOTOR</small><strong>{escape(engine)}</strong></div>'
        '</div></section>',
        unsafe_allow_html=True,
    )


def candidate_stack_html(rows: list[dict[str, str]]) -> str:
    """Human-readable candidate rows built exclusively from recorded values."""
    items = []
    for row in rows:
        tone = row.get("tone", "neutral")
        response = row.get("response", "")
        relation = row.get("relation", "")
        badges = "".join(
            f'<span class="sk-response-badge sk-response-{escape(kind)}">{escape(text)}</span>'
            for kind, text in ((tone, response), ("relation", relation))
            if text
        )
        items.append(
            f'<div class="sk-candidate-row sk-candidate-{escape(tone)}">'
            f'<div><div class="sk-candidate-symbol">{escape(row.get("symbol", "—"))}</div>'
            f'<div class="sk-response-badges">{badges}</div>'
            f'<div class="sk-candidate-meta">{escape(row.get("localization", "—"))}'
            f'<span>{"Network importance" if locale() == "en" else "Ağ önemi"}: {escape(row.get("importance", "—"))}</span></div></div>'
            f'<div class="sk-candidate-delta">{escape(row.get("delta", "—"))}</div>'
            '</div>'
        )
    return '<div class="sk-candidate-stack">' + "".join(items) + '</div>'


def render_empty_workspace() -> None:
    """Render an instructive, product-level first screen without touching analysis state."""
    if locale() == "en":
        st.markdown(
            '<section class="sk-onboarding"><div class="sk-onboarding-copy">'
            '<div class="sk-onboarding-kicker"><span></span>NETWORK BIOLOGY WORKSPACE</div>'
            '<h1>Choose a target.<br><em>See the network response.</em></h1>'
            '<p>Sophiark attenuates target edges and presents relative gains and losses of network importance in a tissue-aware biological context.</p>'
            '<div class="sk-onboarding-cue"><b>01</b> Use <strong>Start analysis</strong> in the upper-right corner <span>→</span></div>'
            '</div><div class="sk-onboarding-preview">'
            '<div class="sk-preview-head"><span>LIVE ANALYSIS FLOW</span><i>Ready</i></div>'
            '<div class="sk-preview-target"><small>EXAMPLE TARGET</small><strong>CFTR</strong><span>Homo sapiens · Lung · Classic</span></div>'
            '<div class="sk-preview-flow"><div><b>1</b><span>Target and tissue<small>Choose the analysis context</small></span></div>'
            '<div><b>2</b><span>Perturbation<small>The engine calculates the network response</small></span></div>'
            '<div><b>3</b><span>Response balance<small>Review gains and losses</small></span></div>'
            '</div></div></section>', unsafe_allow_html=True,
        )
        st.markdown(
            '<div class="sk-capability-grid">'
            '<article><span>01</span><b>Perturbation analysis</b><p>Compare network redistribution after attenuating one or more targets.</p></article>'
            '<article><span>02</span><b>Directed evidence</b><p>Read activation or inhibition records separately from the network response and within their source limits.</p></article>'
            '<article><span>03</span><b>Traceable results</b><p>Gene aliases, cellular location, recorded values, and exports remain in one workspace.</p></article>'
            '</div><div class="sk-scientific-note"><b>Analysis engines</b><br><b>Classic</b> calculates the primary undirected PPI edge-attenuation response. <b>Directed · BETA</b> adds available direction records as a separate signed-evidence view. <b>Evidence · BETA</b> presents evidence-oriented context while retaining the STRING network backbone. <b>Compare</b> places Classic and Directed outputs side by side for the same input.</div>',
            unsafe_allow_html=True,
        )
        return
    st.markdown(
        '<section class="sk-onboarding">'
        '<div class="sk-onboarding-copy">'
        '<div class="sk-onboarding-kicker"><span></span>NETWORK BIOLOGY WORKSPACE</div>'
        '<h1>Bir hedef seçin.<br><em>Ağın yanıtını görün.</em></h1>'
        '<p>Sophiark, hedef kenarlarını zayıflatır; dokuya özgü ağda artan ve azalan '
        'göreli önemi, biyolojik bağlamıyla birlikte görünür kılar.</p>'
        '<div class="sk-onboarding-cue"><b>01</b> Sağ üstte <strong>Analize başla</strong> '
        '<span>→</span></div>'
        '</div>'
        '<div class="sk-onboarding-preview">'
        '<div class="sk-preview-head"><span>CANLI ANALİZ AKIŞI</span><i>Hazır</i></div>'
        '<div class="sk-preview-target"><small>ÖRNEK HEDEF</small><strong>CFTR</strong>'
        '<span>Homo sapiens · Lung · Classic</span></div>'
        '<div class="sk-preview-flow">'
        '<div><b>1</b><span>Hedef ve doku<small>Analiz bağlamını seçin</small></span></div>'
        '<div><b>2</b><span>Pertürbasyon<small>Motor ağı hesaplar</small></span></div>'
        '<div><b>3</b><span>Yanıt dengesi<small>Artış ve kayıpları okuyun</small></span></div>'
        '</div></div></section>',
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div class="sk-capability-grid">'
        '<article><span>01</span><b>Pertürbasyon analizi</b><p>Bir veya daha fazla hedefin ağ üzerindeki yeniden dağılımını karşılaştırın.</p></article>'
        '<article><span>02</span><b>Yönlü kanıt</b><p>Aktivasyon/baskılama kaydını ağ yanıtından ayrı ve kaynak sınırlarıyla okuyun.</p></article>'
        '<article><span>03</span><b>İzlenebilir sonuç</b><p>Gen aliası, hücresel konum, ham değer ve dışa aktarma aynı çalışma alanında.</p></article>'
        '</div>',
        unsafe_allow_html=True,
    )


def highlight_card_html(*, title: str, body: str, meta: str = "", tone: str = "blue") -> str:
    """Shared level-one analytical highlight surface."""
    return (
        f'<div class="sk-highlight-card sk-highlight-{escape(tone)}">'
        f'<div class="sk-highlight-title">{escape(title)}</div>'
        f'<div class="sk-highlight-body">{body}</div>'
        + (f'<div class="sk-highlight-meta">{escape(meta)}</div>' if meta else "")
        + '</div>'
    )


def render_context_bar(*, species: str, tissue: str, beta: bool = False) -> None:
    beta_badge = render_badge("BETA", "beta") if beta else ""
    st.markdown(f'<div class="sk-context"><b>Aktif analiz</b>{render_badge(species, "kanit_guclu")}'
                f'{render_badge(tissue, "kanit_zayif")}{beta_badge}</div>', unsafe_allow_html=True)


def render_context_header(*, target: str, tissue: str, species: str, attenuation: float | None = None,
                          candidate_limit: int | None = None, engine: str | None = None,
                          ready: bool = False) -> None:
    """Compact, result-first analysis identity without exposing implementation details."""
    metadata = [species, tissue]
    if engine:
        metadata.append(engine)
    if attenuation is not None:
        metadata.append((f"Remaining edge weight: {attenuation * 100:.2f}%") if locale() == "en" else f"Kalan kenar ağırlığı: %{attenuation * 100:.2f}")
    if candidate_limit is not None:
        metadata.append(f"{candidate_limit} {'target' if locale() == 'en' else 'hedef'}")
    ready_html = f'<span class="sk-ready-dot">● {"Graph ready" if locale() == "en" else "Graf hazır"}</span>' if ready else ""
    st.markdown(
        '<div class="sk-context-header">'
        f'<div><h1>{target}</h1><div class="sk-section-label">{t("protein_network_perturbation")}</div>'
        f'<div class="sk-context-meta">{" · ".join(metadata)}</div></div>{ready_html}</div>',
        unsafe_allow_html=True,
    )


def render_section_header(title: str, description: str = "", *, label: str = "") -> None:
    label_html = f'<div class="sk-section-label">{label}</div>' if label else ""
    st.markdown(
        f'<div class="sk-section">{label_html}<h2>{title}</h2><p>{description}</p></div>',
        unsafe_allow_html=True,
    )


def render_scientific_note(text: str) -> None:
    st.markdown(f'<div class="sk-scientific-note">{text}</div>', unsafe_allow_html=True)


def render_landing(*, proteins: int, edges: int, enrichment_sources: int, species_count: int) -> str | None:
    if locale() == "en":
        logo = "assets/logo.png"
        try:
            st.image(logo, width=210)
        except Exception:
            pass
        st.markdown('<h1 class="sk-title">Sophiark</h1>', unsafe_allow_html=True)
        st.markdown('<p class="sk-subtitle">Tissue-aware protein-network perturbation analysis.</p>', unsafe_allow_html=True)
        st.markdown(
            '<div class="sk-scientific-note"><b>Choose an analysis engine</b><br>'
            '<b>Classic</b> is the primary undirected PPI edge-attenuation analysis. '
            '<b>Directed · BETA</b> separately applies available direction records. '
            '<b>Evidence · BETA</b> presents the evidence-oriented mode without changing the STRING backbone. '
            '<b>Compare</b> contrasts Classic and Directed results for the same input. '
            'Engine choice changes the analysis view; inspect its availability and provenance notes.</div>',
            unsafe_allow_html=True,
        )
        cards = [
            (":material/bolt:", "Perturbation Analysis", "Attenuate one or more targets and inspect network redistribution and critical changes."),
            (":material/search:", "Target Stress Search", "Find candidates that create the strongest network stress on the target without directly perturbing it."),
            (":material/monitoring:", "Compensation Analysis · BETA", "Inspect potential compensation candidates whose network role increases after perturbation."),
        ]
        selected = None
        tabs = st.tabs([item[1] for item in cards])
        for tab, (icon, title, description) in zip(tabs, cards):
            with tab:
                st.markdown(f'<div style="background:{BG_SURFACE};border:1px solid {BORDER};border-radius:10px;padding:18px;"><h3>{icon} {title}</h3><p>{description}</p></div>', unsafe_allow_html=True)
                if st.button("Start", key=f"landing_{title}", width="stretch"):
                    selected = title.replace(" · BETA", "")
        with st.expander("Network context", expanded=False):
            metrics = st.columns(4)
            for col, value, label in zip(metrics, (proteins, edges, enrichment_sources, species_count), ("Proteins", "Interactions", "Enrichment sources", "Supported species")):
                with col:
                    st.metric(label, f"{value:,}")
        render_quick_start()
        return selected

    st.markdown(
        '<div class="sk-scientific-note"><b>Analiz motorunu seçin</b><br>'
        '<b>Classic</b> ana yönsüz PPI kenar-zayıflatma analizidir. '
        '<b>Directed · BETA</b> mevcut yön kayıtlarını ayrı uygular. '
        '<b>Evidence · BETA</b> STRING omurgasını değiştirmeden kanıt odaklı görünümü sunar. '
        '<b>Compare</b> aynı girdi için Classic ve Directed sonuçlarını karşılaştırır. '
        'Motor seçimi analiz görünümünü değiştirir; kullanılabilirlik ve provenance notlarını inceleyin.</div>',
        unsafe_allow_html=True,
    )
    logo = "assets/logo.png"
    try:
        st.image(logo, width=210)
    except Exception:
        pass
    st.markdown('<h1 class="sk-title">Sophiark</h1>', unsafe_allow_html=True)
    st.markdown('<p class="sk-subtitle">Dokuya duyarlı protein ağı perturbasyon analizi.</p>', unsafe_allow_html=True)
    cards = [(":material/bolt:", "Perturbasyon Analizi", "Bir veya birden fazla geni baskılayın; ağdaki yayılımı ve kritik değişimleri görün."),
             (":material/search:", "Target Stress Search", "Hedefi doğrudan baskılamadan, üzerinde en güçlü ağ stresini oluşturan adayları bulun."),
             (":material/monitoring:", "Compensation Analysis · BETA", "Perturbasyon sonrası ağda rolü artan olası telafi adaylarını inceleyin.")]
    selected = None
    tabs = st.tabs([item[1] for item in cards])
    for tab, (icon, title, description) in zip(tabs, cards):
        with tab:
            st.markdown(f'<div style="background:{BG_SURFACE};border:1px solid {BORDER};border-radius:10px;padding:18px;"><h3>{icon} {title}</h3><p>{description}</p></div>', unsafe_allow_html=True)
            if st.button("Başla", key=f"landing_{title}", width="stretch"):
                selected = title.replace(" · BETA", "")
    with st.expander("Ağ bağlamı", expanded=False):
        metrics = st.columns(4)
        for col, value, label in zip(metrics, (proteins, edges, enrichment_sources, species_count),
                                     ("Protein", "Etkileşim", "Enrichment kaynağı", "Desteklenen tür")):
            with col:
                st.metric(label, f"{value:,}")
    render_quick_start()
    return selected
