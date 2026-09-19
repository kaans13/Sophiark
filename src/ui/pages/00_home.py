"""Intent-led landing page backed by real local dataset counts."""

from __future__ import annotations

from pathlib import Path

import streamlit as st


from src.ui.app_constants import HUMAN_ACTIVE_EDGE_COUNT, HUMAN_ACTIVE_PROTEIN_COUNT


logo = Path("assets/logo.png")
if logo.exists():
    st.image(str(logo), width=180)
st.title("Sophiark")
st.markdown("### Ne yapmak istiyorsunuz?")
st.caption("Dokuya duyarlı protein ağlarında müdahale etkisini simüleyin, hedef adaylarını keşfedin ve sonuçları karşılaştırın.")

left, middle, right = st.columns(3)
with left:
    with st.container(border=True):
        st.markdown("#### :material/bolt: Bir müdahaleyi simüle et")
        st.write("Tek gen, çoklu gen, hastalık hedefi veya hücresel lokalizasyon perturbasyonunun ağ etkisini inceleyin.")
        if st.button("Analize başla", type="primary", width="stretch"):
            st.switch_page("src/ui/pages/01_analiz.py")
with middle:
    with st.container(border=True):
        st.markdown("#### :material/search: Bir hedefi etkileyecek adayları bul")
        st.write("Target Stress Search ile adayları gerçek perturbasyon simülasyonları üzerinden doğrulayın.")
        if st.button("Keşfe geç", width="stretch"):
            st.switch_page("src/ui/pages/02_kesif.py")
with right:
    with st.container(border=True):
        st.markdown("#### :material/monitoring: Ağın yanıtını incele")
        st.write("Compensation Analysis · BETA ile perturbasyon sonrası rol kazanan hesaplamalı adayları görün.")
        if st.button("Yanıt analizini aç", width="stretch"):
            st.switch_page("src/ui/pages/02_kesif.py")

metrics = st.columns(4)
metrics[0].metric("Protein", f"{HUMAN_ACTIVE_PROTEIN_COUNT:,}")
metrics[1].metric("STRING etkileşimi", f"{HUMAN_ACTIVE_EDGE_COUNT:,}")
metrics[2].metric("Enrichment ailesi", "2")
metrics[3].metric("Desteklenen tür", "2")
st.caption("Protein ve etkileşim sayıları yerel insan STRING snapshot'ında combined score ≥ 500 için doğrulanmış metadata'dır.")

with st.expander("Hızlı başlangıç", expanded=True):
    st.markdown(
        """
1. Üstteki bağlam çubuğundan **tür** ve **doku** seçin.
2. **Analiz** sayfasında CFTR, TP53 veya GLP1R gibi bir hedef bulun.
3. Simülasyonu başlatın; ağ hazırlığı ve hesaplama ilerlemesini ekranda izleyin.
4. Önce **Öne Çıkan Bulgu** ve özet metrikleri okuyun; sonra 2B/3B ağ, kritik genler, KEGG/GO ve teknik kanıta inin.
5. Tam tabloyu CSV veya Excel olarak **Rapor ve Export** sekmesinden indirin.
        """
    )
