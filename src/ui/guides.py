"""Clear Turkish in-product usage guidance for non-technical users."""

from __future__ import annotations

import streamlit as st
from src.ui.localization import locale


def render_quick_start() -> None:
    if locale() == "en":
        with st.expander("How to use Sophiark", expanded=True):
            st.markdown("""
1. **Choose the organism.** Human and Mouse use separate data sources and symbol mappings.
2. **Choose a tissue.** A tissue selection builds the stored-expression-filtered network; use `None` to inspect the unfiltered network.
3. **Choose the target.** Use a gene, a mapped disease gene set, or an available cellular-localization input.
4. **Choose the analysis.** Use Perturbation Analysis for target-edge attenuation or Target Stress Search for indirect candidate effects.
5. **Read results in layers.** Start with the summary and network map, then inspect gene cards and raw data as needed.
            """)
        return
    with st.expander("Sophiark nasıl kullanılır?", expanded=True):
        st.markdown("""
1. **Türü seçin.** İnsan ve fare ağları ayrı veri kaynakları ve sembol eşlemeleri kullanır.
2. **Doku seçin.** Doku seçimi, ağın ifade göstermeyen proteinlerden arındırılmış sürümünü oluşturur. Emin değilseniz önce `None` ile genel ağı inceleyin.
3. **Hedefinizi seçin.** Tek gen, birden fazla gen, hastalıkla ilişkili gen seti veya hücresel lokalizasyon kullanabilirsiniz.
4. **Analiz türünü belirleyin.** Bir hedefi baskılamanın sonucunu görmek için Perturbasyon Analizi; hedefe dolaylı etki edecek adayları bulmak için Target Stress Search kullanın.
5. **Sonuçları katmanlı inceleyin.** Önce üstteki özet ve ağ haritasına, sonra gen kartlarına; yalnızca gerektiğinde ham verilere inin.
        """)


def render_mode_guide(mode: str) -> None:
    if mode == "Perturbasyon Analizi" and locale() == "en":
        title = "When to use"
        description = "Use this analysis to examine which regions of the selected biological network are affected by perturbation of one or more targets, where network importance shifts, and which redistribution patterns may emerge."
        steps = [
            "For a single-target analysis, examine the selected gene's direct and indirect system effects on the network.",
            "For a multi-target analysis, compare shared and divergent network responses created by jointly suppressed targets.",
            "Use redistribution candidates to prioritize genes whose relative network importance increases or decreases after perturbation.",
            "Use default parameters for the first analysis; change thresholds and advanced settings only for controlled comparisons.",
            "Treat results as computational hypotheses for follow-up experiments and literature review, not as biological causality or treatment recommendations.",
        ]
        with st.expander(title, expanded=False):
            st.write(description)
            for step in steps:
                st.markdown(f"- {step}")
        return
    guides = {
        "Perturbasyon Analizi": ("Ne zaman kullanılır?", "Bir veya birden fazla hedefe uygulanan pertürbasyonun seçili biyolojik ağda hangi bölgeleri etkilediğini, ağ öneminin hangi genlere kaydığını ve olası yeniden-dağılım örüntülerini incelemek için kullanın.",
                                  ["Tek hedef analizinde, seçilen genin ağ üzerindeki doğrudan ve dolaylı sistem etkisini inceleyin.", "Çoklu hedef analizinde, birlikte baskılanan hedeflerin oluşturduğu ortak veya ayrışan ağ yanıtlarını karşılaştırın.", "Yeniden-dağılım adaylarını, pertürbasyon sonrasında ağ içindeki göreli önemi artan veya azalan genleri önceliklendirmek için kullanın.", "İlk analizde varsayılan parametreleri kullanın; eşik ve gelişmiş ayarları yalnızca kontrollü karşılaştırmalar için değiştirin.", "Sonuçları biyolojik nedensellik veya tedavi önerisi olarak değil, ileri deney ve literatür incelemesi için hesaplamalı hipotezler olarak değerlendirin."]),
        "Target Stress Search": ("Ne zaman kullanılır?", "Bir hedefi doğrudan perturb etmeden, hangi aday müdahalelerin hedef üzerinde daha yüksek öngörülen ağ stresi oluşturduğunu araştırmak istediğinizde kullanın.",
                                ["Tam olarak bir hedef gen seçin.", "İlk çalıştırmada aday limitini 5 bırakın.", "Doğrudan Etki sonucu yoksa bu negatif kanıt değildir; yalnızca yüklü yönlü kanıt bulunmadığını gösterir."]),
        "Compensation Analysis": ("Ne zaman kullanılır? · BETA", "Bir perturbasyonun ardından hangi genlerin ağ içindeki göreli rolünün artabileceğini keşfetmek istediğinizde kullanın.",
                                  ["Tam olarak bir hedef gen seçin.", "Adayları doğrulanmış biyolojik kompanzasyon olarak yorumlamayın.", "Öne çıkan adayları Target Stress Search veya ayrı forward simülasyonlarla çapraz inceleyin."]),
    }
    title, description, steps = guides[mode]
    with st.expander(title, expanded=False):
        st.write(description)
        for step in steps:
            st.markdown(f"- {step}")


def render_results_guide() -> None:
    with st.expander("Sonuçları nasıl okumalıyım?", expanded=False):
        st.markdown("""
- **Özet metrikleri:** Ağdaki genel etkiyi hızlıca gösterir; tek başına biyolojik doğrulama değildir.
- **Ağ haritası:** Etkinin hangi düğümler ve bağlantılar çevresinde yoğunlaştığını görmenizi sağlar.
- **Gen kartları:** Önce biyolojik yorumu okuyun; sonra Ağdaki Önemi, Geçiş Merkeziliği ve doku bilgisini karşılaştırın.
- **Yolak analizi:** Stresli genlerin hangi biyolojik süreçlerde toplandığını gösterir. Ağ etkisini doğrudan yolak aktivasyonu veya inhibisyonu olarak yorumlamayın.
- **Ham veri / export:** Tekrar analiz, raporlama veya bağımsız istatistiksel inceleme için kullanın.
        """)
