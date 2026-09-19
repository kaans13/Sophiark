"""
bio_interpreter.py — Kural Tabanlı Biyolojik Yorumlama Motoru (KTBYM)
Sophiark'ın KEGG/GO zenginleştirme sonuçlarını deterministik, açıklanabilir
kurallarla yorumlayan bir "expert system" katmanıdır.

TASARIM NOTU:
Bu bir LLM DEĞİLDİR ve öyle sunulmamalıdır. Amacı; enrichment sonuçlarını
önceden tanımlı biyolojik bilgi tabanı + şablon varyasyonu + gen-bağlamı
kuralları ile insan-okunur bir rapora dönüştürmektir. Deterministik olması
(aynı girdi -> aynı çıktı) bir zayıflık değil, TÜBİTAK/bilimsel rapor
bağlamında bir AVANTAJDIR: hallüsinasyon riski yoktur, her cümle
kaynak gösterilebilir bilgiye dayanır, ve çıktı denetlenebilir/testedilebilir.

Rapor veya sunumlarda bu modülden bahsederken önerilen ifade:
"Kural tabanlı, gen-bağlamı duyarlı bir biyolojik yorumlama katmanı"
"Expert-system tarzı deterministik yorumlama motoru"
KULLANILMAMASI GEREKEN ifade: "LLM", "yapay zeka destekli yorumlama"
(gerçek bir dil modeli çağrılmıyorsa).
"""

import hashlib
from collections import Counter
import pandas as pd

# ============================================================================
# BİLGİ TABANI — Yolakların biyolojik anlamları ve Sophiark bağlamındaki yorumları
# ============================================================================

KNOWLEDGE_BASE = {
    # ---- CFTR / İyon taşınımı ekseni ----
    "Pancreatic secretion": {
        "description": "Pankreasın sindirim enzimleri ve bikarbonat salgılamasını düzenleyen yolaktır.",
        "impact": "Hedef genin susturulması, pankreas salgı fonksiyonlarını etkileyebilir. Bu, malabsorbsiyon (emilim bozukluğu) ve sindirim sorunları riskini artırabilir.",
        "recommendation": "Pankreas enzim replasman tedavisi veya alternatif salgı yollarını hedefleyen kombinasyon terapileri değerlendirilebilir.",
    },
    "Protein digestion and absorption": {
        "description": "Proteinlerin amino asitlere parçalanması ve bağırsaktan emilimini kapsar.",
        "impact": "Protein sindiriminde aksamaya işaret eder. Uzun süreli tedavilerde beslenme desteği gerekebilir.",
        "recommendation": "Proteaz inhibitörleriyle kombinasyon veya diyet takviyeleri düşünülebilir.",
    },
    "Chloride transport (GO:0006821)": {
        "description": "Klorür iyonlarının hücre zarından taşınmasını sağlayan biyolojik süreçtir.",
        "impact": "Klorür taşınımındaki bozulma, epitel dokularda iyon dengesizliğine yol açabilir. Bu, kistik fibrozis ve benzeri kanalopatilerin temel mekanizmasıdır.",
        "recommendation": "Alternatif klorür kanalları (ANO1, BEST2, SLC26 ailesi) hedeflenerek iyon dengesi yeniden sağlanabilir.",
    },
    "Bile secretion": {
        "description": "Karaciğerin safra asitlerini üretmesi ve safra kesesinden salgılamasını kapsar.",
        "impact": "Safra salgısındaki değişiklik, yağda çözünen vitaminlerin emilimini etkileyebilir ve hepatobiliyer yan etkilere yol açabilir.",
        "recommendation": "Karaciğer fonksiyon testleriyle takip önerilir. Ursodeoksikolik asit gibi safra asidi düzenleyicileri değerlendirilebilir.",
    },
    "Ion transport (GO:0006811)": {
        "description": "İyonların hücre zarından geçişini sağlayan temel biyolojik süreçtir.",
        "impact": "Geniş çaplı iyon transportu etkilenimi, hücrenin elektriksel dengesini ve hacim regülasyonunu bozabilir. Bu, özellikle sinir ve kas dokularında kritiktir.",
        "recommendation": "Etkilenen spesifik iyon kanallarını belirlemek için hedefli analiz önerilir.",
    },
    "Inorganic anion transport (GO:0015698)": {
        "description": "Klorür, bikarbonat, sülfat gibi inorganik anyonların taşınmasını kapsar.",
        "impact": "Anyon taşınımındaki yaygın etkilenim, hücre içi pH dengesini bozabilir ve metabolik asidoz riskini artırabilir.",
        "recommendation": "pH düzenleyici mekanizmaların (karbonik anhidrazlar, SLC4 ailesi) kompanzasyon kapasitesi değerlendirilmelidir.",
    },
    "Regulation of cellular pH (GO:0030641)": {
        "description": "Hücre içi asit-baz dengesinin korunmasını sağlayan düzenleyici mekanizmalardır.",
        "impact": "pH regülasyonundaki bozulma, enzim aktivitelerini ve protein fonksiyonlarını olumsuz etkileyebilir. Hücre stres yanıtını tetikleyebilir.",
        "recommendation": "Karbonik anhidraz inhibitörleri veya SLC9 ailesi modülatörleriyle pH homeostazı desteklenebilir.",
    },
    "Water transport (GO:0006833)": {
        "description": "Suyun hücre zarından aquaporinler aracılığıyla taşınmasını sağlar.",
        "impact": "Su transportunun etkilenmesi, epitel dokularda sıvı dengesizliğine yol açabilir. Dehidratasyon veya ödem riski oluşturabilir.",
        "recommendation": "Aquaporin ekspresyon profili izlenmeli; gerekirse sıvı replasman tedavisi planlanmalıdır.",
    },
    "Collecting duct acid secretion": {
        "description": "Böbrek toplayıcı kanallarında asit sekresyonunu düzenler.",
        "impact": "Böbrek fonksiyonlarında asit-baz dengesi bozulabilir. Renal tübüler asidoz riski değerlendirilmelidir.",
        "recommendation": "Böbrek fonksiyon testleri ve idrar pH takibi önerilir.",
    },
    "Nitrogen metabolism": {
        "description": "Azot içeren bileşiklerin metabolizmasını kapsar.",
        "impact": "Azot metabolizmasındaki değişiklik, amino asit dengesini ve üre döngüsünü etkileyebilir.",
        "recommendation": "Kan üre azotu (BUN) ve amino asit profili takip edilmelidir.",
    },

    # ---- GLP1R / inkretin-glisemi ekseni ----
    "Insulin secretion": {
        "description": "Pankreas beta hücrelerinden insülin salgılanmasını düzenleyen sinyal yolağıdır.",
        "impact": "İnsülin salgı mekanizmasının etkilenmesi, glisemik kontrolde bozulmaya ve postprandiyal glukoz dalgalanmalarına yol açabilir.",
        "recommendation": "Alternatif inkretin reseptörleri (GIPR) veya sülfonilüre sınıfı ajanlarla kombinasyon değerlendirilebilir.",
    },
    "GnRH secretion": {
        "description": "Hipotalamo-hipofiz eksenindeki nöroendokrin salgı düzenlenmesini kapsar.",
        "impact": "Bu eksenin etkilenmesi, dolaylı olarak metabolik ve üreme fonksiyonları ile ilişkili sinyal yollarını etkileyebilir.",
        "recommendation": "Endokrin panel ile eksen bütünlüğünün izlenmesi önerilir.",
    },
    "cAMP signaling pathway": {
        "description": "G-protein bağlı reseptörler üzerinden hücre içi ikincil habercilik yolağıdır (GLP1R dahil birçok reseptörün ortak efektörü).",
        "impact": "cAMP sinyalindeki bozulma, hedef genin bağlı olduğu reseptör ailesinin downstream etkilerini geniş çapta etkileyebilir.",
        "recommendation": "PKA/Epac alt yolaklarının ayrı ayrı değerlendirilmesi, etkinin hangi kolda yoğunlaştığını netleştirebilir.",
    },
    "Neuroactive ligand-receptor interaction": {
        "description": "Nörotransmitter ve peptit hormonların reseptörleriyle etkileşimini kapsayan geniş yolaktır.",
        "impact": "Bu yolağın etkilenmesi, hedef genin merkezi sinir sistemi veya periferik sinyal ağındaki geniş bağlantısına işaret eder.",
        "recommendation": "Doku-spesifik ekspresyon verisiyle (Human Protein Atlas) etkinin lokalize olduğu doku belirlenmelidir.",
    },
    "Cushing syndrome": {
        "description": "Kortikotropin ve steroid hormon sinyalizasyonuyla ilişkili yolaktır (KEGG hastalık haritası).",
        "impact": "Bu yolağın enrichment'ta çıkması, hedef genin steroid/stres hormon eksenleriyle dolaylı etkileşimine işaret eder.",
        "recommendation": "Kortizol ve ilişkili hormon düzeylerinin izlenmesi önerilir.",
    },

    # ---- PSEN1 / Notch-amiloid ekseni ----
    "Notch signaling pathway": {
        "description": "Hücre farklılaşması, gelişim ve hücreler arası iletişimi düzenleyen evrimsel olarak korunmuş sinyal yolağıdır.",
        "impact": "Notch yolağının etkilenmesi, hücre kaderi kararlarında ve doku homeostazında bozulmaya yol açabilir. PSEN1 gibi genler için bu, gama-sekretaz kompleksi aktivitesinin merkezi rolünü yansıtır.",
        "recommendation": "Gama-sekretaz aktivitesinin dolaylı ölçümü (Notch intracellular domain düzeyi) ile doğrulama önerilir.",
    },
    "Alzheimer disease": {
        "description": "Amiloid-beta birikimi ve nörodejenerasyon ile ilişkili KEGG hastalık yolağıdır.",
        "impact": "Bu yolağın güçlü şekilde çıkması, hedef genin amiloidojenik işlem ağındaki merkezi konumunu doğrular; PSEN1 için beklenen bir bulgudur.",
        "recommendation": "APP işleme ürünlerinin (Aβ40/Aβ42 oranı) simülasyon dışı doğrulanması, bulgunun biyolojik geçerliliğini güçlendirir.",
    },
    "Amyloid-beta metabolic process (GO:0050435)": {
        "description": "Amiloid öncül proteininin (APP) işlenmesi ve amiloid-beta peptidlerinin oluşum/temizlenme sürecini kapsar.",
        "impact": "Bu sürecin etkilenmesi, hedef genin nöropatolojik amiloid birikimi ile doğrudan mekanistik bağlantısını gösterir.",
        "recommendation": "Simülasyon sonucu, PSEN1 vaka çalışmalarındaki bilinen literatür bulgularıyla karşılaştırılarak raporlanmalıdır.",
    },
    "Wnt signaling pathway": {
        "description": "Hücre proliferasyonu, polarite ve kader belirlenmesinde rol oynayan, Notch ile çapraz etkileşimli sinyal ağıdır.",
        "impact": "Wnt yolağının etkilenmesi, PSEN1 ile ilişkili gama-sekretaz kompleksinin Notch dışı substratları üzerindeki dolaylı etkisini yansıtabilir.",
        "recommendation": "Çapraz yolak etkileşiminin (Notch-Wnt crosstalk) ayrı bir analizle doğrulanması önerilir.",
    },

    # ---- Genel terimler için fallback şablonları ----
    "_default_transport": {
        "description": "Moleküllerin hücre zarından taşınmasını sağlayan bir süreçtir.",
        "impact": "Bu transport sürecinin etkilenmesi, ilgili moleküllerin hücresel dengesini bozabilir.",
        "recommendation": "Etkilenen spesifik taşıyıcı proteinlerin belirlenmesi ve hedeflenmesi önerilir.",
    },
    "_default_signaling": {
        "description": "Hücresel sinyal iletimini düzenleyen bir yolaktır.",
        "impact": "Sinyal yolağının etkilenmesi, hücresel yanıt mekanizmalarını değiştirebilir.",
        "recommendation": "Alternatif sinyal yolaklarının kompanzasyon kapasitesi değerlendirilmelidir.",
    },
    "_default_metabolism": {
        "description": "Metabolik bir süreci düzenler.",
        "impact": "Metabolik dengenin bozulması, enerji homeostazını etkileyebilir.",
        "recommendation": "Metabolik panel testleriyle takip önerilir.",
    },
    "_default_development": {
        "description": "Hücre veya doku gelişimi/farklılaşmasıyla ilişkili bir süreçtir.",
        "impact": "Gelişimsel programın etkilenmesi, doku homeostazı ve hücre kaderi kararlarında sapmaya yol açabilir.",
        "recommendation": "Etkilenen dokunun gelişimsel bağlamda ayrı değerlendirilmesi önerilir.",
    },
    "_default_generic": {
        "description": "Hücresel fonksiyonla ilişkili bir biyolojik süreçtir.",
        "impact": "Bu sürecin etkilenmesi, ilişkili hücresel mekanizmalarda değişikliğe yol açabilir.",
        "recommendation": "Yolağın hedef genle mekanistik ilişkisinin literatür taramasıyla doğrulanması önerilir.",
    },
}

# Anahtar kelime -> fallback kategori ağırlıklı eşleştirme (kaba `in` yerine)
_FALLBACK_KEYWORDS = {
    "_default_transport": ["transport", "channel", "pump", "efflux", "influx", "secretion"],
    "_default_signaling": ["signal", "pathway", "cascade", "receptor", "transduction"],
    "_default_metabolism": ["metabol", "synthesis", "degradation", "catabol", "anabol"],
    "_default_development": ["development", "differentiation", "morphogenesis", "proliferation"],
}

# ============================================================================
# DETERMİNİSTİK ÇEŞİTLİLİK — aynı girdi hep aynı çıktıyı üretir, ama farklı
# girdiler (gen adı, term adı, p-değeri) farklı cümle varyantları seçer.
# Bu, sabit kalıp doldurma izlenimini azaltır; gerçek rastgelelik değildir
# (rapor tekrarlanabilirliği korunur).
# ============================================================================

def _seeded_choice(options: list, *seed_parts: str) -> str:
    """Girdi parçalarından deterministik bir hash türetip options içinden seçer."""
    key = "|".join(str(s) for s in seed_parts)
    h = int(hashlib.md5(key.encode("utf-8")).hexdigest(), 16)
    return options[h % len(options)]


def _get_significance_phrase(p_adj: float, seed: str = "") -> str:
    if p_adj < 1e-10:
        options = [
            "son derece güçlü ve istatistiksel olarak çok yüksek anlamlılıkla",
            "olağanüstü düzeyde yüksek istatistiksel anlamlılıkla",
            "ihmal edilebilir yanlış-pozitif riskiyle, çok güçlü şekilde",
        ]
    elif p_adj < 1e-5:
        options = [
            "yüksek istatistiksel anlamlılıkla",
            "güçlü bir istatistiksel sinyalle",
            "düşük yanlış-keşif oranıyla, belirgin şekilde",
        ]
    elif p_adj < 0.01:
        options = [
            "anlamlı şekilde",
            "istatistiksel olarak tutarlı biçimde",
            "kabul edilebilir güvenilirlikte",
        ]
    elif p_adj < 0.05:
        options = [
            "sınırda anlamlılıkla",
            "eşiğe yakın, ihtiyatla yorumlanması gereken bir anlamlılıkla",
            "zayıf ama istatistiksel eşiği geçen bir sinyalle",
        ]
    else:
        options = [
            "düşük anlamlılıkla (trend düzeyinde)",
            "istatistiksel eşiğin altında, sadece eğilim olarak",
        ]
    return _seeded_choice(options, seed, f"{p_adj:.3e}")


def _get_impact_level(overlap_ratio: float, seed: str = "") -> str:
    if overlap_ratio > 0.3:
        options = ["ciddi şekilde", "yaygın ve derin biçimde", "ağı geniş çapta etkileyerek"]
    elif overlap_ratio > 0.15:
        options = ["belirgin biçimde", "gözle görülür ölçüde", "orta-yüksek şiddette"]
    elif overlap_ratio > 0.05:
        options = ["kısmen", "sınırlı ama tespit edilebilir biçimde", "ölçülü düzeyde"]
    else:
        options = ["hafif düzeyde", "marjinal ölçüde", "düşük şiddette"]
    return _seeded_choice(options, seed, f"{overlap_ratio:.3f}")


# ============================================================================
# BAĞLAM YÖNETİCİSİ — Hedef gen türüne göre yorum
# ============================================================================

def _get_target_context(hinterland_score: float, category: str, hedef_gen: str = "") -> str:
    if category and "İmparator" in category:
        opts = [
            (
                f"⚠️ Hedef gen, ağın en kritik düğümlerinden biridir (İmparator Hub, skor: {hinterland_score:.1f}). "
                "Susturulması zincirleme etki yaratmış ve birçok hayati yolakta bozulmaya yol açmıştır. "
                "Bu genin doğrudan hedeflenmesi yüksek yan etki riski taşır; kombinasyon terapileri düşünülmelidir."
            ),
            (
                f"⚠️ {hedef_gen or 'Hedef gen'}, ağ topolojisinde İmparator Hub sınıfındadır (skor: {hinterland_score:.1f}). "
                "Bu düzeydeki merkezilik, tek gen inhibisyonunun ağ genelinde öngörülemeyen (off-target benzeri) "
                "etkiler doğurma riskini artırır."
            ),
        ]
        return _seeded_choice(opts, hedef_gen, category)
    elif category and "Stratejik" in category:
        opts = [
            (
                f"Hedef gen, ağın önemli bir dağıtıcısıdır (Stratejik Dağıtıcı, skor: {hinterland_score:.1f}). "
                "Susturulması, bağlı olduğu yolaklarda belirgin etkilere neden olmuştur."
            ),
            (
                f"{hedef_gen or 'Hedef gen'} ağda Stratejik Dağıtıcı rolündedir (skor: {hinterland_score:.1f}); "
                "İmparator Hub'lara göre daha lokalize ama yine de çok-yolaklı bir etki profili beklenir."
            ),
        ]
        return _seeded_choice(opts, hedef_gen, category)
    elif hinterland_score > 50:
        opts = [
            (
                f"Hedef gen, ağda aktif bir geçiş noktasıdır (skor: {hinterland_score:.1f}). "
                "Susturulması belirli yolaklarda etki yaratmış, ancak hücre kısmen kompanse edebilmiştir."
            ),
            (
                f"{hedef_gen or 'Hedef gen'} orta-yüksek merkezilik düzeyindedir (skor: {hinterland_score:.1f}). "
                "Ağ, susturmanın etkisini bir miktar alternatif yollarla dengeleyebilmiş görünüyor."
            ),
        ]
        return _seeded_choice(opts, hedef_gen, "mid")
    else:
        opts = [
            (
                f"Hedef gen, ağda daha lokal bir role sahiptir (skor: {hinterland_score:.1f}). "
                "Susturulması sınırlı sistemik etki yaratmıştır."
            ),
            (
                f"{hedef_gen or 'Hedef gen'} düşük merkezilik gösterir (skor: {hinterland_score:.1f}); "
                "bu, periferik konumunun beklenen bir sonucu olabilir."
            ),
        ]
        return _seeded_choice(opts, hedef_gen, "low")


def _match_knowledge_base(term: str):
    """Terimi bilgi tabanıyla eşleştirir; tam eşleşme yoksa ağırlıklı fallback seçer."""
    term_lower = term.lower()

    # 1) Doğrudan/kısmi eşleşme
    for key, value in KNOWLEDGE_BASE.items():
        if key.startswith("_default_"):
            continue
        if key.lower() in term_lower or term_lower in key.lower():
            return value, key

    # 2) Ağırlıklı anahtar kelime skorlama (kaba `in` yerine)
    scores = {cat: 0 for cat in _FALLBACK_KEYWORDS}
    for cat, keywords in _FALLBACK_KEYWORDS.items():
        for kw in keywords:
            if kw in term_lower:
                scores[cat] += 1

    best_cat = max(scores, key=scores.get)
    if scores[best_cat] > 0:
        return KNOWLEDGE_BASE[best_cat], best_cat

    return KNOWLEDGE_BASE["_default_generic"], "_default_generic"


# ============================================================================
# ANA YORUMLAMA FONKSİYONU
# ============================================================================

def interpret_pathways(
    enrichment_df: pd.DataFrame,
    hedef_gen: str = "",
    hinterland_score: float = 0,
    kategori: str = "",
    signal_loss: float = 0,
    max_pathways: int = 5,
) -> str:
    """
    KEGG/GO zenginleştirme sonuçlarını kural tabanlı olarak yorumlar ve
    Sophiark bağlamında biyolojik bir metin üretir.

    NOT: Bu fonksiyon bir dil modeli ÇAĞIRMAZ. Tüm cümleler önceden
    tanımlanmış bilgi tabanı + deterministik şablon seçiminden gelir.
    """
    if enrichment_df is None or enrichment_df.empty:
        return "Bu simülasyon için anlamlı bir yolak zenginleştirmesi bulunamadı."

    if "Adjusted P-value" in enrichment_df.columns:
        sig_df = enrichment_df[enrichment_df["Adjusted P-value"] < 0.05].head(max_pathways)
    else:
        sig_df = enrichment_df.head(max_pathways)

    if sig_df.empty:
        return "Hiçbir yolak istatistiksel anlamlılık eşiğini (p-adj < 0.05) geçemedi."

    parts = []
    hedef_str = hedef_gen if hedef_gen else "hedef gen"

    parts.append("## 🧬 Sophiark Biyolojik Yorum Raporu\n")
    parts.append(f"**Hedef:** {hedef_str}  \n")
    parts.append(f"**Sinyal Kaybı:** %{signal_loss:.4f}  \n\n")
    parts.append(
        "_Bu rapor kural tabanlı bir yorumlama motoru tarafından, önceden "
        "tanımlanmış biyolojik bilgi tabanı kullanılarak üretilmiştir. "
        "Bir dil modeli çağrısı içermez._\n\n"
    )

    if kategori:
        parts.append(_get_target_context(hinterland_score, kategori, hedef_gen) + "\n\n")

    parts.append(f"### 📊 Etkilenen Başlıca Yolaklar ({len(sig_df)} adet)\n\n")

    for i, (_, row) in enumerate(sig_df.iterrows(), 1):
        term = row.get("Term", row.get("term", "Bilinmeyen Yolak"))
        p_adj = row.get("Adjusted P-value", row.get("Adjusted P-value", 1.0))
        overlap = row.get("Overlap", row.get("overlap", "0/0"))

        try:
            if isinstance(overlap, str) and "/" in overlap:
                hit, total = overlap.split("/")
                overlap_ratio = int(hit) / int(total) if int(total) > 0 else 0
                overlap_str = f"{hit}/{total}"
            else:
                overlap_ratio = 0.05
                overlap_str = str(overlap)
        except (ValueError, ZeroDivisionError):
            overlap_ratio = 0.05
            overlap_str = str(overlap)

        kb_entry, matched_key = _match_knowledge_base(term)

        sig_phrase = _get_significance_phrase(p_adj, seed=f"{hedef_gen}:{term}")
        impact_level = _get_impact_level(overlap_ratio, seed=f"{hedef_gen}:{term}")

        parts.append(f"#### {i}. {term}\n")
        parts.append(f"**Anlamlılık:** {sig_phrase} (düzeltilmiş p-değeri: {p_adj:.2e})\n\n")
        parts.append(f"**Açıklama:** {kb_entry['description']}\n\n")
        parts.append(
            f"**Sophiark Yorumu:** Bu yolak, hedef genin susturulmasından {impact_level} "
            f"etkilenmiştir ({overlap_str} gen overlap). {kb_entry['impact']}\n\n"
        )
        parts.append(f"**Öneri:** {kb_entry['recommendation']}\n\n")
        if matched_key.startswith("_default_"):
            parts.append(
                "_Not: Bu yolak bilgi tabanında doğrudan tanımlı değildir; "
                "kategori bazlı genel yorum uygulanmıştır. Manuel literatür "
                "doğrulaması önerilir._\n\n"
            )
        parts.append("---\n\n")

    parts.append("### 🔬 Genel Değerlendirme\n\n")

    if signal_loss > 1.0:
        opts = [
            (
                f"⚠️ Sinyal kaybı %{signal_loss:.2f} gibi yüksek bir seviyede. "
                "Bu, hedef genin hücre için kritik olduğunu ve susturulmasının ciddi sistemik etkiler "
                "yaratacağını göstermektedir. Doğrudan hedefleme yerine, doz kontrollü inhibisyon veya "
                "kombinasyon terapileri değerlendirilmelidir.\n\n"
            ),
            (
                f"⚠️ %{signal_loss:.2f} düzeyindeki sinyal kaybı, ağ genelinde geniş çaplı bir bozulmaya "
                "işaret ediyor. Bu düzeyde bir etki, hedefin tek başına inhibe edilmesinin klinik olarak "
                "riskli olabileceğini düşündürür.\n\n"
            ),
        ]
        parts.append(_seeded_choice(opts, hedef_gen, "high_loss"))
    elif signal_loss > 0.1:
        opts = [
            (
                f"Sinyal kaybı %{signal_loss:.2f} düzeyindedir. Bu, hedef genin susturulmasının belirli "
                "yolaklarda etki yarattığını, ancak hücrenin kısmi kompanzasyon mekanizmaları devreye "
                "sokabildiğini göstermektedir.\n\n"
            ),
            (
                f"%{signal_loss:.2f} sinyal kaybı orta düzeydedir; etkilenen yolaklar tespit edilebilir "
                "seviyede ama ağın tamamını felç edecek boyutta değildir.\n\n"
            ),
        ]
        parts.append(_seeded_choice(opts, hedef_gen, "mid_loss"))
    else:
        opts = [
            (
                f"Sinyal kaybı %{signal_loss:.4f} gibi düşük bir seviyededir. Bu, hedef genin susturulmasının "
                "sınırlı etki yarattığını ve hücrenin alternatif yollarla bu kaybı büyük ölçüde telafi "
                "edebildiğini göstermektedir. Bu hedef tek başına etkili bir ilaç hedefi olmayabilir; "
                "kombinasyon terapileri düşünülmelidir.\n\n"
            ),
            (
                f"%{signal_loss:.4f} düzeyindeki düşük sinyal kaybı, ağın bu genin kaybına karşı güçlü "
                "bir dayanıklılık (robustness) sergilediğini gösteriyor.\n\n"
            ),
        ]
        parts.append(_seeded_choice(opts, hedef_gen, "low_loss"))

    parts.append(
        "---\n"
        "*Bu rapor, Sophiark Kural Tabanlı Biyolojik Yorumlama Motoru (KTBYM) tarafından "
        "deterministik kurallarla oluşturulmuştur; bir dil modeli çağrısı içermez. "
        "Klinik kararlar için lütfen ıslak laboratuvar validasyonu yapınız.*\n"
    )

    return "".join(parts)


# ============================================================================
# KISA ÖZET FONKSİYONU (Sidebar veya hızlı görünüm için)
# ============================================================================

def quick_summary(
    enrichment_df: pd.DataFrame,
    hedef_gen: str = "",
    signal_loss: float = 0,
) -> str:
    """Zenginleştirme sonuçlarının kısa, tek cümlelik özetini döndürür."""
    if enrichment_df is None or enrichment_df.empty:
        return "Anlamlı yolak bulunamadı."

    if "Adjusted P-value" in enrichment_df.columns:
        sig_df = enrichment_df[enrichment_df["Adjusted P-value"] < 0.05]
    else:
        sig_df = enrichment_df.head(3)

    if sig_df.empty:
        return "Hiçbir yolak anlamlılık eşiğini geçemedi."

    hedef_str = hedef_gen if hedef_gen else "Hedef gen"
    top_term = sig_df.iloc[0].get("Term", sig_df.iloc[0].get("term", "bir yolak"))
    count = len(sig_df)

    if signal_loss > 1.0:
        return (
            f"{hedef_str} susturulması ciddi sistemik etki yarattı (%{signal_loss:.2f} kayıp). "
            f"{count} yolak anlamlı şekilde etkilendi, en güçlü sinyal: {top_term}."
        )
    elif signal_loss > 0.1:
        return (
            f"{hedef_str} susturulması {count} yolakta anlamlı etki yarattı (%{signal_loss:.2f} kayıp). "
            f"Başlıca etkilenen: {top_term}."
        )
    else:
        return (
            f"{hedef_str} susturulması sınırlı etki yarattı (%{signal_loss:.4f} kayıp). "
            f"{count} yolak sınırda anlamlı, en güçlü: {top_term}."
        )
        


# ============================================================================
# EŞİK SABİTLERİ — literatür-motivasyonlu, ampirik olarak seçilmiş.
# Optimize edilmiş hiperparametreler DEĞİLDİR; TÜBİTAK/bilimsel rapor
# bağlamında bu şekilde etiketlenmelidir.
# ============================================================================

ERKEN_BOOST_ESIK_ORANI   = 0.5   # max boost'un %50'sine SF>=0.50'de ulaşılırsa "erken"
ERKEN_SF_ESIGI           = 0.50  # dozların "ilk yarısı" sınırı (SF >= bu değer)
TUTARLI_SUPHELI_ORANI    = 0.5   # dozların en az yarısında tekrar eden gen = tutarlı
IC50_HASSAS_ESIK         = 0.5   # IC50 > bu değer ise "hassas sistem"


# ============================================================================
# AŞAMA 1: Erken Kompanzasyon Tespiti
# ============================================================================

def _detect_early_compensation(dose_df: pd.DataFrame) -> dict:
    """
    Paralog boost eğrisinin ne kadar erken (düşük dozda) yükseldiğini tespit eder.

    Mantık: max boost değerinin %50'sine, dozların ilk yarısında (SF >= 0.50)
    ulaşılıyorsa bu 'erken kompanzasyon' — kötü işaret (hedef kolayca yedeklenir).
    """
    df = dose_df.sort_values("Survival_Fraction", ascending=False).reset_index(drop=True)
    max_boost = df["Paralog_Boost_Factor"].max()

    if max_boost <= 0:
        return {"durum": "YOK", "esik_sf": None, "max_boost": 0.0}

    esik_deger = max_boost * ERKEN_BOOST_ESIK_ORANI

    # SF azalan sırada (0.90 -> 0.01) ilk eşiği aşan noktayı bul
    esigi_asan = df[df["Paralog_Boost_Factor"] >= esik_deger]
    if esigi_asan.empty:
        return {"durum": "YOK", "esik_sf": None, "max_boost": round(max_boost, 6)}

    ilk_asim_sf = esigi_asan["Survival_Fraction"].max()  # en yüksek SF'de (en erken) aşan

    if ilk_asim_sf >= ERKEN_SF_ESIGI:
        durum = "ERKEN"
    else:
        durum = "GEC"

    return {
        "durum": durum,
        "esik_sf": round(float(ilk_asim_sf), 4),
        "max_boost": round(float(max_boost), 6),
    }


# ============================================================================
# AŞAMA 2: Tutarlı Şüpheli Gen Tespiti
# ============================================================================

def _detect_consistent_suspects(dose_df: pd.DataFrame) -> list:
    """
    Congestion_Top_Stressed sütununda, dozların en az yarısında tekrar eden
    genleri bulur. Frekansa göre azalan sırada döner.

    NOT: Şu anki veri sadece isim taşıyor (ΔBC büyüklüğü değil). İleride
    run_pharmacological_dose_response içinde gen-adı + ΔBC çiftleri
    saklanırsa, burada ağırlıklı bir skor (Σ ΔBC / doz_sayısı) kullanılabilir.
    """
    toplam_doz = len(dose_df)
    if toplam_doz == 0:
        return []

    sayac: Counter = Counter()
    for stresli_str in dose_df["Congestion_Top_Stressed"].fillna(""):
        if not stresli_str.strip():
            continue
        genler = [g.strip() for g in stresli_str.split(",") if g.strip()]
        sayac.update(genler)

    esik_sayi = toplam_doz * TUTARLI_SUPHELI_ORANI

    tutarlilar = [
        {"gene": gen, "frekans": sayi, "oran": round(sayi / toplam_doz, 3)}
        for gen, sayi in sayac.items()
        if sayi >= esik_sayi
    ]
    tutarlilar.sort(key=lambda x: x["frekans"], reverse=True)
    return tutarlilar


# ============================================================================
# AŞAMA 3: Hedef Yetersizliği Tespiti
# ============================================================================

def _detect_target_insufficiency(dose_df: pd.DataFrame) -> dict:
    """IC50_Survival_Frac değerini yorumlar."""
    if "IC50_Survival_Frac" not in dose_df.columns or dose_df.empty:
        return {"durum": "BILINMIYOR", "ic50": None}

    ic50_raw = dose_df["IC50_Survival_Frac"].iloc[0]

    if isinstance(ic50_raw, str) and ic50_raw == "Resistant":
        return {"durum": "YETERSIZ", "ic50": "Resistant"}

    try:
        ic50_val = float(ic50_raw)
    except (TypeError, ValueError):
        return {"durum": "BILINMIYOR", "ic50": None}

    if ic50_val > IC50_HASSAS_ESIK:
        durum = "HASSAS"
    else:
        durum = "ORTA"

    return {"durum": durum, "ic50": round(ic50_val, 4)}


# ============================================================================
# AŞAMA 4: Birleştirme — Kural Tabanlı Puanlama
# ============================================================================

def _compute_priority_level(erken: dict, supheli: list, yetersizlik: dict) -> tuple:
    """Üç göstergeyi tek bir öncelik seviyesine indirger (ağırlıksız, kural tabanlı)."""
    puan = 0
    gerekceler = []

    if erken["durum"] == "ERKEN":
        puan += 1
        gerekceler.append("erken_kompanzasyon")

    if len(supheli) >= 2:
        puan += 1
        gerekceler.append("coklu_tutarli_supheli")

    if yetersizlik["durum"] == "YETERSIZ":
        puan += 1
        gerekceler.append("hedef_yetersizligi")

    if puan >= 2:
        seviye = "YUKSEK_ONCELIK"
    elif puan == 1:
        seviye = "ORTA_ONCELIK"
    else:
        seviye = "DUSUK_ONCELIK"

    return seviye, puan, gerekceler


# ============================================================================
# METİN ÜRETİMİ — deterministik, seeded_choice ile varyasyonlu
# (bio_interpreter.py'deki interpret_pathways ile aynı üslup)
# ============================================================================

def _generate_priority_text(
    seviye: str,
    hedef_gen: str,
    erken: dict,
    supheli: list,
    yetersizlik: dict,
) -> str:
    hedef_str = hedef_gen if hedef_gen else "hedef gen"
    parts = []

    parts.append("## 🧬 Sophiark Direnç Deneyi Önceliklendirme Raporu\n\n")
    parts.append(f"**Hedef:** {hedef_str}  \n")
    parts.append(f"**Öncelik Seviyesi:** {seviye.replace('_', ' ')}  \n\n")

    parts.append(
        "_Bu rapor, doz-yanıt simülasyonunun kural tabanlı yorumudur. "
        "Bir dil modeli çağrısı içermez ve gerçek deneysel direnç verisiyle "
        "DOĞRULANMAMIŞTIR._\n\n"
    )

    parts.append("### 📊 Gerekçe\n\n")

    # Gösterge 1 metni
    if erken["durum"] == "ERKEN":
        opts = [
            f"Paralog/kompanzasyon aktivitesi, hedef genin bağlantıları henüz hafif "
            f"zayıflatılmışken (kalan aktivite ≈ {erken['esik_sf']}) belirgin şekilde "
            f"artıyor. Bu, ağ modelinde erken bir yedekleme paterni olduğunu gösteriyor.",
            f"Kalan aktivite fraksiyonu {erken['esik_sf']} civarındayken paralog boost "
            f"zaten yükseliyor — sistemin erken safhada telafi mekanizmalarını devreye "
            f"soktuğuna işaret eden bir topolojik bulgu.",
        ]
        parts.append("- " + _seeded_choice(opts, hedef_gen, "erken") + "\n\n")
    elif erken["durum"] == "GEC":
        opts = [
            "Paralog/kompanzasyon aktivitesi yalnızca yüksek inhibisyon seviyelerinde "
            "belirginleşiyor; erken safhada anlamlı bir yedekleme sinyali gözlenmiyor.",
            "Kompanzasyon mekanizması geç safhada devreye giriyor — bu, düşük/orta "
            "dozlarda hedefin nispeten yedeksiz kaldığına işaret eder.",
        ]
        parts.append("- " + _seeded_choice(opts, hedef_gen, "gec") + "\n\n")
    else:
        parts.append("- Paralog boost sinyali gözlenmedi (bu hedefin BLAST eşiğini aşan güçlü paralogu olmayabilir).\n\n")

    # Gösterge 2 metni
    if supheli:
        gen_listesi = ", ".join(f"**{g['gene']}** ({g['oran']*100:.0f}% dozlarda)" for g in supheli[:5])
        opts = [
            f"Doz serisi boyunca tutarlı şekilde en çok yük altına giren genler: {gen_listesi}. "
            f"Bu genler, olası kompanzasyon yolunun topolojik adaylarıdır.",
            f"Tekrarlanan sıkışıklık paterni şu genlerde gözlendi: {gen_listesi}. "
            f"Deneysel direnç taramasında öncelikli incelenmesi önerilir.",
        ]
        parts.append("- " + _seeded_choice(opts, hedef_gen, "supheli") + "\n\n")
    else:
        parts.append("- Dozlar arasında tutarlı şekilde tekrar eden bir sıkışıklık geni tespit edilmedi.\n\n")

    # Gösterge 3 metni
    if yetersizlik["durum"] == "YETERSIZ":
        opts = [
            "En yüksek inhibisyon seviyesinde bile sinyal kaybı %50 eşiğine ulaşamadı "
            "(IC50: Resistant). Bu, tek başına bu hedefin susturulmasının yetersiz "
            "kalabileceğine dair güçlü bir topolojik sinyaldir.",
            "Model, hiçbir doz seviyesinde hedefi yarı-maksimum düzeyde bastıramadı "
            "(IC50: Resistant) — bu durum ek hedef arayışını gerekçelendiriyor.",
        ]
        parts.append("- " + _seeded_choice(opts, hedef_gen, "yetersiz") + "\n\n")
    elif yetersizlik["durum"] == "HASSAS":
        parts.append(f"- IC50 (kalan aktivite noktası) ≈ {yetersizlik['ic50']} — sistem görece hassas, düşük inhibisyon yeterli olabilir.\n\n")
    elif yetersizlik["durum"] == "ORTA":
        parts.append(f"- IC50 (kalan aktivite noktası) ≈ {yetersizlik['ic50']} — yarı-maksimum etki için orta-yüksek inhibisyon gerekiyor.\n\n")

    parts.append("### 🔬 Önerilen Sonraki Adım\n\n")
    if seviye == "YUKSEK_ONCELIK":
        oneri_opts = [
            "Bu hedef için deneysel direnç taraması (hücre hattı knockout/knockdown + "
            "qPCR veya RNA-Seq ile yukarıda listelenen genlerin ifade takibi) öncelikli "
            "olarak planlanmalıdır.",
            "Yukarıda işaretlenen genlerin knockdown/knockout deneyleriyle doğrulanması, "
            "olası kombinasyon terapisi hedeflerini erken safhada belirlemeye yardımcı olabilir.",
        ]
        parts.append(_seeded_choice(oneri_opts, hedef_gen, "oneri_yuksek") + "\n\n")
    elif seviye == "ORTA_ONCELIK":
        parts.append(
            "Bazı kompanzasyon sinyalleri mevcut ancak tek bir gösterge baskın değil. "
            "İkincil öncelikli bir deneysel doğrulama adayı olarak değerlendirilebilir.\n\n"
        )
    else:
        parts.append(
            "Mevcut simülasyon verisinde belirgin bir erken/tutarlı kompanzasyon paterni "
            "gözlenmedi. Bu, direnç riskinin düşük olduğu anlamına gelmez — yalnızca bu "
            "modelin şu anki göstergelerinde öne çıkan bir sinyal bulunmadığını gösterir.\n\n"
        )

    parts.append(
        "---\n"
        "*⚠️ ÖNEMLİ SINIRLAMA: Bu sıralama, ağ topolojisi ve Hill/paralog matematik "
        "modeline dayanır; gerçek biyokimyasal kinetik veya farmakolojik direnç "
        "verisiyle DOĞRULANMAMIŞTIR. Bir tahmin değil, deneysel önceliklendirme "
        "hipotezidir. Klinik/deneysel kararlar için ıslak laboratuvar doğrulaması "
        "gereklidir.*\n"
    )

    return "".join(parts)


# ============================================================================
# ANA GİRİŞ NOKTASI
# ============================================================================

def interpret_resistance_priority(
    dose_df: pd.DataFrame,
    hedef_gen: str = "",
) -> dict:
    """
    Doz-yanıt simülasyonu çıktısını kural tabanlı olarak yorumlar ve
    direnç deneyi önceliklendirme raporu üretir.

    NOT: Bu fonksiyon bir dil modeli ÇAĞIRMAZ ve bir tahmin motoru DEĞİLDİR.
    Çıktısı, deneysel doğrulama için bir başlangıç noktası önerisidir.

    Parametreler
    ------------
    dose_df : run_pharmacological_dose_response() çıktısı
    hedef_gen : hedef genin adı/ENSP kimliği

    Döner
    -----
    dict:
        {
          "seviye": "YUKSEK_ONCELIK" | "ORTA_ONCELIK" | "DUSUK_ONCELIK",
          "puan": int,
          "gerekceler": list[str],
          "gosterge_1_erken_kompanzasyon": dict,
          "gosterge_2_tutarli_supheliler": list[dict],
          "gosterge_3_yetersizlik": dict,
          "metin": str,  # insan-okunur, Markdown formatlı rapor
        }
    """
    if dose_df is None or dose_df.empty:
        return {
            "seviye": "BILINMIYOR",
            "puan": 0,
            "gerekceler": [],
            "gosterge_1_erken_kompanzasyon": {},
            "gosterge_2_tutarli_supheliler": [],
            "gosterge_3_yetersizlik": {},
            "metin": "Doz-yanıt verisi bulunamadı; direnç önceliklendirmesi yapılamıyor.",
        }

    gerekli_sutunlar = {"Survival_Fraction", "Paralog_Boost_Factor",
                         "Congestion_Top_Stressed", "IC50_Survival_Frac"}
    eksik = gerekli_sutunlar - set(dose_df.columns)
    if eksik:
        return {
            "seviye": "BILINMIYOR",
            "puan": 0,
            "gerekceler": [],
            "gosterge_1_erken_kompanzasyon": {},
            "gosterge_2_tutarli_supheliler": [],
            "gosterge_3_yetersizlik": {},
            "metin": f"Eksik sütun(lar) nedeniyle rapor üretilemedi: {eksik}",
        }

    erken       = _detect_early_compensation(dose_df)
    supheliler  = _detect_consistent_suspects(dose_df)
    yetersizlik = _detect_target_insufficiency(dose_df)

    seviye, puan, gerekceler = _compute_priority_level(erken, supheliler, yetersizlik)

    metin = _generate_priority_text(seviye, hedef_gen, erken, supheliler, yetersizlik)

    return {
        "seviye": seviye,
        "puan": puan,
        "gerekceler": gerekceler,
        "gosterge_1_erken_kompanzasyon": erken,
        "gosterge_2_tutarli_supheliler": supheliler,
        "gosterge_3_yetersizlik": yetersizlik,
        "metin": metin,
    }


# ============================================================================
# HIZLI ÖZET (sidebar / kısa görünüm için — quick_summary ile aynı örüntü)
# ============================================================================

def quick_resistance_summary(dose_df: pd.DataFrame, hedef_gen: str = "") -> str:
    """Direnç önceliklendirmesinin tek cümlelik özetini döndürür."""
    sonuc = interpret_resistance_priority(dose_df, hedef_gen)
    hedef_str = hedef_gen if hedef_gen else "Hedef gen"
    seviye = sonuc["seviye"]

    if seviye == "BILINMIYOR":
        return f"{hedef_str} için direnç önceliklendirmesi yapılamadı (veri eksik)."

    n_supheli = len(sonuc["gosterge_2_tutarli_supheliler"])
    seviye_okunabilir = seviye.replace("_", " ").title()

    if n_supheli > 0:
        ilk_gen = sonuc["gosterge_2_tutarli_supheliler"][0]["gene"]
        return (
            f"{hedef_str}: {seviye_okunabilir} — {n_supheli} tutarlı şüpheli gen "
            f"tespit edildi (en belirgini: {ilk_gen}). Deneysel doğrulama önerilir."
        )
    return f"{hedef_str}: {seviye_okunabilir} — belirgin tutarlı şüpheli gen tespit edilmedi."
