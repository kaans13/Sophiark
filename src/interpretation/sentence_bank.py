"""Curated Turkish sentences for the deterministic interpretation layer."""

from __future__ import annotations


SENTENCES = {
    "positive_delta": "{symbol}, hedef perturbasyonu sonrasında {delta} düzeyinde pozitif göreli PageRank yeniden dağılımı gösterdi.",
    "negative_delta": "{symbol}, hedef perturbasyonu sonrasında {delta} düzeyinde göreli PageRank kaybı gösterdi.",
    "baseline_context": "Ağdaki Önemi ({hinterland}) ve Geçiş Merkeziliği ({bc}) başlangıç topolojik bağlamını tanımlar; perturbasyon yanıtının kendisi değildir.",
    "baseline_missing": "Başlangıç topoloji alanları bu aday için mevcut değil; yorum ΔPageRank ile sınırlıdır.",
    "annotation_context": "Mevcut anotasyonlar {themes} bağlamıyla eşleşir; bu eşleşme mekanizma veya biyolojik durum doğrulaması değildir.",
    "annotation_missing": "Bu aday için kullanılabilir işlevsel anotasyon bulunmadığından ek mekanistik yorum yapılmadı.",
    "low_evidence": "Mevcut yerel anotasyonlar bu aday için güçlü bir işlevsel yorum oluşturmak için yetersiz. Adayın önemi şu anda esas olarak perturbasyon sonrası topolojik değişimden kaynaklanmaktadır.",
    "family_context": "{clusters} içindeki tekrar, aynı aile/kompleks üyelerinin topolojik olarak birlikte görünmesine bağlı olabilir; üyeler bağımsız kanıt olarak toplanmamalıdır.",
    "system_localized": "Hedef-yerel grafik verimliliği değişimi global değişimden daha büyük olduğundan, yapısal etkinin hedef çevresinde yoğunlaşmasıyla uyumludur.",
    "system_neutral": "Global ve hedef-yerel grafik verimliliği değerleri yalnızca ağ yapısının değişimini tanımlar; hücresel veya metabolik verimlilik ölçümü değildir.",
    "validation_transcript_protein": "Önerilen doğrulama: bağımsız transkript ve protein ölçümüyle adayın biyolojik yanıtı değerlendirilmelidir.",
    "validation_localization": "Önerilen doğrulama: adayın ilgili doku/hücredeki lokalizasyonu ve hedefle mekânsal birlikteliği deneysel olarak incelenmelidir.",
    "validation_perturbation": "Önerilen doğrulama: kontrollü hedef perturbasyonu altında zamana bağlı fonksiyonel test uygulanmalıdır.",
    "summary_primary": "{target} perturbasyonu sonrasında pozitif topolojik yeniden dağılımın baskın ekseni {theme} çevresinde yoğunlaşmaktadır.",
    "summary_support": "Bu ekseni, {genes} gibi farklı aile veya alt grup bağlamlarından gelen başlıca adaylar temsil etmektedir.",
    "summary_secondary": "İkincil eksen olarak {theme}, {genes} adaylarıyla görünmektedir.",
    "summary_family_caveat": "{families} içinde birden çok üyenin bulunması nedeniyle, bu üyeler bağımsız biyolojik kanıtlar olarak toplanmamalıdır.",
    "summary_tissue": "Bu okuma, {tissue} için kullanılan ağ bağlamına aittir ve dokuya özgü deneysel doğrulamanın yerini tutmaz.",
    "validation_theme": "Önerilen doğrulama: {theme} ekseni için uygun hücresel fonksiyon testi, kontrollü hedef perturbasyonu ile birlikte değerlendirilmelidir.",
}

# These tokens must not occur as a positive biological assertion in composed
# output. Tests inspect the composer output against this list.
FORBIDDEN_ASSERTION_TOKENS = (
    "activated", "upregulated", "inhibited", "compensated", "rescued",
    "aktive oldu", "aktive edildi", "yukarı düzenlendi", "baskılandı",
)
