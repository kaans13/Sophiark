# Sophiark Manuel Uygulama Kabul Testi

Her testte beklenmeyen teknik hata, boş beyaz ekran, ham exception metni veya tür/doku bağlamıyla uyuşmayan sonuç **başarısızlık** sayılır.

## 1. Açılış ve navigasyon

1. Proje klasöründe `streamlit run app.py` komutunu çalıştırın.
2. Ana Sayfa, Analiz, Keşif, Karşılaştırma ve Ayarlar bağlantılarını sırayla açın.
3. Beklenen: Sayfalar 5–10 saniye içinde kullanılabilir hale gelir; başlıklar ve üst context bar görünür.

## 2. İnsan/fare ve doku context'i

1. Üstte **Değiştir** üzerinden İnsan → Fare geçişi yapın.
2. Fare doku listesini açın, bir doku seçip uygulayın; sonra İnsan'a dönün.
3. Beklenen: Context bar doğru tür/dokuyu yazar; eski türün grafı yeni türde kullanılmaz.

## 3. İnsan CFTR forward simulation

1. İnsan ve `Lung` seçin.
2. Analiz → Forward Simulation → Gen; CFTR'yi seçin.
3. Simülasyonu başlatın.
4. Beklenen: İlerleme durumu görünür; sonuçta Öne Çıkan Bulgu, dört özet metrik ve sonuç sekmeleri oluşur. Lung grafı genel ağdan küçük olmalıdır.

## 4. İnsan TP53 ve çoklu gen

1. TP53 ile tek hedef analizi çalıştırın.
2. Ardından TP53 + MDM2 gibi iki hedef seçip tekrar çalıştırın.
3. Beklenen: Hedefler raporda ayrı izlenebilir; bir genin seçilmesi diğerini sessizce düşürmez.

## 5. Fare forward simulation

1. Fare'yi seçin; Cftr veya bilinen başka bir fare genini çalıştırın.
2. Beklenen: Sonuçlarda gen sembolü birincil, `ENSMUSP...` kimliği ikincil alanda görünür; insan motoru/ENSP kimliği karışmaz.

## 6. Hastalık ve lokalizasyon hedefi

1. İnsan modunda hedef yöntemini Hastalık yapın ve `cystic fibrosis` arayın.
2. Bir hastalık seçip hedefleri getirin.
3. Ayrıca Hücresel lokalizasyon modunda `Mitochondrion` seçin.
4. Beklenen: API yoksa dostça uyarı görünür; uygulama çökmez. Lokalizasyon modu hedef gen zorunluluğu olmadan çalışabilir.

## 7. 2B ve 3B ağ

1. Tamamlanmış bir forward sonucu içinde Ağ sekmesini açın.
2. 2B ve 3B sekmelerini sırayla açın; düğümlerin üzerine gelin.
3. Beklenen: Etiketler kapalıyken haritada ENSP/ENSMUSP metni basılmaz. Hover'da alias ve ikincil protein kimliği görünür. 3B grafik döndürülebilir.

## 8. KEGG/GO, doku ifadesi ve lokalizasyon

1. Biyolojik Bağlam sekmesini açın.
2. Beklenen: KEGG/GO otomatik sonuç veya açık erişim bildirimi; hedef doku ifade grafiği ve varsa hücresel dağılım görünür.
3. Enrichment sonucu pathway aktivasyonu/inhibisyonu olarak yazılmamalıdır.

## 9. Doz-yanıt ve Threshold Sweep

1. İleri Analizler → Doz-Yanıt'ı çalıştırın.
2. İnsan modunda Threshold Sweep'i çalıştırın.
3. Beklenen: İlerleme durumu, tablo/grafik ve export oluşur. Fare modunda sweep'in desteklenmediği açıkça belirtilir.

## 10. Propagation Trace

1. İleri Analizler → Propagation Trace → **Etkinin yayılımını incele**.
2. Layer ve evidence filtrelerini deneyin.
3. Beklenen: Yapısal STRING kenarı yönsüz; yönlü kanıt okla; conflict kesikli olarak gösterilir. Sonuç biyolojik nedensellik olarak sunulmaz.

## 11. Target Stress Search

1. Keşif → Target Stress Search; TP53 seçin, aday limitini 5 bırakın.
2. Doğrudan, ağ aracılı ve birleşik modlardan en az ikisini deneyin.
3. Beklenen: Aday tablosu, ΔPageRank, BC proxy/full alanı, ağ mesafesi ve yönlü kanıt birlikte görünür; boş evidence negatif evidence sayılmaz.

## 12. Compensation Analysis

1. Keşif → Compensation Analysis · BETA; bir hedef seçin.
2. Beklenen: Sonuç tablosu oluşur ve computational prediction uyarısı kalır.

## 13. Tissue Differential

1. Analiz → Dokular Arası Karşılaştırma · BETA.
2. Bir hedef ve en az iki doku seçin.
3. Beklenen: Doku bazlı ağ boyutu, stressed/critical sayısı, efficiency ve Jaccard bileşenleri ayrı gösterilir; gizli tek skor yoktur.

## 14. Karşılaştırma ve export

1. Hem insan hem fare forward analizi sonrası Karşılaştırma sayfasını açın.
2. Rapor ve Export sekmesinden CSV ve Excel indirin.
3. Beklenen: Ortolog tablosu gerçek iki raporu kullanır; Excel dosyası açılır ve sembol/kimlik sütunları ayrıdır.

## 15. State, hata ve yeniden üretilebilirlik

1. Analiz sonucundan Ayarlar'a gidip geri dönün.
2. Doku veya türü değiştirin.
3. Beklenen: Aynı bağlamda sonuç korunur; bağlam değişince eski mutable graf temizlenir. Kullanıcıya traceback/ham exception gösterilmez.
