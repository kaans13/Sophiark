# Sophiark — Perturbation Propagation Trace Manuel Kontrol Listesi

Bu ekran **predicted network propagation** gÃ¶sterir; sinyal iletimi veya
nedensel biyolojik mekanizma kanÄ±tÄ± deÄŸildir.

## Test 1 — Human CFTR

1. Ä°nsan ve uygun bir doku seÃ§in; CFTR ile perturbasyon Ã§alÄ±ÅŸtÄ±rÄ±n.
2. **Etkinin YayÄ±lÄ±mÄ±nÄ± Ä°ncele** dÃ¼ÄŸmesine basÄ±n.
3. Layer 1, Layer 1-2 ve Layer 1-3+ filtrelerini sÄ±rayla deneyin.
4. Significant only filtresinde node sayÄ±sÄ±nÄ±n artmadÄ±ÄŸÄ±nÄ± ve node tablosunda q-value alanÄ±nÄ± kontrol edin.
5. Directed edge varsa oka; yoksa structural-only kenarlara bakÄ±n.

Beklenen: Hedef Layer 0'da, doÄŸrudan komÅŸular Layer 1'de gÃ¶rÃ¼lÃ¼r. Ok olmayan kenarlar yÃ¶nlÃ¼ gibi sunulmaz.

BaÅŸarÄ±sÄ±zlÄ±k: target kayÄ±p, layer sayÄ±larÄ± tutarsÄ±z veya structural edge'lerde yÃ¶n oku.

## Test 2 — Human TP53

1. TP53 ile aynÄ± akÄ±ÅŸÄ± Ã§alÄ±ÅŸtÄ±rÄ±n.
2. Hub aÄŸÄ±nda node limitinin grafiÄŸi binlerce node ile doldurmadÄ±ÄŸÄ±nÄ± kontrol edin.
3. Routes tablosunda her rotanÄ±n **Predicted Network Propagation Route** olarak, nedensel mekanizma iddiasÄ± olmadan gÃ¶rÃ¼ndÃ¼ÄŸÃ¼nÃ¼ kontrol edin.

Beklenen: Bounded node listesi, en fazla ayarlanabilir route sayÄ±sÄ±, q-value/ΔPR alanlarÄ±.

BaÅŸarÄ±sÄ±zlÄ±k: “TP53 X'i aktive eder/inhibe eder” gibi evidence'sÄ±z dil veya sÄ±nÄ±rsÄ±z rota Ã¼retimi.

## Test 3 — Mouse Cftr

1. Fareyi seÃ§in ve Cftr perturbasyonu Ã§alÄ±ÅŸtÄ±rÄ±n.
2. Trace'i aÃ§Ä±n; node ID'lerinin ENSMUSP, sembollerin mouse sembolÃ¼ olduÄŸunu doÄŸrulayÄ±n.
3. Evidence ayrÄ±ntÄ±larÄ±nda mouse OmniPath snapshot bilgisini kontrol edin.

Beklenen: Ä°nsan evidence'Ä± sessizce kullanÄ±lmaz; metadata organism=10090 snapshotÄ±nÄ± gÃ¶sterir.

BaÅŸarÄ±sÄ±zlÄ±k: Human node IDs, boÅŸ snapshotta “no directed evidence” iddiasÄ± veya species karÄ±ÅŸÄ±klÄ±ÄŸÄ±.

## Test 4 — Evidence conflict

1. Evidence tablosunda `conflicting_evidence` kaydÄ± bulunan bir rota seÃ§in.
2. Grafikte kesikli `conflict` iÅŸaretini ve tabloda stimulation/inhibition alanlarÄ±nÄ± kontrol edin.

Beklenen: Conflict tek bir stimulation/inhibition sonucuna indirgenmez.

BaÅŸarÄ±sÄ±zlÄ±k: Ã‡eliÅŸkili kanÄ±tÄ±n tek yÃ¶n veya tek etki olarak sunulmasÄ±.

## Test 5 — No directed evidence

1. Directed evidence only filtresi kapalÄ±yken structural-only rota bulun.
2. AynÄ± filtreyi aÃ§Ä±n.

Beklenen: Structural rota ilk gÃ¶rÃ¼nÃ¼mde kalÄ±r; directed-only filtresinde kaybolabilir. Bu negatif biyolojik kanÄ±t deÄŸildir.

BaÅŸarÄ±sÄ±zlÄ±k: Evidence bulunmayan kenara activation/inhibition oku verilmesi.

## Test 6 — Severity

1. AynÄ± hedef/doku ile Mild ve Near-complete senaryolarÄ±nda ayrÄ± simÃ¼lasyon Ã§alÄ±ÅŸtÄ±rÄ±n.
2. Trace node/route tablolarÄ±nÄ± ve export metadata'daki attenuation severity'yi karÅŸÄ±laÅŸtÄ±rÄ±n.

Beklenen: Trace, her senaryonun Ã¶lÃ§Ã¼lmÃ¼ÅŸ sonucu Ã¼zerinden ayrÄ± Ã¼retilir; parametre metadata'da gÃ¶rÃ¼nÃ¼r.

BaÅŸarÄ±sÄ±zlÄ±k: Ä°ki simÃ¼lasyonun metadata'sÄ±nda aynÄ± attenuation yazmasÄ± veya eski trace'in yeni sonuÃ§ta kalmasÄ±.

## Test 7 — Reproducibility

1. AynÄ± target, tissue, severity ve random seed ile simÃ¼lasyonu iki kez Ã§alÄ±ÅŸtÄ±rÄ±n.
2. Export edilen `nodes.csv`, `routes.csv` ve `metadata.json` dosyalarÄ±nÄ± karÅŸÄ±laÅŸtÄ±rÄ±n.

Beklenen: AynÄ± graph/data snapshot ve seed ile node/route sÄ±ralamasÄ± aynÄ±dÄ±r; timestamp farklÄ± olabilir.

BaÅŸarÄ±sÄ±zlÄ±k: Parametreler aynÄ±yken aÃ§Ä±klanamayan route veya node farkÄ±.
