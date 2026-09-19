# HPA ve Classic CORUM legacy provenance

Bu audit yalnız mevcut artifact'ları okur ve yeniden üretimleri
`data/derived/staging/legacy-provenance` altına yazar. Production veri veya
scientific engine dosyası değiştirilmez.

## HPA

Runtime `data/raw/hinterland_core.db:tissue_expression` tablosunu okur.
`src/data_loader.py` seçilen tissue için pozitif kayıtları alır ve
`protein_id` bazında `MAX(expression_level)` uygular.

Tabloda iki fiziksel ingestion bloğu vardır:

1. Rowid 1–652.103: 13.366 ENSP, 63 tissue ve 651.125 farklı key içeren
   discrete legacy blok.
2. Rowid 652.104–1.433.543: 19.536 ENSP, 40 tissue ve 781.440 farklı key
   içeren TPM blok.

Git commit `ed096a5` içindeki `db_patcher.py`, ilk bloğun pipeline'ını açıklar:
HPA v22 `normal_tissue.tsv`, `High/Medium/Low/Not detected → 3/2/1/0.1`,
gene+tissue cell-type MAX, MyGene `symbol,alias → ensembl.protein`, eşleşmeyeni
düşürme ve tissue title-case. Committed source 1.194.479 satırdan mapping öncesi
655.411 gene+tissue satırı üretir. Production blok 652.103 satırdır. Historical
MyGene response'u saklanmadığı için kalan mapping seçimi raw kaynaktan bugün
deterministik yeniden üretilemez.

TPM zinciri semantik olarak tamdır: `rna_tissue_hpa.tsv` içindeki 806.480 kayıt
`hpa_pivot.pkl` ile exact eşittir. ENSG→ENSP mapping sonrası 781.440 kayıt,
production tablosunun ikinci bloğuyla key/value olarak 781.440/781.440 eşittir.

İki blokta 267.803 ortak key vardır. Runtime MAX sonucu 48.045 key'de legacy,
212.786 key'de TPM, 6.972 key'de eşitlik tarafından belirlenir. Legacy-only
383.322, TPM-only 513.637 key vardır. TPM bloğunda doğal olarak `0.1`, `1`, `2`
veya `3` değerine sahip 44.102 kayıt bulunduğundan yalnız sayıya bakarak source
attribution yapılmaz; lineage fiziksel row block ve source equality ile kurulur.

Karar: `LEGACY_PARTIALLY_REPRODUCIBLE` ve
`SCIENTIFIC_MIGRATION_REQUIRED`. Mevcut composite tablo frozen runtime artifact
olarak korunmalıdır. Yeni sürümler clean HPA builder ile staging'de üretilmeli;
geçiş ayrı scientific/numerical validation gerektirir.

## Classic CORUM

Classic runtime artifact `data/processed/complex_members.pkl` dosyasıdır.
`src/config.py` bunu `dict[symbol, list[complex_name]]` olarak yükler;
`src/biology_logic.py` yalnız hedef ve komşu complex-name setlerinin kesişip
kesişmediğine bakar.

Production pickle'ın gerçek dönüşümü mevcut `data/raw/coreComplexes.txt`
üzerinden semantik olarak tam yeniden üretildi:

1. Organism filtresi uygulanmamış; Human, Mouse, Rat ve diğer tüm kayıtlar
   alınmış.
2. `subunits(Gene name)` primary alanı semicolon ile ayrılmış.
3. Bütün symbol değerleri uppercase yapılmış.
4. Rat `Nephrin MAGI` satırındaki boş member, string conversion nedeniyle
   `NAN` anahtarına dönüşmüş.

Bu pipeline 4.241/4.241 production anahtarını ve her anahtarın complex-name
membership setini tam üretir. Human-only clean builder'a göre görülen 807
production-only anahtarın 799'u non-human uppercase üyelerden, 8'i Human symbol
case normalization'dan açıklanır. `NAN`, 799 non-human grubunun içindedir.
Bu 807 anahtarın 748'i current canonical local symbol listesinde, 762'si daha
geniş local mapping tablosunda bulunur; 45'i geniş mapping'de yoktur. Bu
karakterizasyon lineage'i değiştirmez: anahtarlar alias expansion ile değil,
raw non-human/case-normalized primary alanlardan gelmiştir.
Alias/synonym expansion gerekmez ve üretim scriptinin kendisi repository'de
bulunmamıştır; ancak source ve transformation semantiği içerikten tam ve
deterministik olarak rekonstrükte edilmiştir.

Staging CORUM artifact ile CFTR, TGFBR1 ve ANO2/Lung controlled injection
parity'sinde maksimum Classic fark `7.21e-11`, Directed fark `1.38e-10` oldu.

Karar: `LEGACY_REPRODUCIBLE`. Existing Classic için legacy builder semantiği
korunabilir. Human-only clean CORUM builder ayrı Unified/Evidence context yolu
olarak kalmalıdır; Classic'e geçirilmesi ayrı scientific migration validation
gerektirir.
