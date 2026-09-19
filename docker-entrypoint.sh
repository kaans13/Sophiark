#!/bin/sh
set -eu

mkdir -p /app/outputs/reports /app/outputs/cache /app/outputs/history
mkdir -p /app/outputs_mouse/reports /app/outputs_mouse/cache

if [ ! -f /app/outputs/reports/genom_analiz_sonuc_v5.csv ]; then
  cp /app/runtime_seed/outputs/reports/genom_analiz_sonuc_v5.csv /app/outputs/reports/
fi

if [ ! -f /app/outputs_mouse/reports/mouse_genom_analiz_sonuc_v5.csv ]; then
  cp /app/runtime_seed/outputs_mouse/reports/mouse_genom_analiz_sonuc_v5.csv /app/outputs_mouse/reports/
fi

if [ ! -f /app/outputs/cache/ensp_symbol_cache.pkl ]; then
  cp /app/runtime_seed/outputs/cache/ensp_symbol_cache.pkl /app/outputs/cache/
fi

if [ ! -f /app/outputs/cache/kegg_cache.pkl ]; then
  cp /app/runtime_seed/outputs/cache/kegg_cache.pkl /app/outputs/cache/
fi

exec streamlit run app.py \
  --server.address=0.0.0.0 \
  --server.port="${PORT:-8501}" \
  --server.headless=true
