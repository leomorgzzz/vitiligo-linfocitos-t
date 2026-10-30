#!/bin/bash
# Descargar las matrices de Cell Ranger y metadata de GSE298871 desde GEO
# correr desde la carpeta inicial del repositorio: bash scripts/01_descargar.sh

set -euo pipefail

GEO="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298871"
RAW="data/raw"

mkdir -p "$RAW/matrices"

# matrices de las 17 muestras comprimidas y metadata de cada muestra
wget -c -P "$RAW" "$GEO/suppl/GSE298871_RAW.tar"
wget -c -P "$RAW" "$GEO/matrix/GSE298871_series_matrix.txt.gz"

# extraer sin descomprimir los .gz de adentro, el programa scanpy los lee de esta manera
tar -xf "$RAW/GSE298871_RAW.tar" -C "$RAW/matrices"

echo "archivos extraidos: $(ls "$RAW/matrices" | wc -l) (deberían ser 51)"
