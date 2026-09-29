# Prueba técnica — Bioinformática senior

Análisis de 1000 Genomas (30x, GRCh38): ancestría, variantes clínicas, farmacogenómica y nutrigenómica.

## Requisitos
- macOS (probado en Apple M2) o Linux
- 16 GB de RAM y ~30 GB libres
- Miniforge (conda)

## Instalación
```bash
bash scripts/00_setup.sh          # crea el entorno y descarga PharmCAT
conda activate bioinfo-prueba
bash scripts/00_check.sh          # verifica herramientas y acceso a los datos
```

`environment.yml` lista las herramientas; `environment.lock.yml` fija las versiones exactas instaladas.
En Mac con chip Apple el entorno se instala en modo Intel (`CONDA_SUBDIR=osx-64`) porque ADMIXTURE no tiene binario ARM.

## Estructura
```
scripts/     pasos numerados (00_setup, 01_descarga, 02_qc, ...)
config/      regiones, rutas y semilla
results/     tablas y figuras
tools/       PharmCAT (descargado por 00_setup.sh)
```
