#!/usr/bin/env bash
# 00_setup.sh — Crea el entorno conda y fija las versiones exactas.
set -euo pipefail

# 1) Rosetta 2: permite correr programas Intel en un Mac con chip Apple.
if ! /usr/bin/pgrep -q oahd; then
  echo ">> Instalando Rosetta 2..."
  softwareupdate --install-rosetta --agree-to-license
fi

# 2) Crear el entorno en modo Intel (osx-64).
echo ">> Creando el entorno bioinfo-prueba..."
CONDA_SUBDIR=osx-64 conda env create -f environment.yml

# 3) Dejar el entorno marcado como osx-64 para futuras instalaciones.
conda run -n bioinfo-prueba conda config --env --set subdir osx-64

# 4) Congelar versiones exactas
conda env export -n bioinfo-prueba --no-builds > environment.lock.yml
echo ">> Versiones guardadas en environment.lock.yml"

# 5) Descargar PharmCAT (no está en conda; se usa el paquete oficial "pipeline").
#    Requisitos oficiales: Java >= 17, Python >= 3.10.14, bcftools/htslib >= 1.18.
mkdir -p tools
PHARMCAT_VERSION="3.4.0"
echo ">> Descargando PharmCAT ${PHARMCAT_VERSION}..."
curl -fL -o tools/pharmcat-pipeline.tar.gz \
  "https://github.com/PharmGKB/PharmCAT/releases/download/v${PHARMCAT_VERSION}/pharmcat-pipeline-${PHARMCAT_VERSION}.tar.gz"
mkdir -p tools/pharmcat
tar -xzf tools/pharmcat-pipeline.tar.gz -C tools/pharmcat
conda run -n bioinfo-prueba pip install -r tools/pharmcat/requirements.txt
rm tools/pharmcat-pipeline.tar.gz

echo "Creado en el env, ahora activarlo"