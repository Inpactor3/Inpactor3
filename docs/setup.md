# Guía de instalación y arranque en Visual Studio Code

Proyecto Integrador II · Universidad de Caldas · agosto–diciembre 2026.

Este documento condensa la guía del curso para dejar listo el entorno de trabajo con:

1. Inpactor v1 (referencia histórica, C + MPI).
2. Inpactor2 (el software que Inpactor3 extiende).
3. Repositorio Inpactor3 con la estructura de la propuesta.

## 0. Antes de empezar

### 0.1 Sistema operativo

| Tu equipo | Qué usar |
|---|---|
| Linux (Ubuntu 22.04/24.04) | Directo |
| Windows 10/11 | WSL2 con Ubuntu |
| macOS Intel | `Inpactor2_mac.yml` |
| macOS Apple Silicon (M1–M4) | No soportado; usar Linux, WSL2 o Google Colab |

### 0.2 WSL2 (solo Windows)

En PowerShell como administrador:

```powershell
wsl --install -d Ubuntu
```

Reiniciar el equipo y abrir Ubuntu para completar el usuario.

## 1. Herramientas base

- Git
- Miniconda / Anaconda
- Visual Studio Code + extensiones: Python, Jupyter, Remote - WSL (si aplica)
- Compiladores: `build-essential`, `gcc`, `make` (para Inpactor v1)

```bash
sudo apt update && sudo apt install -y git build-essential wget curl
```

Instalar Miniconda:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda3
source $HOME/miniconda3/etc/profile.d/conda.sh
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
```

## 2. Inpactor v1 (referencia)

Ver README de [Inpactor](https://github.com/simonorozcoarias/Inpactor). Requiere OpenMPI 1.8.8 compilado manualmente, `blast-legacy`, Emboss, Wise2, Mafft, LTR-FINDER. Compilar con `mpicc Inpactor.c -o Inpactor`.

## 3. Inpactor2

```bash
git clone https://github.com/simonorozcoarias/Inpactor2.git
cd Inpactor2
conda env create -f Inpactor2.yml
conda activate Inpactor2
python3 Inpactor2.py -f Testing/toy_genome.fasta -o Testing/toy_execution -a no
```

Verifica que la salida coincida con `Inpactor2_library_successfull_run.fasta` y `Inpactor2_predictions_successfull_run.tab`.

## 4. Inpactor3 (este repositorio)

### 4.1 Inicializar y subir

```bash
cd Inpactor3
git init
git add .
git commit -m "chore: initial scaffold"
git branch -M main
git remote add origin https://github.com/Inpactor3/Inpactor3.git
git push -u origin main
```

### 4.2 Entorno de desarrollo

Se creará un `environment.yml` propio (Python 3.10+, PyTorch o TensorFlow según decisión de equipo) durante la semana 3.

## 5. Entregable semana 2

Según Lineamientos Generales v4:

- [ ] Repositorio creado y compartido con docentes.
- [ ] README y estructura inicial.
- [ ] Bitácora del equipo iniciada en `docs/bitacora.md`.
- [ ] Evidencia de instalación de Inpactor2 (captura de la prueba con `toy_genome.fasta`).
