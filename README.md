
# Inpactor3

Clasificador generalizable de elementos transponibles (plantas, hongos, humano) que sustituye la red `Inpactor2_Class` dentro del pipeline de Inpactor2.

**Proyecto Integrador II** · Universidad de Caldas · agosto–diciembre 2026

## Contexto

Inpactor3 es un reemplazo de la red neuronal de clasificación (`Inpactor2_Class`) del pipeline Inpactor2, con el objetivo de generalizar la detección/clasificación de LTR-retrotransposones a distintos reinos (plantas, hongos, humano), no solo plantas.

Referencias del linaje del proyecto:

- Inpactor v1 (C + MPI) — referencia histórica: https://github.com/simonorozcoarias/Inpactor
- Inpactor2 (Deep Learning, Python) — software que se extiende: https://github.com/simonorozcoarias/Inpactor2

## Estructura del repositorio

```
Inpactor3/
├── data/         # Datasets (no versionar datos pesados; usar DVC o enlaces)
├── notebooks/    # Exploración, EDA, prototipos de arquitecturas
├── src/          # Código fuente del clasificador (entrenamiento, inferencia)
├── models/       # Pesos entrenados / checkpoints
├── tests/        # Pruebas unitarias e integración
├── scripts/      # Utilidades CLI, preprocesamiento, evaluación
├── docs/         # Documentación técnica, bitácora, entregables
└── README.md
```

## Estado

Semana 2 — repositorio inicializado. Próximo entregable según Lineamientos Generales v4.

## Entorno de desarrollo

Ver [`docs/setup.md`](docs/setup.md) para la guía completa de instalación (Linux/WSL2/macOS Intel) siguiendo la guía del curso.

## Licencia

Por definir (se alineará con Inpactor2, GPL-3.0).
