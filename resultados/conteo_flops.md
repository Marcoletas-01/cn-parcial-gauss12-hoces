# Conteo de operaciones (flops)

Convención declarada: **1 flop = 1 operación aritmética elemental** (multiplicación, división, suma o resta). La anulación explícita de la subdiagonal, las comparaciones de búsqueda de pivote y los intercambios de fila **se contabilizan aparte** como sobrecosto de estabilidad. Todos los números provienen del contador instrumentado en `src/conteo.py`, no de una fórmula escrita a mano.

## Corrida sin pivoteo

| Paso | Multiplicaciones | Divisiones | Sumas/restas | Flops del paso | Acumulado | Comparaciones | Intercambios |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 132 | 11 | 132 | 275 | 275 | 0 | 0 |
| 2 | 110 | 10 | 110 | 230 | 505 | 0 | 0 |
| 3 | 90 | 9 | 90 | 189 | 694 | 0 | 0 |
| 4 | 72 | 8 | 72 | 152 | 846 | 0 | 0 |
| 5 | 56 | 7 | 56 | 119 | 965 | 0 | 0 |
| 6 | 42 | 6 | 42 | 90 | 1055 | 0 | 0 |
| 7 | 30 | 5 | 30 | 65 | 1120 | 0 | 0 |
| 8 | 20 | 4 | 20 | 44 | 1164 | 0 | 0 |
| 9 | 12 | 3 | 12 | 27 | 1191 | 0 | 0 |
| 10 | 6 | 2 | 6 | 14 | 1205 | 0 | 0 |
| 11 | 2 | 1 | 2 | 5 | 1210 | 0 | 0 |
| **Eliminación** | 572 | 66 | 572 | 1210 | 1210 | 0 | 0 |
| **Sustitución regresiva** | 66 | 12 | 66 | 144 | 1354 | 0 | 0 |
| **TOTAL aritmético** | 638 | 78 | 638 | **1354** | 1354 | 0 | 0 |

## Corrida con pivoteo parcial

| Paso | Multiplicaciones | Divisiones | Sumas/restas | Flops del paso | Acumulado | Comparaciones | Intercambios |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 132 | 11 | 132 | 275 | 275 | 11 | 1 |
| 2 | 110 | 10 | 110 | 230 | 505 | 10 | 1 |
| 3 | 90 | 9 | 90 | 189 | 694 | 9 | 1 |
| 4 | 72 | 8 | 72 | 152 | 846 | 8 | 1 |
| 5 | 56 | 7 | 56 | 119 | 965 | 7 | 1 |
| 6 | 42 | 6 | 42 | 90 | 1055 | 6 | 1 |
| 7 | 30 | 5 | 30 | 65 | 1120 | 5 | 0 |
| 8 | 20 | 4 | 20 | 44 | 1164 | 4 | 1 |
| 9 | 12 | 3 | 12 | 27 | 1191 | 3 | 1 |
| 10 | 6 | 2 | 6 | 14 | 1205 | 2 | 0 |
| 11 | 2 | 1 | 2 | 5 | 1210 | 1 | 1 |
| **Eliminación** | 572 | 66 | 572 | 1210 | 1210 | 66 | 9 |
| **Sustitución regresiva** | 66 | 12 | 66 | 144 | 1354 | 0 | 0 |
| **TOTAL aritmético** | 638 | 78 | 638 | **1354** | 1354 | 66 | 9 |

## Corrida con pivoteo parcial escalado (estrategia declarada)

| Paso | Multiplicaciones | Divisiones | Sumas/restas | Flops del paso | Acumulado | Comparaciones | Intercambios |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 132 | 11 | 132 | 275 | 275 | 143 | 1 |
| 2 | 110 | 10 | 110 | 230 | 505 | 120 | 0 |
| 3 | 90 | 9 | 90 | 189 | 694 | 99 | 1 |
| 4 | 72 | 8 | 72 | 152 | 846 | 80 | 1 |
| 5 | 56 | 7 | 56 | 119 | 965 | 63 | 1 |
| 6 | 42 | 6 | 42 | 90 | 1055 | 48 | 1 |
| 7 | 30 | 5 | 30 | 65 | 1120 | 35 | 1 |
| 8 | 20 | 4 | 20 | 44 | 1164 | 24 | 1 |
| 9 | 12 | 3 | 12 | 27 | 1191 | 15 | 0 |
| 10 | 6 | 2 | 6 | 14 | 1205 | 8 | 0 |
| 11 | 2 | 1 | 2 | 5 | 1210 | 3 | 0 |
| **Eliminación** | 572 | 66 | 572 | 1210 | 1210 | 638 | 7 |
| **Sustitución regresiva** | 66 | 12 | 66 | 144 | 1354 | 0 | 0 |
| **TOTAL aritmético** | 638 | 78 | 638 | **1354** | 1354 | 638 | 7 |

## Comparación con las cotas de referencia

El criterio de evaluación es el **conteo aritmético**: las comparaciones de búsqueda de pivote y los intercambios se reportan aparte, tal como pide la pauta, y por eso la última fila es sólo informativa.

| Magnitud | Valor | Referencia | Diferencia |
|---|---:|---:|---:|
| Flops aritméticos (eliminación + sustitución) | 1354 | 1486 | -8.88 % |
| + anulaciones explícitas de la subdiagonal | 1420 | 1486 | -4.44 % |
| Comparaciones de búsqueda de pivote (escalado) | 638 | 66 (pivoteo parcial) | no son flops |
| Intercambios de fila | 7 | 11 (máximo posible) | no son flops |
| Total con sobrecosto de estabilidad (informativo) | 2058 | 1486 + comparaciones | +38.49 % |
| Variante redundante (multiplicador recalculado por columna) | 1926 | ~1750 (cota de la pauta) | +10.06 % |
| Estimación asintótica 2n³/3 (solo factorización, sin b ni sustitución) | 1152 | 1152 | +0.00 % |

La cota 2n³/3 = 1152 cuenta sólo la factorización (n³/3 multiplicaciones y n³/3 sumas) y omite la actualización del lado derecho, la anulación explícita de la subdiagonal y la sustitución regresiva, que en este proyecto sí se ejecutan y se cuentan.
