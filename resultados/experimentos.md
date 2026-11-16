# Experimentos comparativos (rama `pruebas`)

Comparación de variantes de pivoteo y de la variante redundante. Los dígitos correctos se cuentan contra la solución decimal exacta, verificada con residuo nulo.

| Variante | Flops | Comparaciones | Intercambios | Mult. máximo | Crecimiento | Bits perdidos (prom.) | ‖r‖∞ | Dígitos (mín.) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Sin pivoteo | 1354 | 0 | 0 | 2.7526e+13 | 2641.2243 | 2.93 | 1.7718e-07 | 7 |
| Pivoteo parcial | 1354 | 66 | 9 | 6.1988e-01 | 1.0000 | 1.06 | 5.3730e-10 | 9 |
| Pivoteo parcial escalado | 1354 | 638 | 7 | 2.4656e+09 | 1.7832 | 0.78 | 9.7126e-10 | 11 |
| Pivoteo parcial escalado de escalas fijas | 1354 | 198 | 7 | 1.8812e+07 | 1.8208 | 0.79 | 5.4256e-10 | 10 |
| Escalado + multiplicador recalculado por columna | 1926 | 638 | 7 | 2.4656e+09 | 1.7832 | 0.78 | 9.7126e-10 | 11 |

La columna de flops es aritmética (1 flop = 1 operación elemental). Las comparaciones y los intercambios son sobrecosto de estabilidad y se informan aparte.
