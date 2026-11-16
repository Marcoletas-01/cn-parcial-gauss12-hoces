"""Experimentos comparativos de estrategias de pivoteo (rama `pruebas`).

Este script **no** forma parte de la entrega evaluada (rama `main`): sirve
para sustentar con mediciones propias las afirmaciones del informe sobre el
compromiso entre costo operacional y exactitud de mantisa.

Compara las variantes de pivoteo disponibles (`ninguno`, `parcial`, `escalado`
y `escalado_fijo`, que usa escalas de fila fijas al estilo LAPACK) más la
variante redundante que recalcula el multiplicador en cada columna.

Uso:

    python src/experimentos.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from conteo import Contador  # noqa: E402
from gauss import (  # noqa: E402
    a_flotante,
    b_flotante,
    eliminar,
    leer_matriz,
    leer_vector,
)
from verificacion import (  # noqa: E402
    detectar_solucion_exacta,
    digitos_correctos,
    referencia_decimal,
    residuo_extendido,
)

DATOS = RAIZ / "datos"
RESULTADOS = RAIZ / "resultados"


def main() -> None:
    """Ejecuta el comparativo y escribe `resultados/experimentos.md`."""
    RESULTADOS.mkdir(parents=True, exist_ok=True)
    matriz_cad = leer_matriz(DATOS / "matriz_a.csv")
    vector_cad = leer_vector(DATOS / "vector_b.csv")
    A = a_flotante(matriz_cad)
    b = b_flotante(vector_cad)

    x_referencia, _ = referencia_decimal(matriz_cad, vector_cad, 80)
    deteccion = detectar_solucion_exacta(matriz_cad, vector_cad, x_referencia)
    x_verdadera = deteccion[0] if deteccion else x_referencia

    variantes = [
        ("Sin pivoteo", dict(pivote="ninguno")),
        ("Pivoteo parcial", dict(pivote="parcial")),
        ("Pivoteo parcial escalado", dict(pivote="escalado")),
        ("Pivoteo parcial escalado de escalas fijas", dict(pivote="escalado_fijo")),
        (
            "Escalado + multiplicador recalculado por columna",
            dict(pivote="escalado", recalcular_multiplicadores=True),
        ),
    ]

    lineas = [
        "# Experimentos comparativos (rama `pruebas`)",
        "",
        "Comparación de variantes de pivoteo y de la variante redundante. "
        "Los dígitos correctos se cuentan contra la solución decimal exacta, "
        "verificada con residuo nulo.",
        "",
        "| Variante | Flops | Comparaciones | Intercambios | Mult. máximo | "
        "Crecimiento | Bits perdidos (prom.) | ‖r‖∞ | Dígitos (mín.) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]

    for nombre, opciones in variantes:
        corrida = eliminar(A, b, contador=Contador(), **opciones)
        residuo = residuo_extendido(matriz_cad, vector_cad, corrida.x, 80)
        digitos = digitos_correctos(corrida.x, x_verdadera)
        minimo = min(fila["digitos"] for fila in digitos)

        lineas.append(
            f"| {nombre} | {corrida.contador.flops_totales} | "
            f"{corrida.contador.comparaciones_totales} | {corrida.n_intercambios} | "
            f"{corrida.max_multiplicador:.4e} | {corrida.crecimiento:.4f} | "
            f"{corrida.bits_cancelacion_promedio:.2f} | "
            f"{float(residuo.norma_inf):.4e} | {minimo} |"
        )
    lineas.append("")
    lineas.append(
        "La columna de flops es aritmética (1 flop = 1 operación elemental). "
        "Las comparaciones y los intercambios son sobrecosto de estabilidad y "
        "se informan aparte."
    )
    (RESULTADOS / "experimentos.md").write_text(
        "\n".join(lineas) + "\n", encoding="utf-8"
    )
    print("\n".join(lineas))


if __name__ == "__main__":
    main()
