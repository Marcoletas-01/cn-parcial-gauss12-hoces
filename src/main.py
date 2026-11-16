"""Punto de entrada único del proyecto.

Uso:

    python src/main.py

Reproduce todo el contenido de `resultados/` sin intervención manual:
`solucion.csv`, `traza_eliminacion.txt`, `conteo_flops.md` y `residuo.txt`.
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
    analisis_casi_dependencia,
    condicionamiento,
    control_biblioteca,
    detectar_solucion_exacta,
    digitos_correctos,
    norma_infinito,
    referencia_decimal,
    refinamiento_iterativo,
    residuo_extendido,
)

DATOS = RAIZ / "datos"
RESULTADOS = RAIZ / "resultados"

DIGITOS_EXTENDIDOS = 80
COTA_MINIMA = 1486
TEORICO_2N3_SOBRE_3 = 2 * 12 ** 3 / 3


def _formatear(valor: float, digitos: int = 10) -> str:
    """Formato de tabla para un flotante."""
    return f"{valor:.{digitos}e}"


def _formatear_exacto(valor: float) -> str:
    """Formato con los 17 dígitos significativos de la doble precisión."""
    return f"{valor:.17g}"


def main() -> None:
    """Ejecuta las tres corridas y regenera todos los resultados."""
    RESULTADOS.mkdir(parents=True, exist_ok=True)

    matriz_cad = leer_matriz(DATOS / "matriz_a.csv")
    vector_cad = leer_vector(DATOS / "vector_b.csv")
    A = a_flotante(matriz_cad)
    b = b_flotante(vector_cad)
    n = len(b)

    corridas = {}
    for modo in ("ninguno", "parcial", "escalado"):
        corridas[modo] = eliminar(A, b, pivote=modo, contador=Contador())

    # Corrida de control: multiplicadores recalculados por columna.
    corrida_redundante = eliminar(
        A, b, pivote="escalado", recalcular_multiplicadores=True, contador=Contador()
    )

    # Referencia en precisión extendida.
    x_referencia, residuo_referencia = referencia_decimal(
        matriz_cad, vector_cad, DIGITOS_EXTENDIDOS
    )
    # Solución verdadera: decimal finito que reproduce b sin residuo.
    deteccion = detectar_solucion_exacta(matriz_cad, vector_cad, x_referencia)
    if deteccion is not None:
        x_verdadera, decimales_verdadera = deteccion
    else:
        x_verdadera, decimales_verdadera = x_referencia, 0

    # Verificaciones de cada corrida.
    analisis = {}
    for modo, corrida in corridas.items():
        residuo = residuo_extendido(matriz_cad, vector_cad, corrida.x, DIGITOS_EXTENDIDOS)
        x_refinado, historial = refinamiento_iterativo(
            matriz_cad,
            vector_cad,
            corrida.l,
            corrida.u,
            corrida.perm,
            corrida.x,
            iteraciones=3,
            digitos=DIGITOS_EXTENDIDOS,
        )
        digitos = digitos_correctos(corrida.x, x_verdadera)
        digitos_refinados = digitos_correctos(x_refinado, x_verdadera)
        condicion = condicionamiento(A, corrida.l, corrida.u, corrida.perm)
        analisis[modo] = {
            "corrida": corrida,
            "residuo": residuo,
            "refinamiento": historial,
            "x_refinado": x_refinado,
            "digitos": digitos,
            "digitos_refinados": digitos_refinados,
            "condicionamiento": condicion,
            "residuo_doble": norma_infinito(
                [b[i] - sum(A[i][j] * corrida.x[j] for j in range(n)) for i in range(n)]
            ),
        }

    control = control_biblioteca(A, b)
    dependencia = analisis_casi_dependencia(A)

    _escribir_solucion(corridas, analisis, x_referencia, x_verdadera)
    _escribir_traza(corridas)
    _escribir_conteo(corridas, corrida_redundante)
    _escribir_residuo(
        corridas, analisis, control, residuo_referencia, dependencia, x_verdadera
    )

    _resumen(corridas, analisis)


# ---------------------------------------------------------------------- #
# Escritura de resultados
# ---------------------------------------------------------------------- #
def _escribir_solucion(corridas, analisis, x_referencia, x_verdadera) -> None:
    """Genera `resultados/solucion.csv`."""
    ruta = RESULTADOS / "solucion.csv"
    lineas = [
        "indice,"
        "x_sin_pivote,digitos_sin_pivote,"
        "x_parcial,digitos_parcial,"
        "x_escalado,digitos_escalado,"
        "x_referencia_80_digitos,x_verdadera,error_relativo_escalado"
    ]
    for i in range(len(corridas["escalado"].x)):
        fila = [str(i + 1)]
        for modo in ("ninguno", "parcial", "escalado"):
            fila.append(_formatear_exacto(corridas[modo].x[i]))
            fila.append(str(analisis[modo]["digitos"][i]["digitos"]))
        fila.append(str(x_referencia[i]))
        fila.append(str(x_verdadera[i]))
        fila.append(f"{analisis['escalado']['digitos'][i]['error_relativo']:.6e}")
        lineas.append(",".join(fila))
    ruta.write_text("\n".join(lineas) + "\n", encoding="utf-8")


def _escribir_traza(corridas) -> None:
    """Genera `resultados/traza_eliminacion.txt`."""
    titulos = {
        "ninguno": "Corrida SIN pivoteo",
        "parcial": "Corrida con pivoteo PARCIAL",
        "escalado": "Corrida con pivoteo PARCIAL ESCALADO",
    }
    bloques = []
    for modo, corrida in corridas.items():
        bloques.append("=" * 78)
        bloques.append(f"{titulos[modo]}  (matriz tras cada uno de los 11 pasos)")
        bloques.append("=" * 78)
        bloques.append(
            f"Multiplicador máximo: {corrida.max_multiplicador:.6e}   "
            f"Intercambios de fila: {corrida.n_intercambios}   "
            f"Factor de crecimiento: {corrida.crecimiento:.6e}"
        )
        for indice, matriz in enumerate(corrida.traza, start=1):
            bloques.append("")
            bloques.append(f"Paso {indice} (k = {indice})")
            bloques.append("-" * 78)
            for numero, fila in enumerate(matriz, start=1):
                bloques.append(
                    f"F{numero:02d} " + " ".join(f"{valor: .12e}" for valor in fila)
                )
        bloques.append("")
    (RESULTADOS / "traza_eliminacion.txt").write_text(
        "\n".join(bloques) + "\n", encoding="utf-8"
    )


def _escribir_conteo(corridas, corrida_redundante) -> None:
    """Genera `resultados/conteo_flops.md`."""
    titulos = {
        "ninguno": "Corrida sin pivoteo",
        "parcial": "Corrida con pivoteo parcial",
        "escalado": "Corrida con pivoteo parcial escalado (estrategia declarada)",
    }
    lineas: list[str] = []
    lineas.append("# Conteo de operaciones (flops)")
    lineas.append("")
    lineas.append(
        "Convención declarada: **1 flop = 1 operación aritmética elemental** "
        "(multiplicación, división, suma o resta). La anulación explícita de la "
        "subdiagonal, las comparaciones de búsqueda de pivote y los intercambios "
        "de fila **se contabilizan aparte** como sobrecosto de estabilidad. "
        "Todos los números provienen del contador instrumentado en `src/conteo.py`, "
        "no de una fórmula escrita a mano."
    )
    lineas.append("")

    for modo, corrida in corridas.items():
        contador = corrida.contador
        lineas.append(f"## {titulos[modo]}")
        lineas.append("")
        lineas.append(
            "| Paso | Multiplicaciones | Divisiones | Sumas/restas | Flops del paso | "
            "Acumulado | Comparaciones | Intercambios |"
        )
        lineas.append("|---:|---:|---:|---:|---:|---:|---:|---:|")
        acumulado = 0
        for paso in contador.pasos:
            acumulado += paso.flops
            lineas.append(
                f"| {paso.paso} | {paso.multiplications} | {paso.divisiones} | "
                f"{paso.sumas_restas} | {paso.flops} | {acumulado} | "
                f"{paso.comparaciones} | {paso.intercambios} |"
            )
        eliminacion = contador.eliminacion
        sustitucion = contador.sustitucion
        lineas.append(
            f"| **Eliminación** | {eliminacion.multiplications} | "
            f"{eliminacion.divisiones} | {eliminacion.sumas_restas} | "
            f"{eliminacion.flops} | {eliminacion.flops} | "
            f"{eliminacion.comparaciones} | {eliminacion.intercambios} |"
        )
        lineas.append(
            f"| **Sustitución regresiva** | {sustitucion.multiplications} | "
            f"{sustitucion.divisiones} | {sustitucion.sumas_restas} | "
            f"{sustitucion.flops} | "
            f"{contador.flops_eliminacion + contador.flops_sustitucion} | 0 | 0 |"
        )
        lineas.append(
            f"| **TOTAL aritmético** | "
            f"{eliminacion.multiplications + sustitucion.multiplications} | "
            f"{eliminacion.divisiones + sustitucion.divisiones} | "
            f"{eliminacion.sumas_restas + sustitucion.sumas_restas} | "
            f"**{contador.flops_totales}** | {contador.flops_totales} | "
            f"{contador.comparaciones_totales} | {contador.intercambios_totales} |"
        )
        lineas.append("")

    escalado = corridas["escalado"].contador
    aritmetico = escalado.flops_totales
    con_anulaciones = aritmetico + escalado.asignaciones_totales
    total_estabilidad = escalado.costo_total_con_estabilidad
    lineas.append("## Comparación con las cotas de referencia")
    lineas.append("")
    lineas.append(
        "El criterio de evaluación es el **conteo aritmético**: las comparaciones "
        "de búsqueda de pivote y los intercambios se reportan aparte, tal como "
        "pide la pauta, y por eso la última fila es sólo informativa."
    )
    lineas.append("")
    lineas.append("| Magnitud | Valor | Referencia | Diferencia |")
    lineas.append("|---|---:|---:|---:|")
    lineas.append(
        f"| Flops aritméticos (eliminación + sustitución) | {aritmetico} | "
        f"{COTA_MINIMA} | {(aritmetico - COTA_MINIMA) / COTA_MINIMA * 100:+.2f} % |"
    )
    lineas.append(
        f"| + anulaciones explícitas de la subdiagonal | {con_anulaciones} | "
        f"{COTA_MINIMA} | {(con_anulaciones - COTA_MINIMA) / COTA_MINIMA * 100:+.2f} % |"
    )
    lineas.append(
        f"| Comparaciones de búsqueda de pivote (escalado) | "
        f"{escalado.comparaciones_totales} | 66 (pivoteo parcial) | no son flops |"
    )
    lineas.append(
        f"| Intercambios de fila | {escalado.intercambios_totales} | 11 (máximo posible) "
        f"| no son flops |"
    )
    lineas.append(
        f"| Total con sobrecosto de estabilidad (informativo) | {total_estabilidad} | "
        f"{COTA_MINIMA} + comparaciones | "
        f"{(total_estabilidad - COTA_MINIMA) / COTA_MINIMA * 100:+.2f} % |"
    )
    lineas.append(
        f"| Variante redundante (multiplicador recalculado por columna) | "
        f"{corrida_redundante.contador.flops_totales} | ~1750 (cota de la pauta) | "
        f"{(corrida_redundante.contador.flops_totales - 1750) / 1750 * 100:+.2f} % |"
    )
    lineas.append(
        f"| Estimación asintótica 2n³/3 (solo factorización, sin b ni sustitución) | "
        f"{TEORICO_2N3_SOBRE_3:.0f} | 1152 | "
        f"{(TEORICO_2N3_SOBRE_3 - 1152) / 1152 * 100:+.2f} % |"
    )
    lineas.append("")
    lineas.append(
        "La cota 2n³/3 = 1152 cuenta sólo la factorización (n³/3 multiplicaciones y "
        "n³/3 sumas) y omite la actualización del lado derecho, la anulación "
        "explícita de la subdiagonal y la sustitución regresiva, que en este "
        "proyecto sí se ejecutan y se cuentan."
    )
    (RESULTADOS / "conteo_flops.md").write_text(
        "\n".join(lineas) + "\n", encoding="utf-8"
    )


def _escribir_residuo(
    corridas, analisis, control, residuo_referencia, dependencia, x_verdadera
) -> None:
    """Genera `resultados/residuo.txt`."""
    titulos = {
        "ninguno": "SIN pivoteo",
        "parcial": "pivoteo PARCIAL",
        "escalado": "pivoteo PARCIAL ESCALADO",
    }
    lineas: list[str] = []
    lineas.append("VERIFICACIÓN DEL SISTEMA 12x12")
    lineas.append("=" * 78)
    lineas.append(
        "Residuo r = b - A x evaluado en precisión extendida (decimal, 80 dígitos); "
        "x se convierte a Decimal de forma exacta desde la doble precisión."
    )
    lineas.append("")
    for modo, corrida in corridas.items():
        datos = analisis[modo]
        residuo = datos["residuo"]
        lineas.append("-" * 78)
        lineas.append(f"Corrida: {titulos[modo]}")
        lineas.append("-" * 78)
        lineas.append("  i   r(i)")
        for i, valor in enumerate(residuo.componentes, start=1):
            lineas.append(f"{i:3d}   {valor:.6e}")
        lineas.append(f"  ||r||_inf (extendido)   = {residuo.norma_inf:.6e}")
        lineas.append(f"  ||r||_inf / ||b||_inf   = {residuo.norma_inf_relativo:.6e}")
        lineas.append(
            f"  ||r||_inf (doble precisión) = {datos['residuo_doble']:.6e}"
        )
        lineas.append("")
        lineas.append("  Refinamiento iterativo (residuo extendido + corrección con L, U)")
        lineas.append(
            "   iter        ||r||_inf        ||r||_inf/||b||_inf      ||d||_inf        paso"
        )
        for paso in datos["refinamiento"]:
            lineas.append(
                f"   {paso.iteracion:4d}   {paso.norma_residuo:.6e}   "
                f"{paso.norma_residuo_relativa:.6e}   {paso.correccion_inf:.6e}   "
                f"{paso.paso:.6e}"
            )
        digitos = [fila["digitos"] for fila in datos["digitos"]]
        digitos_ref = [fila["digitos"] for fila in datos["digitos_refinados"]]
        lineas.append("")
        lineas.append(
            f"  Dígitos correctos (sin refinar): {digitos}  -> mínimo {min(digitos)}"
        )
        lineas.append(
            f"  Dígitos correctos (refinado)   : {digitos_ref}  -> mínimo {min(digitos_ref)}"
        )
        lineas.append("")
        lineas.append("  Bits de mantisa perdidos")
        lineas.append(
            f"    factor de crecimiento            = {datos['corrida'].crecimiento:.6e}"
            f"  => log2 = {math.log2(datos['corrida'].crecimiento):.2f} bits"
        )
        lineas.append(
            f"    bits perdidos por cancelación     = {datos['corrida'].bits_cancelacion_promedio:.2f}"
            f" en promedio (máximo {datos['corrida'].bits_cancelacion_maximo:.2f})"
            f" sobre {datos['corrida'].actualizaciones} actualizaciones"
        )
        lineas.append("")

    condicion = analisis["escalado"]["condicionamiento"]
    lineas.append("-" * 78)
    lineas.append("CONDICIONAMIENTO")
    lineas.append("-" * 78)
    lineas.append(f"  sigma_max (SVD, referencia)    = {condicion['sigma_max']:.6e}")
    lineas.append(f"  sigma_min (SVD, referencia)    = {condicion['sigma_min']:.6e}")
    lineas.append(f"  kappa_2 (SVD, referencia)     = {condicion['kappa_2']:.6e}")
    lineas.append(
        f"  sigma_max (factorización propia) = {condicion['sigma_max_propio']:.6e}"
    )
    lineas.append(
        f"  sigma_min (factorización propia) = {condicion['sigma_min_propio']:.6e}"
    )
    lineas.append(
        f"  kappa_2 propio (σmax/σmin)    = {condicion['kappa_2_propio']:.6e}"
    )
    lineas.append(f"  eps maquina (IEEE 754)       = {condicion['epsilon']:.6e}")
    lineas.append(
        f"  cota relativa kappa_2 * eps   = {condicion['cota_relativa']:.6e} "
        f"(=> {condicion['digitos_cota']:.2f} dígitos garantizados por la cota)"
    )
    lineas.append(
        f"  residuo de la referencia (Decimal, 80 dígitos) = {residuo_referencia:.6e}"
    )
    lineas.append("")
    lineas.append("  Solución verdadera (decimal finito, verificada con residuo nulo):")
    for i, valor in enumerate(x_verdadera, start=1):
        lineas.append(f"    x({i:2d}) = {valor}")
    lineas.append("")

    lineas.append("-" * 78)
    lineas.append("ANÁLISIS DE CASI DEPENDENCIA LINEAL (filas casi combinación de otras)")
    lineas.append("-" * 78)
    lineas.append(
        "Para cada terna de filas se ajustó por mínimos cuadrados la fila "
        "objetivo como combinación lineal de dos filas; se listan los diez "
        "mejores ajustes (menor residuo relativo)."
    )
    lineas.append("")
    lineas.append(
        "  fila objetivo   filas combinadas      alfa           beta"
        "       error relativo"
    )
    for fila in dependencia:
        lineas.append(
            f"  {fila['fila_objetivo']:>8}      "
            f"{fila['fila_i']:>3} y {fila['fila_j']:<3}   "
            f"{fila['alfa']:>12.6e}  {fila['beta']:>12.6e}   "
            f"{fila['error_relativo']:.6e}"
        )
    lineas.append("")

    lineas.append("-" * 78)
    lineas.append("CONTROL CON BIBLOTECA (no reemplaza la implementación propia)")
    lineas.append("-" * 78)
    escalado = corridas["escalado"].x
    x_numpy = control["x_numpy"]
    x_scipy = control["x_scipy"]
    lineas.append("   i   x_propio (escalado)        x_numpy                diferencia rel.")
    for i in range(len(escalado)):
        dif = abs(escalado[i] - x_numpy[i]) / abs(x_numpy[i])
        lineas.append(
            f"  {i + 1:3d}   {escalado[i]: .17g}   {x_numpy[i]: .17g}   {dif:.3e}"
        )
    if x_scipy is not None:
        dif_max = max(
            abs(escalado[i] - x_scipy[i]) / abs(x_scipy[i]) for i in range(len(escalado))
        )
        lineas.append(
            f"  Diferencia relativa máxima frente a scipy.linalg.solve: {dif_max:.3e}"
        )
    (RESULTADOS / "residuo.txt").write_text("\n".join(lineas) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------- #
# Resumen en pantalla
# ---------------------------------------------------------------------- #
def _resumen(corridas, analisis) -> None:
    """Imprime el resumen de la ejecución."""
    print("=" * 78)
    print("ELIMINACIÓN GAUSSIANA 12x12 - RESUMEN")
    print("=" * 78)
    print(
        f"{'Estrategia':<28}{'Flops':>8}{'Mult.máx':>14}{'Crecim.':>12}"
        f"{'||r||inf':>13}{'dígitos':>9}"
    )
    for modo in ("ninguno", "parcial", "escalado"):
        corrida = corridas[modo]
        datos = analisis[modo]
        digitos = min(fila["digitos"] for fila in datos["digitos"])
        print(
            f"{modo:<28}{corrida.contador.flops_totales:>8}"
            f"{corrida.max_multiplicador:>14.4e}{corrida.crecimiento:>12.4e}"
            f"{float(datos['residuo'].norma_inf):>13.3e}{digitos:>9d}"
        )
    condicion = analisis["escalado"]["condicionamiento"]
    print("-" * 78)
    print(f"kappa_2 = {condicion['kappa_2']:.6e}   cota kappa_2*eps = {condicion['cota_relativa']:.3e}")
    print(f"Archivos generados en: {RESULTADOS}")


if __name__ == "__main__":
    main()
