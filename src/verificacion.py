"""Verificación independiente del resultado: residuo en precisión extendida,
refinamiento iterativo, condicionamiento y control con biblioteca.

La precisión extendida se logra con `decimal` (biblioteca estándar), que
trabaja con base 10 exacta: los valores de `datos/matriz_a.csv` y
`datos/vector_b.csv` se leen como cadenas y se convierten a `Decimal` sin
perder ni un dígito, a diferencia de lo que ocurre al pasarlos por `float`.

Ninguna rutina de biblioteca resuelve el sistema: `numpy` y `scipy` sólo se
usan como control y para estimar σ_min (SVD), y el refinamiento reutiliza los
factores L y U calculados por `gauss.py`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from decimal import Decimal, localcontext

import numpy as np

from gauss import resolver_con_factores


# ---------------------------------------------------------------------- #
# Utilidades
# ---------------------------------------------------------------------- #
def norma_infinito(v: list[float]) -> float:
    """Norma infinito (máximo absoluto) de un vector."""
    return max((abs(valor) for valor in v), default=0.0)


def a_decimal(cadenas: list[list[str]], digitos: int) -> list[list[Decimal]]:
    """Convierte la matriz de cadenas a `Decimal` con el contexto dado."""
    with localcontext() as contexto:
        contexto.prec = digitos
        return [[Decimal(valor) for valor in fila] for fila in cadenas]


def b_decimal(cadenas: list[str], digitos: int) -> list[Decimal]:
    """Convierte el vector de cadenas a `Decimal` con el contexto dado."""
    with localcontext() as contexto:
        contexto.prec = digitos
        return [Decimal(valor) for valor in cadenas]


def float_a_decimal_exacto(valor: float, digitos: int) -> Decimal:
    """Representación decimal exacta (o redondeada a `digitos`) de un float."""
    with localcontext() as contexto:
        contexto.prec = digitos
        return Decimal(valor)


# ---------------------------------------------------------------------- #
# Referencia en precisión extendida
# ---------------------------------------------------------------------- #
def referencia_decimal(
    A_cad: list[list[str]], b_cad: list[str], digitos: int = 80
) -> tuple[list[Decimal], Decimal]:
    """Resuelve el sistema en precisión extendida con `Decimal`.

    Sirve como solución de referencia para contar los dígitos correctos de
    la solución en doble precisión. Devuelve la solución y la norma infinito
    del residuo en precisión extendida.
    """
    with localcontext() as contexto:
        contexto.prec = digitos
        A = a_decimal(A_cad, digitos)
        b = b_decimal(b_cad, digitos)
        n = len(b)

        for k in range(n - 1):
            # Pivoteo parcial en precisión extendida (sólo para la referencia).
            fila = max(range(k, n), key=lambda i: abs(A[i][k]))
            if fila != k:
                A[k], A[fila] = A[fila], A[k]
                b[k], b[fila] = b[fila], b[k]
            pivote = A[k][k]
            for i in range(k + 1, n):
                m = A[i][k] / pivote
                A[i][k] = Decimal(0)
                for j in range(k + 1, n):
                    A[i][j] -= m * A[k][j]
                b[i] -= m * b[k]

        x = [Decimal(0)] * n
        for i in range(n - 1, -1, -1):
            acumulado = b[i]
            for j in range(i + 1, n):
                acumulado -= A[i][j] * x[j]
            x[i] = acumulado / A[i][i]

        residuo = [
            b[i] - sum(A[i][j] * x[j] for j in range(n)) for i in range(n)
        ]
        return x, max(abs(valor) for valor in residuo)


# ---------------------------------------------------------------------- #
# Residuo en precisión extendida
# ---------------------------------------------------------------------- #
@dataclass
class Residuo:
    """Residuo `r = b - A x` evaluado en precisión extendida."""

    componentes: list[Decimal]
    norma_inf: Decimal
    norma_inf_relativo: Decimal
    digitos: int


def residuo_extendido(
    A_cad: list[list[str]],
    b_cad: list[str],
    x: list[float],
    digitos: int = 80,
) -> Residuo:
    """Calcula `r = b - A x` en precisión extendida (decimal, 80 dígitos).

    `x` proviene de la corrida en doble precisión; se convierte a `Decimal`
    de forma exacta (representación binaria del flotante) para que el residuo
    mida el error real de la solución en doble precisión y no un redondeo
    artificial.
    """
    with localcontext() as contexto:
        contexto.prec = digitos
        A = a_decimal(A_cad, digitos)
        b = b_decimal(b_cad, digitos)
        n = len(b)
        x_dec = [float_a_decimal_exacto(valor, digitos) for valor in x]

        componentes = [
            b[i] - sum(A[i][j] * x_dec[j] for j in range(n)) for i in range(n)
        ]
        norma = max(abs(valor) for valor in componentes)
        norma_b = max(abs(valor) for valor in b)
        return Residuo(
            componentes=componentes,
            norma_inf=norma,
            norma_inf_relativo=(norma / norma_b) if norma_b else norma,
            digitos=digitos,
        )


# ---------------------------------------------------------------------- #
# Refinamiento iterativo
# ---------------------------------------------------------------------- #
@dataclass
class PasoRefinamiento:
    """Estado del refinamiento iterativo en una iteración dada."""

    iteracion: int
    norma_residuo: float
    norma_residuo_relativa: float
    correccion_inf: float
    paso: float


def refinamiento_iterativo(
    A_cad: list[list[str]],
    b_cad: list[str],
    L: list[list[float]],
    U: list[list[float]],
    perm: list[int],
    x_inicial: list[float],
    iteraciones: int = 3,
    digitos: int = 80,
) -> tuple[list[float], list[PasoRefinamiento]]:
    """Refinamiento iterativo con residuo en precisión extendida.

    En cada vuelta se calcula el residuo `r = b - A x` en `Decimal` (80
    dígitos), se convierte a doble precisión y se resuelve `A d = r` con los
    **factores ya calculados** (sin volver a triangularizar), actualizando
    `x <- x + d`.
    """
    x = list(x_inicial)
    historial: list[PasoRefinamiento] = []
    for it in range(1, iteraciones + 1):
        residuo = residuo_extendido(A_cad, b_cad, x, digitos)
        r_float = [float(valor) for valor in residuo.componentes]
        d = resolver_con_factores(L, U, perm, r_float)
        correccion = max(abs(valor) for valor in d)
        nuevo = [x[i] + d[i] for i in range(len(x))]
        paso = max(abs(nuevo[i] - x[i]) for i in range(len(x)))
        x = nuevo
        residuo_final = residuo_extendido(A_cad, b_cad, x, digitos)
        historial.append(
            PasoRefinamiento(
                iteracion=it,
                norma_residuo=float(residuo_final.norma_inf),
                norma_residuo_relativa=float(residuo_final.norma_inf_relativo),
                correccion_inf=correccion,
                paso=paso,
            )
        )
    return x, historial


# ---------------------------------------------------------------------- #
# Condicionamiento
# ---------------------------------------------------------------------- #
def estimar_valores_singulares(
    A: list[list[float]],
    L: list[list[float]],
    U: list[list[float]],
    perm: list[int],
    iteraciones: int = 200,
) -> tuple[float, float]:
    """Estima σ_max y σ_min usando **sólo** la factorización propia.

    * σ_max por iteración de la potencia sobre `A^T A` (multiplicaciones por
      la transpuesta mediante la factorización LU).
    * σ_min por iteración inversa sobre `A^T A`, es decir aplicando
      `A^-T A^-1`, que es el método clásico para el valor singular más
      pequeño.

    No se forma la inversa de A: sólo se resuelven sistemas triangulares.
    """
    n = len(A)
    generador = np.random.default_rng(20261116)

    # σ_max: iteración de la potencia sobre A^T A con cociente de Rayleigh.
    # Se usan productos matriz-vector explícitos (no resoluciones).
    v = list(generador.normal(size=n))
    norma = math.sqrt(sum(valor * valor for valor in v))
    v = [valor / norma for valor in v]
    lambda_max = 0.0
    for _ in range(iteraciones):
        w = [sum(A[j][i] * v[j] for j in range(n)) for i in range(n)]  # A^T v
        z = [sum(A[j][i] * w[i] for i in range(n)) for j in range(n)]  # A A^T v
        lambda_max = sum(v[i] * z[i] for i in range(n)) / sum(
            valor * valor for valor in v
        )
        norma = math.sqrt(sum(valor * valor for valor in z))
        if norma == 0.0:
            break
        v = [valor / norma for valor in z]
    sigma_max = math.sqrt(lambda_max)

    # σ_min: iteración inversa sobre A^T A (aplica A^-1 y luego A^-T).
    w = list(generador.normal(size=n))
    norma = math.sqrt(sum(valor * valor for valor in w))
    w = [valor / norma for valor in w]
    mu = 0.0
    for _ in range(iteraciones):
        y = resolver_con_factores(L, U, perm, w)  # A^-1 w
        z = _resolver_transpuesto(L, U, perm, y)  # A^-T A^-1 w
        mu = sum(w[i] * z[i] for i in range(n)) / sum(valor * valor for valor in w)
        norma = math.sqrt(sum(valor * valor for valor in z))
        if norma == 0.0:
            break
        w = [valor / norma for valor in z]
    sigma_min = 1.0 / math.sqrt(mu) if mu > 0 else float("inf")
    return sigma_max, sigma_min


def _resolver_transpuesto(
    L: list[list[float]], U: list[list[float]], perm: list[int], c: list[float]
) -> list[float]:
    """Resuelve `A^T w = c` reutilizando la factorización `P A = L U`."""
    n = len(L)
    # A^T = U^T L^T P  con  P la permutación acumulada.
    # U^T es triangular inferior: sustitución progresiva.
    v = [0.0] * n
    for i in range(n):
        acumulado = c[i]
        for j in range(i):
            acumulado -= U[j][i] * v[j]
        v[i] = acumulado / U[i][i]
    # L^T es triangular superior unitaria: sustitución regresiva.
    u = [0.0] * n
    for i in range(n - 1, -1, -1):
        acumulado = v[i]
        for j in range(i + 1, n):
            acumulado -= L[j][i] * u[j]
        u[i] = acumulado
    # w = P^T u
    w = [0.0] * n
    for k in range(n):
        w[perm[k]] = u[k]
    return w


def condicionamiento(
    A: list[list[float]], L: list[list[float]], U: list[list[float]], perm: list[int]
) -> dict:
    """Estima el condicionamiento del sistema.

    * `kappa_2`: razón entre el mayor y el menor valor singular (SVD de numpy,
      usada como herramienta de análisis, no como resolutora).
    * `kappa_2_propio`: razón estimada con los valores singulares obtenidos
      por iteración de la potencia y iteración inversa sobre los factores LU
      calculados por `gauss.py`.
    """
    matriz = np.array(A, dtype=float)
    valores = np.linalg.svd(matriz, compute_uv=False)
    sigma_max, sigma_min = estimar_valores_singulares(A, L, U, perm)
    kappa2 = float(valores[0] / valores[-1])
    kappa2_propio = sigma_max / sigma_min
    eps = float(np.finfo(float).eps)
    return {
        "sigma_max": float(valores[0]),
        "sigma_min": float(valores[-1]),
        "kappa_2": kappa2,
        "sigma_max_propio": sigma_max,
        "sigma_min_propio": sigma_min,
        "kappa_2_propio": kappa2_propio,
        "epsilon": eps,
        "cota_relativa": kappa2 * eps,
        "digitos_cota": max(0.0, -math.log10(kappa2 * eps)),
    }


# ---------------------------------------------------------------------- #
# Control con biblioteca y conteo de dígitos
# ---------------------------------------------------------------------- #
def control_biblioteca(A: list[list[float]], b: list[float]) -> dict:
    """Contrasta el resultado propio con `numpy` y `scipy` (solo control)."""
    matriz = np.array(A, dtype=float)
    vector = np.array(b, dtype=float)
    x_numpy = np.linalg.solve(matriz, vector)
    resultado = {"x_numpy": x_numpy}
    try:
        import scipy.linalg

        x_scipy = scipy.linalg.solve(matriz, vector)
        resultado["x_scipy"] = x_scipy
    except Exception:  # pragma: no cover - scipy es dependencia declarada
        resultado["x_scipy"] = None
    return resultado


def analisis_casi_dependencia(A: list[list[float]]) -> list[dict]:
    """Detecta filas casi linealmente dependientes.

    Para cada par de filas `(i, j)` se busca la combinación lineal
    `alfa * fila_i + beta * fila_j` que mejor aproxima una tercera fila `k`
    (mínimos cuadrados normales resueltos en aritmética propia) y se mide el
    residuo relativo. Un residuo pequeño revela una casi dependencia lineal
    oculta: la matriz no es singular, pero está mal condicionada.
    """
    n = len(A)
    resultados: list[dict] = []
    for k in range(n):
        for i in range(n):
            if i == k:
                continue
            for j in range(i + 1, n):
                if j == k:
                    continue
                ri, rj, rk = A[i], A[j], A[k]
                g11 = sum(valor * valor for valor in ri)
                g12 = sum(ri[t] * rj[t] for t in range(n))
                g22 = sum(valor * valor for valor in rj)
                c1 = sum(ri[t] * rk[t] for t in range(n))
                c2 = sum(rj[t] * rk[t] for t in range(n))
                determinante = g11 * g22 - g12 * g12
                if determinante == 0.0:
                    continue
                alfa = (c1 * g22 - c2 * g12) / determinante
                beta = (g11 * c2 - g12 * c1) / determinante
                norma = math.sqrt(sum(valor * valor for valor in rk))
                if norma == 0.0:
                    continue
                error = math.sqrt(
                    sum(
                        (rk[t] - alfa * ri[t] - beta * rj[t]) ** 2 for t in range(n)
                    )
                ) / norma
                resultados.append(
                    {
                        "fila_objetivo": k + 1,
                        "fila_i": i + 1,
                        "fila_j": j + 1,
                        "alfa": alfa,
                        "beta": beta,
                        "error_relativo": error,
                    }
                )
    resultados.sort(key=lambda fila: fila["error_relativo"])
    return resultados[:10]


def detectar_solucion_exacta(
    A_cad: list[list[str]],
    b_cad: list[str],
    x_referencia: list[Decimal],
    decimales: tuple[int, ...] = tuple(range(4, 17)),
    digitos: int = 80,
) -> tuple[list[Decimal], int] | None:
    """Busca un redondeo de la referencia que reproduce `b` de forma exacta.

    `b` fue construido con aritmética racional exacta a partir de la
    solución verdadera; en consecuencia existe un decimal finito `x*` tal que
    `A x* = b` sin error residual en aritmética decimal. Redondeando la
    solución de referencia a distintos números de decimales se puede
    encontrar ese `x*` y certificarlo con un residuo exactamente nulo.
    """
    with localcontext() as contexto:
        contexto.prec = digitos
        A = a_decimal(A_cad, digitos)
        b = b_decimal(b_cad, digitos)
        n = len(b)
        for decimales_candidatos in decimales:
            x = [
                valor.quantize(Decimal(1).scaleb(-decimales_candidatos))
                for valor in x_referencia
            ]
            residuo = max(
                abs(b[i] - sum(A[i][j] * x[j] for j in range(n))) for i in range(n)
            )
            if residuo == 0:
                return x, decimales_candidatos
    return None


def digitos_correctos(
    x: list[float], x_referencia: list[Decimal]
) -> list[dict]:
    """Estima los dígitos significativos correctos de cada incógnita.

    Para cada componente se calcula el error relativo frente a la solución de
    referencia en precisión extendida; el número de dígitos correctos es
    `floor(-log10(error relativo))`, redondeado hacia abajo, que es la
    definición usual de cifras significativas correctas.
    """
    filas: list[dict] = []
    for i, (valor, referencia) in enumerate(zip(x, x_referencia), start=1):
        diferencia = abs(Decimal(valor) - referencia)
        error_relativo = diferencia / abs(referencia) if referencia != 0 else diferencia
        if error_relativo == 0:
            digitos = 17
        else:
            digitos = int(math.floor(-math.log10(float(error_relativo))))
        filas.append(
            {
                "indice": i,
                "valor": valor,
                "referencia": referencia,
                "error_absoluto": float(diferencia),
                "error_relativo": float(error_relativo),
                "digitos": digitos,
            }
        )
    return filas
