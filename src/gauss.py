"""Eliminación gaussiana instrumentada para el sistema denso de 12 incógnitas.

Este módulo implementa la eliminación gaussiana en Python puro (listas de
flotantes). No se usa ninguna rutina de biblioteca para resolver el sistema:
numpy y scipy aparecen únicamente como *control* en `verificacion.py`.

Características del algoritmo:

* Una **única factorización** LU por corrida; el sistema nunca se resuelve
  calculando la inversa de A.
* El multiplicador `m(i,k) = a(i,k)/a(k,k)` se calcula **una sola vez por
  fila** y se reutiliza en todas las actualizaciones de esa fila
  (columnas `k+1..n` y el lado derecho `b`).
* El barrido del paso `k` recorre las columnas desde `k` en adelante: la
  subdiagonal ya anulada no se vuelve a tocar.
* Búsqueda de pivote en tres variantes: `"ninguno"`, `"parcial"` y
  `"escalado"` (pivoteo parcial escalado por filas).
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field

from conteo import Contador

MODOS_PIVOTE = ("ninguno", "parcial", "escalado", "escalado_fijo")


# ---------------------------------------------------------------------- #
# Lectura de datos
# ---------------------------------------------------------------------- #
def leer_matriz(ruta) -> list[list[str]]:
    """Lee `datos/matriz_a.csv` y devuelve las 144 cadenas sin modificar."""
    filas: list[list[str]] = []
    with open(ruta, "r", encoding="utf-8", newline="") as archivo:
        for registro in csv.reader(archivo):
            valores = [v.strip() for v in registro if v.strip() != ""]
            if valores:
                filas.append(valores)
    if not filas or len(filas) != len(filas[0]):
        raise ValueError("matriz_a.csv debe ser una matriz cuadrada de 12x12")
    return filas


def leer_vector(ruta) -> list[str]:
    """Lee `datos/vector_b.csv` y devuelve las 12 cadenas sin modificar."""
    valores: list[str] = []
    with open(ruta, "r", encoding="utf-8", newline="") as archivo:
        for registro in csv.reader(archivo):
            for campo in registro:
                campo = campo.strip()
                if campo != "":
                    valores.append(campo)
    return valores


def a_flotante(cadenas: list[list[str]]) -> list[list[float]]:
    """Convierte la matriz de cadenas a flotantes de doble precisión."""
    return [[float(valor) for valor in fila] for fila in cadenas]


def b_flotante(cadenas: list[str]) -> list[float]:
    """Convierte el vector de cadenas a flotantes de doble precisión."""
    return [float(valor) for valor in cadenas]


# ---------------------------------------------------------------------- #
# Resultado de una corrida
# ---------------------------------------------------------------------- #
@dataclass
class Resultado:
    """Todo lo que produce una corrida de eliminación."""

    x: list[float]
    u: list[list[float]]
    l: list[list[float]]
    perm: list[int]
    y: list[float]
    multiplicadores: list[tuple[int, int, float]]
    max_multiplicador: float
    max_abs_inicial: float
    max_abs_durante: float
    crecimiento: float
    traza: list[list[list[float]]]
    n_intercambios: int
    contador: Contador
    pivote: str
    recalcular_multiplicadores: bool = False
    bits_cancelacion_promedio: float = 0.0
    bits_cancelacion_maximo: float = 0.0
    actualizaciones: int = 0

    @property
    def factor_crecimiento(self) -> float:
        """Crecimiento de los elementos respecto del máximo inicial."""
        return self.crecimiento


@dataclass
class _Pivote:
    """Resultado interno de la búsqueda de pivote de un paso."""

    fila: int
    comparaciones: int = 0
    intercambiada: bool = False
    valor: float = 0.0
    escalas: list[float] = field(default_factory=list)


# ---------------------------------------------------------------------- #
# Núcleo del algoritmo
# ---------------------------------------------------------------------- #
def eliminar(
    A: list[list[float]],
    b: list[float],
    pivote: str = "escalado",
    recalcular_multiplicadores: bool = False,
    contador: Contador | None = None,
) -> Resultado:
    """Triangulariza A con eliminación gaussiana y resuelve el sistema.

    Parámetros
    ----------
    A : matriz cuadrada densa (lista de listas de flotantes).
    b : lado derecho.
    pivote : `"ninguno"`, `"parcial"` o `"escalado"`.
    recalcular_multiplicadores : si es `True`, el multiplicador se recalcula
        en cada columna en lugar de reutilizarse (solo para demostrar el costo
        del trabajo redundante; el resultado numérico es idéntico).
    contador : contador de operaciones a instrumentar.

    Devuelve
    -------
    Resultado con la solución, los factores L y U, la traza por paso y las
    estadísticas de estabilidad.
    """
    if pivote not in MODOS_PIVOTE:
        raise ValueError(f"modo de pivote desconocido: {pivote}")

    contador = contador if contador is not None else Contador()
    n = len(A)
    M = [fila[:] for fila in A]
    rhs = b[:]
    L = [[0.0] * n for _ in range(n)]
    perm = list(range(n))
    multiplicadores: list[tuple[int, int, float]] = []
    traza: list[list[list[float]]] = []

    max_abs_inicial = max(abs(valor) for fila in M for valor in fila)
    max_abs_durante = max_abs_inicial
    max_multiplicador = 0.0
    n_intercambios = 0
    cancelacion = _Cancelacion()
    escalas_fijas = [max(abs(valor) for valor in fila) for fila in A]

    for k in range(n - 1):
        paso = contador.iniciar_paso(k + 1)
        elegido = _buscar_pivote(M, k, n, pivote, contador, escalas_fijas)

        if elegido.fila != k:
            M[k], M[elegido.fila] = M[elegido.fila], M[k]
            rhs[k], rhs[elegido.fila] = rhs[elegido.fila], rhs[k]
            for j in range(k):
                L[k][j], L[elegido.fila][j] = L[elegido.fila][j], L[k][j]
            perm[k], perm[elegido.fila] = perm[elegido.fila], perm[k]
            contador.intercambiar_filas(1)
            n_intercambios += 1

        paso.pivote = k
        L[k][k] = 1.0

        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]
            contador.dividir()
            L[i][k] = m
            multiplicadores.append((i, k, m))
            if abs(m) > max_multiplicador:
                max_multiplicador = abs(m)
            if abs(m) > paso.multiplicador_maximo:
                paso.multiplicador_maximo = abs(m)

            if recalcular_multiplicadores:
                # Trabajo redundante: el multiplicador se recalcula por
                # columna. Los operandos son los mismos, luego el resultado
                # numérico no cambia; solo crece el costo.
                for j in range(k + 1, n):
                    m_j = M[i][k] / M[k][k]
                    contador.dividir()
                    antes = M[i][j]
                    M[i][j] = antes - m_j * M[k][j]
                    contador.multiplicar()
                    contador.sumar_restar()
                    _registrar_cancelacion(antes, m_j * M[k][j], M[i][j], cancelacion)
                m_b = M[i][k] / M[k][k]
                contador.dividir()
                antes = rhs[i]
                rhs[i] = antes - m_b * rhs[k]
                contador.multiplicar()
                contador.sumar_restar()
                _registrar_cancelacion(antes, m_b * rhs[k], rhs[i], cancelacion)
            else:
                for j in range(k + 1, n):
                    antes = M[i][j]
                    M[i][j] = antes - m * M[k][j]
                    contador.multiplicar()
                    contador.sumar_restar()
                    _registrar_cancelacion(antes, m * M[k][j], M[i][j], cancelacion)
                antes = rhs[i]
                rhs[i] = antes - m * rhs[k]
                contador.multiplicar()
                contador.sumar_restar()
                _registrar_cancelacion(antes, m * rhs[k], rhs[i], cancelacion)

            M[i][k] = 0.0
            contador.asignar()

        traza.append([fila[:] for fila in M])
        max_actual = max(abs(valor) for fila in M for valor in fila)
        if max_actual > max_abs_durante:
            max_abs_durante = max_actual
        contador.cerrar_paso()

    # La última fila no elimina nada: sólo aporta la unidad de la diagonal
    # inferior (no se cuenta como paso de eliminación).
    L[n - 1][n - 1] = 1.0

    # Sustitución regresiva.
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        paso_s = contador.iniciar_paso_sustitucion(i + 1)
        acumulado = rhs[i]
        for j in range(i + 1, n):
            acumulado = acumulado - M[i][j] * x[j]
            contador.multiplicar()
            contador.sumar_restar()
        x[i] = acumulado / M[i][i]
        contador.dividir()
        contador.cerrar_paso()

    return Resultado(
        x=x,
        u=M,
        l=L,
        perm=perm,
        y=rhs,
        multiplicadores=multiplicadores,
        max_multiplicador=max_multiplicador,
        max_abs_inicial=max_abs_inicial,
        max_abs_durante=max_abs_durante,
        crecimiento=max_abs_durante / max_abs_inicial,
        traza=traza,
        n_intercambios=n_intercambios,
        contador=contador,
        pivote=pivote,
        recalcular_multiplicadores=recalcular_multiplicadores,
        bits_cancelacion_promedio=cancelacion.promedio,
        bits_cancelacion_maximo=cancelacion.maximo,
        actualizaciones=cancelacion.total,
    )


class _Cancelacion:
    """Instrumentación del error de redondeo por cancelación catastrófica.

    En cada actualización `nuevo = antes - producto` los bits de mantisa
    que se pierden se estiman con
    `log2((|antes| + |producto|) / |nuevo|)`: es la amplificación del error
    de redondeo relativa debida a la resta. Es instrumentación de análisis y
    no forma parte del conteo de flops.
    """

    def __init__(self) -> None:
        self.suma = 0.0
        self.maximo = 0.0
        self.total = 0

    def registrar(self, antes: float, producto: float, nuevo: float) -> None:
        """Acumula los bits perdidos por la cancelación de una resta."""
        if nuevo == 0.0:
            bits = 53.0
        else:
            bits = math.log2((abs(antes) + abs(producto)) / abs(nuevo))
            if bits < 0.0:
                bits = 0.0
        self.suma += bits
        self.total += 1
        if bits > self.maximo:
            self.maximo = bits

    @property
    def promedio(self) -> float:
        """Bits perdidos en promedio por actualización."""
        return self.suma / self.total if self.total else 0.0


def _registrar_cancelacion(
    antes: float, producto: float, nuevo: float, cancelacion: _Cancelacion
) -> None:
    """Envoltorio de `._Cancelacion.registrar` para las actualizaciones."""
    cancelacion.registrar(antes, producto, nuevo)


def _buscar_pivote(
    M: list[list[float]],
    k: int,
    n: int,
    modo: str,
    contador: Contador,
    escalas_fijas: list[float] | None = None,
) -> _Pivote:
    """Busca el pivote del paso `k` según el modo pedido."""
    if modo == "ninguno":
        return _Pivote(fila=k)

    if modo == "parcial":
        mejor = abs(M[k][k])
        elegido = k
        comparaciones = 0
        for i in range(k + 1, n):
            valor = abs(M[i][k])
            comparaciones += 1
            if valor > mejor:
                mejor = valor
                elegido = i
        contador.comparar(comparaciones)
        return _Pivote(fila=elegido, comparaciones=comparaciones, valor=mejor)

    if modo == "escalado_fijo":
        # Escalas de fila calculadas una sola vez sobre la matriz original:
        # el pivoteo escalado "clásico" (LAPACK), más barato en comparaciones.
        escalas = escalas_fijas or [max(abs(M[i][j]) for j in range(n)) for i in range(n)]
        escala_inicial = escalas[k]
        mejor = abs(M[k][k]) / escala_inicial if escala_inicial else 0.0
        elegido = k
        comparaciones = 0
        for i in range(k + 1, n):
            escala = escalas[i]
            valor = abs(M[i][k]) / escala if escala else 0.0
            comparaciones += 1
            if valor > mejor:
                mejor = valor
                elegido = i
        if k == 0:
            # El cálculo de las escalas fijas (máximo absoluto de cada fila)
            # se hace una sola vez: n filas con n-1 comparaciones cada una.
            comparaciones += n * (n - 1)
        contador.comparar(comparaciones)
        return _Pivote(fila=elegido, comparaciones=comparaciones, valor=mejor, escalas=escalas)

    # Pivoteo parcial escalado: se normaliza cada fila candidata por su
    # escala (máximo absoluto de la fila en el bloque activo) y se elige el
    # mayor cociente |a(i,k)| / escala_i.
    escalas = [max(abs(M[i][j]) for j in range(k, n)) for i in range(k, n)]
    escala_inicial = escalas[0]
    mejor = abs(M[k][k]) / escala_inicial if escala_inicial else 0.0
    elegido = k
    comparaciones = 0
    for i in range(k + 1, n):
        escala = escalas[i - k]
        valor = abs(M[i][k]) / escala if escala else 0.0
        comparaciones += 1
        if valor > mejor:
            mejor = valor
            elegido = i
    # El cálculo de las escalas es una reducción por máximo: se contabilizan
    # como comparaciones de apoyo, no como flops.
    comparaciones += sum(max(n - k - 1, 0) for _ in range(n - k))
    contador.comparar(comparaciones)
    return _Pivote(fila=elegido, comparaciones=comparaciones, valor=mejor, escalas=escalas)


def resolver_con_factores(
    L: list[list[float]],
    U: list[list[float]],
    perm: list[int],
    b: list[float],
) -> list[float]:
    """Resuelve A x = b reutilizando una factorización LU ya calculada.

    Se usa en el refinamiento iterativo: no se vuelve a triangularizar el
    sistema, de modo que el costo de esta resolución no altera el conteo de
    la corrida principal.
    """
    n = len(L)
    y = [b[perm[k]] for k in range(n)]
    for i in range(1, n):
        acumulado = y[i]
        for j in range(i):
            acumulado -= L[i][j] * y[j]
        y[i] = acumulado
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        acumulado = y[i]
        for j in range(i + 1, n):
            acumulado -= U[i][j] * x[j]
        x[i] = acumulado / U[i][i]
    return x
