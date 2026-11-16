"""Instrumentación del costo operacional de la eliminación gaussiana.

Convención de conteo declarada (se repite en `informe.md` y en
`resultados/conteo_flops.md`):

* **1 flop = 1 operación aritmética elemental**: multiplicación, división,
  suma o resta entre dos números en punto flotante (doble precisión IEEE 754).
* Las **asignaciones** (en particular la anulación explícita de la subdiagonal
  con `a[i][k] = 0.0`) **no** son operaciones aritméticas: se contabilizan
  aparte para no inflar artificialmente el conteo.
* Las **comparaciones** de la búsqueda de pivote y los **intercambios** de fila
  tampoco son flops: se contabilizan aparte como *sobrecosto de estabilidad*
  exigido por la pauta.
* La resta `a - b` cuenta como 1 suma/resta; la actualización
  `a[i][j] -= m * a[k][j]` cuenta como 1 multiplicación + 1 resta.

El contador no usa ninguna fórmula escrita a mano: cada operación del
algoritmo llama explícitamente a un método del contador, de modo que el total
publicado proviene del código y no de un cálculo teórico.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Paso:
    """Detalle de las operaciones realizadas en un paso k de la eliminación."""

    paso: int
    multiplications: int = 0
    divisiones: int = 0
    sumas_restas: int = 0
    asignaciones: int = 0
    comparaciones: int = 0
    intercambios: int = 0
    filas_intercambiadas: int = 0
    pivote: int = -1
    multiplicador_maximo: float = 0.0

    @property
    def flops(self) -> int:
        """Operaciones aritméticas elementales del paso."""
        return self.multiplications + self.divisiones + self.sumas_restas

    @property
    def total_con_apoyo(self) -> int:
        """Flops del paso más el trabajo de apoyo (asignaciones, comparaciones)."""
        return self.flops + self.asignaciones + self.comparaciones


@dataclass
class Contador:
    """Contador de operaciones de la eliminación gaussiana."""

    pasos: list[Paso] = field(default_factory=list)
    pasos_sustitucion: list[Paso] = field(default_factory=list)
    _paso: Paso | None = None

    # ------------------------------------------------------------------ #
    # Control del paso en curso
    # ------------------------------------------------------------------ #
    def iniciar_paso(self, numero: int) -> Paso:
        """Abre un nuevo paso de eliminación y lo devuelve."""
        paso = Paso(paso=numero)
        self.pasos.append(paso)
        self._paso = paso
        return paso

    def iniciar_paso_sustitucion(self, numero: int) -> Paso:
        """Abre un nuevo paso de la sustitución regresiva."""
        paso = Paso(paso=numero)
        self.pasos_sustitucion.append(paso)
        self._paso = paso
        return paso

    def cerrar_paso(self) -> None:
        """Cierra el paso en curso."""
        self._paso = None

    # ------------------------------------------------------------------ #
    # Operaciones elementales
    # ------------------------------------------------------------------ #
    def multiplicar(self, veces: int = 1) -> None:
        """Registra `veces` multiplicaciones."""
        if self._paso is not None:
            self._paso.multiplications += veces

    def dividir(self, veces: int = 1) -> None:
        """Registra `veces` divisiones."""
        if self._paso is not None:
            self._paso.divisiones += veces

    def sumar_restar(self, veces: int = 1) -> None:
        """Registra `veces` sumas o restas."""
        if self._paso is not None:
            self._paso.sumas_restas += veces

    def asignar(self, veces: int = 1) -> None:
        """Registra `veces` asignaciones (no son flops)."""
        if self._paso is not None:
            self._paso.asignaciones += veces

    def comparar(self, veces: int = 1) -> None:
        """Registra `veces` comparaciones de la búsqueda de pivote."""
        if self._paso is not None:
            self._paso.comparaciones += veces

    def intercambiar_filas(self, filas: int = 1) -> None:
        """Registra el intercambio de `filas` filas (conteo de estabilidad)."""
        if self._paso is not None:
            self._paso.intercambios += 1
            self._paso.filas_intercambiadas += filas

    # ------------------------------------------------------------------ #
    # Totales
    # ------------------------------------------------------------------ #
    @property
    def eliminacion(self) -> Paso:
        """Totales de la fase de eliminación."""
        return self._sumar(self.pasos)

    @property
    def sustitucion(self) -> Paso:
        """Totales de la fase de sustitución regresiva."""
        return self._sumar(self.pasos_sustitucion)

    @property
    def flops_eliminacion(self) -> int:
        """Flops aritméticos de la eliminación."""
        return self.eliminacion.flops

    @property
    def flops_sustitucion(self) -> int:
        """Flops aritméticos de la sustitución regresiva."""
        return self.sustitucion.flops

    @property
    def flops_totales(self) -> int:
        """Flops aritméticos de eliminación más sustitución."""
        return self.flops_eliminacion + self.flops_sustitucion

    @property
    def asignaciones_totales(self) -> int:
        """Asignaciones (anulaciones explícitas de la subdiagonal)."""
        return self.eliminacion.asignaciones + self.sustitucion.asignaciones

    @property
    def comparaciones_totales(self) -> int:
        """Comparaciones hechas en las búsquedas de pivote."""
        return self.eliminacion.comparaciones

    @property
    def intercambios_totales(self) -> int:
        """Intercambios de fila realizados."""
        return self.eliminacion.intercambios

    @property
    def costo_total_con_estabilidad(self) -> int:
        """Flops aritméticos más el trabajo de apoyo del pivoteo."""
        return (
            self.flops_totales
            + self.asignaciones_totales
            + self.comparaciones_totales
        )

    @staticmethod
    def _sumar(pasos: list[Paso]) -> Paso:
        total = Paso(paso=-1)
        for paso in pasos:
            total.multiplications += paso.multiplications
            total.divisiones += paso.divisiones
            total.sumas_restas += paso.sumas_restas
            total.asignaciones += paso.asignaciones
            total.comparaciones += paso.comparaciones
            total.intercambios += paso.intercambios
            total.filas_intercambiadas += paso.filas_intercambiadas
        return total
