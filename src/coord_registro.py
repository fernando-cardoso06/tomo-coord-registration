"""Registro de coordenadas de nódulos entre aquisições tomográficas.

Problema: um conjunto de nódulos foi localizado manualmente (em pixel/fatia)
num exame de referência (ex: um phantom fixo, escaneado repetidas vezes).
Cada nova aquisição tem sua própria origem, espaçamento e — no caso de um
objeto rotacionado entre escaneamentos — um ângulo de rotação diferente em
torno do eixo Z. O objetivo é levar as posições conhecidas do exame de
referência para o sistema de coordenadas de uma nova aquisição.

A estratégia adotada aqui não usa registro de imagem completo (ex: ICP ou
otimização de mútua informação); em vez disso, usa dois pontos de calibração
indicados manualmente na nova aquisição (o centro do objeto e um ponto de
referência angular) para estimar a rotação e a translação em Z necessárias.
É uma solução mais simples e mais barata computacionalmente, adequada
quando o objeto é rígido e sofre apenas rotação em torno de Z + deslocamento
em Z entre aquisições — o caso de um phantom de calibração fixo em um
suporte giratório.
"""

from dataclasses import dataclass, field

import numpy as np


@dataclass
class ReferenciaCalibracao:
    """Parâmetros do exame de referência usados para converter pixel/fatia
    em coordenadas físicas (mm), e para servir de base à comparação angular
    com novas aquisições.
    """

    x_ref: float
    y_ref: float
    s_ref: float          # fatia (slice) de referência para o centro
    z1_ref: float          # posição Z (mm) da última fatia da referência
    z0_ref: float          # posição Z (mm) da primeira fatia da referência
    spacing_ref: tuple      # (dx, dy, dz) em mm/pixel
    theta_ref: float        # ângulo de referência (rad) no exame original
    theta_z_ref: float      # posição Z (mm) da fatia usada para medir theta_ref

    def slice_para_z(self, fatia: float) -> float:
        """Converte índice de fatia em posição Z (mm), usando a referência."""
        return self.z1_ref - (self.s_ref - fatia) * self.spacing_ref[2]

    def pixel_para_mm(self, valor_em_pixels: float) -> float:
        """Converte uma distância em pixels para mm no plano XY."""
        return valor_em_pixels * self.spacing_ref[0]

    def cartesiano(self, x: float, y: float, s: float) -> np.ndarray:
        """Converte (x, y, fatia) em coordenadas cartesianas (mm) relativas
        ao centro do exame de referência."""
        dx = self.pixel_para_mm(x - self.x_ref)
        # y cresce no sentido contrário ao convencional em matplotlib
        dy = self.pixel_para_mm(y + self.y_ref)
        z = self.slice_para_z(s)
        return np.array([dx, dy, z])

    def cilindrico(self, x: float, y: float, s: float) -> np.ndarray:
        """Converte (x, y, fatia) em coordenadas cilíndricas (rho, theta, Z)
        relativas ao centro do exame de referência."""
        dx, dy, z = self.cartesiano(x, y, s)
        rho = np.hypot(dx, dy)
        theta = np.arctan2(dy, dx)
        return np.array([rho, theta, z])


def pontos_para_cilindricas(pontos_pixel_fatia: np.ndarray, ref: ReferenciaCalibracao) -> np.ndarray:
    """Converte um array (N, 3) de pontos (x, y, fatia) em coordenadas
    cilíndricas (rho, theta, Z), usando os parâmetros da referência."""
    return np.array([ref.cilindrico(x, y, s) for x, y, s in pontos_pixel_fatia])


def registrar_em_nova_aquisicao(
    pontos_cilindricos_ref: np.ndarray,
    angulo_medido: float,
    delta_z: float,
    theta_ref: float,
    closest_z_fn,
) -> np.ndarray:
    """Aplica a rotação (d_theta) e a translação em Z (delta_z) estimadas por
    calibração manual, levando pontos conhecidos do exame de referência para
    o sistema de coordenadas de uma nova aquisição.

    Parameters
    ----------
    pontos_cilindricos_ref : np.ndarray
        Array (N, 3) com (rho, theta, Z) dos pontos no exame de referência.
    angulo_medido : float
        Ângulo (rad) medido na nova aquisição a partir de dois cliques de
        calibração (centro do objeto + ponto de referência angular).
    delta_z : float
        Deslocamento em Z (mm) entre a referência e a nova aquisição,
        estimado a partir dos limites conhecidos do volume.
    theta_ref : float
        Ângulo de referência (rad) no exame original, usado como base da
        rotação relativa.
    closest_z_fn : Callable[[float], float]
        Função que ajusta uma posição Z (mm) para a fatia mais próxima
        disponível na nova aquisição.

    Returns
    -------
    np.ndarray
        Array (N, 3) com (rho, theta, Z) dos pontos já registrados no
        sistema de coordenadas da nova aquisição.
    """
    d_theta = angulo_medido - theta_ref
    return np.array([
        [rho, theta + d_theta, closest_z_fn(z + delta_z)]
        for rho, theta, z in pontos_cilindricos_ref
    ])
