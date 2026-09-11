"""Geração de um volume 3D sintético para demonstração pública do pipeline.

Os dados reais deste projeto são tomografias de um phantom de calibração,
armazenadas em um cluster institucional com acesso restrito (dados de
posicionamento de exame, sujeitos a política de acesso do departamento).
Para permitir que o pipeline seja executado e revisado publicamente sem
depender desse acesso, este módulo gera um volume sintético com estrutura
equivalente: um "phantom" simples (um cilindro com algumas esferas em
posições conhecidas, simulando nódulos) mais ruído gaussiano, salvo como
array numpy — sem passar por nenhum arquivo DICOM real.

Isso não é uma simulação fisicamente precisa de um exame de tomografia; é
uma estrutura mínima suficiente para exercitar o carregamento, a conversão
de escala e o registro de coordenadas do pipeline real.
"""

import numpy as np


def gerar_volume_sintetico(
    shape: tuple = (200, 256, 256),
    n_nodulos: int = 6,
    seed: int = 42,
) -> tuple:
    """Gera um volume 3D sintético simulando um phantom cilíndrico com
    nódulos esféricos, mais ruído — para uso em demonstração do pipeline.

    Parameters
    ----------
    shape : tuple
        Dimensões do volume (fatias, altura, largura).
    n_nodulos : int
        Número de esferas ("nódulos") a posicionar aleatoriamente dentro
        do cilindro do phantom.
    seed : int
        Semente do gerador aleatório, para reprodutibilidade.

    Returns
    -------
    volume : np.ndarray
        Volume 3D sintético em unidades arbitrárias (análogas a HU).
    posicoes_nodulos : np.ndarray
        Array (n_nodulos, 3) com as posições (x, y, fatia) de cada nódulo
        sintético em pixels/índice de fatia — equivalente ao que seria
        obtido por marcação manual num exame real.
    """
    rng = np.random.default_rng(seed)
    n_fatias, altura, largura = shape

    volume = np.full(shape, -1000.0, dtype=np.float32)  # fundo ~ ar em HU

    cy, cx = altura / 2, largura / 2
    raio_phantom = min(altura, largura) * 0.35

    yy, xx = np.mgrid[0:altura, 0:largura]
    dist_centro = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)
    mascara_phantom = dist_centro < raio_phantom

    for fatia in range(n_fatias):
        volume[fatia][mascara_phantom] = 0.0  # material do phantom ~ água em HU

    posicoes_nodulos = []
    for _ in range(n_nodulos):
        raio_n = rng.uniform(0, raio_phantom * 0.8)
        angulo_n = rng.uniform(0, 2 * np.pi)
        x = cx + raio_n * np.cos(angulo_n)
        y = cy + raio_n * np.sin(angulo_n)
        fatia = rng.uniform(n_fatias * 0.15, n_fatias * 0.85)
        raio_esfera = rng.uniform(3, 6)

        f0, f1 = max(0, int(fatia - raio_esfera)), min(n_fatias, int(fatia + raio_esfera))
        for f in range(f0, f1):
            dz = f - fatia
            r_max_fatia = np.sqrt(max(raio_esfera ** 2 - dz ** 2, 0))
            dist_nodulo = np.sqrt((yy - y) ** 2 + (xx - x) ** 2)
            volume[f][dist_nodulo < r_max_fatia] = 200.0  # nódulo ~ tecido denso

        posicoes_nodulos.append([x, y, fatia])

    volume += rng.normal(0, 15, size=shape).astype(np.float32)  # ruído do sensor

    return volume, np.array(posicoes_nodulos)
