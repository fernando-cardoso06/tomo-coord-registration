"""Carregamento de séries DICOM em volume 3D.

Lê os headers de todos os arquivos de um diretório em paralelo (sem os
pixels, que são mais pesados), ordena as fatias pela posição Z real no
espaço do paciente e só então lê os pixels — também em paralelo — na ordem
correta. Essa separação em duas fases é o que permite paralelizar a leitura
com segurança: a ordem final depende de metadado (ImagePositionPatient), não
da ordem de conclusão das threads.
"""

import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pydicom


def carregar_serie_dicom(caminho_diretorio: str, max_workers: int = 8):
    """Carrega uma série DICOM de um diretório e retorna um volume 3D em HU.

    Parameters
    ----------
    caminho_diretorio : str
        Diretório contendo os arquivos DICOM de uma única série/exame.
    max_workers : int
        Número de threads usadas na leitura paralela de headers e pixels.

    Returns
    -------
    volume_3d : np.ndarray
        Volume 3D (fatias, altura, largura) já convertido para unidades
        Hounsfield (HU), ordenado pela posição Z real no espaço do paciente.
    referencia : pydicom.Dataset
        Header DICOM da primeira fatia, usado como referência de metadados
        (spacing, orientação, origem) para o volume completo.
    z_positions : list[float]
        Posição Z (mm) de cada fatia, na mesma ordem do volume.

    Raises
    ------
    ValueError
        Se nenhum arquivo DICOM válido for encontrado no diretório.
    """

    def ler_header(caminho_arquivo):
        try:
            ds = pydicom.dcmread(caminho_arquivo, stop_before_pixels=True)
            if hasattr(ds, "ImagePositionPatient"):
                return (caminho_arquivo, ds)
        except Exception:
            pass
        return None

    candidatos = [
        os.path.join(caminho_diretorio, f)
        for f in os.listdir(caminho_diretorio)
        if os.path.isfile(os.path.join(caminho_diretorio, f))
        and (f.lower().endswith(".dcm") or "." not in f)
    ]

    # Fase 1: leitura paralela apenas dos headers (leve), para descobrir a
    # ordem real das fatias antes de pagar o custo de ler os pixels.
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        resultados_headers = list(filter(None, executor.map(ler_header, candidatos)))

    if not resultados_headers:
        raise ValueError("Nenhum arquivo DICOM válido encontrado.")

    # Ordenação pela posição Z real no espaço do paciente — o nome do
    # arquivo não é uma fonte confiável de ordem entre aquisições.
    resultados_headers.sort(key=lambda x: float(x[1].ImagePositionPatient[2]))

    caminhos_ordenados = [caminho for caminho, _ in resultados_headers]
    referencia = resultados_headers[0][1]
    z_positions = [float(header.ImagePositionPatient[2]) for _, header in resultados_headers]

    # Fase 2: leitura paralela dos pixels, já na ordem correta de cada fatia.
    def ler_pixels(args):
        indice, caminho = args
        ds = pydicom.dcmread(caminho)
        return indice, ds.pixel_array

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futuros = {
            executor.submit(ler_pixels, (i, c)): i
            for i, c in enumerate(caminhos_ordenados)
        }
        slices = [None] * len(caminhos_ordenados)
        for futuro in as_completed(futuros):
            indice, array = futuro.result()
            slices[indice] = array

    # Empilha as fatias e converte de valores brutos do sensor para
    # unidades Hounsfield usando o slope/intercept do header de referência.
    volume_3d = np.stack(slices).astype(np.float32)

    slope = float(getattr(referencia, "RescaleSlope", 1))
    intercept = float(getattr(referencia, "RescaleIntercept", 0))
    volume_3d *= slope
    volume_3d += intercept

    return volume_3d, referencia, z_positions
