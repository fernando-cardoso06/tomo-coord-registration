# Registro de coordenadas entre aquisições tomográficas

Parte do meu projeto de iniciação científica no IFUSP (Uso de Aprendizado
de Máquina em Áreas Interdisciplinares Relacionadas à Saúde), este repositório 
isola o pipeline de carregamento de imagens DICOM e de registro de coordenadas
usado para localizar nódulos de um phantom de calibração de forma consistente
entre diferentes aquisições tomográficas.

## O problema

Um phantom de calibração com nódulos esféricos em posições conhecidas
é escaneado repetidas vezes, sendo girado em seu suporte entre algumas das
aquisições. Cada scan novo tem sua própria origem espacial e também,
devido à rotação, um deslocamento angular em relação às posições já
conhecidas dos nódulos. O metadado de posicionamento do equipamento nem
sempre é confiável o suficiente para recuperar esse deslocamento sozinho.

Este pipeline resolve isso com uma calibração leve: dois pontos marcados
manualmente numa nova aquisição (o centro do phantom e um ponto de
referência angular) são suficientes para estimar a rotação e a translação
em Z entre exames, e então levar as posições conhecidas dos nódulos para o
sistema de coordenadas da nova aquisição.

## Estrutura

```
src/
  dicom_io.py         # carregamento paralelo de séries DICOM em volume 3D
  coord_registro.py   # transformação e registro de coordenadas entre exames
  dados_sinteticos.py # geração de volume sintético para demonstração pública
notebooks/
  demo_registro_coordenadas.ipynb  # pipeline completo, executável, sobre dado sintético
```

## Sobre os dados

As tomografias reais deste projeto ficam em um cluster institucional com
acesso restrito e não podem ser publicadas aqui. Para que o
pipeline seja executável e revisável publicamente, o notebook de
demonstração roda sobre um volume sintético: um phantom cilíndrico
simples com nódulos esféricos e ruído, gerado em memória por
`dados_sinteticos.py`. A lógica de carregamento e registro é exatamente a
mesma usada com dados reais; apenas a origem do volume muda.

## Resultado no pipeline real

Otimizando o carregamento DICOM com leitura paralela de headers e pixels
(com uso da função `dicom_io.carregar_serie_dicom`), o tempo de processamento
de um exame de referência com 2.915 imagens caiu de ~170s para ~105s, o que 
expressa uma redução de ~38% em comparação com uma versão sequencial
do mesmo carregamento.

## Rodando a demonstração

```bash
pip install -r requirements.txt
jupyter notebook notebooks/demo_registro_coordenadas.ipynb
```

O notebook gera o volume sintético, simula uma nova aquisição rotacionada,
aplica o registro de coordenadas e mostra visualmente a posição dos
nódulos antes e depois, sem depender de nenhum dado externo.
