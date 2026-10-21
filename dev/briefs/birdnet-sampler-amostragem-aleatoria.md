# Brief: modos de amostragem aleatória no birdnet-validation-sampler

Para uma sessão do Claude Code aberta em
`~/GitHub/birdnet-validation-sampler` (repositório
https://github.com/diogoprov/birdnet-validation-sampler, licença MIT,
DOI 10.5281/zenodo.22696210). Converse com o Diogo em português; código e
comentários seguem o estilo do repositório (português, tidyverse, `|>`).
Antes de mexer no código, invoque as skills
`/andrej-karpathy-skills:karpathy-guidelines` e, antes do commit,
`/mattpocock-skills:code-review`.

## Problema

Hoje `extrair_clipes_validacao.R` escolhe os clipes pelas *N* predições de
maior confidence por local × espécie × data/noite × hora (bloco
`group_by(...) |> arrange(desc(confidence)) |> slice_head(n = config$n_top)`,
por volta da linha 375), depois aplica `meta_total` ou `max_por_especie`.
Isso serve bem a dois usos atuais e **deve continuar como está e como
padrão**:

- retreino do classificador (clipes positivos);
- calibração score → P(verdadeiro positivo) no Estágio 2 (`calibracao.R`):
  como o modelo é condicional ao score, selecionar pelo score não enviesa a
  curva de calibração.

Mas não serve para estimar a **taxa de falsos positivos** em modelos de
ocupação com falsos positivos, como os do pacote R crywolf
(https://provetelab.org/crywolf/). Esses modelos supõem que as ocorrências
checadas foram escolhidas **independentemente do score**. Clipes com os
maiores scores têm menos falsos positivos do que uma amostra aleatória, e a
taxa estimada sai baixa demais. O README do crywolf já avisa isso e diz para
"draw the clips at random"; falta o sampler oferecer essa opção.

## O que implementar

Um parâmetro novo em `config`, por exemplo
`modo_amostragem = c("top", "aleatorio_deteccoes", "aleatorio_gravacoes")`,
com `"top"` como padrão (comportamento atual idêntico, inclusive os mesmos
clipes com a mesma semente).

### 1. `"aleatorio_deteccoes"` (desenho de confirmação de observações, Chambert et al. 2015)

- Unidade de sorteio: a **ocasião**, não a predição. Ocasião = a unidade que
  vai entrar no modelo de ocupação como uma ocasião automatizada; no
  crywolf costuma ser a gravação (arquivo de áudio) ou o minuto. Deixe
  configurável (`unidade_ocasiao = "arquivo"` ou `"minuto"`, este último a
  partir de `file_offset`).
- Candidatas: ocasiões com **pelo menos uma predição da espécie com
  confidence ≥ `limiar`**, onde `limiar` é o mesmo limiar que será usado na
  análise (parâmetro novo; não confundir com `conf_minima`, que hoje só
  filtra candidatas no modo `"top"`).
- Sorteio aleatório simples das ocasiões candidatas, com semente fixa,
  dentro de estratos que **não dependam do score** (local, espécie, noite,
  hora são aceitáveis; `meta_total` e `balancear_por_hora` podem continuar
  valendo, desde que a escolha dentro do estrato seja aleatória e não por
  confidence).
- Para cada ocasião sorteada, cortar **todas** as predições da espécie com
  confidence ≥ `limiar` daquela ocasião, não só uma. O desenho supõe que
  todas as observações de uma ocasião checada foram checadas; o crywolf
  (`fp_aggregate()`) avisa quando uma ocasião tem mais detecções do que
  clipes revisados.

### 2. `"aleatorio_gravacoes"` (desenho de múltiplos métodos, Miller et al. 2011)

- O método verificado precisa ser independente do classificador: sorteio
  aleatório de **gravações (ou minutos) entre todas as gravadas**, com ou
  sem predição, para anotação por ouvido. É o que foi feito com o conjunto
  de teste do artigo do curiango (≈1% das gravações, proporcional ao esforço
  de cada local).
- Isso exige a lista de todos os arquivos de áudio, não só as selection
  tables (que listam apenas predições). Use `list.files()` sobre os
  diretórios de áudio (já há `remap_de`/`remap_para` para os caminhos) ou
  aceite uma tabela de gravações.
- O clipe é a gravação (ou o minuto) inteira, não uma janela de 3 s.
- Fração ou número por local configurável (por exemplo
  `fracao_gravacoes = 0.01`, proporcional ao esforço por padrão).

### Manifesto e pós-validação

Mantenha os nomes de coluna atuais do `manifesto_clipes.csv` e do
`validacao_por_clipe.csv`, porque `crywolf::fp_read_validation()` lê
`local`, `species_code`, `audio_origem`, `inicio_clipe`, `fim_clipe`,
`confidence` e `veredito` ("positivo"/"negativo"). Acrescente:

- `modo_amostragem`;
- `ocasiao_id` (arquivo, ou arquivo + minuto);
- `n_deteccoes_ocasiao` (predições da espécie ≥ `limiar` na ocasião), para
  conferir que todas foram cortadas;
- `prob_inclusao` (probabilidade de a ocasião ser sorteada no estrato),
  útil se um dia for preciso ponderar.

No modo `"aleatorio_gravacoes"`, `confidence` pode ficar `NA` e
`species_code` deve indicar a espécie a anotar.

`pos_validacao.R` calcula precisão por faixa de confidence e o limiar para
uma precisão-alvo. Nos modos aleatórios, essas tabelas continuam válidas
(a amostra aleatória dá uma estimativa não enviesada da precisão no
limiar), mas confira que nada ali presume o `rank` do modo `"top"`.

## Como verificar

1. Com `modo_amostragem = "top"` e a mesma semente, o manifesto sai
   **idêntico** ao atual (compare com um manifesto gerado antes da mudança,
   em `simular = TRUE`).
2. Com `"aleatorio_deteccoes"`, a distribuição de confidence dos clipes
   sorteados deve ser parecida com a de todas as predições ≥ `limiar`
   (compare quantis); com `"top"`, ela fica concentrada no topo.
3. Toda ocasião sorteada tem todos os seus clipes no manifesto
   (`n_deteccoes_ocasiao` = clipes da ocasião).
4. Rode em `simular = TRUE` num local pequeno antes de cortar áudio.
5. Atualize o README (tabela de parâmetros e uma seção curta explicando
   quando usar cada modo, com a ressalva de que `"top"` não serve para
   estimar taxas de falsos positivos) e o CITATION.cff/NEWS se houver
   release. Referências citáveis: Miller et al. 2011, Ecology 92:1422–1428,
   doi:10.1890/10-1396.1; Chambert et al. 2015, Ecology 96:332–339,
   doi:10.1890/14-1507.1. Não invente outras.

## Fora do escopo

- Não mude o comportamento padrão nem os nomes de arquivos de saída.
- Não implemente o ajuste do modelo de ocupação: isso é do crywolf
  (`fp_read_validation()` + `fp_aggregate()` + `fp_data()` + `fp_occu()`).
