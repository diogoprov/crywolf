# Proposta: calibração com covariável ambiental ou índice acústico (parte 3)

Escrita em 28/10/2026 pela sessão do crywolf, em resposta ao pedido da sessão
do sheepdog (partes 1 e 2 em diogoprov/sheepdog#55). **Aprovada pelo Diogo
em 28/10/2026**, com as ressalvas da sessão do sheepdog marcadas abaixo. O
lado do crywolf entra depois do artigo da MEE. Contexto: Larissa e Liliana sugeriram
sortear os clipes revisados ao longo de um gradiente ambiental e de um índice
acústico, a partir do resultado do fpcalib (p10 variando com a covariável de
ocupação).

## Resumo

1. Covariável contínua na calibração: o sheepdog calcula P(TP) em cada
   predição e sorteia as versões (`occasion_draws.parquet`, que já existe no
   1.2). O crywolf só lê. Esquema 1.3 com um único acréscimo:
   `calibration$covariates`. Com covariável contínua, o export **não** traz
   `calibration_draws.parquet`.
2. Modelo: análise de sensibilidade 2 × 2 (calibração com e sem x;
   `fp ~ 1` e `fp ~ x`), como a sazonal do Curiango. A recomendação de qual
   usar como principal espera um cenário do fpcalib.
3. Pesos 1/`inclusion_prob`: não implementar. Regra: o que estratifica o
   sorteio entra na calibração (faixa como chave ou x como covariável).

## 1. Esquema 1.3

### Por que o sheepdog calcula P(TP), e não o crywolf

Três caminhos possíveis:

- **Grade score × covariável** em `calibration_draws.parquet`: a grade em x
  teria de cobrir todos os valores observados (a covariável é da ocasião,
  com tantos valores quanto ocasiões), o que já é P(TP) por predição, com
  erro de interpolação a mais.
- **Exportar os coeficientes** para o crywolf calcular P(TP): o crywolf teria
  de remontar a matriz de delineamento da fórmula do brms (padronização,
  splines, níveis novos dos efeitos aleatórios, `new_levels`). Duplica o
  modelo do sheepdog, e cada fórmula nova quebra o contrato.
- **O sheepdog sorteia as versões** com `occasion_draws()`, calculando P(TP)
  de cada predição com o modelo que ele mesmo ajustou. O crywolf já lê esse
  arquivo desde o 2.1.0. **Proposto.**

### O que muda no esquema (1.3, só acréscimo)

- `schema.json`: `calibration$covariates`, os nomes das covariáveis
  contínuas da fórmula da calibração. Cada uma é coluna de `occasions.csv`
  (e está em `occasion_covariates`). Vazio ou ausente = calibração só com
  score e chaves, como no 1.2.
- Com `calibration$covariates` não vazio:
  - `occasion_draws.parquet` é obrigatório (com regra, semente, `p_min` e
    número de draws em `schema.json`, como no 1.2);
  - `calibration_draws.parquet` **não** é escrito. Motivo: o crywolf 2.1.0
    monta as versões a partir dessa grade sempre que a regra ou a semente
    pedidas não batem com as gravadas; com a grade presente, ele ignoraria a
    covariável sem avisar. Sem a grade, o 2.1.0 para com erro ("The export
    has no calibration draws"), o que é seguro.
- Faixa como chave (`tp ~ logit_conf + (1 | temp_bin)`) continua como no
  1.1/1.2, com grade; não usa `calibration$covariates`.
- **Ressalva do sheepdog (28/10/2026):** hoje o `occasion_draws()` trabalha
  na pasta exportada, a partir da grade, e a pasta não guarda o ajuste. Com
  covariável contínua, as versões são sorteadas dentro do
  `export_occupancy_data()`, que tem o ajuste em mãos: o export ganha
  `n_versions`, `rule`, `p_min` e `seed` e grava o `occasion_draws.parquet`
  direto. O `occasion_draws()` sobre a pasta dá erro claro nesse caso. O
  `schema.json` também registra a padronização da covariável, se a fórmula
  usar a versão padronizada.

Fica de reserva, só se for preciso sortear com outras regras ou sementes no
crywolf: `prediction_draws.parquet` (draw × predição, com P(TP)). Não
proposto agora: é grande (predições × draws) e o `occasion_draws()` cobre o
uso atual.

### Divisão

**sheepdog** (a sessão do sheepdog se ofereceu para o export e o teste de
contrato):

- `fit_calibration()` com termos contínuos e `calibration_ppc(group = )`
  pela faixa da covariável, antes de exportar;
- export 1.3: `calibration$covariates`, `occasion_draws.parquet`
  obrigatório e sem `calibration_draws.parquet` quando há covariável
  contínua;
- fixture `export-1.3` (por exemplo `tp ~ logit_conf + temp`) e teste de
  contrato: todas as ocasiões em `occasion_draws`, `y` NA só nas ocasiões
  sem valor de covariável ou chave, `schema.json` coerente;
- aviso no export quando uma variável usada nos estratos do
  `random_detections` não está na fórmula da calibração (item 3). **Adiado
  (ressalva do sheepdog):** o export não sabe quais foram os estratos, porque
  o manifesto e a validação não guardam o argumento `strata`. Registrá-lo
  mexe no formato do manifesto, fixo desde a v1.2.0. Até a 1.3, valem a
  mensagem do `select_clips()` e a vinheta `validation-effort`, que ganha
  uma frase sobre sítio com modelo global e alocação não proporcional (o
  `random_detections` sempre estratifica por sítio × espécie).

**crywolf** (pequeno; entra quando o Diogo decidir, já que o foco agora é o
artigo da MEE):

- `fp_draws_data()`: com `calibration$covariates` não vazio, usa só
  `occasion_draws.parquet`, sem exigir `calibration_draws.parquet`; se a
  regra, a semente, o `p_min` ou o número de draws pedidos não estiverem
  gravados, erro dizendo para exportar de novo com eles
  (`sheepdog::export_occupancy_data(n_versions = , rule = , seed = )`);
- `fp_read_sheepdog()`: mostrar as covariáveis da calibração no `print()`;
- testes com a fixture 1.3; `dev/check-sheepdog-fixture.R` passa a aceitar
  o nome da fixture (hoje fixo em `export-1.2`);
- NEWS e help; parágrafo na vinheta sobre a análise 2 × 2 (item 2).

## 2. Mesma covariável na calibração e no modelo de ocupação

O que já se sabe:

- fpcalib: quando p10 varia com a covariável de ocupação, `fp ~ x` corrige o
  viés de `occ_x` (resultado relatado pela sessão do sheepdog; não conferi
  os números aqui).
- Curiango (26/10/2026, `roadmap-integration.md`): com a calibração sem mês,
  `fp_temp` = +0,17; com a calibração por mês, −0,04. O `fp ~ x` acusa o
  erro da calibração ao longo de x.

Proposta para o `fp_occu_mc()`, sem código novo (ele já aceita `fp`):

| | `fp ~ 1` | `fp ~ x` |
|---|---|---|
| calibração sem x | o fluxo atual | `fp_x` mede o erro da calibração ao longo de x |
| calibração com x | o que a parte 3 permite | `fp_x` perto de 0 se a calibração com x ajusta |

- Reportar como análise de sensibilidade, como a sazonal do Curiango, e
  conferir antes o `calibration_ppc()` por faixa.
- Ressalvas, as mesmas do mês: P(TP | score, x) absorve a prevalência do
  canto ao longo de x (a hipótese ainda não testada que a Fase 2 deve
  conferir primeiro), e os minutos revisados entram na calibração e na
  verificação do `fp_data()`.
- Para decidir qual linha é a principal: um cenário do fpcalib com p10
  variando com x, calibração com e sem x, e `fp ~ 1` e `fp ~ x`, medindo
  viés e cobertura de `occ_x`. Sem isso, não há base para recomendar.
- O caminho com princípio continua sendo a Fase 2 (x nas distribuições f1 e
  f0 dos scores), depois do `fp_occu_dyn()`.

## 3. Pesos 1/`inclusion_prob`

- O sorteio estratificado por faixa tem probabilidade de inclusão exata. Para
  um modelo condicional P(TP | score, ...), um sorteio que depende só de
  variáveis do modelo é ignorável: se a faixa (ou x) está na calibração, não
  há o que ponderar.
- Se a faixa não está na calibração e o TP varia com ela, o ajuste sem pesos
  estima P(TP | score) na mistura de faixas da amostra, não na das
  predições. Os pesos corrigem a estimativa pontual, mas com pesos o brms dá
  uma pseudo-posterior: os pesos multiplicam a log-verossimilhança, e a
  incerteza não é a do desenho. É justamente a incerteza que o
  `fp_occu_mc()` propaga.
- Proposta: não implementar pesos. O que estratifica o sorteio entra na
  calibração; o `select_clips()` já lembra disso, e o export pode avisar
  (item 1, lado do sheepdog).
- O mesmo raciocínio vale para os estratos que já existem: se o sítio
  estratifica o `random_detections` com alocação não proporcional e a
  calibração escolhida não tem sítio, o modelo global é uma média na mistura
  de sítios da amostra. No Curiango isso não pesa: os minutos anotados foram
  sorteados em proporção ao esforço de cada gravador, e o BirdNET ficou com
  o modelo global.
