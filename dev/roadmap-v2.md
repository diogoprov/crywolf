# Roadmap do crywolf 2.0.0

Escrito em 22/10/2026, a partir do 1.6.0.9000 (sete vinhetas novas,
sheepdog via R-universe). Decisões do Diogo em 22/10/2026 marcadas com
**decidido**; o resto é plano, sujeito a revisão a cada marco.

## O que o 2.0.0 significa

- **Decidido:** o 2.0.0 sai **antes** do aceite na MEE. Não espera o artigo.
- **Decidido:** tentar o **CRAN via `instantiate`**. Isto substitui a
  decisão de 15/10/2026 de só reavaliar o CRAN depois do aceite.
- Três critérios para o 2.0.0:
  1. **API estável**: argumentos, nomes de saída e classes congelados;
     mudanças incompatíveis só aqui, anunciadas no 1.9.
  2. **Qualidade automatizada**: a validação manual feita a cada mudança no
     Stan (draws idênticos com opções desligadas, log_lik R = Stan,
     recuperação) vira teste no CI.
  3. **Pronto para CRAN**: `R CMD check --as-cran` limpo com o backend do
     `instantiate`.

Fora do escopo do 2.0.0 (decidido em 14/10/2026): multiestação/dinâmico,
multiespécies, efeitos aleatórios em todos os submodelos. O efeito
aleatório na detecção (abaixo) é uma exceção pontual, decidida pelo Diogo
por ser comum em ecoacústica. **Decidido em 22/10/2026:** o modelo
dinâmico (multiestação) e o multiespécies voltam ao plano, depois do
2.0.0 e nessa ordem de prioridade (seções "Depois do 2.0.0" no fim).

## Ordem dos marcos

**Decidido:** qualidade antes de ciência (1.7 → 1.8 → 1.9 → 2.0), porque o
teste de regressão protege as mudanças seguintes.

## 1.7.0: infraestrutura e qualidade

| Item | Esforço | Pronto quando |
|---|---|---|
| Teste de regressão da posterior: resumo com semente fixa e 1 thread em 3–4 desenhos (Miller + gravador, Chambert, `fp = NULL`, Markov), comparado com snapshot | P | `devtools::test()` falha se uma média posterior mudar além da tolerância; roda no CI |
| Validações de `dev/validate-*.R` como testes (identidade com opções desligadas; log_lik R vs Stan em todas as combinações) | M | as combinações rodam em `tests/` com `skip_on_cran()` |
| Hash do Stan gravado nos `.rds` pré-computados + `dev/precompute-all.R` | P | a vinheta avisa quando o `.rds` veio de outro Stan; um script refaz tudo |
| Reprodutibilidade em paralelo: gerar os conjuntos de `fp_calibrate()` no processo principal | P | `parallel = TRUE` e `FALSE` dão os mesmos dados com a mesma semente |
| Refatorar `site_loglik()`/`partial_loglik()` (≈50 argumentos posicionais) com tuplas | M | **feito em 22/10/2026**; exige CmdStan ≥ 2.40 (decisão do Diogo); `log_lik` idêntico ao 1.6.0 em sete desenhos e draws idênticos em seis |
| Cache do CmdStan no GitHub Actions | P | R-CMD-check abaixo de ~6 min |
| `print.fp_data`: separar sítios verificados de ocasiões confirmadas (Chambert) | P | duas linhas distintas na saída |

## 1.8.0: ciência

| Item | Esforço | Pronto quando |
|---|---|---|
| **Efeito aleatório de sítio na detecção verdadeira** (**decidido**): `fp_occu(random_det = TRUE)`, intercepto aleatório em logit(p11) por sítio (na curiango, por sítio-noite), análogo ao `random_fp`; desligado por padrão | M–G | draws idênticos com a opção desligada; log_lik R = Stan; recuperação em simulação; SBC uniforme; `simulate()`, `fp_ppc()` e `fp_sites()` cientes da opção |
| Exemplo na curiango (**decidido**): ajustar com `random_det` (com e sem Markov) e ver se o PPC, hoje com p = 0, melhora | P (computação) | resultados pré-computados e seção na vinheta `markov` ou vinheta própria |
| Limiar 0,2 da curiango com `fp ~ temp + rh` (com `fp ~ 1` o efeito da temperatura troca de sinal) | P (computação ~30 min) | resultado na vinheta `threshold` |
| `fp_threshold_sweep()` derivando MaxF1/MaxPrecision de um conjunto de validação | P–M | **feito na 1.7.0** (branch do Codex `feat/threshold-selection`, mesclado e ajustado) |
| LOO com efeitos de sítio (`random_fp`, `random_det`): `log_lik` integrado sobre o efeito do sítio | M | k de Pareto aceitável sem `reloo` no exemplo da vinheta `fp_variation` |
| Desenho de Chambert no caminho rápido (GLM por sítio) | M | só com dado real grande que justifique; medir antes |

Para o `random_det`, conferir a identificabilidade: com `random_fp` e
`random_det` ligados juntos, os dois efeitos por sítio competem pelas mesmas
detecções. Medir antes de permitir a combinação.

**Andamento em 23/10/2026:**

- Limiar 0,2 com `fp ~ temp + rh`: **rodado** (`vignettes/threshold-fp-results.rds`).
  `occ_temp` vai de +0,16 (0,04; 0,28) com `fp ~ 1` para +0,09 (−0,02;
  0,21), e `fp_temp` = +0,27; em 0,4 nada muda. O FP ligado à temperatura
  explica parte da troca de sinal, não toda. Falta o texto na vinheta
  `threshold`.
- Protótipo do `random_det` (`dev/random-det/`), 20 réplicas: recupera
  `sd_det_site` sem custo quando não há heterogeneidade; o `random_fp` sozinho
  absorve heterogeneidade de detecção como se fosse de FP (sd_fp 1,43 com
  verdade 0); os dois juntos recuperam os dois sd e o declive, mas com
  R-hat > 1,01 em 30–60% dos ajustes. **Decidido pelo Diogo:** implementar,
  permitindo a combinação, com aviso de mixing.
- `random_det` **implementado no pacote** (Stan com quatro chamadas do
  `reduce_sum` e zeros como dados; draws idênticos com a opção desligada em
  quatro desenhos; `log_lik` R = Stan ≤ 5e-14 nos caminhos GLM, vetorizado,
  Chambert e combinado; `loo()` integra em grade de Gauss-Hermite 2D com os
  dois efeitos, 30 nós por efeito, diferença de 6e-4 contra 60). Saída
  `sd_det_site`, prior `fp_priors(det_site_sd = 1)`.
- **Feito em 24/10/2026:** 100 réplicas dos cenários (1600 ajustes),
  SBC com `random_det` (200 réplicas, ranks uniformes), curiango com
  `random_det` com e sem Markov (LOO prefere os dois termos, 113, EP 22),
  vinheta `site_effects` com os gráficos, texto do limiar 0,2 na vinheta
  `threshold`. O PPC da curiango deixa de acusar a variação entre noites,
  mas isso é esperado com um efeito por noite (o `fp_ppc()` replica cada
  sítio com o próprio efeito); a evidência é o LOO.
- **Resta da 1.8:** só o desenho de Chambert no caminho GLM (opcional).
  Próximo: congelar os nomes da saída (`help("crywolf-output")`), agora que
  o `sd_det_site` existe.

## 1.9.0: congelamento da API

Sem o pacote lifecycle (decidido em 23/10/2026): os avisos de transição
são mensagens do próprio crywolf.

- **Decidido:** `threads_per_chain = "auto"` passa a ser o padrão no 2.0.
  **Feito na 1.9 (23/10/2026):**
  - O help de `fp_occu()` explica as opções: `"auto"` (rápido, draws não
    reproduzíveis exatamente pela semente, mesma posterior), `1`
    (reprodução exata), ou um número.
  - Mensagem, uma vez por sessão, para quem chama `fp_occu()` sem passar
    `threads_per_chain`. É uma condição com classe
    (`crywolf_threads_default`): os ajustes repetidos do próprio pacote
    (`.quiet_threads()`: calibração, SBC, reloo, varredura de limiar) a
    descartam sem gastá-la.
  - Testes, vinhetas (chunks ao vivo e o código mostrado do curiango; o
    código mostrado do fp_variation não tem `seed` e ficou como estava) e
    os scripts `dev/precompute-*` que rodavam com 1
    thread fixam `threads_per_chain = 1`; os que já usavam mais threads
    (markov, threshold, threshold-fp) ficaram como estavam.
  - Teto de 2 núcleos sob `R CMD check` (`_R_CHECK_LIMIT_CORES_` definido e
    diferente de `"false"`): threads ≤ 2 e `parallel_chains` ≤
    `2 %/% threads`. Menos cadeias em paralelo não mudam os draws.
  - **Feito (26/10/2026, no 2.0):** padrão `"auto"`; mensagem da 1.9 e
    conselho depois de ajustes lentos apagados (com `.quiet_threads()`);
    `fp_calibrate()`, `fp_calibrate_fit()` e `fp_sbc()` passam
    `threads_per_chain = 1` (sobrescrevível em `fit_args`). Antes: trocar o
    padrão, apagar a mensagem e decidir o padrão
    dos ajustes internos repetidos. `fp_calibrate()` e `fp_sbc()` rodam um
    ajuste por worker (o help do `fp_calibrate()` recomenda
    `threads_per_chain = 1`, "the default"); devem continuar com 1,
    passando-o explicitamente, para não multiplicar núcleos por workers.
- **Decidido (23/10/2026):** efeitos aleatórios continuam como argumentos
  lógicos separados (`random_recorder`, `random_fp` e, quando entrar, o
  `random_det`), sem uma interface única do tipo `random = ~ ...`. Motivo:
  o escopo é estreito (três efeitos fixos, cada um com prior e saída
  própria: `sd_recorder`, `sd_fp_site`, `sd_det_site`), e argumentos lógicos
  são mais fáceis de documentar e validar.
- **Decidido (23/10/2026):** `fp_data_wide()` fica na API do 2.0, como
  ponte para quem vem do `unmarked::occuFP()` com matrizes sítio × ocasião
  (`vignette("unmarked")`).
- Formato legado do sampler (`validacao_por_clipe.csv`): lido em todo o
  2.x; remoção no 3.0 **anunciada** no help de `fp_read_validation()` e no
  NEWS (23/10/2026), sem aviso em tempo de execução.
- Nomes da saída (`occ_*`, `det_*`, `fp_*`, `verified_*`, `p10_bar`, ...):
  congelar em `help("crywolf-output")`. **Feito na 1.9.0 (24/10/2026)**,
  com o `sd_det_site` já incluído; seção "Stability" na página. Ficam para
  o 2.0: congelar argumentos e classes (critério 1) e estender a garantia
  aos resultados de outras funções (`fp_read_sheepdog()`, `loo()` etc.).
  **Feito (26/10/2026):** `help("crywolf-api")` com o que fica congelado
  (argumentos: nomes, ordem e padrões; classes; componentes e colunas da
  seção Value) e a lista de funções estáveis; `fp_occu_mc()`,
  `fp_draws_data()` e `fp_to_unmarked()` marcadas como experimentais;
  `tests/testthat/test-api.R` grava num snapshot a assinatura das funções
  estáveis e dos métodos (sem as reexportações).
- **Release 1.9.0 (24/10/2026, decisão do Diogo):** um único release com a
  ciência da 1.8 e a API da 1.9; não houve tag 1.8.0.
- Exigências mínimas: R ≥ 4.2; CmdStan ≥ 2.40 desde a 1.7.0 (subiu de 2.34
  com o refactor com tuplas, por decisão do Diogo em 22/10/2026).

## Fora do código: manuscrito e documentação (pedido em 23/10/2026)

Ordem pedida pelo Diogo, depois de fechar o `random_det`:

1. Revisar o manuscrito no Google Docs à luz das simulações e benchmarks
   recentes (LOO integrado, `random_fp`/`random_det`, limiar 0,2, threads)
   e ver o que melhorar no estudo de simulação do artigo (compêndio
   fpcalib).
2. Avançar no texto: Discussão, Introdução, Abstract, material
   suplementar, cover letter e o pacote de submissão para a MEE.
3. Traduzir para pt-BR as principais vinhetas do crywolf, no formato do
   sheepdog (`articles/pt/`, como em
   https://provetelab.org/sheepdog/articles/pt/sheepdog.html).

## 2.0.0: lançamento e CRAN

- Backend via `instantiate`: **conferido em 22/10/2026**
  (`dev/cran-instantiate-check.md`, patch em
  `dev/instantiate-prototype.patch`). `stan_threads` funciona, os draws são
  idênticos, e o `--as-cran` sem CmdStan dá só a NOTE esperada. Falta rodar
  o check nas três plataformas no CI.
- `R CMD check --as-cran` limpo nas três plataformas; exemplos e vinhetas
  dentro dos limites de tempo e de núcleos do CRAN (as vinhetas longas já
  são pré-computadas).
- **Feito (26/10/2026, PR #5, CI verde nas três plataformas):** patch
  aplicado, com quatro correções da revisão: sem executável instalado,
  `fp_compile()` vai direto ao cache (o `instantiate` compilaria um sem
  threads na biblioteca); o `install.libs.R` não compila na instalação
  temporária do `devtools::load_all()` (`devtools_install_`) e uma falha de
  compilação não derruba a instalação; o CI ganha o R-universe do Stan
  (o pak não lê `Additional_repositories`), remove o `cmdstanr` nos jobs
  sem CmdStan e roda os jobs com CmdStan sem `--as-cran` (o WARNING de
  flags `-Wno-*` do CmdStan). Draws idênticos ao `main`.
- NEWS do 2.0.0 com a lista de mudanças incompatíveis; CITATION.cff
  atualizado; a citação do artigo da MEE entra quando ele for aceito
  (2.0.x).
- Release no GitHub, DOI de versão no Zenodo, build no R-universe e
  submissão ao CRAN.
- **Lançado em 26/10/2026:** PR #6 (merge 69483a6), tag v2.0.0, DOI de
  versão 10.5281/zenodo.23677193. `cran-comments.md` pronto; a submissão
  ao CRAN fica com o Diogo.

## 2.1.0 (27/10/2026)

Release só da integração com o sheepdog, experimental
(`roadmap-integration.md`): leitura dos esquemas 1.1 e 1.2 do export
(ffd4ffc, d5ddc78), PR #7 (merge f4610c8), DOI de versão
10.5281/zenodo.23692218. Interface estável sem mudanças.

**Próximo release (decidido em 02/11/2026):** esperar o aceite do 2.1.0 no
CRAN (submetido em 01/11/2026, em inspeção manual). Depois do aceite,
trazer o `cran-comments.md` do branch `cran-2.1`, lançar o 2.2.0 no GitHub
e no Zenodo com o `fp_occu_dyn()` e a vinheta `dynamic`, e submeter o 2.2.0
ao CRAN uns 1–2 meses depois do aceite do 2.1.0 (a política do CRAN pede
atualizações espaçadas), talvez já com o começo do multiespécies.

**Atualização (05/11/2026):** o 2.1.0 segue em `incoming/newbies` no CRAN, e o
Diogo decidiu não esperar o aceite para o release no GitHub e no Zenodo. O 2.2.0
sai quando o ajuste de 16 espécies do FrogNet (`crywolf-ms-fu`) chegar e for
conferido. O envio do 2.2.0 ao CRAN continua depois do aceite do 2.1.0 e só com
pedido do Diogo. Se o CRAN pedir correções no 2.1.0, a resposta sai do branch de
envio, não da `main`.

## Riscos

- O `instantiate` pode não suportar threading do jeito que o pacote usa;
  nesse caso, decidir entre CRAN sem threads e ficar no R-universe.
- `random_det` com poucos ocasiões por sítio pode misturar mal, como o
  `random_fp`; medir R-hat e ESS em simulação antes de documentar.
- Mudar o padrão de threads altera resultados de quem reexecuta análises
  antigas com a mesma semente; o NEWS e a mensagem do 1.9 precisam dizer
  isso com clareza.

## Depois do 2.0.0, prioridade 1: modelo dinâmico (multiestação)

**Decidido em 22/10/2026:** vem antes do multiespécies, por exigir menos
esforço e reaproveitar mais do que já existe. Função nova,
`fp_occu_dyn()`; a API do `fp_occu()` congelada no 2.0 não muda.
**Decidido em 27/10/2026:** vem também antes da Fase 2 da integração
(`fp_occu_score()`, ver `roadmap-integration.md`).

O período primário (temporada) é um argumento do `fp_occu_dyn()`, não uma
coluna do esquema do sheepdog: ele sai de `night` ou `recording_date`,
que as ocasiões do export já têm (inclusive as sem predição), ou de uma
covariável de ocasião. Se a calibração usar o mesmo agrupamento (por
exemplo `calibration_keys = "month"`), a chave já vem como texto nas
ocasiões.

### O modelo

O estado de cada sítio muda entre temporadas como uma cadeia de Markov de
dois estados: ocupação inicial $\psi_1$, colonização $\gamma_{it}$ e
extinção $\epsilon_{it}$ (ou persistência $\phi = 1 - \epsilon$), com
covariáveis de sítio e de temporada. Dentro de cada temporada, a
observação é a do `fp_occu()`: desenhos de Miller e de Chambert,
verificação, gravador, covariáveis de falso positivo, termo de Markov
dentro da temporada.

A verossimilhança marginaliza os estados pelo algoritmo forward de um
modelo oculto de Markov. As emissões são as log-verossimilhanças dos ramos
ocupado e desocupado de cada sítio-temporada, que o Stan atual já calcula
por sítio. Como a transição muda com as covariáveis, o forward é escrito à
mão (dois estados, `log_sum_exp`), em vez do `hmm_marginal()`, que exige
transição fixa. O `reduce_sum` continua sobre sítios, cada fatia com sua
sequência de temporadas.

### Passos

| Passo | Pronto quando |
|---|---|
| Dados: coluna de temporada no `fp_data()` (automatizados, verificados e confirmados); sítios sem dados numa temporada contribuem só com a transição | `print()` e os avisos de desenho mostram sítios por temporada |
| Stan: emissões por sítio-temporada reaproveitando `site_loglik()`; forward com transições por sítio-temporada | com uma temporada, log-densidade igual à do `fp_occu()` (~1e-12); verossimilhança em R igual ao `log_lik` do Stan |
| Fórmulas `psi1`, `col` e `ext`, com nomes fixados antes do primeiro release que as contenha | saídas `psi1_*`, `col_*`, `ext_*` documentadas na página de nomes da saída |
| Estados latentes: forward-filter backward-sample em R para ocupação por temporada | `fp_sites()` por temporada; `simulate()` e `fp_ppc()` com dinâmica |
| Quantidades derivadas: ocupação por temporada, rotatividade, ocupação de equilíbrio quando a transição é constante | na saída de `summary()` |
| Identificação com falsos positivos: um falso positivo isolado parece colonização seguida de extinção | simulação mostrando quanta verificação por temporada é preciso para recuperar colonização e extinção; resultado na vinheta |
| Validação | recuperação em simulação, SBC uniforme, calibração por `fp_calibrate()` com desenho dinâmico |
| Vinheta | curiango com as noites de cada sítio como temporadas (uso noite a noite), comparada por LOO com o modelo de sítio-noites independentes |

**Passos 1 e 2 feitos em 01/11/2026.** `fp_data(season = )`;
`fp_occu_dyn()` interno (não exportado até o passo 3 fixar os nomes), com
colonização e extinção só com intercepto. As emissões ficaram em
`src/stan/fp_emissions.stanfunctions`, incluído pelo `fp_occu.stan` (draws
idênticos aos de antes) e pelo `fp_occu_dyn.stan`. Cada sítio-temporada é
uma unidade, numerada `(sítio - 1) * nseason + temporada`, e os dados de
observação das unidades saem do `.build_stan_data()` atual. Os efeitos
aleatórios de sítio são por sítio. O Stan já recebe `Xcol` e `Xext` por
unidade, então o passo 3 só mexe em R. Com uma temporada, a log-densidade
é igual à do `fp_occu()` a ~3e-14 (`dev/check-dyn-one-season.R`), e o
`log_lik` também (`test-fp_occu_dyn.R`). Com várias temporadas, a
verossimilhança em R (`.dyn_loglik()`) é igual ao `log_lik` do Stan.
PR #8.

**01/11/2026, decidido com o Diogo (pedido do sheepdog):** temporadas que
são datas (`Date`, ou texto `YYYY-MM-DD`/`YYYYMMDD`, como `night` e
`recording_date` do export) são preenchidas pelo `fp_data()`. As datas que
faltam entre a primeira e a última entram como temporadas sem dados, no
passo do menor intervalo. Assim cada transição vale um passo (uma noite),
mesmo em noites em que nenhum sítio gravou. Gravações em dias alternados
por sítio já eram tratadas pela grade sítio × temporada. Foi descartada a
alternativa de transição dependente do intervalo (φ^Δt, offset log Δt):
ela complica as fórmulas sem ganho de interpretação. Número e factor
continuam como vêm (consecutivas = uma transição), com aviso se o
espaçamento for desigual. O esquema do sheepdog não muda.

**Passo 3 feito em 01/11/2026 (PR 2):** `fp_occu_dyn(psi1, col, ext, ...)`
exportado como experimental. As covariáveis de temporada entram por
`fp_data(season_covs = )`, por temporada ou por sítio e temporada, e
valem também para `det` e `fp`. As covariáveis de uma transição são as da
temporada em que ela entra. Os nomes são `psi1_*`, `col_*` e `ext_*`
(`crywolf-output`). Há `print()`, `summary()` e `as_draws_df()`. As
quantidades `p11_ref`, `p10_ref`, `p10_bar` e `p_verified` são as do
`fp_occu()`, com `p10_bar` calculado sobre sítio-temporadas.

**Passo 4 feito em 01/11/2026 (PR 3):** FFBS em R (`.dyn_ffbs()`), com a
probabilidade suavizada de cada sítio-temporada pelo forward-backward e uma
amostra conjunta dos estados por draw. `fp_sites()` com modelo dinâmico dá
uma linha por sítio-temporada (`psi` pelas covariáveis, `occupied` dados
todos os dados do sítio). `fp_ppc()` soma detecções por sítio-temporada,
com os estados do FFBS. `simulate()` sorteia estados novos por
colonização e extinção. Simulador `simulate_fp_dyn()` experimental. Em
casos pequenos, o FFBS bate com a enumeração das sequências de estados
(1e-10) e as frequências amostradas com as probabilidades exatas; um
ajuste curto (150 sítios × 4 temporadas) recupera todos os parâmetros
dentro dos intervalos de 95%. PR #8 com merge no main em 01/11/2026.

**Passo 5 feito em 01/11/2026 (PR 4):** quantidades derivadas no
`generated quantities` do `fp_occu_dyn.stan`, por draw: `psi_mean[t]`
(ocupação esperada média sobre sítios), `psi_fs[t]` (ocupação realizada,
com FFBS no Stan), `turnover[t]` (t ≥ 2) e `psi_eq` (só com `col` e `ext`
só com intercepto). Ficaram no Stan, e não em R, para não guardar matrizes
draws × ocasiões e para entrarem no `as_draws_df()` com R-hat e ESS. O
`summary()` ganha a coluna `season`; o `print()` mostra só os parâmetros;
os diagnósticos de convergência continuam só sobre os parâmetros.
Verificado: `psi_mean`, `turnover` e `psi_eq` iguais às fórmulas a partir
dos coeficientes (1e-8), média de `psi_fs` igual à do forward-backward em
R, e `psi_fs` cobrindo a ocupação realizada simulada.

**Divisão dos passos seguintes em PRs (aprovada em 01/11/2026):**
PR 3, passo 4 (FFBS, `fp_sites()` por temporada, `simulate()`/`fp_ppc()` e
simulador dinâmico); PR 4, passo 5; PR 5, passos 6 e 7 (rodadas na AWS);
PR 6, passo 8. O escopo de cada um é confirmado com o Diogo antes.

**Estimativa para os passos 6 e 7 na AWS (01/11/2026).** No Mac, com 4
cadeias de 1000+1000 em 1 núcleo, um ajuste leva 32 s com 100 sítios × 5
temporadas × 4 ocasiões e 87 s com 200 × 10 × 4, recuperando os valores
simulados. Com o fator de 1,14 da c7a medido no fpcalib, ~7.300 ajustes
(identificação ~4.500, recuperação ~1.200, SBC ~1.000, `fp_calibrate()`
~600) somam ~170 núcleo-hora. Isso dá ~1 h por lote numa c7a.48xlarge
spot, ~US$ 3 a 5 no total, com teto de US$ 10 a 15 se o desenho crescer.

**Passos 6 e 7 rodados em 01/11/2026 (PR 5), sem o `fp_calibrate()`.**
Scripts em `dev/dyn-sim/` e tabelas em `dev/dyn-sim/results.md`. Foram
4.450 tarefas numa c7a.24xlarge spot (94 workers), em ~2,8 h e com uma
interrupção retomada. Nenhum dos 6.700 ajustes falhou. Na identificação e
na recuperação, três dos 5.700 ajustes tiveram R-hat ≥ 1,01 ou
divergências.

- Identificação (passo 6), 100 sítios, p11 = 0,6, 50 réplicas por célula.
  Com `fp = ~ 1`, colonização, extinção, `psi1`, `p10_bar` e `psi_fs` têm
  viés < 0,02 e cobertura de 0,88 a 1,00 em todas as 45 células. Isso vale
  inclusive sem nenhuma sítio-temporada verificada. Neste desenho, a
  separação entre p11 e p10 e o prior do `fp_gap` bastam para identificar
  o modelo. Verificar 40% das sítio-temporadas estreita o intervalo de
  `p10_bar` em 20 a 25% e o de colonização em 5 a 12%. O modelo sem falsos
  positivos superestima a colonização em +0,06, +0,15 e +0,28 com p10 de
  0,02, 0,05 e 0,10, e a cobertura cai para 0 a 0,82. A extinção é bem
  menos afetada (viés < 0,035). O desenho não testou p11 baixo, em que a
  verificação deve pesar mais.
- Recuperação com covariáveis (passo 7), 300 réplicas por desenho.
  Coberturas de 0,92 a 0,98 para todos os coeficientes nos quatro
  desenhos (100 ou 200 sítios × 5 ou 10 temporadas). Os vieses dos
  declives são pequenos diante do RMSE e caem com o tamanho do desenho.
- SBC (passo 7), 1.000 réplicas. Ranks uniformes (qui-quadrado, 10
  classes) para 8 dos 10 parâmetros. O intercepto da extinção tem p =
  0,003, com excesso no meio do histograma (posterior um pouco larga
  demais, isto é, conservadora). O intercepto de fp tem p = 0,048. Cinco
  ajustes tiveram divergências.

**Passo 8 feito em 01–02/11/2026 (PR #14, merge 7d2586f; vinheta em PR
seguinte).** LOO não serve: o log_lik do `fp_occu_dyn()` é por sítio e a
curiango tem 7 sítios. A comparação com sítio-noites independentes foi por
validação cruzada em blocos de 7 noites (10 dobras iguais para todos os
sítios; primeira e última semana sempre no treino, para a grade de noites
ser a mesma), em `dev/dyn-kfold/`, numa c7a.24xlarge spot (~2 h). Modelos:
`col`/`ext ~ temp_d + rh_d` (médias diárias entre sítios, que existem em
todas as noites) contra `occ ~ temp_d + rh_d`, ambos com `det ~ bio + aci`,
`fp ~ temp + rh` e `random_recorder`. Resultado: elpd dinâmico −
independente = +191,5 (EP 29,0, provavelmente um pouco pequeno), melhor nas
10 dobras e nos 7 sítios. O ajuste independente da dobra 8 ficou preso numa
moda (R-hat 1,54) e foi refeito com outra semente; o elpd quase não mudou.
Ajuste com todos os dados (`dev/precompute-dynamic.R`, 35 min em 8
núcleos): uma noite em uso deixa de ser usada com probabilidade 0,23
(sequências de ~4,4 noites), uma noite sem uso passa a ser usada com 0,045;
noites úmidas iniciam menos sequências. O PPC dá p = 0, com discrepância
parecida com a do modelo independente: a falta de ajuste está dentro das
noites (minutos independentes), não na dinâmica. `markov` com o dinâmico
não foi ajustado. Vinheta `dynamic` com os resultados pré-computados.

**Validação concluída (05/11/2026).** Critérios de saída de experimental em
`dev/dyn-val/README.md`, definidos antes dos resultados.

| Critério | Situação |
|---|---|
| 1, convergência | não atendido como escrito: 9 células com p11 = 0,15 sem verificação; `random_det` com sd = 1 precisa de 2000 + 2000. O Diogo considera pontual |
| 2, recuperação | atendido em 45 de 45 células; cobertura mínima 0,906, viés máximo 0,024 |
| 3 a 7, limites, calibração, SBC, `random_det`, `markov` | atendidos, com os limites descritos na vinheta |
| 8, segundo conjunto de dados | atendido com os pikas de Goldman et al. (2023), em `dev/dyn-pika/`: `dyn` 178,7 (EP 42,6) acima de `ind` por LOO de parcelas |
| 9, um release sem mudança de interface | pendente: 2.2.0, depois os nomes entram em `tests/testthat/test-api.R` |
| 10, documentação dos limites | atendido |

Na curiango, a validação cruzada em blocos pôs `dyn_mk` em primeiro: 206,9
(EP 33,4) acima de `ind_mk` e 481,3 (EP 120,1) acima de `dyn`. Isso respondeu ao
risco abaixo sobre dinâmica e efeito de noite: `random_det` melhora `ind`, mas
menos que `markov`. A vinheta `dynamic` tem a DAG do modelo e figuras da validação,
do k-fold e dos pikas (PR #22).

**Depois do 2.2.0:**
1. Feito: `fp_occu_dyn()` e `simulate_fp_dyn()` saíram de experimental e entraram
   no snapshot da API (critério 9); chega aos usuários no release seguinte.
2. `fp_calibrate()` com desenhos dinâmicos, previsto no passo de validação e não
   feito. Hoje a função só simula uma temporada, e a vinheta avisa disso.
3. Opcional: dados de lobos de Miller et al. (2013), pedidos à Montana FWP em
   05/11/2026, como comparação com temporadas independentes.

### Riscos

- Tempo: a observação domina o custo hoje, e ela não muda; o forward
  acrescenta pouco. Medir mesmo assim.
- Com poucas temporadas e pouca verificação por temporada, colonização e
  extinção podem ficar mal identificadas diante dos falsos positivos.
- Na curiango, a dinâmica noite a noite e o efeito aleatório de noite na
  detecção (1.8.0) disputam a mesma variação; comparar os dois.

## Depois do 2.0.0, prioridade 2: modelo multiespécies (MSO)

**Decidido em 22/10/2026:** fica para um 2.x ou para o 3.0.0, sem entrar no
caminho do 2.0. A API do `fp_occu()` congelada no 2.0 não muda; o modelo
multiespécies é uma função nova.

**Decidido em 02/11/2026:** a Fase A vem agora, antes da Fase 2 da
integração (`fp_occu_score()`, `roadmap-integration.md`). Do lado do
sheepdog (sheepdog#74): fixture `export-1.3-ms` com três espécies (uma
rara, só em dois sítios) e calibração `(logit_conf | species) + (1 |
site:species)`. Do lado do crywolf: copiar a fixture para os testes e para
o `dev/check-sheepdog-fixture.R`, e ler todas as espécies juntas
(`fp_read_sheepdog(species = NULL)` → `fp_data()` em formato longo; hoje o
`fp_draws_data()` para com mais de uma espécie). Dados reais para a
vinheta (decidido em 02/11/2026, `roadmap-integration.md` do sheepdog): o
FrogNet (DiverSons, 7 sítios, 2021–2025), num export fora dos repositórios
dos pacotes, depois do depósito dos dados. O VicFrogNET foi descartado por
não publicar a lista de gravações por sítio (sem ela não há
não-detecções).

### Fase A (2.x, depois do dinâmico): comunidade com falsos positivos por espécie

**Ordem decidida em 05/11/2026:** um efeito por unidade em logit(p10) só entra no
pacote depois de `u`, falso positivo por sítio × espécie e sazonalidade (PRs 3 a
5), com recuperação por simulação e a escolha entre `sigma_t` comum ou por espécie.

`fp_occu_ms()`: cada espécie tem ocupação, detecção verdadeira e taxa de
falso positivo próprias, com os coeficientes sorteados de distribuições da
comunidade (hierárquico, parametrização não centrada). O falso positivo de
cada espécie é independente das outras: o classificador erra ao marcar a
espécie, não ao trocá-la por outra.

| Passo | Pronto quando |
|---|---|
| Dados: `fp_data()` com coluna de espécie (formato longo), verificação por espécie | `print()` e os avisos de desenho funcionam por espécie |
| Stan: z marginalizado por par sítio × espécie; `reduce_sum` sobre os pares; reaproveitar a verossimilhança por sítio do `fp_occu()` | com uma espécie, a log-densidade é igual à do `fp_occu()` (diferença ~1e-12) |
| Hierarquia: médias e desvios da comunidade para occ, det e fp; gap por espécie (`alpha_1s - gamma_0s > 0`) | recuperação em simulação com 5, 20 e 50 espécies; SBC uniforme |
| Saídas: riqueza por sítio, ocupação condicional por espécie e sítio, `p10` por espécie | `fp_sites()`, `fp_ppc()`, `simulate()` e `predict()` com espécie |
| Custo e mixing: medir antes de documentar | tempo e R-hat em 20 espécies num conjunto do tamanho da curiango; documentar `threads_per_chain` |
| Vinheta com dados simulados e, se houver, dados reais de várias espécies classificadas pelo BirdNET | resultados pré-computados como nas vinhetas atuais |

Riscos: o tempo cresce mais ou menos linearmente com o número de espécies;
os desvios da comunidade podem misturar mal com poucas espécies (como o
`random_fp`); a verificação por espécie pode ser rara demais para as
espécies raras, e o empréstimo de informação da comunidade é justamente o
que se quer medir.

### Estado da Fase A em 14/11/2026

Os cinco PRs de `dev/proposal-mso.md` estão no main, com o
`fp_occu_ms()` ainda marcado como em desenvolvimento.

| Passo | Estado |
|---|---|
| Dados por espécie (PR 1) e Stan com hierarquia (PR 2, #17) | feito; igualdade com `fp_occu()` espécie a espécie |
| Declives por espécie (PR 3, #30) | feito |
| Saídas: `fp_sites()`, `fp_ppc()`, `simulate()`, `fp_richness()`, `simulate_fp_ms()` (PR 4, #31); `predict()` | feito (`predict()` em 14/11/2026) |
| Falso positivo por grupo × espécie (`fp_group`) e por unidade (`random_fp`) (#32) | feito |
| Simulações e SBC na AWS, vinheta `multispecies` (PR 5, #34) | feito com o desenho reduzido de 11/11 (SBC com 200 réplicas, não 1.000); SBC sem sinal de má calibração |
| Parametrização centrada (`centered`, `centered_occ`) (#34, #36) | feito; 50 espécies: convergência 10% (padrão), 75% (`centered`), 100% (os dois) |
| Custo num conjunto do tamanho da curiango; `threads_per_chain` | feito em 17/11/2026 (#44): 8 threads por cadeia, 9,5 h em 32 núcleos; 42 divergências, a explicar com o ajuste `full-sdrec` (#43) |
| Dados reais (FrogNet) | à espera do depósito no Zenodo |

### Atualização de 17/11/2026

- **O FrogNet deixou de ser critério para sair de "em desenvolvimento"**
  (decidido pelo Diogo em 17/11/2026). O passo 4 vira uma versão posterior
  da vinheta, depois do depósito no Zenodo. O passo 5 passa a depender só
  dos passos 2 e 3.
- **Passo 3**: feito no #44 (vinheta e `?fp_occu_ms`). Falta um parágrafo:
  se as 42 divergências vêm da fronteira de `sd_recorder` em 0 (verdadeiro
  0 na simulação). O ajuste `full-sdrec`, com `sd_recorder = 0.5`, roda na
  AWS (run `crywolf-ms-cost`).
- **Passo 2**: a grade (`dev/ms-sim/run_sims.R`, 23 células, 20 réplicas)
  roda na AWS em on-demand: células 1–18 (run `crywolf-ms`), 21–22
  (`crywolf-ms-long`) e 23 (`crywolf-ms-c23`). A célula 19 (`ms`, 20
  espécies × 300 sítios) ficou de fora por custo (16/11/2026). Com a grade,
  a regra do padrão sai do `summarise.R` e entra na vinheta.

### Atualização de 19/11/2026: multiespécies fora de "em desenvolvimento"

- **Passo 2 feito.** A grade terminou (440 ajustes, 20 por célula). O Diogo
  escolheu o padrão `"auto"` para `centered` e `centered_occ`: os dois
  `TRUE` com ≥ 20 espécies e ≥ 100 sítios, `FALSE` caso contrário. Com 30
  sítios ou 5 espécies, não centrar converge mais; com 20+ espécies e 100+
  sítios, centrar os dois converge 100% e é o mais rápido. Designs fora da
  grade (10 espécies, mais de 10 ocasiões) não foram testados.
- **Passo 3 fechado sem o `full-sdrec`.** O ajuste com `sd_recorder = 0.5`
  bateu no limite de 30 h (cadeias de 16, 22 e 26 h; a quarta em 80% da
  amostragem). O Diogo decidiu documentar sem rerodar: a origem das 42
  divergências fica em aberto na vinheta.
- **Passo 5 feito.** `fp_occu_ms()`, `simulate_fp_ms()`, `fp_richness()`,
  as classes `fp_occu_ms` e `fp_sim_ms` e o argumento `species` do
  `fp_data()` entram na interface estável a partir da 2.3.0. Release só a
  pedido.
- **Merge do #46 em 19/11/2026** (bcac82c no main). A versão de
  desenvolvimento segue 2.2.0.9000. A 2.3.0 sai só a pedido do Diogo.

### O que falta no multiespécies depois de 19/11/2026

- **Release 2.3.0**: só a pedido. A NEWS já traz as três entradas do
  multiespécies.
- **FrogNet na vinheta** (passo 4): depois do depósito no Zenodo, sem
  data.
- **Integração com o sheepdog**: a sessão do sheepdog respondeu em
  19/11/2026. O ρ implausível na curiango vinha da mistura de f1, e a PPC
  ainda reprova a autocorrelação. O export-1.3-ms continua válido; o
  esquema atual é o 1.4. Plano: (1) exemplo ponta a ponta com
  `fp_occu_ms()` no export-1.3-ms e (2) `fp_occu_mc()` aceitando dados
  multiespécies, **feitos no #48 (19/11/2026)**; (3) `fp_occu_score()` com
  uma espécie: o estágio 1 do sheepdog só com o ct01+ct02 da curiango (pesos
  1/π₁, ct03 só como sensibilidade) foi refeito em 19/11/2026, seguido do
  estágio 2 com réplicas do bootstrap (20/11), do estágio 1 por faixas de
  score (21/11) e da cobertura de ψ (22/11). Dos critérios para levar o v1 ao
  crywolf, o 2 passou; falta o 3, a rodada aleatória do FrogNet (plano em
  `dev/fase2/PLANO-v1-markov.md`, no sheepdog; atualizado em 23/11/2026);
  (4) `fp_review_priority()`, ligado ao `select_directed()` do sheepdog,
  **feito no #62 (23/11/2026)**, experimental: `size` = q(1 − q) de
  `fp_sites()`, sobre `fp_occu()` e `fp_occu_ms()`. O sheepdog aceita a saída
  como `units` e avisa antes de mudar a chave das unidades, a coluna de
  tamanho ou o esquema do export.
- **Pontos fracos da grade**: o efeito do gravador só foi testado com 10
  ocasiões. Designs com 10 espécies ou mais ocasiões ficaram de fora.

### Próximos passos do multiespécies (proposta de 14/11/2026, a decidir)

1. **Feito em 14/11/2026.** **`predict()` para `fp_occu_ms`**: ocupação e detecção por espécie em
   novos dados, como o `predict.fp_occu()`. Pronto quando: com uma espécie,
   igual ao `predict()` do `fp_occu()`.
2. **Feito em 19/11/2026 (#46): `"auto"` pela grade.** **Padrão da parametrização**: decidir se `centered` e `centered_occ`
   continuam `FALSE` por padrão, viram `TRUE`, ou passam a depender do
   número de ocasiões e de sítios. Hoje só há dois pontos de evidência:
   50 espécies em 100 sítios (os dois ajudam) e 2 espécies em 40 sítios
   (o β centrado piora). Uma grade pequena de simulações (sítios ×
   espécies) daria a regra; custo na AWS a estimar com piloto.
3. **Feito em 17/11/2026 (#44); `full-sdrec` documentado sem rerodar.** **Custo e mixing com dados do tamanho real**: 20 espécies num conjunto
   do tamanho da curiango (ou do FrogNet), com `threads_per_chain`
   medido; documentar na vinheta e em `?fp_occu_ms`.
4. **FrogNet na vinheta**, depois do depósito: fs como modelo principal e
   o misto como sensibilidade (decidido na sessão FrogNet em 11/11/2026),
   com `centered = TRUE`.
5. **Feito em 19/11/2026 (#46), sem o passo 4.** **Sair de "em desenvolvimento"**: com 1 a 4 feitos, propor a API como
   estável numa release (2.3.0?). Release só a pedido.
6. **Integração com o sheepdog (Fase 2, `fp_occu_score()`)**: vem depois
   da Fase A (decidido em 02/11/2026). O lado do sheepdog relatou em
   13/11/2026 resultados implausíveis na curiango (ρ pode absorver a
   dependência dos scores); esperar o diagnóstico de lá.
7. **Opções adiadas das perguntas de 02/11/2026**: correlação entre os
   interceptos dos submodelos (LKJ) e efeito do gravador por espécie. Só
   com um caso real que peça.

### Fase B (3.0.0): modelo conjunto (JSDM) com fatores latentes

Acrescentar fatores latentes na ocupação, para capturar a coocorrência
residual entre espécies (o "J" dos modelos conjuntos de distribuição de
espécies).

- Referência de implementação: o jsdmstan (Fiona Seaton, NERC-CEH,
  https://nerc-ceh.github.io/jsdmstan/, GPL >= 3, Zenodo
  10.5281/zenodo.13356516) ajusta GLLVM e MGLMM em Stan, mas sem modelo de
  observação: não tem ocupação com detecção imperfeita, falsos positivos
  nem verificação, e usa rstan (o crywolf usa cmdstanr). Serve de
  referência para a parametrização dos fatores latentes, não como
  dependência. A licença é compatível com a do crywolf; reaproveitar código
  exige crédito explícito. Conferir antes como o jsdmstan identifica as
  cargas dos fatores (restrições de sinal e de triângulo).
- Pronto quando: recuperação das correlações residuais em simulação; o
  modelo da Fase A é o caso com zero fatores (mesma log-densidade);
  comparação por LOO entre 0, 1 e 2 fatores numa vinheta.
- Por que é 3.0.0: muda as saídas e a interface multiespécies de forma
  incompatível com a Fase A, se a Fase A já estiver publicada.

**Avaliação de 25/11/2026 (com o MSO estável e o RN validado):**
- **Distância:** de 2 a 3 semanas de trabalho, no ritmo do RN e do
  multiespécies. Hoje não há código nem protótipo.
- **O que já serve:**
  - A estrutura do `fp_occu_ms()`. As unidades sítio × espécie são
    independentes dado o sítio. Os fatores entram como λⱼ′zᵢ no preditor da
    ocupação, e dado zᵢ a verossimilhança continua fatorando por unidade. A
    parte de detecção, falsos positivos e verificação fica intacta.
  - A integração por Gauss–Hermite do `loo()`, já usada para os efeitos de
    sítio. O `log_lik` passa a ser agrupado por sítio, e não por unidade,
    porque zᵢ é comum às espécies do sítio.
  - A receita de validação: piloto local, SBC e grade na AWS.
- **Riscos, em ordem:**
  1. Identificação das cargas: triângulo inferior com diagonal positiva.
     Mesmo assim, cada cadeia pode convergir para um sinal diferente dos
     fatores, e o R-hat explode. Ver antes como o jsdmstan resolve isso.
  2. Coocorrência aparente por detecção correlacionada, e não por ocupação
     correlacionada, com falsos positivos. A pergunta da simulação é se as
     correlações residuais são recuperáveis com o nível de verificação dos
     dados reais. O RN mostrou que a resposta pode ser "parcialmente".
  3. Custo: nsítios × K parâmetros a mais. No FrogNet (5 M linhas), o
     ajuste já é pesado, mas o custo extra fica só na ocupação.
  4. Saídas novas: a matriz de correlação residual ΛΛ′, `print`,
     `fp_sites()`, `simulate()`, `fp_ppc()` e uma vinheta com 0, 1 e 2
     fatores comparados por LOO.
- **Versão:** o critério "zero fatores = a mesma log-densidade da Fase A"
  permite fazer o JSDM como um argumento `n_factors = 0` do `fp_occu_ms()`.
  Assim nada quebra, e ele pode sair numa 2.x, e não necessariamente na 3.0.0.
- **Primeiro passo:** protótipo de viabilidade em `dev/`, de 1 a 2 dias,
  como a etapa 1 do RN. O foco são os riscos 1 e 2.
- **Protótipo feito (25/11/2026, branch `jsdm-proto`, `dev/jsdm/README.md`,
  sem merge):**
  - Risco 1 resolvido: diagonal livre com sinal fixado depois da amostragem,
    como no jsdmstan. A diagonal positiva misturou mal.
  - Risco 2 resolvido para a detecção correlacionada: um efeito de sítio na
    detecção, comum às espécies, remove o fator espúrio.
  - O padrão de ΛΛ′ é recuperado com 100 sítios, mas a magnitude não.
  - Risco novo, a confusão entre espécies no classificador. Ela vira
    associação espúria em ΛΛ′ já com confusão moderada (0,15 por ocasião),
    e o efeito de detecção não corrige.
  - O `fp_occu()` e o `fp_occu_ms()` absorvem a confusão moderada. Com
    confusão forte (0,30), a ocupação sai enviesada para cima, em 0,6 a 1,7
    logit. A premissa agora está documentada nos dois.
- **Decidido pelo Diogo em 25/11/2026: JSDM arquivado.**
  - O JSDM só faria sentido com um termo de confusão. Declarar os pares a
    priori não basta, porque só os dados verificados separam confusão de
    co-ocorrência.
  - Ele é o modelo mais frágil à confusão e recupera só o padrão de ΛΛ′, com
    custo alto. Para o público do crywolf, o retorno não compensa.
  - A branch `jsdm-proto` fica como registro.
- **Condição para reabrir:** a validação registrar a espécie verdadeira dos
  falsos positivos, para que a confusão seja estimada como uma matriz de
  confusão. O sheepdog já permite declarar essa espécie no fluxo de
  validação, mas é opcional, e o `fp_read_sheepdog()` ainda não lê o campo.
  Seria um modelo de classificação errada entre espécies (ver "Fora do
  plano"), um projeto de pesquisa maior que o JSDM.
- **Alternativa em avaliação:** exportar draws posteriores dos estados de
  ocupação do `fp_occu_ms()` para JSDMs externos (Hmsc, gllvm), por
  imputação múltipla. Viabilidade em `dev/jsdm-export/`, com um conjunto de
  dados:
  - Funciona e é barato.
  - No gllvm, o padrão de ΛΛ′ melhora em relação às detecções brutas
    (correlação com a verdade de 0,91 contra 0,82).
  - A força das associações sai atenuada (|ΛΛ′| médio de 0,28, com verdade
    de 0,69). Isso é esperado, porque o `fp_occu_ms()` sorteia os estados das
    espécies de forma independente.
  - O Hmsc ficou inconclusivo e precisa de configuração própria.
  - Próximo passo, se for adiante: uma grade local só com gllvm e, depois,
    uma função experimental de exportação com essas limitações documentadas.

### Fora do plano

Confusão entre espécies (o classificador marca A quando quem canta é B):
exige uma matriz de classificação e estado latente multinomial, com
verificação desenhada para identificá-la. É problema de pesquisa, não de
implementação; não entra no crywolf.


## Depois das tarefas atuais: Royle–Nichols com falsos positivos (pedido de 19/11/2026)

**Origem:** pedido da sessão FrogNet, com os dados ACAD (Panwar et al. 2026,
privados). O `random_det` absorve a variação de atividade (abundância ou taxa
de canto) entre sítios. A pergunta do artigo é outra: se os falsos negativos
do BirdNET têm estrutura de sítio além da atividade.

**Modelo proposto:**
- N_i ~ Poisson(λ_i) e p11_i = 1 − (1 − r_i)^N_i, com um efeito aleatório de
  sítio em r além de λ.
- Mantém p10, as verificadas e a confirmação de Chambert.
- Fica como opção do `fp_occu()`, sem dinâmico nem multiespécies.

**Decidido pelo Diogo em 19/11/2026:** a etapa 1 entra no roadmap e só
começa depois das tarefas em andamento. A implementação, as etapas 2 a 5,
fica para depois da etapa 1.

**Andamento (24/11/2026):**
- **Etapa 1, feita (#65):** favorável com 200 sítios; parcialmente favorável
  no tamanho do ACAD (21 sítios). Resultados em `dev/rn-fp/README.md`.
- **Etapa 2, feita (#66):** o Diogo decidiu por uma função própria,
  `fp_occu_rn()`, e não por uma opção do `fp_occu()`. Falsos positivos
  aditivos nos sítios ocupados. Entram a confirmação de Chambert e os efeitos
  de sítio em r, log λ e p10, além do efeito por gravador. Ficam de fora as
  covariáveis de fp por ocasião e o `markov`.
- **Grade do desenho ACAD realista, pedida pelo FrogNet (#67):** só
  confirmação. O posterior de σ_r acompanha o prior, então a análise do ACAD
  deve mostrar σ_r sob mais de um prior. Detalhes em
  `dev/rn-fp/acad-real/README.md`.
- **Etapa 3, feita (#68):** `loo()`, `simulate()`, `fp_sites()`, `fp_ppc()`
  e `fp_loo_weights()` para `fp_occu_rn`. Todas usam o posterior de N de cada
  sítio, calculado em R e testado contra o `log_lik` do Stan. `predict()`,
  `coef()`, `fp_effects()` e `plot()` ficam para depois, se forem pedidos.
- **Etapa 4, feita (25/11/2026):** SBC (2 × 500) e grade geral (72 células ×
  100) na AWS, por cerca de US$10. O SBC está calibrado; só `beta[2]` no
  desenho A ficou no limite. O poder para σr depende sobretudo de sítios e
  ocasiões, e pouco do esquema de anotação. Com 100 sítios, 1000 + 1000
  iterações dão R-hat de 1,01–1,03. Detalhes em `dev/rn-fp/sim/README.md`.
- **Etapa 5, feita (25/11/2026):** vinheta `royle-nichols` (lê os
  resultados pré-computados das etapas 1 e 4), nomes das saídas em
  `crywolf-output`, e `fp_occu_rn()` e `simulate_fp_rn()` passam à interface
  estável a partir da release seguinte à 2.3.0 (número a decidir).

**Etapa 1, o ponto de decisão (1–2 dias):**
- Simular dados RN com falsos positivos e ajustar um protótipo em Stan fora do
  pacote, em `dev/`.
- Pergunta: um efeito de sítio em r, somado à variação em λ, é recuperável
  sem anotação completa?
- Se não for, o modelo não responde à pergunta. O caminho é então um GLMM dos
  verdadeiros positivos do BirdNET com a atividade anotada da ACAD como offset,
  fora do crywolf.

**Etapas seguintes, só se a etapa 1 for favorável (estimativa: 1,5–2,5
semanas no total, extrapolada do `fp_occu_dyn()` e do multiespécies):**
2. Programa Stan e interface em R (2–3 dias). Marginaliza N_i de 0 a N_max,
   o que multiplica o custo da verossimilhança por cerca de N_max. É preciso
   decidir como r11 das verificadas depende de N.
3. Saídas (`fp_ppc()`, `loo()`, `fp_sites()`, `simulate()`, `print()`) e
   testes, incluindo a igualdade entre a verossimilhança em R e a do Stan (2
   dias).
4. SBC e grade de simulações (2–3 dias de relógio). Provavelmente na AWS,
   com o custo confirmado antes.
5. Vinheta, documentação, review e check (1–2 dias).

A versão de contagem (minutos com detecção por manhã, com mistura de minutos
verdadeiros e falsos) somaria cerca de uma semana.

## 2.3.0 (preparada em 19/11/2026)

Multiespécies estável (#46), `fp_occu_mc()` multiespécies (#48), correções do
CRAN (#49), `curiango$verified$test_set` (#50) e a correção do `fp_ppc()` com
confirmação (#51). O Diogo pediu o release em 19/11/2026, mas a tag espera a
resposta da Tessa sobre a autoria (20/11/2026): se ela ficar, reverter a saída
dela no `DESCRIPTION`, no `CITATION.cff`, no NEWS e no `cran-comments.md`
antes da tag. Ajustar `date-released` ao dia do release.
