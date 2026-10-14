# Esboço: vignettes/crywolf.Rmd, "Introduction to crywolf"

Formato de referência: "Introduction to spOccupancy" (Doser lab, artigo
`modelfitting`). Vinheta longa de tutorial; títulos e prosa em inglês; sem
travessão na prosa; descrição do modelo bayesiano segundo o BARG (Kruschke
2021) e o checklist de Davis & Kay (2023).

Legenda: **[chunk]** = bloco de código; **[fig]** = figura logo após o chunk;
`[CITAR: ...]` = referência ainda não verificada (não entra sem conferir).

---

## 0. Cabeçalho e setup

- YAML: `title: "Introduction to crywolf"` (decidido), `html_vignette`,
  `VignetteIndexEntry{Introduction to crywolf}`.
- [chunk] `knitr::opts_chunk$set(..., eval = crywolf:::.has_cmdstan())`,
  `fig.width`, `fig.height`, `out.width`.
- [chunk] `set.seed(2026)` logo no início, como no spOccupancy.

## 1. Introduction

- O problema em dois parágrafos: classificadores automáticos (ex.: BirdNET)
  geram falsos negativos e falsos positivos; ocupação padrão só trata os
  primeiros; falsos positivos inflam psi e podem enviesar declives.
- O que o crywolf faz: modelo de múltiplos métodos (site confirmation) de
  Miller et al. (2011), em Stan via cmdstanr, com taxa de FP constante ou
  função de covariáveis de sítio, intercepto aleatório de gravador, LOO por
  sítio, simulação e calibração.
- Mensagem central (Provete et al., in preparation): quando o FP acompanha um
  gradiente correlacionado com a covariável de ocupação, a taxa constante
  perde a cobertura do declive e verificar mais sítios não resolve; modelar o
  FP resolve em parte (`fp = ~ x`) ou quase totalmente (`fp = ~ w`).
- Lista dos modelos/funções cobertos, com links para a referência:
  `fp_data()`, `fp_data_wide()`, `fp_occu()`, `fp_priors()`, métodos,
  `simulate_fp()`, `fp_calibrate()`.
- Installation: chunk `eval = FALSE` com cmdstanr, `install_cmdstan()`,
  `remotes::install_github("diogoprov/crywolf")`; `library(crywolf)`.

## 2. Example data

- 2.1 **The multiple-methods design**: o que é uma ocasião automática e uma
  verificada; verificação em subconjunto de sítios; detecções verificadas são
  certas.
- 2.2 **Simulated data**: [chunk] `simulate_fp(nsite = 300, nocc = 20,
  psi = 0.5, b1 = -0.8, p11 = 0.6, p10 = 0.05, conf_effect = 1.5,
  conf_cor = 0.7, frac_verified = 0.3, seed = 2026)`; explicar `x`
  (covariável de ocupação, padronizada) e `w` (gradiente que confunde; não
  observado na prática, presente aqui como controle positivo).
  - [chunk] `sim`, `head(sim$data$auto)`, `head(sim$data$sites)`,
    `head(sim$data$verified)`, `sim$truth`.
  - [fig] p10 verdadeiro por sítio contra `x`, cor por `w` (ggplot2, já em
    Suggests). Mostra visualmente a confusão estruturada.
- 2.3 **Your own data**: formato longo de `fp_data()` (três tabelas,
  colunas exigidas, `recorder`) com exemplo mínimo de brinquedo que roda sem
  Stan; menção a `fp_data_wide()` para matrizes estilo unmarked (link para a
  futura vinheta `unmarked`).
- (Futuro: dados reais do curiango, quando liberados.)

## 3. The multiple-methods false-positive occupancy model

### Basic model description

Equações em display (MathJax), vetores em negrito, cada símbolo em prosa:

- Processo ecológico:
  $z_i \sim \text{Bernoulli}(\psi_i)$, $\text{logit}(\psi_i) = \mathbf{x}_i^\top \boldsymbol{\beta}$.
- Método automático:
  $y^{A}_{ik} \mid z_i \sim \text{Bernoulli}(z_i\, p_{11,ik} + (1 - z_i)\, p_{10,i})$,
  $\text{logit}(p_{11,ik}) = \mathbf{v}_{ik}^\top \boldsymbol{\alpha} + u_{r[ik]}$,
  $\text{logit}(p_{10,i}) = (\alpha_1 - \Delta) + \mathbf{w}_i^\top \boldsymbol{\gamma}$.
- Método verificado: $y^{V}_{im} \mid z_i \sim \text{Bernoulli}(z_i\, r_{11})$,
  $\text{logit}(r_{11}) = \delta_0$.
- Intercepto aleatório: $u_r \sim N(0, \sigma_u^2)$.
- Parágrafo sobre $\Delta$ (`fp_gap`): $\Delta > 0$ garante
  $p_{11} > p_{10}$ na referência (covariáveis em zero) e elimina a troca de
  rótulos; por isso centrar as covariáveis.
- Parágrafo sobre o papel das verificações: ancoram $p_{10}$; sítios com
  detecção verificada têm $z_i = 1$ conhecido.
- Marginalização de $z_i$ (fórmula da verossimilhança por sítio) em uma
  frase + equação; detalhes ficam para a vinheta `stan_model`.

### Priors

- Tabela com padrões de `fp_priors()`: interceptos N(0, 1,5); inclinações
  N(0, 2,5); $\Delta \sim N(2, 1{,}5)$ truncada em 0; $\sigma_u \sim$
  half-N(0, 2). Justificativa (escala logito, covariáveis padronizadas).
- [chunk] `fp_priors()`; `fp_priors(slope_sd = 1)` como exemplo de mudança.
- BARG: dizer explicitamente que as prioris são declaradas e mostrar como
  fazer análise de sensibilidade (ajustar com prioris alternativas; chunk
  curto ou só texto, para não pesar).

## 4. Fitting the constant-rate model with fp_occu()

- Argumentos: `occ`, `det`, `fp`, `random_recorder`, `priors`, `chains`,
  `iter_warmup`, `iter_sampling`, `adapt_delta`, `seed`.
- [chunk] `m_const <- fp_occu(sim$data, occ = ~ x, det = ~ 1, fp = ~ 1,
  random_recorder = FALSE, seed = 1)`; `m_const`.
- Interpretação da saída impressa, linha por linha (nomes `occ_*`, `det_*`,
  `fp_*`, `fp_gap`, `p_verified`, `p11_ref`, `p10_ref`, `p10_bar`).

## 5. Convergence and model adequacy

- R-hat (< 1,01), ESS bulk/tail (> 400), divergências e treedepth: o que
  `fp_occu()` já reporta e quando avisa.
- [chunk] `summary(m_const)`; comentar as colunas `rhat`, `ess_bulk`,
  `ess_tail`.
- [fig] traceplots e densidades com bayesplot sobre `as_draws_df()`
  (`mcmc_trace()`, `mcmc_dens_overlay()`, `mcmc_rhat()`).
- **Posterior predictive check** com `fp_ppc()`: [chunk]
  `fp_ppc(m_const, seed = 1)`; explicar a discrepância (Freeman-Tukey) por
  sítio e o p-valor bayesiano; [fig] `bayesplot::ppc_bars(ppc$y, ppc$yrep)`
  e dispersão de `fit_rep` contra `fit_obs`. Ressalva: p-valor perto de 0,5
  não prova que a estrutura do FP está certa (o modelo constante pode
  passar no PPC e ainda ter o declive enviesado), por isso a calibração.
- Interpretação dos resultados só depois dessas checagens.

## 6. Modelling the false-positive rate

- 6.1 Por que o declive sai enviesado: com o FP acompanhando `w`
  (correlacionado com `x`), o modelo constante atribui FP ao processo de
  ocupação. Comparar `coef(m_const, "occ")` com o valor verdadeiro (−0,8).
- 6.2 `fp = ~ x`: [chunk] `m_x`; remédio prático, sempre disponível;
  resultados do estudo de simulação (viés reduzido em 34–62%, cobertura
  parcial; melhora com mais verificação).
- 6.3 `fp = ~ w`: [chunk] `m_w`; controle positivo; quando há uma covariável
  medida do que engana o classificador (ex.: abundância de espécie
  parecida). Mencionar proxies imperfeitos (correlação 0,7: viés reduzido em
  17–37%).
- [fig] estimativas do declive de ocupação (média e IC 95%) dos três modelos,
  com linha no valor verdadeiro.
- Nota: um único ajuste não diz se os ICs estão calibrados; link para a
  vinheta `calibration`.

## 7. Model comparison with LOO

- `log_lik` por sítio com z marginalizado, então `loo()` é
  leave-one-site-out.
- [chunk] `loo::loo_compare(list(constant = loo(m_const), fp_x = loo(m_x),
  fp_w = loo(m_w)))`; interpretar `elpd_diff` e `se_diff`; checar Pareto k.
- Ressalva: o LOO avalia predição, não viés do declive; o modelo com melhor
  elpd não é necessariamente o que recupera o declive.

## 8. Prediction

- [chunk] `predict(m_x, newdata, process = "occ")` e `process = "fp"`
  ao longo de um gradiente de `x`.
- [fig] psi e p10 previstos contra `x`, com faixa de IC, para `m_const` e
  `m_x` lado a lado.
- `summary = FALSE` para obter a matriz de draws.

## 9. Detection covariates and recorder effects

- Covariáveis de ocasião (`bio`, `aci`) e de sítio em `det`; intercepto
  aleatório de gravador.
- Simular um segundo conjunto com `a_bio`, `a_aci`, `sd_recorder > 0`.
- [chunk] `fp_occu(sim2$data, occ = ~ x, det = ~ bio + aci, fp = ~ x,
  random_recorder = TRUE)` (ver decisão 2 sobre rodar ou não).
- Observação sobre custo: covariáveis por ocasião desfazem o colapso
  binomial; o ajuste fica mais lento.

## 10. Next steps

- Links para as vinhetas `calibration`, `simulation`, `unmarked`,
  `stan_model` (só quando existirem).

## 11. References

Autor-ano; só entram referências conferidas.

- Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C., Bailey,
  L. L., & Weir, L. A. (2011). *Ecology*, 92, 1422–1428.
  doi:10.1890/10-1396.1 (já no pacote)
- Chambert, T., Miller, D. A. W., & Nichols, J. D. (2015). *Ecology*, 96,
  332–339. doi:10.1890/14-1507.1 (já no pacote; citar só se a seção de
  "outros desenhos" entrar)
- Kruschke (2021), *Nat. Hum. Behav.*, doi:10.1038/s41562-021-01177-7
  (BARG; provavelmente não citado no texto, só seguido)
- `[CITAR: Clare et al. 2021]`: conciliação com o resultado de efeitos
  ≤ 0,5; conferir referência completa
- `[CITAR: Vehtari, Gelman & Gabry 2017, PSIS-LOO]`: conferir
- `[CITAR: Royle & Link 2006, falsos positivos em ocupação]`: conferir
- `[CITAR: Doser et al. 2022, spOccupancy]`: só se mencionado; conferir
- `[CITAR: BirdNET, Kahl et al. 2021]`: conferir
- Provete et al., in preparation

---

## Decisões (14/10/2026)

1. PPC: feito com a função `fp_ppc()` (R/fp_ppc.R), calculada em R a partir
   dos draws; o Stan não mudou.
2. Computação no site: ainda em aberto (ver conversa).
3. Gráficos de diagnóstico: bayesplot (em Suggests).
4. Título: "Introduction to crywolf".
