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

Fora do escopo (decidido em 14/10/2026): multiestação/dinâmico,
multiespécies, efeitos aleatórios em todos os submodelos. O efeito
aleatório na detecção (abaixo) é uma exceção pontual, decidida pelo Diogo
por ser comum em ecoacústica.

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
| Refatorar `site_loglik()`/`partial_loglik()` (≈50 argumentos posicionais) com tuplas | M | draws idênticos ao 1.6.0 (diferença 0) em todos os desenhos do teste de regressão |
| Cache do CmdStan no GitHub Actions | P | R-CMD-check abaixo de ~6 min |
| `print.fp_data`: separar sítios verificados de ocasiões confirmadas (Chambert) | P | duas linhas distintas na saída |

## 1.8.0: ciência

| Item | Esforço | Pronto quando |
|---|---|---|
| **Efeito aleatório de sítio na detecção verdadeira** (**decidido**): `fp_occu(random_det = TRUE)`, intercepto aleatório em logit(p11) por sítio (na curiango, por sítio-noite), análogo ao `random_fp`; desligado por padrão | M–G | draws idênticos com a opção desligada; log_lik R = Stan; recuperação em simulação; SBC uniforme; `simulate()`, `fp_ppc()` e `fp_sites()` cientes da opção |
| Exemplo na curiango (**decidido**): ajustar com `random_det` (com e sem Markov) e ver se o PPC, hoje com p = 0, melhora | P (computação) | resultados pré-computados e seção na vinheta `markov` ou vinheta própria |
| Limiar 0,2 da curiango com `fp ~ temp + rh` (com `fp ~ 1` o efeito da temperatura troca de sinal) | P (computação ~30 min) | resultado na vinheta `threshold` |
| `fp_threshold_sweep()` derivando MaxF1/MaxPrecision de um conjunto de validação | P–M | limiares calculados dentro da função, como opção |
| LOO com efeitos de sítio (`random_fp`, `random_det`): `log_lik` integrado sobre o efeito do sítio | M | k de Pareto aceitável sem `reloo` no exemplo da vinheta `fp_variation` |
| Desenho de Chambert no caminho rápido (GLM por sítio) | M | só com dado real grande que justifique; medir antes |

Para o `random_det`, conferir a identificabilidade: com `random_fp` e
`random_det` ligados juntos, os dois efeitos por sítio competem pelas mesmas
detecções. Medir antes de permitir a combinação.

## 1.9.0: congelamento da API

Depreciações com `lifecycle::deprecate_warn()` no 1.9 e remoção no 2.0.

- **Decidido:** `threads_per_chain = "auto"` passa a ser o padrão no 2.0.
  - O help de `fp_occu()` explica as opções: `"auto"` (rápido, draws não
    reproduzíveis exatamente pela semente, mesma posterior), `1`
    (reprodução exata), ou um número.
  - No 1.9: mensagem avisando a mudança de padrão para quem não passar o
    argumento.
  - Testes, vinhetas e scripts `dev/precompute-*` passam a fixar
    `threads_per_chain = 1` onde a reprodução importa.
  - **Atenção ao CRAN:** exemplos, testes e vinhetas não podem usar mais de
    2 núcleos. O `"auto"` precisa respeitar `_R_CHECK_LIMIT_CORES_` (ou ter
    teto de 2 sob `R CMD check`), contando cadeias × threads.
- Efeitos aleatórios com argumentos paralelos (`random_recorder`,
  `random_fp`, `random_det`): manter assim, dado o escopo estreito.
- `fp_data_wide()`: manter, como ponte do unmarked.
- Formato legado do sampler (`validacao_por_clipe.csv`): manter no 2.x e
  anunciar a remoção para o 3.0.
- Nomes da saída (`occ_*`, `det_*`, `fp_*`, `verified_*`, `p10_bar`, ...):
  congelar e documentar numa página única.
- Exigências mínimas (R ≥ 4.2, CmdStan ≥ 2.34): subir só se o refactor com
  tuplas ou o `instantiate` exigir.

## 2.0.0: lançamento e CRAN

- Backend via `instantiate`: **conferido em 22/10/2026**
  (`dev/cran-instantiate-check.md`, patch em
  `dev/instantiate-prototype.patch`). `stan_threads` funciona, os draws são
  idênticos, e o `--as-cran` sem CmdStan dá só a NOTE esperada. Falta rodar
  o check nas três plataformas no CI.
- `R CMD check --as-cran` limpo nas três plataformas; exemplos e vinhetas
  dentro dos limites de tempo e de núcleos do CRAN (as vinhetas longas já
  são pré-computadas).
- NEWS do 2.0.0 com a lista de mudanças incompatíveis; CITATION.cff
  atualizado; a citação do artigo da MEE entra quando ele for aceito
  (2.0.x).
- Release no GitHub, DOI de versão no Zenodo, build no R-universe e
  submissão ao CRAN.

## Riscos

- O `instantiate` pode não suportar threading do jeito que o pacote usa;
  nesse caso, decidir entre CRAN sem threads e ficar no R-universe.
- `random_det` com poucos ocasiões por sítio pode misturar mal, como o
  `random_fp`; medir R-hat e ESS em simulação antes de documentar.
- Mudar o padrão de threads altera resultados de quem reexecuta análises
  antigas com a mesma semente; o NEWS e a mensagem do 1.9 precisam dizer
  isso com clareza.
