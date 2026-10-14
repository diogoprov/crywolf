# Como montar o repositório do pacote (roteiro para o MacBook)

Este diretório `dev/` fica fora do pacote construído (está no `.Rbuildignore`).

## 1. Nome

O pacote já foi renomeado para `crywolf` (14/10/2026). O `rename_package.sh`
da raiz pode ser apagado. A memória do projeto para o Claude Code está em
`CLAUDE.md`, na raiz.

## 2. Criar o repositório

```bash
git init -b main
git add .
git commit -m "Initial package scaffold"
gh repo create diogoprov/crywolf --public --source=. --push
```

## 3. Primeiros comandos no R

```r
# uma vez: CmdStan
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
# cmdstanr::install_cmdstan()

# LICENSE.md (GPL-3) já criado com usethis::use_gpl_license(3, include_future = FALSE)
devtools::document()          # regenera man/ e NAMESPACE
devtools::test()              # 69 testes; os que usam Stan levam ~1 min
devtools::check()             # a vinheta principal ajusta 3 modelos (~5 min)
pkgdown::build_site()         # opcional, local
```

Para publicar o site do pkgdown: em Settings > Pages, escolher a branch
`gh-pages` depois que o workflow `pkgdown.yaml` rodar pela primeira vez.

## 4. O que já foi testado (14/10/2026, Linux, R 4.3.3, CmdStan 2.37.0, cmdstanr 0.9.0)

- `R CMD check --no-manual`: 0 erros, 0 notas, 1 aviso que vem só do
  ambiente da sessão (o locale en_US.UTF-8 não existia lá). Inclui os 69
  testes e as duas vinhetas.
- O modelo Stan novo dá exatamente a mesma log-densidade que a versão com
  laço por ocasião (diferença < 1e-12 em 5 pontos da posterior). As
  ocasiões com covariáveis idênticas dentro do sítio viram linhas binomiais,
  o que é exato e deixou um ajuste com 300 sítios × 20 ocasiões cerca de 9
  vezes mais rápido (22 s contra ~200 s, 4 cadeias em 4 núcleos) quando a
  detecção não tem covariáveis por ocasião.
- Recuperação de parâmetros: 30 conjuntos simulados (300 sítios, 10
  ocasiões, confusão estruturada, `fp = ~ w`): cobertura do declive de
  ocupação 0,93 e do efeito da confusão 0,97.
- A vinheta principal reproduz o padrão do artigo nos dados de exemplo: o
  modelo de taxa constante estima o declive em −0,25 (verdade −0,8),
  `fp = ~ x` em −0,56 e `fp = ~ w` em −0,78; o LOO prefere `fp = ~ w`.
- `fp_verification_curve()` em paralelo (furrr, 2 workers) roda de ponta a
  ponta num desenho mínimo.

## 5. Pendências

- Ordem dos autores do meio (Liliana, Tessa) no `DESCRIPTION` e no
  `CITATION.cff`. Diogo primeiro e Larissa por último já estão fixados.
- Depois do release v0.1.0 do compêndio fpcalib, citar o DOI dele no README.
- Conferir se `simulate_fp(seed = s)` reproduz exatamente os conjuntos do
  compêndio: a ordem das chamadas ao gerador é a mesma, mas isso ainda não
  foi comparado diretamente.
- Ideias para versões seguintes: o desenho de confirmação de observações
  (observation confirmation) de Chambert et al. (2015); `plot()` para `fp_calibration`;
  ajuste a dados reais do curiango como terceira vinheta.
