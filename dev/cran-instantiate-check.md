# Conferência: crywolf no CRAN via instantiate (22/10/2026)

Protótipo feito numa cópia do 1.6.0.9000 (commit 34ca3ec), sem mexer no
pacote. As mudanças estão em `dev/instantiate-prototype.patch`, refeito em
22/10/2026 sobre a 1.7.0.9000 (commit 30e8dcd; Stan com tuplas, CmdStan
≥ 2.40), que se aplica limpo ao `HEAD` com `git apply`. Ambiente: macOS (M5 Pro), R 4.6,
CmdStan 2.40, cmdstanr 0.9.0.9002, instantiate 0.2.4 (CRAN, publicado em
07/10/2026).

## O que o protótipo muda

- `inst/stan/fp_occu.stan` passa para `src/stan/`; `instantiate::stan_package_configure()`
  gera `src/install.libs.R`, `src/Makevars*` e `cleanup*`.
- `src/install.libs.R` compila com `cpp_options = list(stan_threads = TRUE)`.
- DESCRIPTION: `instantiate (>= 0.2.4)` em Imports; `cmdstanr` sai de
  Imports e vai para Suggests (continua em `Additional_repositories`);
  `Remotes` removido (o CRAN acusa campo desconhecido).
- `fp_compile()`: usa o executável compilado na instalação
  (`instantiate::stan_package_model()`); se ele não existir ou não tiver
  threads, compila no cache do usuário como hoje. Esse caminho é o dos
  binários do CRAN, que saem sem o modelo.
- `fp_stan_file()` aponta para `bin/stan/fp_occu.stan` do pacote instalado.
- `.has_cmdstan()` confere também se o cmdstanr está instalado.
- Vinheta `unmarked`: os chunks que ajustam modelos passam a exigir CmdStan
  (era a única vinheta que falhava sem ele).

## Guarda de versão na instalação (22/10/2026)

O Stan com tuplas exige CmdStan 2.40. `src/install.libs.R` só compila na
instalação se encontrar CmdStan 2.40 ou mais novo; caso contrário instala
assim mesmo e avisa. Testado em macOS:

| Situação | Resultado |
|---|---|
| Instalar com CmdStan 2.40 | modelo compilado na instalação; `fp_compile()` usa esse executável, com threads |
| Instalar com CmdStan 2.39 | instala; mensagem "CmdStan 2.40 or later not found (found 2.39.0)"; sem executável |
| Usar esse pacote ainda com 2.39 | `fp_compile()` para com "crywolf needs CmdStan 2.40 or later (found 2.39.0)" |
| O mesmo pacote depois de passar para 2.40 | compila no cache do usuário no primeiro uso, com threads |

## Resultados

| Teste | Resultado |
|---|---|
| Instalar do código-fonte com CmdStan | modelo compilado na instalação (10 s), `stan_threads` ligado |
| Ajuste com o modelo instalado, 1 thread | draws idênticos ao crywolf 1.6.0 (diferença 0) |
| Ajuste com 2 threads | roda |
| Instalar sem CmdStan (como no CRAN) | instala; só copia o `.stan` e avisa |
| Pacote instalado sem modelo, depois com CmdStan | compila no cache com threads; draws idênticos (diferença 0) |
| `R CMD check --as-cran`, sem CmdStan, `_R_CHECK_FORCE_SUGGESTS_=false` | 1 NOTE, 55 s; 212 testes passam, 10 de Stan pulados; vinhetas e exemplos OK |
| `R CMD check --as-cran`, com CmdStan | 1 NOTE + 1 WARNING local (abaixo) |
| Tarball | 3,7 MB (limite do CRAN: 5 MB) |

A NOTE é a esperada: submissão nova, versão de desenvolvimento, e
`cmdstanr`/`sheepdog` em Suggests fora do CRAN, com a disponibilidade
confirmada pelo próprio check via `Additional_repositories`.

O WARNING "compilation flags used" (`-Wno-*`) vem dos makefiles do CmdStan
e só aparece quando o CmdStan compila o modelo durante o check. O CRAN não
tem CmdStan, então ele não ocorre lá.

## Precedente

Sete pacotes no CRAN usam o instantiate (sccomp, hdbayes, footBayes,
clmstan, bvarnet, BSCB, agregR). O footBayes 2.1.0 usa o mesmo esquema
(instantiate em Imports, cmdstanr em Suggests com
`Additional_repositories`, testes com `skip_on_cran()` e
`stan_cmdstan_exists()`) e passa sem problemas em todas as plataformas do
CRAN. Nenhum deles compila com `stan_threads`; o crywolf seria o primeiro,
mas isso só importa onde há CmdStan, ou seja, nunca nas máquinas do CRAN.

## O que ainda não foi testado

- Windows e Linux: o protótipo só rodou no macOS. Antes da submissão,
  rodar o check sem CmdStan no GitHub Actions nas três plataformas (o
  footBayes indica que funciona, mas não com este pacote).
- O padrão `threads_per_chain = "auto"` do 2.0: sob `R CMD check` precisa
  respeitar o limite de 2 núcleos do CRAN. Hoje os testes de Stan são
  pulados no CRAN, então o risco está nos exemplos e vinhetas que rodem
  Stan (nenhum roda sem CmdStan).
- Revisão humana do CRAN: o pacote depende de software externo (CmdStan)
  para a função principal; os pacotes acima foram aceitos nessas
  condições, mas o revisor pode pedir ajustes na documentação.

## Conclusão

Viável. O protótipo passa no `--as-cran` com uma NOTE aceitável e mantém os
draws idênticos. O trabalho do 2.0.0 é aplicar o patch, rodar o check nas
três plataformas no CI, e tratar o padrão de threads.
