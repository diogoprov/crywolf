# Roadmap de integração sheepdog + crywolf

O plano conjunto está no repositório do sheepdog:
https://github.com/diogoprov/sheepdog/blob/main/dev/roadmap-integration.md
(PR https://github.com/diogoprov/sheepdog/pull/35). Este arquivo só aponta
para lá, para não haver duas versões.

O que toca o crywolf, em ordem de prioridade:

- **Fase 0, contrato de dados:** `fp_read_sheepdog()` lê o diretório
  escrito por `sheepdog::export_occupancy_data()` (esquema
  `sheepdog-occupancy`, v1) e recusa ou avisa quando a validação vem de
  `sampling_mode = "top"`. **Feito em 24/10/2026** (lado do sheepdog no
  PR diogoprov/sheepdog#36): descarta os clipes "top" com aviso
  (`keep_top = TRUE` os mantém), devolve as tabelas no formato de
  `fp_aggregate(by = c("site", "occasion_id"))` e `fp_data()`; teste de
  contrato com a exportação real do sheepdog, que roda no CI quando o
  r-universe publicar o sheepdog 3.2.1.9000.
- **Fase 1, propagação Monte Carlo:** `fp_occu_mc()` ajusta o modelo nos
  históricos de `sheepdog::history_draws()` e junta as posteriores.
  **Feito em 24/10/2026 (crywolf 1.9.0.9000):** `fp_occu_mc()` (mistura
  igual; `summary()` com variância dentro e entre ajustes; com M = 1, draws
  idênticos ao `fp_occu()`) e `fp_draws_data()`, que monta as versões a
  partir de `calibration_draws.parquet` (regra `"sample"` por padrão,
  decisão do Diogo; `"threshold"` como alternativa). Falta: escolher M
  (estabilidade com M = 10, 20, 50) e a comparação na curiango.
  **25/10/2026:** uma primeira rodada (M = 50, `fp = NULL` e `fp ~ 1`) foi
  descartada: o export tinha um artefato (score do minuto repetido em
  cada linha, ~7 detecções por minuto) e a calibração era por janela, não
  por minuto (corrigido em sheepdog#39). Com o export corrigido, a regra
  `"sample"` dá 7,9% dos minutos com detecção (Custom) e `"threshold"`
  1,3%. Os 21,4% das gravações anotadas não servem de referência para o
  total (meses concentrados).
  **Rodada no export corrigido (25/10/2026, `dev/mc-curiango/results.rds`):**
  M = 50, 7,7% dos minutos com detecção por versão. 99 de 100 ajustes
  convergiram (o 50º de `fp ~ temp + rh`, R-hat 1,53, fica fora: valores
  de M = 30). M = 10 já basta. Ocupação, `fp ~ 1`: temperatura −0,07
  (−0,24; 0,11), umidade −0,64, lua 0,15; `fp ~ temp + rh`: temperatura
  −0,15 (−0,30; 0,01), umidade −0,54, lua 0,14, e `fp_temp` = +0,17 (0,09;
  0,24). Limiar 0,4 com `fp ~ 1`: −0,25, −0,78, 0,14. Leitura: 30–70% da
  variância dos efeitos de ocupação vem da calibração (o limiar fixo a
  esconde); o FP volta a crescer com a temperatura, como no limiar 0,2, e
  modelá-lo aproxima o efeito da temperatura do valor do limiar 0,4.
  Ideia: `fp_occu_mc()` poderia avisar ou excluir ajustes que não
  convergem.
- **Fase 2, modelo conjunto score–ocupação:** `fp_occu_score()`, função
  nova, com as distribuições dos scores dos verdadeiros e dos falsos
  positivos, o termo de Markov e o desenho de verificação. **[decidir]**
  onde entra na fila depois do 2.0.0: antes do `fp_occu_dyn()`, depois dele
  ou dentro dele.
- **Fase 3:** `fp_review_priority()` para a camada dirigida da validação.
- **Fases 4 e 5:** atividade vocal e comunidade, alinhadas com a Fase A do
  multiespécies do [roadmap v2](roadmap-v2.md). A confusão entre espécies
  continua fora do plano.

Nada disso muda a API congelada no 2.0.0: tudo entra em funções novas.
