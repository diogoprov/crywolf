# Roadmap de integração sheepdog + crywolf

O plano conjunto está no repositório do sheepdog:
https://github.com/diogoprov/sheepdog/blob/main/dev/roadmap-integration.md
(PR https://github.com/diogoprov/sheepdog/pull/35). Este arquivo só aponta
para lá, para não haver duas versões.

O que toca o crywolf, em ordem de prioridade:

- **Fase 0, contrato de dados:** `fp_read_sheepdog()` lê o diretório
  escrito por `sheepdog::export_occupancy_data()` (esquema
  `sheepdog-occupancy`, v1) e recusa ou avisa quando a validação vem de
  `sampling_mode = "top"`.
- **Fase 1, propagação Monte Carlo:** `fp_occu_mc()` ajusta o modelo nos
  históricos de `sheepdog::history_draws()` e junta as posteriores.
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
