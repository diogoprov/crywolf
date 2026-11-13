# Simulations of the multispecies model: results

Written by `dev/ms-sim/summarise.R`.

## Fits

| nspecies | scenario | model | n | converged | with_divergences | max_rhat | min_ess | median_hours |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  5 | base | ms | 50 | 0.98 | 0.02 | 1.01 |  922 | 0.7100 |
|  5 | rare | ms | 50 | 0.90 | 0.10 | 1.01 |  902 | 0.6780 |
|  5 | rare | single | 50 | 1.00 | 0.00 | 1.00 | 1390 | 0.0177 |
|  5 | unverified | ms | 50 | 0.96 | 0.04 | 1.01 |  962 | 0.5920 |
|  5 | unverified | single | 50 | 1.00 | 0.00 | 1.00 | 1260 | 0.0132 |
| 10 | fs | fs | 50 | 0.94 | 0.00 | 1.02 |  456 | 1.8400 |
| 10 | fs | mixed | 50 | 0.92 | 0.02 | 1.02 |  255 | 2.0000 |
| 10 | mixed | fs | 50 | 1.00 | 0.00 | 1.01 |  502 | 1.8500 |
| 10 | mixed | mixed | 50 | 0.42 | 0.02 | 1.03 |  122 | 2.0300 |
| 20 | base | ms | 50 | 0.64 | 0.00 | 1.02 |  300 | 4.3100 |
| 20 | rare | ms | 50 | 0.68 | 0.00 | 1.02 |  430 | 4.6500 |
| 20 | rare | single | 50 | 0.96 | 0.02 | 1.02 |  346 | 0.0200 |
| 20 | unverified | ms | 50 | 0.50 | 0.00 | 1.02 |  336 | 4.3000 |
| 20 | unverified | single | 50 | 1.00 | 0.00 | 1.01 |  812 | 0.0169 |

## Recovery: community parameters

| nspecies | scenario | param | n | bias | rmse | coverage | width |
| --- | --- | --- | --- | --- | --- | --- | --- |
|  5 | base | det_(Intercept) | 50 |  0.010300 | 0.2260 | 0.94 | 1.2000 |
|  5 | base | det_bio | 50 |  0.013500 | 0.0569 | 0.92 | 0.1970 |
|  5 | base | fp_(Intercept) | 50 | -0.061900 | 0.4540 | 1.00 | 2.5000 |
|  5 | base | fp_gap | 50 |  0.072100 | 0.3940 | 1.00 | 2.1900 |
|  5 | base | occ_(Intercept) | 50 |  0.123000 | 0.4080 | 0.94 | 1.8200 |
|  5 | base | occ_x | 50 | -0.026300 | 0.2420 | 0.98 | 1.2700 |
|  5 | base | sd_det_(Intercept) | 50 |  0.124000 | 0.2090 | 0.98 | 1.0200 |
|  5 | base | sd_log_fp_gap | 50 |  0.139000 | 0.2130 | 0.88 | 0.8230 |
|  5 | base | sd_occ_(Intercept) | 50 |  0.000991 | 0.3090 | 0.98 | 1.3400 |
|  5 | base | sd_occ_x | 50 |  0.077600 | 0.2580 | 0.96 | 1.1500 |
|  5 | base | verified_(Intercept) | 50 |  0.002660 | 0.1390 | 0.98 | 0.5990 |
|  5 | rare | det_(Intercept) | 50 | -0.007600 | 0.2190 | 0.92 | 1.1700 |
|  5 | rare | det_bio | 50 | -0.000243 | 0.0565 | 0.94 | 0.2160 |
|  5 | rare | fp_(Intercept) | 50 | -0.046300 | 0.3900 | 1.00 | 2.1800 |
|  5 | rare | fp_gap | 50 |  0.038700 | 0.3450 | 0.96 | 1.8500 |
|  5 | rare | occ_(Intercept) | 50 | -0.339000 | 0.5060 | 0.96 | 2.3400 |
|  5 | rare | occ_x | 50 |  0.052600 | 0.2560 | 0.96 | 1.3000 |
|  5 | rare | sd_det_(Intercept) | 50 |  0.085400 | 0.2160 | 0.98 | 1.0200 |
|  5 | rare | sd_log_fp_gap | 50 |  0.053500 | 0.1900 | 0.94 | 0.7000 |
|  5 | rare | sd_occ_(Intercept) | 50 |  0.361000 | 0.4300 | 0.92 | 1.5400 |
|  5 | rare | sd_occ_x | 50 |  0.076400 | 0.2340 | 0.96 | 1.1900 |
|  5 | rare | verified_(Intercept) | 50 | -0.019600 | 0.1750 | 0.96 | 0.6750 |
|  5 | unverified | det_(Intercept) | 50 | -0.012100 | 0.2370 | 0.96 | 1.1900 |
|  5 | unverified | det_bio | 50 |  0.001640 | 0.0636 | 0.90 | 0.2080 |
|  5 | unverified | fp_(Intercept) | 50 | -0.017800 | 0.2820 | 1.00 | 2.3500 |
|  5 | unverified | fp_gap | 50 |  0.005690 | 0.2840 | 1.00 | 2.0400 |
|  5 | unverified | occ_(Intercept) | 50 | -0.078600 | 0.5080 | 0.90 | 1.8500 |
|  5 | unverified | occ_x | 50 | -0.077300 | 0.3340 | 0.96 | 1.4000 |
|  5 | unverified | sd_det_(Intercept) | 50 |  0.111000 | 0.2670 | 0.92 | 1.0100 |
|  5 | unverified | sd_log_fp_gap | 50 |  0.102000 | 0.1970 | 0.90 | 0.7660 |
|  5 | unverified | sd_occ_(Intercept) | 50 |  0.008250 | 0.2880 | 1.00 | 1.3700 |
|  5 | unverified | sd_occ_x | 50 |  0.156000 | 0.2740 | 0.92 | 1.2200 |
|  5 | unverified | verified_(Intercept) | 50 |  0.009140 | 0.1600 | 1.00 | 0.7080 |
| 20 | base | det_(Intercept) | 50 | -0.002900 | 0.1430 | 0.86 | 0.4820 |
| 20 | base | det_bio | 50 |  0.003400 | 0.0232 | 0.96 | 0.0988 |
| 20 | base | fp_(Intercept) | 50 |  0.001810 | 0.2330 | 0.94 | 0.9980 |
| 20 | base | fp_gap | 50 | -0.004710 | 0.2010 | 0.92 | 0.8840 |
| 20 | base | occ_(Intercept) | 50 | -0.003130 | 0.2060 | 1.00 | 0.9470 |
| 20 | base | occ_x | 50 |  0.005840 | 0.1240 | 0.94 | 0.5330 |
| 20 | base | sd_det_(Intercept) | 50 |  0.016600 | 0.0824 | 0.98 | 0.3740 |
| 20 | base | sd_log_fp_gap | 50 |  0.021100 | 0.0798 | 0.90 | 0.2530 |
| 20 | base | sd_occ_(Intercept) | 50 |  0.041100 | 0.1920 | 0.90 | 0.7250 |
| 20 | base | sd_occ_x | 50 |  0.018700 | 0.1080 | 0.98 | 0.4690 |
| 20 | base | verified_(Intercept) | 50 | -0.013000 | 0.0608 | 0.98 | 0.3060 |
| 20 | rare | det_(Intercept) | 50 | -0.005450 | 0.1150 | 0.94 | 0.4840 |
| 20 | rare | det_bio | 50 |  0.005390 | 0.0228 | 1.00 | 0.1000 |
| 20 | rare | fp_(Intercept) | 50 |  0.029300 | 0.1960 | 1.00 | 0.9770 |
| 20 | rare | fp_gap | 50 | -0.034800 | 0.1560 | 0.98 | 0.8600 |
| 20 | rare | occ_(Intercept) | 50 | -0.095000 | 0.2350 | 1.00 | 1.0700 |
| 20 | rare | occ_x | 50 |  0.021000 | 0.1360 | 0.92 | 0.5110 |
| 20 | rare | sd_det_(Intercept) | 50 |  0.017800 | 0.0804 | 1.00 | 0.3800 |
| 20 | rare | sd_log_fp_gap | 50 |  0.018500 | 0.0592 | 0.96 | 0.2520 |
| 20 | rare | sd_occ_(Intercept) | 50 |  0.181000 | 0.2460 | 0.88 | 0.8110 |
| 20 | rare | sd_occ_x | 50 | -0.014800 | 0.1260 | 0.90 | 0.4630 |
| 20 | rare | verified_(Intercept) | 50 |  0.007910 | 0.0605 | 0.98 | 0.3020 |
| 20 | unverified | det_(Intercept) | 50 | -0.023200 | 0.1150 | 0.98 | 0.4960 |
| 20 | unverified | det_bio | 50 |  0.003130 | 0.0228 | 0.96 | 0.0975 |
| 20 | unverified | fp_(Intercept) | 50 | -0.010900 | 0.2260 | 0.96 | 0.9730 |
| 20 | unverified | fp_gap | 50 | -0.012300 | 0.1910 | 0.94 | 0.8440 |
| 20 | unverified | occ_(Intercept) | 50 |  0.052600 | 0.2380 | 0.90 | 0.9130 |
| 20 | unverified | occ_x | 50 | -0.012100 | 0.1410 | 0.90 | 0.5280 |
| 20 | unverified | sd_det_(Intercept) | 50 |  0.037800 | 0.0887 | 0.98 | 0.3880 |
| 20 | unverified | sd_log_fp_gap | 50 |  0.006890 | 0.0699 | 0.92 | 0.2440 |
| 20 | unverified | sd_occ_(Intercept) | 50 | -0.008800 | 0.1560 | 0.94 | 0.7020 |
| 20 | unverified | sd_occ_x | 50 |  0.009020 | 0.1200 | 0.94 | 0.4690 |
| 20 | unverified | verified_(Intercept) | 50 | -0.007570 | 0.0681 | 0.98 | 0.3020 |

## Recovery: species values, pooled over species

| nspecies | scenario | param | n | bias | rmse | coverage | width |
| --- | --- | --- | --- | --- | --- | --- | --- |
|  5 | base | det_(Intercept)[sp] |  250 | -0.002940 | 0.1230 | 0.940 | 0.4730 |
|  5 | base | fp_gap[sp] |  250 |  0.214000 | 0.8390 | 0.956 | 2.1500 |
|  5 | base | occ_(Intercept)[sp] |  250 |  0.040300 | 0.2700 | 0.952 | 1.0100 |
|  5 | base | occ_x[sp] |  250 | -0.032700 | 0.2520 | 0.944 | 0.9890 |
|  5 | base | p10_ref[sp] |  250 | -0.000404 | 0.0122 | 0.952 | 0.0446 |
|  5 | base | p11_ref[sp] |  250 | -0.000572 | 0.0292 | 0.940 | 0.1110 |
|  5 | rare | det_(Intercept)[sp] |  250 | -0.015500 | 0.1960 | 0.920 | 0.6050 |
|  5 | rare | fp_gap[sp] |  250 |  0.075600 | 0.7080 | 0.932 | 1.4600 |
|  5 | rare | occ_(Intercept)[sp] |  250 |  0.032400 | 0.3150 | 0.940 | 1.2200 |
|  5 | rare | occ_x[sp] |  250 |  0.049600 | 0.3020 | 0.928 | 1.0900 |
|  5 | rare | p10_ref[sp] |  250 |  0.000198 | 0.0125 | 0.944 | 0.0403 |
|  5 | rare | p11_ref[sp] |  250 | -0.003310 | 0.0455 | 0.920 | 0.1410 |
|  5 | unverified | det_(Intercept)[sp] |  250 | -0.008310 | 0.1540 | 0.936 | 0.5110 |
|  5 | unverified | fp_gap[sp] |  250 |  0.164000 | 1.0700 | 0.940 | 1.9000 |
|  5 | unverified | occ_(Intercept)[sp] |  250 |  0.011900 | 0.3110 | 0.944 | 1.0900 |
|  5 | unverified | occ_x[sp] |  250 | -0.035800 | 0.2780 | 0.960 | 1.0700 |
|  5 | unverified | p10_ref[sp] |  250 | -0.000310 | 0.0128 | 0.948 | 0.0426 |
|  5 | unverified | p11_ref[sp] |  250 | -0.001870 | 0.0355 | 0.936 | 0.1180 |
| 20 | base | det_(Intercept)[sp] | 1000 | -0.010600 | 0.1310 | 0.955 | 0.4870 |
| 20 | base | fp_gap[sp] | 1000 |  0.019100 | 0.3700 | 0.957 | 1.2200 |
| 20 | base | occ_(Intercept)[sp] | 1000 | -0.006080 | 0.2520 | 0.960 | 1.0200 |
| 20 | base | occ_x[sp] | 1000 | -0.002360 | 0.2450 | 0.951 | 0.9550 |
| 20 | base | p10_ref[sp] | 1000 |  0.000355 | 0.0129 | 0.945 | 0.0443 |
| 20 | base | p11_ref[sp] | 1000 | -0.002390 | 0.0308 | 0.955 | 0.1140 |
| 20 | rare | det_(Intercept)[sp] | 1000 |  0.001820 | 0.1390 | 0.958 | 0.5150 |
| 20 | rare | fp_gap[sp] | 1000 |  0.005220 | 0.3390 | 0.950 | 1.1800 |
| 20 | rare | occ_(Intercept)[sp] | 1000 |  0.004590 | 0.2930 | 0.938 | 1.0600 |
| 20 | rare | occ_x[sp] | 1000 |  0.008140 | 0.2590 | 0.939 | 0.9510 |
| 20 | rare | p10_ref[sp] | 1000 |  0.000403 | 0.0126 | 0.945 | 0.0437 |
| 20 | rare | p11_ref[sp] | 1000 |  0.000582 | 0.0327 | 0.958 | 0.1210 |
| 20 | unverified | det_(Intercept)[sp] | 1000 |  0.003060 | 0.1280 | 0.948 | 0.4810 |
| 20 | unverified | fp_gap[sp] | 1000 |  0.007780 | 0.3450 | 0.941 | 1.1900 |
| 20 | unverified | occ_(Intercept)[sp] | 1000 |  0.012600 | 0.2590 | 0.940 | 0.9980 |
| 20 | unverified | occ_x[sp] | 1000 | -0.008720 | 0.2550 | 0.935 | 0.9440 |
| 20 | unverified | p10_ref[sp] | 1000 |  0.000588 | 0.0126 | 0.948 | 0.0439 |
| 20 | unverified | p11_ref[sp] | 1000 |  0.000813 | 0.0299 | 0.948 | 0.1120 |

## First species: multispecies (ms) vs alone (single)

| nspecies | scenario | model | param | n | bias | rmse | coverage | width |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  5 | rare | ms | det_(Intercept)[sp1] | 50 | -0.039600 | 0.31600 | 0.88 | 1.0700 |
|  5 | rare | ms | fp_gap[sp1] | 50 |  0.216000 | 1.35000 | 0.88 | 2.2300 |
|  5 | rare | ms | occ_(Intercept)[sp1] | 50 |  0.208000 | 0.43100 | 0.96 | 1.8700 |
|  5 | rare | ms | occ_x[sp1] | 50 |  0.069500 | 0.35300 | 0.94 | 1.4200 |
|  5 | rare | ms | p10_ref[sp1] | 50 | -0.001020 | 0.00696 | 0.96 | 0.0286 |
|  5 | rare | ms | p11_ref[sp1] | 50 | -0.008630 | 0.07360 | 0.88 | 0.2440 |
|  5 | rare | single | det_(Intercept)[sp1] | 50 | -0.157000 | 0.48400 | 0.84 | 1.3000 |
|  5 | rare | single | fp_gap[sp1] | 50 | -0.057400 | 0.45400 | 0.90 | 1.5200 |
|  5 | rare | single | occ_(Intercept)[sp1] | 50 |  0.245000 | 0.43700 | 0.96 | 1.9300 |
|  5 | rare | single | occ_x[sp1] | 50 |  0.076500 | 0.43200 | 0.94 | 1.7300 |
|  5 | rare | single | p10_ref[sp1] | 50 | -0.001310 | 0.00683 | 0.96 | 0.0294 |
|  5 | rare | single | p11_ref[sp1] | 50 | -0.034600 | 0.10400 | 0.84 | 0.2780 |
|  5 | unverified | ms | det_(Intercept)[sp1] | 50 |  0.015000 | 0.19600 | 0.92 | 0.5570 |
|  5 | unverified | ms | fp_gap[sp1] | 50 |  0.341000 | 1.32000 | 0.90 | 2.7800 |
|  5 | unverified | ms | occ_(Intercept)[sp1] | 50 |  0.032800 | 0.35100 | 0.94 | 1.1500 |
|  5 | unverified | ms | occ_x[sp1] | 50 | -0.022800 | 0.27200 | 0.98 | 1.0800 |
|  5 | unverified | ms | p10_ref[sp1] | 50 | -0.000646 | 0.01190 | 0.94 | 0.0409 |
|  5 | unverified | ms | p11_ref[sp1] | 50 |  0.003970 | 0.04400 | 0.92 | 0.1280 |
|  5 | unverified | single | det_(Intercept)[sp1] | 50 | -0.000877 | 0.19200 | 0.88 | 0.5840 |
|  5 | unverified | single | fp_gap[sp1] | 50 |  0.035800 | 0.34100 | 0.96 | 1.3700 |
|  5 | unverified | single | occ_(Intercept)[sp1] | 50 |  0.026800 | 0.34600 | 0.92 | 1.2200 |
|  5 | unverified | single | occ_x[sp1] | 50 | -0.030800 | 0.33300 | 0.96 | 1.2200 |
|  5 | unverified | single | p10_ref[sp1] | 50 |  0.000539 | 0.01130 | 0.96 | 0.0421 |
|  5 | unverified | single | p11_ref[sp1] | 50 |  0.000437 | 0.04310 | 0.88 | 0.1330 |
| 20 | rare | ms | det_(Intercept)[sp1] | 50 | -0.002100 | 0.27300 | 0.94 | 1.0900 |
| 20 | rare | ms | fp_gap[sp1] | 50 | -0.048300 | 0.34600 | 0.96 | 1.3100 |
| 20 | rare | ms | occ_(Intercept)[sp1] | 50 |  0.206000 | 0.45400 | 0.94 | 1.7900 |
| 20 | rare | ms | occ_x[sp1] | 50 |  0.108000 | 0.34000 | 0.94 | 1.2900 |
| 20 | rare | ms | p10_ref[sp1] | 50 |  0.000103 | 0.00697 | 0.92 | 0.0272 |
| 20 | rare | ms | p11_ref[sp1] | 50 |  0.000979 | 0.06380 | 0.94 | 0.2530 |
| 20 | rare | single | det_(Intercept)[sp1] | 50 | -0.211000 | 0.50700 | 0.92 | 1.3900 |
| 20 | rare | single | fp_gap[sp1] | 50 | -0.233000 | 0.54400 | 0.92 | 1.6100 |
| 20 | rare | single | occ_(Intercept)[sp1] | 50 |  0.183000 | 0.55100 | 0.96 | 2.4100 |
| 20 | rare | single | occ_x[sp1] | 50 |  0.111000 | 0.47500 | 0.94 | 2.0500 |
| 20 | rare | single | p10_ref[sp1] | 50 | -0.000201 | 0.00731 | 0.92 | 0.0306 |
| 20 | rare | single | p11_ref[sp1] | 50 | -0.040200 | 0.09890 | 0.92 | 0.2920 |
| 20 | unverified | ms | det_(Intercept)[sp1] | 50 | -0.010900 | 0.15900 | 0.88 | 0.5180 |
| 20 | unverified | ms | fp_gap[sp1] | 50 |  0.073100 | 0.44400 | 0.92 | 1.3000 |
| 20 | unverified | ms | occ_(Intercept)[sp1] | 50 |  0.033200 | 0.25700 | 0.96 | 1.0800 |
| 20 | unverified | ms | occ_x[sp1] | 50 | -0.004550 | 0.25200 | 0.96 | 0.9870 |
| 20 | unverified | ms | p10_ref[sp1] | 50 | -0.003790 | 0.01730 | 0.90 | 0.0534 |
| 20 | unverified | ms | p11_ref[sp1] | 50 | -0.001700 | 0.03680 | 0.88 | 0.1200 |
| 20 | unverified | single | det_(Intercept)[sp1] | 50 | -0.030800 | 0.19900 | 0.86 | 0.5770 |
| 20 | unverified | single | fp_gap[sp1] | 50 |  0.028800 | 0.39800 | 0.92 | 1.2800 |
| 20 | unverified | single | occ_(Intercept)[sp1] | 50 |  0.049200 | 0.29000 | 0.96 | 1.2600 |
| 20 | unverified | single | occ_x[sp1] | 50 | -0.084300 | 0.33000 | 0.94 | 1.2500 |
| 20 | unverified | single | p10_ref[sp1] | 50 | -0.002200 | 0.01730 | 0.92 | 0.0577 |
| 20 | unverified | single | p11_ref[sp1] | 50 | -0.006170 | 0.04520 | 0.86 | 0.1320 |

## fs and mixed

| scenario | model | param | n | bias | rmse | coverage | width |
| --- | --- | --- | --- | --- | --- | --- | --- |
| fs | fs | det_(Intercept) |  50 |  0.017100 | 0.1580 | 1.000 | 0.7160 |
| fs | fs | det_(Intercept)[sp] | 500 | -0.008550 | 0.1230 | 0.958 | 0.4860 |
| fs | fs | det_bio |  50 | -0.003040 | 0.0341 | 0.960 | 0.1440 |
| fs | fs | fp_(Intercept) |  50 |  0.007930 | 0.3490 | 0.940 | 1.4700 |
| fs | fs | fp_gap |  50 |  0.009160 | 0.3160 | 0.960 | 1.2900 |
| fs | fs | fp_gap[sp] | 500 | -0.007560 | 0.3530 | 0.950 | 1.3600 |
| fs | fs | occ_(Intercept) |  50 | -0.018000 | 0.3070 | 0.920 | 1.2200 |
| fs | fs | occ_(Intercept)[sp] | 500 | -0.000275 | 0.2880 | 0.930 | 1.0300 |
| fs | fs | occ_x |  50 | -0.006450 | 0.1760 | 0.980 | 0.8430 |
| fs | fs | occ_x[sp] | 500 | -0.024200 | 0.2410 | 0.958 | 0.9890 |
| fs | fs | p10_ref[sp] | 500 |  0.000793 | 0.0144 | 0.960 | 0.0571 |
| fs | fs | p11_ref[sp] | 500 | -0.002060 | 0.0294 | 0.958 | 0.1140 |
| fs | fs | sd_det_(Intercept) |  50 |  0.033400 | 0.1420 | 0.940 | 0.5890 |
| fs | fs | sd_fp_group |  50 | -0.011500 | 0.1310 | 0.940 | 0.4540 |
| fs | fs | sd_log_fp_gap |  50 |  0.021800 | 0.0992 | 0.960 | 0.4110 |
| fs | fs | sd_occ_(Intercept) |  50 | -0.067800 | 0.2490 | 0.960 | 0.9600 |
| fs | fs | sd_occ_x |  50 |  0.072300 | 0.2020 | 0.900 | 0.7540 |
| fs | fs | verified_(Intercept) |  50 |  0.006940 | 0.1260 | 0.880 | 0.4360 |
| fs | mixed | det_(Intercept) |  50 |  0.019700 | 0.1590 | 1.000 | 0.7140 |
| fs | mixed | det_(Intercept)[sp] | 500 | -0.006190 | 0.1240 | 0.962 | 0.4900 |
| fs | mixed | det_bio |  50 | -0.000211 | 0.0349 | 0.960 | 0.1440 |
| fs | mixed | fp_(Intercept) |  50 | -0.002980 | 0.3620 | 0.940 | 1.4800 |
| fs | mixed | fp_gap |  50 |  0.022700 | 0.3280 | 0.960 | 1.3000 |
| fs | mixed | fp_gap[sp] | 500 |  0.011300 | 0.3580 | 0.958 | 1.3900 |
| fs | mixed | occ_(Intercept) |  50 | -0.029700 | 0.3030 | 0.920 | 1.2300 |
| fs | mixed | occ_(Intercept)[sp] | 500 | -0.014800 | 0.2900 | 0.930 | 1.0400 |
| fs | mixed | occ_x |  50 | -0.010600 | 0.1760 | 0.980 | 0.8390 |
| fs | mixed | occ_x[sp] | 500 | -0.027400 | 0.2430 | 0.960 | 0.9970 |
| fs | mixed | p10_ref[sp] | 500 |  0.000569 | 0.0147 | 0.960 | 0.0578 |
| fs | mixed | p11_ref[sp] | 500 | -0.001510 | 0.0296 | 0.962 | 0.1150 |
| fs | mixed | sd_det_(Intercept) |  50 |  0.034100 | 0.1430 | 0.960 | 0.5940 |
| fs | mixed | sd_fp_group |  50 | -0.040800 | 0.1590 | 0.940 | 0.4950 |
| fs | mixed | sd_log_fp_gap |  50 |  0.026300 | 0.1020 | 0.960 | 0.4170 |
| fs | mixed | sd_occ_(Intercept) |  50 | -0.062600 | 0.2510 | 0.940 | 0.9640 |
| fs | mixed | sd_occ_x |  50 |  0.074100 | 0.2040 | 0.900 | 0.7500 |
| fs | mixed | verified_(Intercept) |  50 |  0.015900 | 0.1270 | 0.920 | 0.4410 |
| mixed | fs | det_(Intercept) |  50 | -0.028800 | 0.1550 | 0.980 | 0.7620 |
| mixed | fs | det_(Intercept)[sp] | 500 | -0.007220 | 0.1310 | 0.950 | 0.4810 |
| mixed | fs | det_bio |  50 | -0.005850 | 0.0436 | 0.900 | 0.1400 |
| mixed | fs | fp_(Intercept) |  50 |  0.024400 | 0.3390 | 0.960 | 1.5900 |
| mixed | fs | fp_gap |  50 | -0.053200 | 0.3240 | 0.960 | 1.4000 |
| mixed | fs | fp_gap[sp] | 500 | -0.009290 | 0.4870 | 0.930 | 1.5800 |
| mixed | fs | occ_(Intercept) |  50 |  0.103000 | 0.3040 | 0.920 | 1.3200 |
| mixed | fs | occ_(Intercept)[sp] | 500 |  0.019300 | 0.2870 | 0.948 | 1.0600 |
| mixed | fs | occ_x |  50 | -0.007810 | 0.1870 | 0.940 | 0.8520 |
| mixed | fs | occ_x[sp] | 500 | -0.002060 | 0.2590 | 0.936 | 1.0100 |
| mixed | fs | p10_ref[sp] | 500 |  0.002740 | 0.0211 | 0.932 | 0.0630 |
| mixed | fs | p11_ref[sp] | 500 | -0.001680 | 0.0305 | 0.950 | 0.1120 |
| mixed | fs | sd_det_(Intercept) |  50 |  0.072000 | 0.1650 | 0.920 | 0.6210 |
| mixed | fs | sd_fp_group |  50 |  0.021600 | 0.1260 | 0.920 | 0.4560 |
| mixed | fs | sd_log_fp_gap |  50 |  0.060600 | 0.1430 | 0.880 | 0.4640 |
| mixed | fs | sd_occ_(Intercept) |  50 |  0.008920 | 0.2730 | 0.900 | 1.0100 |
| mixed | fs | sd_occ_x |  50 |  0.074800 | 0.1850 | 0.960 | 0.7690 |
| mixed | fs | verified_(Intercept) |  50 | -0.015300 | 0.1260 | 0.900 | 0.4330 |
| mixed | mixed | det_(Intercept) |  50 | -0.025800 | 0.1560 | 0.980 | 0.7680 |
| mixed | mixed | det_(Intercept)[sp] | 500 | -0.003210 | 0.1320 | 0.952 | 0.4910 |
| mixed | mixed | det_bio |  50 |  0.002210 | 0.0437 | 0.920 | 0.1450 |
| mixed | mixed | fp_(Intercept) |  50 |  0.006200 | 0.3480 | 0.980 | 1.6300 |
| mixed | mixed | fp_gap |  50 | -0.032000 | 0.3330 | 0.960 | 1.4300 |
| mixed | mixed | fp_gap[sp] | 500 |  0.025500 | 0.5040 | 0.942 | 1.6600 |
| mixed | mixed | occ_(Intercept) |  50 |  0.065400 | 0.2990 | 0.920 | 1.3400 |
| mixed | mixed | occ_(Intercept)[sp] | 500 | -0.024800 | 0.2830 | 0.958 | 1.1000 |
| mixed | mixed | occ_x |  50 | -0.021300 | 0.1970 | 0.940 | 0.8630 |
| mixed | mixed | occ_x[sp] | 500 | -0.014400 | 0.2610 | 0.938 | 1.0300 |
| mixed | mixed | p10_ref[sp] | 500 |  0.002750 | 0.0214 | 0.944 | 0.0654 |
| mixed | mixed | p11_ref[sp] | 500 | -0.000755 | 0.0308 | 0.952 | 0.1140 |
| mixed | mixed | sd_det_(Intercept) |  50 |  0.076900 | 0.1680 | 0.940 | 0.6280 |
| mixed | mixed | sd_fp_group |  50 | -0.028500 | 0.1410 | 0.920 | 0.5580 |
| mixed | mixed | sd_fp_site |  50 | -0.005300 | 0.1590 | 1.000 | 0.7220 |
| mixed | mixed | sd_log_fp_gap |  50 |  0.072000 | 0.1520 | 0.860 | 0.4770 |
| mixed | mixed | sd_occ_(Intercept) |  50 |  0.024600 | 0.2780 | 0.920 | 1.0200 |
| mixed | mixed | sd_occ_x |  50 |  0.082700 | 0.1900 | 0.960 | 0.7820 |
| mixed | mixed | verified_(Intercept) |  50 |  0.006790 | 0.1240 | 0.860 | 0.4390 |

## Simulation-based calibration

Fits: 200, failed: 0, with divergences: 8.

| variable | n | chisq_p | mean_rank |
| --- | --- | --- | --- |
| det_(Intercept) | 200 | 0.2830 | 46.6 |
| det_(Intercept)[sp1] | 200 | 0.9600 | 48.7 |
| det_bio | 200 | 0.3420 | 49.6 |
| fp_gap | 200 | 0.7300 | 46.8 |
| occ_(Intercept) | 200 | 0.0690 | 48.9 |
| occ_(Intercept)[sp1] | 200 | 0.5040 | 49.2 |
| occ_x | 200 | 0.4940 | 47.4 |
| p10_ref[sp1] | 200 | 0.4190 | 50.6 |
| sd_det_(Intercept) | 200 | 0.0192 | 49.7 |
| sd_log_fp_gap | 200 | 0.5750 | 49.3 |
| sd_occ_(Intercept) | 200 | 0.7200 | 46.0 |
| sd_occ_x | 200 | 0.4750 | 49.5 |
| verified_(Intercept) | 200 | 0.4100 | 50.5 |

