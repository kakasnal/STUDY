# LIN28A upgrade v3 submission

## Topic

LIN28A에 강하게 결합한 mRNA가 Lin28a knockdown 후 ribosome density 증가를 보이는지 확인하고, ER/membrane/secretory pathway 관련 후보 direct target을 우선순위화합니다.

## Why this is an upgrade

이전 과제의 `CLIP enrichment` vs `ribosome density change` 산점도를 유지하면서, 아래를 추가했습니다.

1. `read-counts.txt` 자동 탐색 및 없을 경우 `featureCounts`로 재생성 시도
2. CPM normalization
3. 여러 CLIP replicate 자동 탐지 및 median CLIP enrichment 계산
4. Spearman correlation + bootstrap 95% confidence interval
5. strong binder vs weak binder Mann-Whitney U test + permutation test
6. threshold sensitivity analysis
7. mouse localization annotation을 이용한 ER/membrane/secretory 후보 direct target ranking
8. 선택적으로 2주차 RPF start-codon profile, 3주차 let-7 CLIP entropy 분석 재실행

## How to run

1. `LIN28A_upgrade_v3_submission.ipynb`를 Colab에서 엽니다.
2. Google Drive에 이전 과제의 `binfo1-datapack1` 폴더 또는 `read-counts.txt`가 있는지 확인합니다.
3. `Runtime -> Run all`을 실행합니다.
4. 마지막 checklist cell에서 `PASS`를 확인합니다.
5. 생성된 `LIN28A_upgrade_v3/figures/`, `LIN28A_upgrade_v3/results/` 폴더와 노트북을 GitHub에 push합니다.

## Required output files

- `results/lin28a_filtered_gene_metrics.csv`
- `results/lin28a_candidate_direct_targets_ranked.csv`
- `results/correlation_bootstrap_summary.csv`
- `results/strong_vs_weak_binding_test.csv`
- `results/threshold_sensitivity_analysis.csv`
- `results/summary.txt`
- `figures/fig1_clip_vs_rden_scatter_binned.png`
- `figures/fig2_rden_by_binding_group_boxplot.png`
- `figures/fig3_threshold_sensitivity.png`

## Submission warning

노트북 파일만 push하고 실행 결과가 없는 상태면 감점 또는 0점이 될 수 있습니다. 반드시 `Run all` 후 output이 남아 있는 `.ipynb`를 저장해서 제출하세요.
