# Stage-3 Laughney label-confirmation policy (pre-registered addendum item 5, fixed BEFORE use)

**Status:** policy document, written 2026-08-21 before `stage3_confirm_labels.py` was run and before
any cross-sample switch statistic exists. Implements `manuscript/13_reliability_positioning.md`,
Addendum 2026-08-21 05:30, items 5 and 6. Every number that the implementation later produces is
PROVISIONAL until its verifier passes. The only label-side inspection that preceded this table was
(i) the vocabulary of the two CellTypist models and which labels actually occur in the per-cell file,
(ii) the curated-label x CellTypist-label cross-tabulation (needed to know which labels exist), and
(iii) per-GSM marker-gene availability in the GEX objects. No switch-side agreement number was computed.

## 0. Inputs (read-only; nothing under /mnt/ssd2 is modified)

| input | path | notes |
|---|---|---|
| curated labels | `/mnt/ssd2/Laugney_Aligned/laughney_ref/laughney_celltypes_curated.parquet` | 29,063 cells, columns `cell` (= full id `<GSM>_<16nt>`), `laughney_celltype` (13 classes, argmax of `score_genes` over `curated_markers.json`, no "unknown" by construction), `score` |
| CellTypist per-cell calls | `/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed/runs/B1_cohort_full/B2_celltypist/percell_labels.csv` | produced by `/mnt/ssd2/Laugney_Aligned/celltypist_work/annotate.py` (CellTypist 1.7.1): GEX `layers['counts']` -> `normalize_total(1e4)` -> `log1p`; **`new_label`** = `Human_Lung_Atlas.pkl` (v2, 2023-05-17) majority-voting label, **`conf`** = `conf_score` of that majority-voting label (per cell), labels with `conf < 0.5` already written as `unknown`; **`imm_label`** = `Immune_All_Low.pkl` (v2, 2022-07-16) majority-voting label, **no confidence recorded**; column `barcode` = full id |
| GEX objects | `.../runs/B1_cohort_full/B2_gex_celltyping/<GSM>_gex_labeled.h5ad` (17 files) | `layers['counts']` integral UMI counts; obs_names = full ids; genes were filtered upstream with `sc.pp.filter_genes(min_cells=3)` per GSM, so a marker gene can be **absent** from a GSM's var (see 3.4) |
| sample table | `scripts/stage3/laughney_samples.tsv` | `dataset_id`, `patient` (17 GSMs = 14 patients), `group` (Normal / StageI / IVprimary / Met) |

**Join rule (binding):** every join is on the FULL cell id `'<GSM>_<16nt barcode>'` (parquet `cell` ==
csv `barcode` == h5ad `obs_names`). Never on the bare 16-nt barcode (1,070 cells collide across GSMs).
All three sources cover exactly the same 29,063 ids (checked before writing this document).

## 1. Which CellTypist call is "the CellTypist call"

The per-cell file carries two model outputs. The pre-registered rule says *"agreement with the CellTypist
call (`percell_labels.csv`, conf >= 0.5)"*. The only call in the file that has a confidence is the
Human_Lung_Atlas majority-voting label (`new_label`/`conf`), so:

* **Primary call = `new_label` (Human_Lung_Atlas), confidence = `conf`.** Confirmation requires
  `conf >= 0.5` (this coincides with the file's own `unknown` threshold; cells with `new_label == unknown`
  therefore never confirm).
* **`imm_label` (Immune_All_Low) is used for exactly one thing:** the Treg-vs-T_cell split, because the
  Human_Lung_Atlas model has no regulatory-T class (its T classes are `CD4 T cells`, `CD8 T cells`,
  `T cells proliferating`). It carries no confidence, so it is never used alone to confirm a cell.
* `imm_label` agreement with the curated class is additionally **reported** per cell (`imm_agree`) and per
  type as a diagnostic; it does NOT enter `confirmed`/`final_label`.

## 2. Immune types: fixed mapping table (curated label <- accepted Human_Lung_Atlas labels)

A cell of curated class C is confirmed iff `conf >= 0.5` AND `new_label` is in the accepted set for C
AND the extra condition (Treg/T_cell only) holds. Labels marked * exist in the model vocabulary but were
not observed in this cohort; they are listed so the table is complete and does not change if the file is
regenerated.

| curated | accepted `new_label` (Human_Lung_Atlas) | extra condition on `imm_label` (Immune_All_Low) |
|---|---|---|
| B_cell | `B cells` | none |
| Plasma_IG | `Plasma cells` | none (the model has no Plasmablast class) |
| Monocyte | `Classical monocytes`, `Non-classical monocytes` | none |
| Macrophage | `Alveolar macrophages`, `Alveolar Mph CCL3+`, `Alveolar Mph MT-positive`*, `Alveolar Mph proliferating`, `Monocyte-derived Mph`, `Interstitial Mph perivascular` | none |
| Dendritic | `DC1`, `DC2`, `Migratory DCs`, `Plasmacytoid DCs` | none |
| NK | `NK cells` | none |
| T_cell | `CD4 T cells`, `CD8 T cells`, `T cells proliferating` | `imm_label` NOT in {`Regulatory T cells`, `Treg(diff)`*} (a cell the immune model calls regulatory is not confirmed as generic T_cell; it is dropped, not moved to Treg) |
| Treg | `CD4 T cells`, `T cells proliferating` | `imm_label` IN {`Regulatory T cells`, `Treg(diff)`*} |
| Mast | `Mast cells` | none |

Explicitly NOT accepted (would be relabelling, which the pre-registration forbids): `Monocyte-derived Mph`
for Monocyte; `DC2` for B_cell or Monocyte; `CD8 T cells` for NK (the cytotoxic-T/NK ambiguity is the
known weak spot of the curated NK class); anything for a curated class outside the table.

## 3. Non-immune types: canonical-marker gate (Epithelial_Tumor, Fibroblast, Endothelial, Pericyte)

CellTypist is ignored for these four classes (Human_Lung_Atlas sub-types and the immune model's
`Epithelial cells` / `Fibroblasts` / `Endothelial cells` catch-alls are reported as diagnostics only).

### 3.1 Marker sets (pre-registered, verbatim)

| curated | markers |
|---|---|
| Epithelial_Tumor | EPCAM, KRT8, KRT18 |
| Fibroblast | COL1A2, DCN, LUM |
| Endothelial | PECAM1, VWF, CLDN5 |
| Pericyte | RGS5, ACTA2, PDGFRB |

### 3.2 Score

Per GSM, from `layers['counts']` (never the stored `.X`): `sc.pp.normalize_total(target_sum=1e4)` ->
`sc.pp.log1p` -> `sc.tl.score_genes(gene_list=<markers present in this GSM's var>, ctrl_size=50,
n_bins=25, random_state=0)` for each of the four sets, scored for EVERY cell of the GSM (so the four
scores are comparable within a cell). scanpy 1.11.1 (`tools/PeakATail/.venv`).

### 3.3 Gate (both conditions, evaluated on the cell's own curated class T)

* **A (argmax):** `score_T` is finite AND `score_T > max(score_U for U != T, U finite)` -- strictly greater.
  An exact tie with any other class fails A. If every other class is NaN (3.4), A reduces to `score_T` finite.
* **B (floor):** `score_T >= P25_T`, where `P25_T = numpy.percentile(x, 25)` (linear interpolation) over
  `x` = the `score_T` values of ALL cells in the cohort (all 17 GSMs pooled) whose curated class is T AND
  that pass A. **"Confirmed" for the percentile is therefore defined as "own-type score is highest"** (A),
  exactly as the addendum says; P25 is computed once, cohort-wide, and recorded in the manifest. If no cell
  of class T passes A anywhere, P25_T is undefined and nobody of class T is confirmed.
* `confirmed = A and B`. Reason codes: `ok`, `own_markers_absent` (3.4), `not_argmax` (incl. ties),
  `below_p25`.

### 3.4 NaN / absent-gene handling

The GEX objects are gene-filtered per GSM (`min_cells=3`), so a marker can be missing from a GSM
(observed: CLDN5 and/or VWF missing in 7 GSMs; COL1A2+DCN+LUM all missing in GSM3516675-Normal).
Rules: a set is scored on the markers that are present (1-3 genes); if **none** of its three markers is
present in a GSM, that set's score is NaN for every cell of that GSM. A NaN score never wins an argmax and
never blocks another class (it is skipped in `max(score_U)`); a cell whose own class scores NaN is not
confirmed (`own_markers_absent`). NaN scores are excluded from the P25 pool. Absent-marker status per
GSM x set is written to the manifest. (A missing gene means < 3 expressing cells in that GSM, i.e. the set
would have scored ~0 anyway; treating it as NaN is the conservative choice.)

## 4. Minimum cell count and what enters the switch test

* Unconfirmed cells are **dropped, never relabelled** (`final_label` empty).
* A (GSM, `final_label`) group with **< 20 confirmed cells is dropped for that GSM** (`retained = False` in
  `per_gsm_retained.tsv`; `stage3_label_inject.py --min-cells 20` applies the same rule at injection).
* Per-GSM retained counts per type (`per_gsm_retained.tsv`) are written and reported before any switch
  statistic is viewed (addendum item 5, last sentence).

## 5. Item 6: Epithelial in Normal GSMs, and pair naming

* For cells whose curated class is `Epithelial_Tumor` in a GSM with `group == Normal`
  (GSM3516666, GSM3516673, GSM3516675, GSM3516676), `final_label = Epithelial`; in every other GSM it
  stays `Epithelial_Tumor`. The marker gate (3.3) is computed on the curated class, so the P25 floor is one
  cohort-wide value for `Epithelial_Tumor` regardless of Normal/tumour origin; the relabel is applied after
  the gate, for reporting and for pair naming.
* **A switch pair is defined on `final_label`.** Therefore `Epithelial` (Normal GSMs) and
  `Epithelial_Tumor` (tumour GSMs) are DIFFERENT labels: a pair such as `Epithelial_vs_T_cell` can only
  arise in Normal GSMs and `Epithelial_Tumor_vs_T_cell` only in tumour GSMs, and the two are never
  pooled when counting replication across patients. Consequently any pair involving `Epithelial` has at
  most 4 patients (LX675, LX682, LX684, LX685) and any pair involving `Epithelial_Tumor` at most 13.

## 6. Outputs of `stage3_confirm_labels.py` (all under `/mnt/ssd0/emaout/peakatail_benchmark/stage3_laughney_v2/labels/`)

| file | content |
|---|---|
| `confirmed_labels.tsv` | one row per cell: `cell` (full id), `gsm`, `patient`, `group`, `curated`, `celltypist` (= `new_label`), `conf`, `imm_label`, `imm_agree`, `confirmed` (True/False), `reason`, `final_label`, and for the four non-immune classes the four marker scores + `p25_threshold` |
| `per_gsm_retained.tsv` | GSM x final_label: `curated_n`, `confirmed_n`, `fraction`, `retained` (>= 20), `patient`, `group` |
| `per_type_summary.tsv` | cohort totals per final_label: curated n, confirmed n, fraction, GSMs retained, patients retained |
| `purity_check.tsv` | per (GSM, final_label) with >= 20 confirmed cells: CP10k median and fraction-expressing of the check genes (MS4A1/CD79A B_cell; CD3E T_cell & Treg; LYZ Monocyte/Macrophage; NKG7/GNLY NK; JCHAIN Plasma_IG; TPSAB1 Mast; EPCAM Epithelial*; COL1A2 Fibroblast; PECAM1 Endothelial; plus RGS5 Pericyte / FOXP3 Treg as extras) computed for ALL check genes in every group, curated set ("before") vs confirmed set ("after") |
| `manifest.json` | input paths + mtimes + md5, scanpy/anndata versions, P25 thresholds, absent-marker table, all counts, `status: PROVISIONAL` |

## 7. Things this policy deliberately does not do

* It does not tune any threshold on the data (0.5 and P25 are the pre-registered values).
* It does not use the curated `score` column.
* It does not re-cluster, re-annotate or re-run CellTypist; it only reads the existing calls.
* It does not change the pre-registered min-cells (20) or the patient definition.
