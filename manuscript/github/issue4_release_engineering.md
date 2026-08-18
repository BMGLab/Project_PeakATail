For a Q1 methods submission (Genome Biology / Genome Research / Bioinformatics tier), reviewers must be able to install the tool in one step and cite an archived release. Current state and the gap:

- [ ] **CI running the test suite** — `.github/workflows/` currently contains only `docs.yml` (docs deploy). The ~800-test pytest suite runs nowhere automatically. Add a workflow: matrix over Python 3.11/3.12, `pip install -e .[test]`, `pytest`, on PRs to `develop`.
- [ ] **PyPI package** — `pip index versions peakatail` → not found. `pyproject.toml` is already package-shaped (v0.2.0); needs a build + upload flow (ideally a tag-triggered release workflow).
- [ ] **bioconda recipe** (after PyPI) — the audience installs via conda.
- [ ] **CITATION.cff** — none exists; the citation block in the README has a placeholder author.
- [ ] **Zenodo DOI** — enable the GitHub–Zenodo integration and cut `v0.3.0` when the benchmark changes land; the paper's Availability section needs the DOI.
- [ ] **README command table** — lists 8 of the 12 actual commands (missing `reannotate`, `switch trend`, `switch combine`, `collapse`).
- [ ] **License/availability statement** consistency check for the journal's software policy (MIT ✓).

None of this blocks analyses; all of it blocks submission. Estimated ~1 week total, parallelizable with the science.
