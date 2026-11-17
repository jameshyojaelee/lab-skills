# Agent instructions

This repo holds analysis skills. For a matching task, read the skill's `SKILL.md` first and follow it:

- Gel images, band intensity, ligation efficiency, % full length -> `skills/gel-quant/SKILL.md`
- qPCR, standard curve, copies per well or per bead -> `skills/qpcr-quant/SKILL.md`
- TapeStation, library QC, fragment size -> `skills/tapestation-qc/SKILL.md`
- Library quantification by qPCR (KAPA) -> `skills/library-quant/SKILL.md`
- Countess cell counts, viability, fluorescent fraction, seeding -> `skills/cell-count/SKILL.md`
- Lab notebook entries -> `skills/lab-notebook-entry/SKILL.md`
- Team chat updates -> `skills/slack-update/SKILL.md`

Reporting rules for every skill: `_shared/CONVENTIONS.md`.
Figure rules for every figure, in any skill or ad hoc: `_shared/FIGURES.md`; draw with `_shared/labstyle.py`
(`labstyle.figure()` to create, `labstyle.save()` to write PDF + PNG + caption and check the rules).

Rules:
- Run the scripts with a config; do not rewrite the analysis inline. If the script cannot do
  what is needed, extend it, add a test, and keep project-specific values in configs.
- Run `python3 tests/test_skills.py` after any script change. Results for the examples must not change
  unless the change is intended and the expected files are updated in the same commit.
- Report with the QC and reporting sections of the SKILL.md: define every metric, list controls
  present and missing, mark bounds from clipping, label claims verified / inferred / estimate.
- Every figure goes through `labstyle.save()`; fix every FIGURE WARNING before handing a figure over, then
  open the PNG and look at it. Captions go in the notebook or chat, not in the image.
- Plain language, no em dashes.
- Keep this repo generic: no project names, sample names, people, database or channel IDs. Those go in
  the private config (`$LAB_SKILLS_CONFIG`, else `~/.config/lab-skills/config.yaml`; see `config/README.md`). Examples use synthetic data.
