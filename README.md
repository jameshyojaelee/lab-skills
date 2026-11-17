# lab-skills

Reusable, instrument-level analysis skills for wet-lab data.
Each skill is a folder with a `SKILL.md` (when to use, inputs, steps, QC, how to report),
a script, a config template, a worked example and notes on past pitfalls. The same files
work for Claude Code, Codex, Gemini and people.

| Skill | What it does |
| --- | --- |
| [`gel-quant`](skills/gel-quant/SKILL.md) | Gel band quantification: band areas, step and overall ligation efficiency, relative amounts, pmol from standards, clipping bounds |
| [`qpcr-quant`](skills/qpcr-quant/SKILL.md) | Absolute qPCR: standard curve, Ct to copies, melt QC, single-bead wells with Poisson check, bulk copies per bead |
| [`tapestation-qc`](skills/tapestation-qc/SKILL.md) | Agilent TapeStation: peak and region size, molarity, dimer fraction, warnings, ladder-sized traces, labelled gel |
| [`cell-count`](skills/cell-count/SKILL.md) | Countess 3 / 3 FL: concentration, viability, fluorescent % with CIs, replicate summary, seeding volumes |
| [`lab-notebook-entry`](skills/lab-notebook-entry/SKILL.md) | Bare-minimum electronic notebook entries (destination from private config) |
| [`slack-update`](skills/slack-update/SKILL.md) | Short team updates from notebook entries, posted only after approval |

All skills follow [`_shared/CONVENTIONS.md`](_shared/CONVENTIONS.md). Destinations such as notebook
databases and chat channels live in a private config outside the repo ([`config/`](config/README.md)).

## Layout

```
skills/<name>/
  SKILL.md              instructions (agents read this first)
  scripts/              runnable code, config-driven
  configs/template.yaml start here for a new experiment
  examples/<expt>/      real config + data + expected output
  references/           method.md (why), lessons.md (pitfalls from past runs)
tests/                  tests on synthetic data with known answers
_shared/CONVENTIONS.md  reporting rules all skills follow
config/                 how private destinations are configured (example only)
```

## Setup

```
pip install -r requirements.txt
python3 tests/test_skills.py        # or: python3 -m pytest -q tests
```

## Use from agents

Link each skill into the tool's skill folder so every tool reads this one copy:

```
for d in ~/.claude/skills ~/.codex/skills ~/.agents/skills; do
  mkdir -p "$d"
  for s in "$PWD"/skills/*; do ln -sfn "$s" "$d/$(basename "$s")"; done
done
```

Agents without skill support: point them at `AGENTS.md`.

## Quick start

```
cp skills/gel-quant/configs/template.yaml /path/to/gel/config.yaml   # edit lanes, bands, metrics
python3 skills/gel-quant/scripts/gel_quant.py /path/to/gel/config.yaml

cp skills/qpcr-quant/configs/template.yaml /path/to/run/config.yaml  # edit plate layout
python3 skills/qpcr-quant/scripts/qpcr_quant.py /path/to/run/config.yaml --export run.xls
```

## Contributing

After a new experiment: add lessons to `references/lessons.md`, add a worked example under
`examples/` with its expected output, and extend `tests/`. Keep experiment-specific values in
configs kept outside this repo; examples here use generic or synthetic data.
