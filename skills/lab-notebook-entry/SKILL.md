---
name: lab-notebook-entry
description: Create or update an electronic lab notebook entry (Notion database, Obsidian or markdown) for bench, computational, animal, meeting, logistics or reading work, in a bare-minimum format. Use when asked to log, write up or update a notebook entry.
---

# Lab notebook entries

The destination (database, people, icons, section names) comes from the private config
(`$LAB_SKILLS_CONFIG`, else `~/.config/lab-skills/config.yaml`; see `config/README.md`). If neither exists, ask where the entry goes.

## Scope: bare minimum

An entry records what to do and what happened.

In the entry: experiment design; master mix and reagent tables (5% overage unless stated); plate or
tube maps; steps, cycling conditions, cleanup volumes; results as recorded, including instrument values
and figures; one line on how results were analysed (name the skill and config).

Not in the entry: interpretation, hypotheses, candidate-cause or decision tables, warnings and cautions,
a deviations section, rationale, comparisons between protocols. Those go to the discussion or meeting page.

## Properties and body

- Title `YYYYMMDD short description`; date = date the work was done.
- Performer as a user mention, not plain text, where the system supports it.
- Status "in progress" until there is an outcome; brief results in plain text, facts and numbers only.
- Body sections from the config (default: Summary, Notes, Endpoint). Multi-day work uses dated
  `Day N (date)` headings with `TO FILL` for days not yet done.
- Link the entry to its task or project record when the system has one.

## Data integrity

Transcribe instrument output as recorded; mark unresolved values `TO FILL`; record reagent lots when
asked; say "undiluted", never "neat".

## Workflow

1. Confirm date, activity type and performer before creating the entry.
2. Create with properties set in the same call; re-read before editing (people edit entries between turns).
3. Edit with small targeted replacements; for many sections, replace the whole body and re-read to
   confirm nothing was dropped. Images in Notion: upload via the file-upload API, then insert.
4. Write only when the user asks in that message.
