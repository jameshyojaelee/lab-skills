---
name: slack-update
description: Draft and, after approval, post a short team chat update (Slack) from a lab notebook entry or analysis result, routed to the right channel. Use when asked to post, share or give an update on lab work.
---

# Chat update from a notebook entry

Channels and member IDs come from the private config (`$LAB_SKILLS_CONFIG`, else
`~/.config/lab-skills/config.yaml`). If a topic has no
channel there, ask which channel to use.

## Rules

- Never post without the user approving the exact draft in this conversation.
- Never edit the notebook from this skill; if the entry needs changes, say so.
- Only numbers that appear in the entry or analysis output; no recomputing or re-rounding.
- Mark plan items with no recorded decision as "(proposed)".
- Plain English, short lines, no em dashes, no emojis.

## Draft format (150 words or fewer)

```
**<Experiment>: update (<MM-DD>)**

**Status:** one line.
**Key numbers:** 2-4 bullets, each with what the number is (metric definition).
**Controls:** what was and was not included, if results are quantitative.
**Plan:** next steps with dates.
**Asks:** @person and the specific request.
```

Leave out empty sections. Show the draft, the target channel and the source link; post after approval;
report the message link. For follow-ups on the same experiment, ask whether to reply in the thread.
