# Figure rules (all skills)

Every figure is made with `_shared/labstyle.py`: `labstyle.figure()` to create it, `labstyle.save()` to
write it. `save()` checks the rules marked (checked) below and prints or raises on violations.

## One message per figure
- A figure answers one question. Two questions, two figures.
- Raw data image (gel, trace, micrograph) and summary plots are separate figures.
- Panel titles are short labels (an arm, a sample, a metric), not findings and not methods.
- No figure-wide title on a single panel (checked). The finding goes in the chat or meeting text.

## Text stays out of the image
- Methods, formulas, caveats, acquisition settings and provenance go in the caption
  (`save(..., caption=...)` writes `<stem>.caption.txt`), never inside the figure (checked: no multi-line
  or long figure text).
- Value labels only where they replace an axis reading, at most one per bar or point, and never on top
  of a data point or another label (checked).
- Legends: only when direct labels do not fit; at most 6 entries (checked). Explain marker styles
  (filled / open) in the caption when there is a legend already.
- Axis labels name the quantity and unit: "Full length (%)", "Live cells (per ml)". Shared axes are
  labelled once.
- Tick labels are single-line and short: "Day 2", "1 µM", "0.6". No decimals that carry no meaning.

## Size and type (from the lab publication theme)
- Sizes (inches): single panel 3.5 x 2.6; two or three panels in a row 7.0 x 2.4-2.6; 2 x 2 7.0 x 4.8;
  max 7.2 x 9 (checked).
- Arial / Helvetica (Liberation Sans as fallback), 7 pt text, 8 pt titles, 6-6.5 pt ticks; nothing under
  6 pt (checked).
- Output: vector PDF with editable text (pdf.fonttype 42) is the master; a 300-dpi PNG is written for
  Notion and Slack previews only.

## Colour
- Stages and amounts: the blue family `STAGES` (dark = furthest along / most).
- Conditions: `COLORS` in fixed order (blue, orange, green, purple, grey). Colour follows the entity: the
  same condition keeps its colour across figures and filters.
- Reference or target lines: dashed grey with a small right-aligned label.
- White background, no grid (light horizontal rules only when reading values matters), no 3D, no
  shadows, no colour-only encoding (add marker shape or a label).

## Chart type by question
| Question | Chart |
| --- | --- |
| Efficiency or yield per condition | grouped bars, target line, value inside bar top |
| Amount per sample, with composition | stacked bars on an absolute scale (not 100%) |
| Change over time | lines with points and CI error bars; no value labels; repeats as separate open markers |
| Distribution / replicates | strip plot with median line |
| Standard curve | points + fit line, efficiency and R² in one small corner note |
| Gel or trace | the image with lane labels, group brackets and size ticks; nothing else |

## Before handing a figure over
- Open the PNG and look at it: overlaps, clipped labels, empty space, legibility at Notion width.
- The caption states: what each mark is, the metric formula, the interval type, n, and what is excluded.
