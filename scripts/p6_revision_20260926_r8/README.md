# Revision support scripts

This directory contains the saved-table additions used by repository version
1.1.0.

- `render_workflow.py` renders current Figure 1 from supplied processed tables.
- `render_current_figures.py` preserves the earlier data-figure renderers and
  maps their Figure 1-8 outputs to current Figure 2-9; Figure S1 and S2 keep
  their names.
- `summarize_rule_comparison.py` reproduces the four descriptive Table 1 rows
  (61, 40, 48, and 35 nonredundant pairs) from saved P2 objects.
- `define_equations.py` regenerates `equations.json`, the MathML definitions
  for the six main and three supplementary equations.
- `native_equations.py` is an optional Windows Word helper. It is outside the
  offline reproduction route and requires `python-docx`, `lxml`, and a local
  Microsoft Office `MML2OMML.XSL` transform.

From the repository root:

```bash
python scripts/p6_revision_20260926_r8/define_equations.py
python scripts/p6_revision_20260926_r8/summarize_rule_comparison.py \
  --project-root . --output-dir results/rule-comparison
python scripts/p6_revision_20260926_r8/render_current_figures.py \
  --project-root . --output-dir results/figures --figures all
```

The equation helper provides `append_equation(document, identifier,
xsl_path=...)`; it does not generate or redistribute a manuscript.
