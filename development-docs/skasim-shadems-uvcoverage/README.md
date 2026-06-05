# shadeMS UV-Coverage Weblog Plot

## Goal

Generate a UV-coverage PNG from each run's Measurement Set using `shadems`, record it in the run manifest, and display it beside the telescope layout in the weblog.

## Verified command shape

The quick verification used:

```bash
shadems <measurement-set> \
  --xaxis u \
  --yaxis v \
  --xcanvas 600 \
  --ycanvas 600 \
  --spread-pix 2 \
  --dir <work-dir> \
  --png <run-id>_uvcoverage.png \
  --title "<run-id> uv coverage" \
  --xlabel u \
  --ylabel v \
  --no-lim-save
```

shadeMS does not expose an equal-aspect CLI switch. A square canvas is used here because U and V share units, but exact equal data aspect would require symmetric explicit U/V limits or a shadeMS renderer change.
