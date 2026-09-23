# Station reference ET reports

Point this at any AgriMet station and a date range and it builds a
self-contained, interactive HTML report of daily and hourly ASCE reference ET.

```console
pip install -e .                  # refet reads its version via importlib.metadata
pip install pandas pytz
python examples/make_report.py examples/stations/faln_2015.toml
```

| File | What it is |
|---|---|
| `make_report.py` | computes the data, derives the facts, renders the template |
| `agrimet.py` | client for the USBR AgriMet/Hydromet web service |
| `report_template.html` | all the chart code, styling and prose |
| `stations/*.toml` | one config per report |

Two configs ship as worked examples:

- **`faln_2015.toml`** — Fallon, NV for 2015 from the CSVs in `tests/data/`,
  validated against the Ref-ET `.out` files beside them.
- **`faln_2024_live.toml`** — the same station fetched live from AgriMet for a
  six-month window, with no reference data.

## Reporting a new station

Copy `faln_2024_live.toml`, change the station id and the dates:

```toml
start = 2024-04-01
end   = 2024-09-30

[source]
type = "agrimet"
station = "ABEI"
cache = "examples/.cache"
```

Latitude, longitude, elevation, timezone and the station's display name are read
from the published AgriMet station list, so they don't belong in the config
unless you want to override them. Responses are cached under `[source] cache`,
so re-running a build doesn't re-hit the service. Date ranges may span any
period, including multiple years; the axes, tick spacing, monthly table and
narrative all adapt.

Two things the lookup can't give you:

- **Anemometer height** defaults to **3 m**, the AgriMet network standard
  ([sensors.html](https://www.usbr.gov/pn/agrimet/aginfo/sensors.html)). It is
  *not* published per station, and it materially changes ETr, so override it in
  `[station] zw` if you know a site differs.
- **Great Plains stations** aren't served by the `pn-bin` endpoints this client
  uses; it raises a clear error rather than returning empty data.

## How the narrative stays honest

This is the part worth understanding before editing anything.

Every number that appears in the page's prose is a **derived fact**, computed in
`derive_facts()`. The template writes them as `{{placeholder}}`; page scripts
read them as `F_.name`. Nothing is rounded or reworded in the HTML, so the page
and the computation cannot drift apart.

The build **fails** rather than producing a quietly wrong report when:

- the template uses a placeholder no fact defines, or
- a derived fact is referenced nowhere, or
- a fact can't be derived from the data at all.

The prose *around* those numbers is hand-authored and stays yours. But claims
that are only true for some stations are wrapped in a conditional section, so
they simply don't appear when the data doesn't support them:

```html
The band is widest at {{summer_open}}.<!--if night_et--> A net
{{night_share}} of ET falls overnight.<!--endif-->
```

Flags are set in `build_payload()`. `validation` drops the whole comparison
panel when no reference file is configured; `night_et` gates the overnight-ET
claim.

**This is not hypothetical.** The overnight claim was originally written by hand
from looking at the Fallon 2015 heatmap, which shows obvious streaks after dark.
Across the full year that sum is actually **−1.75 mm** — net negative, with
2,152 of 3,284 night hours below zero. The gate suppresses the claim for that
record and shows it for Fallon Apr–Sep 2024, where the net share really is 2%.

## Two data notes

**The diurnal heatmap is binned on standard time, not clock time.** AgriMet
timestamps are local with DST, so binning on the clock jumps an hour at each
transition and smears the seasonal sunrise/sunset envelope.

**Hourly agreement is weaker at night.** For Fallon 2015, hourly ETr RMSE against
Ref-ET is ~0.022 mm with a max deviation of ~0.60 mm, and the largest gaps fall
in low-radiation hours — the ones `tests/test_hourly_refet_output.py` excludes
with its `RS <= 1.0` filter. Daylight hours and both daily surfaces agree to
well under 0.05 mm.

## Missing data

AgriMet spells gaps at least five ways: `NO RECORD`, an empty field, `m`, the
`998877` sentinel, and values carrying a trailing quality flag (`-459.69-`,
`57.86^`). `agrimet.clean_value()` handles all of them and range-checks the
result. Hourly rows are *omitted* rather than flagged, so gaps appear as blanks
in the series and grey columns in the heatmap.
