# Station reference ET reports

Point `refet-report` at any AgriMet station and a date range and it builds a
self-contained, interactive HTML report of daily and hourly ASCE reference ET.

```console
pip install refet[report]
refet-report examples/stations/faln_2015.toml
```

The report builder lives in the package (`refet.report`); the AgriMet client
is `refet.io.agrimet`. This directory holds worked example configs:

- **`stations/faln_2015.toml`**: Fallon, NV for 2015 from the CSVs in
  `tests/data/`, validated against the Ref-ET `.out` files beside them.
- **`stations/faln_2024_live.toml`**: the same station fetched live from
  AgriMet for a six-month window, with no reference data.

## Reporting a new station

Copy `faln_2024_live.toml`, change the station id and the dates:

```toml
start = 2024-04-01
end   = 2024-09-30

[source]
type = "agrimet"
station = "ABEI"
cache = ".cache"
```

Latitude, longitude, elevation, timezone and the station's display name are
read from the published AgriMet station list, so they don't belong in the
config unless you want to override them. Relative paths in a config resolve
against the config file's directory. Date ranges may span any period,
including multiple years; the axes, tick spacing, monthly table and narrative
all adapt.

Two things the lookup can't give you:

- **Anemometer height** defaults to **3 m**, the AgriMet network standard
  ([sensors.html](https://www.usbr.gov/pn/agrimet/aginfo/sensors.html)). It
  is *not* published per station, and it materially changes ETr, so override
  it in `[station] zw` if you know a site differs.
- **Great Plains stations** aren't served by the `pn-bin` endpoints this
  client uses; it raises a clear error rather than returning empty data.

## How the narrative stays honest

Every number that appears in the page's prose is a **derived fact**, computed
in `refet.report.facts.derive_facts()`. The template writes them as
`{{placeholder}}`; page scripts read them as `F_.name`. Nothing is rounded or
reworded in the HTML, so the page and the computation cannot drift apart.

The build **fails** rather than producing a quietly wrong report when the
template uses a placeholder no fact defines, a derived fact is referenced
nowhere, or a fact can't be derived from the data at all.

Claims that are only true for some stations are wrapped in a conditional
section, so they simply don't appear when the data doesn't support them:

```html
The band is widest at {{summer_open}}.<!--if night_et--> A net
{{night_share}} of ET falls overnight.<!--endif-->
```

Flags are set in `refet.report.compute.build_payload()`: `validation` drops
the whole comparison panel when no reference file is configured; `night_et`
gates the overnight-ET claim.

## Data notes

**The diurnal heatmap is binned on standard time, not clock time.** AgriMet
timestamps are local with DST, so binning on the clock jumps an hour at each
transition and smears the seasonal sunrise/sunset envelope.

**Missing data**: AgriMet spells gaps at least five ways: `NO RECORD`, an
empty field, `m`, the `998877` sentinel, and values carrying a trailing
quality flag (`-459.69-`, `57.86^`). `refet.io.agrimet.clean_value()` handles
all of them and range-checks the result against per-parameter plausibility
ranges. Hourly rows are *omitted* rather than flagged, so gaps appear as
blanks in the series and grey columns in the heatmap.
