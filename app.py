"""
G4All: interactive explorer for the curated G-quadruplex sequence database.

Run locally:      streamlit run app.py
Deploy:           push to GitHub, point Streamlit Community Cloud at this repo.

This build is schema-robust: it discovers columns at runtime and repairs
encoding/mojibake, so new rows OR new/renamed columns in G4All.csv keep working
without editing this file.
"""

import base64
import inspect
import io
import re

import pandas as pd
import plotly.express as px
import streamlit as st

# ------------------------------------------------------------------------ branding
# The mark is a G-quadruplex: three stacked G-tetrads seen edge-on, the phosphate
# backbone running down the flanks with propeller loops, and two K+ ions in the
# central channel. Kept inline so app.py deploys as a single file.
ICON_SVG = """\
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 96 96" width="96" height="96" role="img" aria-label="G-quadruplex icon">
  <title>G-quadruplex</title>
  <defs>
    <linearGradient id="plate" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0%" stop-color="#1f77b4" stop-opacity="0.85"/>
      <stop offset="100%" stop-color="#9467bd" stop-opacity="0.85"/>
    </linearGradient>
  </defs>

  <!-- phosphate backbone: two strands down the flanks, propeller loops top and bottom -->
  <g fill="none" stroke="#5b6b7c" stroke-width="4" stroke-linecap="round" opacity="0.9">
    <path d="M16 24 L16 72"/>
    <path d="M80 24 L80 72"/>
    <path d="M16 24 C 22 10, 44 6, 56 14"/>
    <path d="M80 72 C 74 86, 52 90, 40 82"/>
  </g>

  <!-- three stacked G-tetrads seen edge-on -->
  <g fill="url(#plate)" stroke="#0f4c7a" stroke-width="1.5" stroke-linejoin="round">
    <polygon points="48,14 82,26 48,38 14,26"/>
    <polygon points="48,38 82,50 48,62 14,50"/>
    <polygon points="48,62 82,74 48,86 14,74"/>
  </g>

  <!-- guanine partition inside each tetrad -->
  <g stroke="#ffffff" stroke-width="1.2" opacity="0.55">
    <path d="M14 26 L82 26 M48 14 L48 38"/>
    <path d="M14 50 L82 50 M48 38 L48 62"/>
    <path d="M14 74 L82 74 M48 62 L48 86"/>
  </g>

  <!-- coordinating K+ ions in the central channel -->
  <g fill="#ff7f0e" stroke="#ffffff" stroke-width="1">
    <circle cx="48" cy="38" r="4.5"/>
    <circle cx="48" cy="62" r="4.5"/>
  </g>
</svg>
"""


def favicon():
    """Small raster of the same mark for the browser tab (Pillow ships with Streamlit)."""
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return "\U0001f537"
    size, k = 128, 128 / 96
    img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    bone = (91, 107, 124, 235)
    d.line([(16 * k, 24 * k), (16 * k, 72 * k)], fill=bone, width=int(4 * k))
    d.line([(80 * k, 24 * k), (80 * k, 72 * k)], fill=bone, width=int(4 * k))
    for top, colour in ((14, (31, 119, 180, 230)),
                        (38, (90, 110, 190, 230)),
                        (62, (148, 103, 189, 230))):
        d.polygon([(48 * k, top * k), (82 * k, (top + 12) * k),
                   (48 * k, (top + 24) * k), (14 * k, (top + 12) * k)],
                  fill=colour, outline=(15, 76, 122, 255))
    for cy in (38, 62):
        r = 4.5 * k
        d.ellipse([48 * k - r, cy * k - r, 48 * k + r, cy * k + r], fill=(255, 127, 14, 255))
    return img


st.set_page_config(page_title="G4All", page_icon=favicon(), layout="wide")

# Streamlit renamed use_container_width to width="stretch"; support both.
FIT = ({"width": "stretch"}
       if "width" in inspect.signature(st.dataframe).parameters
       else {"use_container_width": True})
FIT_PLOT = ({"width": "stretch"}
            if "width" in inspect.signature(st.plotly_chart).parameters
            else {"use_container_width": True})

DATA_PATH = "G4All.csv"  # may actually be xlsx despite the extension; handled below

# Displayed decimals. Experimental Tm is not meaningful past 0.01 °C, and neither
# is a G4Hunter score, so the raw values are rounded for display and export.
TM_DECIMALS = 2
SCORE_DECIMALS = 2
GC_DECIMALS = 1

# Preferred display order. Missing ones are skipped; anything not listed is
# appended, so new columns still show when "Show all columns" is on.
PREFERRED_ORDER = [
    "Type", "Sequence", "Length (nt)", "G4Hunter score", "G4Hmax", "Conclusion",
    "Final Tm", "Number of independent Tm determinations", "Number of Tm determinations",
    "Quadparser state", "GC content (%)", "Total G count",
    "Topology (100 mM KCl)", "Name", "Reference", "Study type", "Origin",
]

# Conclusion labels that count as "forms a G4" at each level of strictness.
STABLE_G4_LABELS = {"g4", "stable g4"}
UNSTABLE_G4_LABELS = {"unstable g4", "unstable"}

# ------------------------------------------------------------------ colour scheme
# Blue = G4, orange = No G4, purple = Unstable, neutral grey = Not sure.
COLOR_MAP = {
    # Conclusion
    "G4": "#1f77b4", "Stable G4": "#1f77b4",
    "No G4": "#ff7f0e",
    "Unstable G4": "#9467bd", "Unstable": "#9467bd",
    "Not sure": "#AEB7C2", "Unsure": "#AEB7C2",
    # Quadparser state
    "positive": "#1f77b4", "negative": "#ff7f0e",
    # Type
    "DNA": "#1f77b4", "RNA": "#9467bd",
}
# Fallback for any category not in the map (new labels, new columns).
PALETTE = ["#1f77b4", "#ff7f0e", "#9467bd", "#AEB7C2",
           "#2ca02c", "#d62728", "#8c564b", "#17becf"]

CAT_ORDER = {
    "Conclusion": ["G4", "Unstable G4", "No G4", "Not sure"],
    "Quadparser state": ["positive", "negative"],
    "Type": ["DNA", "RNA"],
}


def cat_order_for(col, frame):
    """category_orders for plotly: known order first, then any unlisted values."""
    if not col or col not in frame.columns:
        return {}
    known = CAT_ORDER.get(col, [])
    present = frame[col].dropna().astype(str).unique().tolist()
    ordered = [v for v in known if v in present] + [v for v in present if v not in known]
    return {col: ordered}


# ---------------------------------------------------------------------------- loading

def _norm(s) -> str:
    """Lowercased alphanumeric-only key for tolerant column matching."""
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _fix_mojibake(s):
    """Undo the classic UTF-8-decoded-as-Latin-1 double-encoding (Â°C -> °C, etc.).

    Only strings that are valid double-encoded UTF-8 are changed; everything
    else (plain ASCII, already-correct unicode, genuine Latin-1) is returned
    untouched.
    """
    if not isinstance(s, str):
        return s
    try:
        return s.encode("latin-1").decode("utf-8")
    except (UnicodeDecodeError, UnicodeEncodeError):
        return s


def _read_and_clean(path_or_buf) -> pd.DataFrame:
    """Read CSV/xlsx from a path or an uploaded buffer, then repair + coerce."""
    # Sniff the magic number; xlsx (a zip) starts with 'PK'.
    if hasattr(path_or_buf, "read"):
        head = path_or_buf.read(4)
        path_or_buf.seek(0)
    else:
        with open(path_or_buf, "rb") as fh:
            head = fh.read(4)

    if head[:2] == b"PK":
        df = pd.read_excel(path_or_buf, engine="openpyxl")
    else:
        df = None
        for enc in ("utf-8", "utf-8-sig", "cp1252", "latin-1"):
            try:
                if hasattr(path_or_buf, "seek"):
                    path_or_buf.seek(0)
                df = pd.read_csv(path_or_buf, encoding=enc)
                break
            except (UnicodeDecodeError, UnicodeError):
                continue
        if df is None:  # last resort: never crash on a bad byte
            if hasattr(path_or_buf, "seek"):
                path_or_buf.seek(0)
            df = pd.read_csv(path_or_buf, encoding="latin-1", encoding_errors="replace")

    # Repair mojibake in headers and in every string cell.
    df.columns = [_fix_mojibake(c) for c in df.columns]
    obj_cols = df.select_dtypes(include=["object", "string"]).columns
    for c in obj_cols:
        df[c] = df[c].map(_fix_mojibake)

    # Promote numeric-looking text columns (e.g. r (K+) stored as strings).
    for c in obj_cols:
        nonnull = df[c].notna().sum()
        if not nonnull:
            continue
        coerced = pd.to_numeric(df[c], errors="coerce")
        if coerced.notna().sum() >= 0.9 * nonnull:
            df[c] = coerced

    return df


@st.cache_data(show_spinner="Loading G4All…")
def load_data(path: str) -> pd.DataFrame:
    return _read_and_clean(path)


# ---------------------------------------------------------------------- schema helpers

def resolve(df: pd.DataFrame, *keys, prefix=None, contains=None):
    """Return the real column name matching any normalized key / prefix / substring."""
    norms = {_norm(c): c for c in df.columns}
    for k in keys:
        if _norm(k) in norms:
            return norms[_norm(k)]
    if prefix:
        p = _norm(prefix)
        for n, real in norms.items():
            if n.startswith(p):
                return real
    if contains:
        sub = _norm(contains)
        for n, real in norms.items():
            if sub in n:
                return real
    return None


def numeric_columns(df: pd.DataFrame, min_coverage=0.0):
    """Numeric columns that vary, ordered by coverage (drops all-NaN / constant)."""
    out = []
    for c in df.select_dtypes(include="number").columns:
        s = df[c].dropna()
        if s.nunique() > 1 and len(s) >= min_coverage * len(df):
            out.append(c)
    out.sort(key=lambda c: df[c].notna().mean(), reverse=True)
    return out


def categorical_columns(df: pd.DataFrame, exclude=(), max_card=30):
    """Low-cardinality text columns suitable for multiselect filters."""
    out = []
    for c in df.select_dtypes(include=["object", "string"]).columns:
        if c in exclude:
            continue
        n = df[c].nunique(dropna=True)
        if 1 < n <= max_card:
            out.append(c)
    return out


def numeric_slider(df: pd.DataFrame, col, label=None, step=None):
    """Sidebar range slider for a numeric column; returns a boolean mask."""
    if not col or col not in df.columns or df[col].dropna().empty:
        return pd.Series(True, index=df.index)
    lo, hi = float(df[col].min()), float(df[col].max())
    if lo == hi:
        return pd.Series(True, index=df.index)
    name = label or col
    sel_lo, sel_hi = st.sidebar.slider(name, lo, hi, (lo, hi), step=step)
    keep_na = st.sidebar.checkbox(f"…keep rows with no {name}", value=True, key=f"na_{col}")
    in_range = df[col].between(sel_lo, sel_hi)
    return in_range | (df[col].isna() & keep_na)


# ---------------------------------------------------------------------------- data load

try:
    df = load_data(DATA_PATH)
except FileNotFoundError:
    st.title("G4All")
    up = st.file_uploader("G4All.csv not found next to app.py, upload it", type=["csv", "xlsx"])
    if not up:
        st.stop()
    df = _read_and_clean(up)

# Resolve the columns the UI cares about, once, against whatever schema loaded.
SEQ_COL = resolve(df, prefix="Sequence")
TM_COLS = [c for c in df.columns if re.match(r"(?i)^tm\s*\d+\b", str(c))]
COL = {
    "g4hunter": resolve(df, "G4Hunter score", contains="g4hunter"),
    "g4hmax": resolve(df, "G4Hmax"),
    "length": resolve(df, "Length (nt)", contains="length"),
    "gc": resolve(df, "GC content (%)", contains="gccontent"),
    "gcount": resolve(df, "Total G count", contains="totalgcount"),
    "final_tm": resolve(df, "Final Tm (°C)", "Final Tm", contains="finaltm"),
    "type": resolve(df, "Type"),
    "conclusion": resolve(df, "Conclusion"),
    "quadparser": resolve(df, "Quadparser state", contains="quadparser"),
}

def header():
    """Title row with the quadruplex mark."""
    b64 = base64.b64encode(ICON_SVG.encode("utf-8")).decode("ascii")
    st.markdown(
        '<div style="display:flex;align-items:center;gap:0.6rem;margin-bottom:0.2rem">'
        f'<img src="data:image/svg+xml;base64,{b64}" width="58" height="58" alt="G4">'
        '<h1 style="margin:0;font-size:2.6rem;letter-spacing:-0.01em">G4All</h1></div>',
        unsafe_allow_html=True,
    )


header()
st.caption(
    "Interactive explorer for the curated G-quadruplex sequence database. "
    "Search and filter on the left; everything below reacts live."
)

mask = pd.Series(True, index=df.index)

# ------------------------------------------------------------------------ search (top)
st.sidebar.header("Search")
if SEQ_COL:
    query = st.sidebar.text_input(
        "Sequence (substring / motif; regex if ticked)", "",
        placeholder="e.g. GGGTTAGGG",
    ).strip().upper()
    use_regex = st.sidebar.checkbox("Treat as regex (e.g. G{3,})", value=False)
    exact = st.sidebar.checkbox("Exact match", value=False)
    if query:
        seqs = df[SEQ_COL].astype(str).str.upper()
        if exact:
            mask &= seqs == query
        else:
            try:
                mask &= seqs.str.contains(query, na=False, regex=use_regex)
            except re.error:
                st.sidebar.warning("Invalid regex, ignored.")
else:
    st.sidebar.caption("No sequence column found in this file.")

st.sidebar.markdown("---")

# -------------------------------------------------------------------------------- filters
st.sidebar.header("Filters")

# Categorical filters: a curated priority order, then any other low-card column,
# so brand-new categorical columns get filters automatically.
priority_cat = [c for c in (COL["type"], COL["conclusion"], COL["quadparser"],
                            resolve(df, "Study type"), resolve(df, "Origin")) if c]
auto_cat = [c for c in categorical_columns(df, exclude={SEQ_COL} | set(priority_cat)) ]
for col in priority_cat + auto_cat:
    opts = sorted(df[col].dropna().astype(str).unique().tolist())
    chosen = st.sidebar.multiselect(col, opts, default=opts)
    mask &= df[col].astype(str).isin(chosen) | df[col].isna()

st.sidebar.markdown("---")
mask &= numeric_slider(df, COL["g4hunter"], "G4Hunter score", step=0.1)
mask &= numeric_slider(df, COL["length"], "Length (nt)", step=1.0)
mask &= numeric_slider(df, COL["gc"], "GC content (%)", step=1.0)
mask &= numeric_slider(df, COL["final_tm"], "Final Tm (°C)", step=0.5)

# Any other numeric column stays reachable without cluttering the sidebar.
primary_num = {COL["g4hunter"], COL["length"], COL["gc"], COL["final_tm"]}
extra_num = [c for c in numeric_columns(df) if c not in primary_num]
with st.sidebar.expander("More numeric filters"):
    picked = st.multiselect("Add filters for", extra_num, key="extra_num")
for col in picked:
    mask &= numeric_slider(df, col)

fdf = df[mask]

# -------------------------------------------------------------------------------- metrics
c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("Sequences", f"{len(fdf):,}",
          delta=f"{len(fdf) - len(df):,}" if len(fdf) != len(df) else None)

if COL["conclusion"]:
    concl = fdf[COL["conclusion"]].astype(str).str.strip().str.lower()
    n_stable = concl.isin(STABLE_G4_LABELS).sum()
    n_unstable = concl.isin(UNSTABLE_G4_LABELS).sum()
    stable_pct = 100 * n_stable / len(fdf) if len(fdf) else 0
    any_pct = 100 * (n_stable + n_unstable) / len(fdf) if len(fdf) else 0
    c2.metric(
        "Stable G4",
        f"{stable_pct:.0f}%",
        delta=f"{any_pct:.0f}% incl. unstable",
        delta_color="off",
        help=(
            "Share of the sequences currently in view whose experimental verdict is "
            '"G4", i.e. a G4 that folds **and** is stable under the assay conditions. '
            'Sequences labelled "Unstable G4" (G4 observed but low thermal stability) '
            "are excluded from this figure and counted in the second line; "
            '"No G4" and "Not sure" are excluded from both.'
        ),
    )

if COL["type"]:
    dna = (fdf[COL["type"]] == "DNA").sum()
    c3.metric("DNA / RNA", f"{dna:,} / {len(fdf) - dna:,}")
if COL["g4hunter"] and fdf[COL["g4hunter"]].notna().any():
    c4.metric("Median G4Hunter", f"{fdf[COL['g4hunter']].median():.2f}")
if COL["final_tm"] and fdf[COL["final_tm"]].notna().any():
    c5.metric("Median Tm", f"{fdf[COL['final_tm']].median():.1f} °C")

if COL["conclusion"]:
    st.caption(
        f'**Stable G4** counts rows with Conclusion = "G4" ({n_stable:,} of {len(fdf):,} '
        f"in view, {stable_pct:.1f}%). Adding the {n_unstable:,} rows labelled "
        f'"Unstable G4" (G4 formation observed but not thermally stable) gives '
        f"{any_pct:.1f}%. The remainder are \"No G4\" or \"Not sure\". "
        "Use the Conclusion filter on the left to restrict any view to a single verdict."
    )

tab_table, tab_dist, tab_scatter = st.tabs(["📋 Table", "📊 Distributions", "🔬 Relationships"])

# Ordered column list: preferred first (resolved), then the rest.
ordered = []
for key in PREFERRED_ORDER:
    real = resolve(df, key, prefix=key, contains=key)
    if real and real not in ordered:
        ordered.append(real)
ordered += [c for c in df.columns if c not in ordered]
core_cols = [c for c in ordered if c not in TM_COLS]

# --------------------------------------------------------------- display precision / QC
ROUNDING = {}
for _c in TM_COLS + [COL["final_tm"]]:
    if _c:
        ROUNDING[_c] = TM_DECIMALS
for _c in (COL["g4hunter"], COL["g4hmax"], resolve(df, "r (K+)"), resolve(df, "r (Na+)")):
    if _c:
        ROUNDING[_c] = SCORE_DECIMALS
if COL["gc"]:
    ROUNDING[COL["gc"]] = GC_DECIMALS


def column_config(frame):
    """Number formats so nothing is shown to more digits than it is known to."""
    cfg = {}
    for c, d in ROUNDING.items():
        if c in frame.columns and pd.api.types.is_numeric_dtype(frame[c]):
            cfg[c] = st.column_config.NumberColumn(c, format=f"%.{d}f")
    return cfg


def round_for_export(frame):
    return frame.round({c: d for c, d in ROUNDING.items() if c in frame.columns})


# -------------------------------------------------------------------------------- table
with tab_table:
    show_all = st.checkbox(
        f"Show all columns (incl. {len(TM_COLS)} Tm columns)" if TM_COLS else "Show all columns",
        value=False,
    )
    cols = ordered if show_all else core_cols
    st.dataframe(fdf[cols], height=560, column_config=column_config(fdf), **FIT)

    # utf-8-sig so the ° / arrow render correctly if opened in Excel.
    raw = st.checkbox(
        "Export raw unrounded values", value=False,
        help="By default the download carries the same rounding as the table above.",
    )
    export = fdf if raw else round_for_export(fdf)
    csv_bytes = export.to_csv(index=False).encode("utf-8-sig")
    st.download_button(
        "⬇️ Download filtered subset (CSV)",
        data=csv_bytes,
        file_name="G4All_filtered.csv",
        mime="text/csv",
    )

# -------------------------------------------------------------------------------- distributions
num_cols = numeric_columns(fdf) or numeric_columns(df)
cat_for_color = [c for c in (COL["conclusion"], COL["type"], COL["quadparser"]) if c]

with tab_dist:
    if not num_cols:
        st.info("No numeric columns to plot.")
    else:
        cola, colb = st.columns(2)
        with cola:
            xcol = st.selectbox("Variable", num_cols, index=0)
        with colb:
            color_by = st.selectbox("Colour by", cat_for_color + [None], index=0)
        if len(fdf):
            fig = px.histogram(fdf, x=xcol, color=color_by, marginal="box",
                               nbins=50, barmode="overlay", opacity=0.75,
                               color_discrete_map=COLOR_MAP,
                               color_discrete_sequence=PALETTE,
                               category_orders=cat_order_for(color_by, fdf))
            fig.update_layout(height=520, legend_title_text=color_by or "")
            st.plotly_chart(fig, **FIT_PLOT)
        else:
            st.info("No rows match the current filters.")

# -------------------------------------------------------------------------------- relationships
with tab_scatter:
    if len(num_cols) < 2:
        st.info("Need at least two numeric columns for a scatter plot.")
    else:
        colx, coly, colc = st.columns(3)
        with colx:
            xax = st.selectbox("X", num_cols, index=0, key="sx")
        with coly:
            yax = st.selectbox("Y", num_cols, index=min(len(num_cols) - 1, 4), key="sy")
        with colc:
            cby = st.selectbox("Colour", cat_for_color or [None], index=0, key="sc")
        plot_df = fdf.dropna(subset=[xax, yax])
        if len(plot_df):
            hover = [c for c in (resolve(df, "Name"), SEQ_COL, resolve(df, "Reference"))
                     if c]
            fig = px.scatter(plot_df, x=xax, y=yax, color=cby,
                             hover_data=hover, opacity=0.6,
                             color_discrete_map=COLOR_MAP,
                             color_discrete_sequence=PALETTE,
                             category_orders=cat_order_for(cby, plot_df))
            fig.update_layout(height=560)
            st.plotly_chart(fig, **FIT_PLOT)
            st.caption(f"{len(plot_df):,} points (rows missing X or Y are dropped).")
        else:
            st.info("Not enough non-missing data for this pair.")
