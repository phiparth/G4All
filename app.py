"""
G4All — interactive explorer for the curated G-quadruplex sequence database.

Run locally:      streamlit run app.py
Deploy:           push to GitHub, point Streamlit Community Cloud at this repo.

This build is schema-robust: it discovers columns at runtime and repairs
encoding/mojibake, so new rows OR new/renamed columns in G4All.csv keep working
without editing this file.
"""

import io
import re

import pandas as pd
import plotly.express as px
import streamlit as st

st.set_page_config(page_title="G4All", page_icon="🧬", layout="wide")

DATA_PATH = "G4All.csv"  # may actually be xlsx despite the extension; handled below

# Preferred display order. Missing ones are skipped; anything not listed is
# appended, so new columns still show when "Show all columns" is on.
PREFERRED_ORDER = [
    "Type", "Sequence", "G4Hunter score", "G4Hmax", "Conclusion", "Final Tm",
    "Length (nt)", "Quadparser state", "GC content (%)", "Total G count",
    "Topology (100 mM KCl)", "Name", "Reference", "Study type", "Origin",
]

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
    st.title("🧬 G4All")
    up = st.file_uploader("G4All.csv not found next to app.py — upload it", type=["csv", "xlsx"])
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

st.title("🧬 G4All")
st.caption(
    "Interactive explorer for the curated G-quadruplex sequence database. "
    "Filter on the left; everything below reacts live."
)

# -------------------------------------------------------------------------------- filters
st.sidebar.header("Filters")
mask = pd.Series(True, index=df.index)

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

st.sidebar.markdown("---")
if SEQ_COL:
    query = st.sidebar.text_input(
        "Sequence search (substring / motif; regex if ticked)", ""
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
                st.sidebar.warning("Invalid regex — ignored.")

fdf = df[mask]

# -------------------------------------------------------------------------------- metrics
c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("Sequences", f"{len(fdf):,}",
          delta=f"{len(fdf) - len(df):,}" if len(fdf) != len(df) else None)
if COL["conclusion"]:
    g4_frac = (fdf[COL["conclusion"]] == "G4").mean() * 100 if len(fdf) else 0
    c2.metric("Confirmed G4", f"{g4_frac:.0f}%")
if COL["type"]:
    dna = (fdf[COL["type"]] == "DNA").sum()
    c3.metric("DNA / RNA", f"{dna:,} / {len(fdf) - dna:,}")
if COL["g4hunter"] and fdf[COL["g4hunter"]].notna().any():
    c4.metric("Median G4Hunter", f"{fdf[COL['g4hunter']].median():.2f}")
if COL["final_tm"] and fdf[COL["final_tm"]].notna().any():
    c5.metric("Median Tm", f"{fdf[COL['final_tm']].median():.1f} °C")

tab_table, tab_dist, tab_scatter = st.tabs(["📋 Table", "📊 Distributions", "🔬 Relationships"])

# Ordered column list: preferred first (resolved), then the rest.
ordered = []
for key in PREFERRED_ORDER:
    real = resolve(df, key, prefix=key, contains=key)
    if real and real not in ordered:
        ordered.append(real)
ordered += [c for c in df.columns if c not in ordered]
core_cols = [c for c in ordered if c not in TM_COLS]

# -------------------------------------------------------------------------------- table
with tab_table:
    show_all = st.checkbox(
        f"Show all columns (incl. {len(TM_COLS)} Tm columns)" if TM_COLS else "Show all columns",
        value=False,
    )
    cols = ordered if show_all else core_cols
    st.dataframe(fdf[cols], use_container_width=True, height=560)

    # utf-8-sig so the ° / arrow render correctly if opened in Excel.
    csv_bytes = fdf.to_csv(index=False).encode("utf-8-sig")
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
            st.plotly_chart(fig, use_container_width=True)
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
            st.plotly_chart(fig, use_container_width=True)
            st.caption(f"{len(plot_df):,} points (rows missing X or Y are dropped).")
        else:
            st.info("Not enough non-missing data for this pair.")
