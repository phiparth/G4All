"""
G4All — interactive explorer for the curated G-quadruplex sequence database.

Run locally:      streamlit run app.py
Deploy:           push to GitHub, point Streamlit Community Cloud at this repo.
"""

import io

import pandas as pd
import plotly.express as px
import streamlit as st

st.set_page_config(page_title="G4All Explorer", page_icon="🧬", layout="wide")

DATA_PATH = "G4All.csv"  # note: this file is actually xlsx despite the extension

# Columns worth showing by default; the 17 individual Tm columns stay hidden
# behind a toggle so the table doesn't drown in NaNs.
CORE_COLS = [
    "Type", "Sequence (5'→3')", "G4Hunter score", "Conclusion", "Final Tm (°C)",
    "Length (nt)", "Quadparser state", "GC content (%)", "Total G count",
    "Topology (100 mM KCl)", "Name", "Reference", "Study type", "Origin",
]


@st.cache_data(show_spinner="Loading G4All…")
def load_data(path: str) -> pd.DataFrame:
    """Load the dataset, transparently handling the xlsx-as-csv situation."""
    with open(path, "rb") as fh:
        head = fh.read(4)
    # xlsx (and any zip) begins with the 'PK\x03\x04' magic number.
    if head[:2] == b"PK":
        df = pd.read_excel(path, engine="openpyxl")
    else:
        df = pd.read_csv(path)
    return df


def numeric_slider(df: pd.DataFrame, col: str, label: str | None = None, step=None):
    """Build a two-sided range slider for a numeric column; returns a boolean mask."""
    if col not in df.columns or df[col].dropna().empty:
        return pd.Series(True, index=df.index)
    lo, hi = float(df[col].min()), float(df[col].max())
    if lo == hi:
        return pd.Series(True, index=df.index)
    sel_lo, sel_hi = st.sidebar.slider(
        label or col, lo, hi, (lo, hi), step=step
    )
    keep_na = st.sidebar.checkbox(f"…keep rows with no {label or col}", value=True,
                                  key=f"na_{col}")
    in_range = df[col].between(sel_lo, sel_hi)
    return in_range | (df[col].isna() & keep_na)


df = load_data(DATA_PATH)

st.title("🧬 G4All")
st.caption(
    "Interactive explorer for the curated G-quadruplex sequence database. "
    "Filter on the left; everything below reacts live."
)

# ----------------------------------------------------------------------------- filters
st.sidebar.header("Filters")

mask = pd.Series(True, index=df.index)

for col in ["Type", "Conclusion", "Quadparser state", "Study type", "Origin"]:
    if col in df.columns:
        opts = sorted(df[col].dropna().unique().tolist())
        chosen = st.sidebar.multiselect(col, opts, default=opts)
        # rows whose value is chosen, plus NaN rows (kept so nothing silently vanishes)
        mask &= df[col].isin(chosen) | df[col].isna()

st.sidebar.markdown("---")
mask &= numeric_slider(df, "G4Hunter score", "G4Hunter score", step=0.1)
mask &= numeric_slider(df, "Length (nt)", "Length (nt)", step=1.0)
mask &= numeric_slider(df, "GC content (%)", "GC content (%)", step=1.0)
mask &= numeric_slider(df, "Final Tm (°C)", "Final Tm (°C)", step=0.5)

st.sidebar.markdown("---")
query = st.sidebar.text_input("Sequence search (substring / motif, or exact if ticked)", "").strip().upper()
exact = st.sidebar.checkbox("Exact match", value=False)
if query and "Sequence (5'→3')" in df.columns:
    seqs = df["Sequence (5'→3')"].str.upper()
    mask &= (seqs == query) if exact else seqs.str.contains(query, na=False)

fdf = df[mask]

# ----------------------------------------------------------------------------- metrics
c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("Sequences", f"{len(fdf):,}", delta=f"{len(fdf) - len(df):,}" if len(fdf) != len(df) else None)
if "Conclusion" in fdf:
    g4_frac = (fdf["Conclusion"] == "G4").mean() * 100 if len(fdf) else 0
    c2.metric("Confirmed G4", f"{g4_frac:.0f}%")
if "Type" in fdf:
    dna = (fdf["Type"] == "DNA").sum()
    c3.metric("DNA / RNA", f"{dna:,} / {len(fdf) - dna:,}")
if "G4Hunter score" in fdf and fdf["G4Hunter score"].notna().any():
    c4.metric("Median G4Hunter", f"{fdf['G4Hunter score'].median():.2f}")
if "Final Tm (°C)" in fdf and fdf["Final Tm (°C)"].notna().any():
    c5.metric("Median Tm", f"{fdf['Final Tm (°C)'].median():.1f} °C")

tab_table, tab_dist, tab_scatter = st.tabs(["📋 Table", "📊 Distributions", "🔬 Relationships"])

# ----------------------------------------------------------------------------- table
with tab_table:
    show_all = st.checkbox("Show all columns (incl. Tm1–Tm17)", value=False)
    cols = list(df.columns) if show_all else [c for c in CORE_COLS if c in df.columns]
    st.dataframe(fdf[cols], use_container_width=True, height=560)

    csv_bytes = fdf.to_csv(index=False).encode("utf-8")
    st.download_button(
        "⬇️ Download filtered subset (CSV)",
        data=csv_bytes,
        file_name="G4All_filtered.csv",
        mime="text/csv",
    )

# ----------------------------------------------------------------------------- distributions
with tab_dist:
    numeric_cols = [c for c in ["G4Hunter score", "Length (nt)", "GC content (%)",
                                "Total G count", "Final Tm (°C)"] if c in fdf.columns]
    cola, colb = st.columns(2)
    with cola:
        xcol = st.selectbox("Variable", numeric_cols, index=0)
    with colb:
        color_by = st.selectbox(
            "Colour by",
            [c for c in ["Conclusion", "Type", "Quadparser state", None] if c in fdf.columns or c is None],
            index=0,
        )
    if len(fdf):
        fig = px.histogram(
            fdf, x=xcol, color=color_by, marginal="box", nbins=50,
            barmode="overlay", opacity=0.75,
        )
        fig.update_layout(height=520, legend_title_text=color_by or "")
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("No rows match the current filters.")

# ----------------------------------------------------------------------------- relationships
with tab_scatter:
    numeric_cols = [c for c in ["G4Hunter score", "Length (nt)", "GC content (%)",
                                "Total G count", "Final Tm (°C)"] if c in fdf.columns]
    colx, coly, colc = st.columns(3)
    with colx:
        xax = st.selectbox("X", numeric_cols, index=0, key="sx")
    with coly:
        yax = st.selectbox("Y", numeric_cols, index=min(4, len(numeric_cols) - 1), key="sy")
    with colc:
        cby = st.selectbox(
            "Colour", [c for c in ["Conclusion", "Type", "Quadparser state"] if c in fdf.columns],
            index=0, key="sc",
        )
    plot_df = fdf.dropna(subset=[xax, yax])
    if len(plot_df):
        fig = px.scatter(
            plot_df, x=xax, y=yax, color=cby,
            hover_data=[c for c in ["Name", "Sequence (5'→3')", "Reference"] if c in plot_df.columns],
            opacity=0.6,
        )
        fig.update_layout(height=560)
        st.plotly_chart(fig, use_container_width=True)
        st.caption(f"{len(plot_df):,} points (rows missing X or Y are dropped).")
    else:
        st.info("Not enough non-missing data for this pair.")
