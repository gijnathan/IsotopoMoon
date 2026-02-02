import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

fname = "time_evolution_cpp.txt"

# Skip the first 2 preamble rows, row 3 becomes the header
df = pd.read_csv(fname, skiprows=2)

# ---- helpers ----
def num(col):
    return pd.to_numeric(df[col], errors="coerce")

def first_existing(options):
    for c in options:
        if c in df.columns:
            return c
    return None

# ---- choose time column ----
t_col = first_existing(["time (Myr)", "time", "t (Myr)", "t_Myr"])
if t_col is None:
    raise ValueError(f"Couldn't find a time column. Columns are:\n{df.columns.tolist()}")

t = num(t_col)

# ---- map δ34S columns (edit these to match your CSV header names) ----
# These names are guesses based on your earlier "d34S M" example.
col_mantle = first_existing(["d34S M", "d34S_mantle", "mantle d34S", "d34S mantle"])
col_frost  = first_existing(["d34S F", "d34S_frost", "frost d34S"])
col_ss     = first_existing(["d34S SS", "d34S_sulfate", "d34S sulfates", "d34S silsulf"])
col_space  = first_existing(["d34S S", "d34S_space", "space d34S"])

# Optional: "outgassing" and "instantaneous space" (if you output them)
col_outgas = first_existing(["d34S outgassing", "d34S_outgas", "d34S outg", "d34S O"])
col_space_inst = first_existing(["d34S inst space", "d34S_space_inst", "d34S instantaneous space"])

# Optional: "Rd mantle/space" style curves (if present)
col_rd_mantle = first_existing(["d34S Rd mantle", "d34S_Rd_mantle", "Rd mantle d34S"])
col_rd_space  = first_existing(["d34S Rd space", "d34S_Rd_space", "Rd space d34S"])
col_rd_space_inst = first_existing(["d34S Rd inst space", "d34S_Rd_space_inst", "Rd inst space d34S"])

# ---- observed values (set these from your script / constants) ----
# Replace with your actual observed mean + 1σ (or whatever you use).
d34S_obs   = 0.0
d34S_obs_sd = 0.0

# If you have these in the CSV as single-valued columns, you can auto-read them:
obs_col = first_existing(["d34S_obs", "d34S observed", "obs d34S"])
sd_col  = first_existing(["d34Ssd", "d34S_sd", "d34S obs sd", "d34Ssd_obs"])
if obs_col is not None:
    d34S_obs = float(num(obs_col).dropna().iloc[0])
if sd_col is not None:
    d34S_obs_sd = float(num(sd_col).dropna().iloc[0])

# ---- plot ----
plt.figure(figsize=(8, 5))

# observed band + observed line (like your Plotly polygon + horizontal line)
xmin, xmax = -1000, 10000
if d34S_obs_sd != 0.0:
    plt.fill_between(
        [xmin, xmax],
        [d34S_obs - d34S_obs_sd, d34S_obs - d34S_obs_sd],
        [d34S_obs + d34S_obs_sd, d34S_obs + d34S_obs_sd],
        alpha=0.25,
        linewidth=0,
        label="observed ± sd",
    )
plt.plot([xmin, xmax], [d34S_obs, d34S_obs], linestyle="-", linewidth=2, label="observed")

# Rd curves (black solid/dash/dot)
if col_rd_mantle:
    plt.plot(t, num(col_rd_mantle), linestyle="-", linewidth=2, label="Rd mantle")
if col_rd_space:
    plt.plot(t, num(col_rd_space), linestyle="--", linewidth=2, label="Rd space")
if col_rd_space_inst:
    plt.plot(t, num(col_rd_space_inst), linestyle=":", linewidth=2, label="Rd inst space")

# Reservoir curves (mantle/firebrick, outgassing firebrick dashed, frost orange, sulfates steelblue, space grey)
if col_mantle:
    plt.plot(t, num(col_mantle), linestyle="-", linewidth=2, label="mantle")
if col_outgas:
    plt.plot(t, num(col_outgas), linestyle="--", linewidth=2, label="outgassing")
if col_frost:
    plt.plot(t, num(col_frost), linestyle="-", linewidth=2, label="frost")
if col_ss:
    plt.plot(t, num(col_ss), linestyle="-", linewidth=2, label="sulfates/silicates")
if col_space:
    plt.plot(t, num(col_space), linestyle="-", linewidth=2, label="space")
if col_space_inst:
    plt.plot(t, num(col_space_inst), linestyle="--", linewidth=2, label="inst space")

# match your Plotly axis styling
plt.ylim(-100, 500)
plt.yticks(np.arange(-100, 501, 100))
plt.xlabel("Time (Myr)")
plt.ylabel("δ34S_VCDT (‰)")
plt.title("δ34S vs Time")
plt.legend(loc="best")
plt.tight_layout()
plt.show()
