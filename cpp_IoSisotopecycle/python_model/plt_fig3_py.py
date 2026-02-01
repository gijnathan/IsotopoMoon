import pandas as pd
import matplotlib.pyplot as plt

fname = "time_evolution.csv"

# Skip the first 2 preamble rows, row 3 becomes the header
df = pd.read_csv(fname, skiprows=2)

# Convert to numeric (blanks -> NaN)
df["time (Myr)"] = pd.to_numeric(df["time (Myr)"], errors="coerce")
df["d34S M"] = pd.to_numeric(df["d34S M"], errors="coerce")

# Plot mantle d34S vs time
plt.figure()
plt.plot(df["time (Myr)"], df["d34S M"])
plt.xlabel("Time (Myr)")
plt.ylabel("δ34S Mantle (‰)")
plt.title("Mantle δ34S vs Time")
plt.show()
