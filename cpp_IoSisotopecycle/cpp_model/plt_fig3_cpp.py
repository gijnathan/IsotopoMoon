import pandas as pd
import matplotlib.pyplot as plt

# Your C++ output file
fname = "time_evolution_cpp.txt"

# Skip the first 2 preamble rows, row 3 becomes the header
df = pd.read_csv(fname, skiprows=2)

# Ensure columns are numeric (blanks become NaN automatically)
df["time (Myr)"] = pd.to_numeric(df["time (Myr)"], errors="coerce")
df["d34S M"] = pd.to_numeric(df["d34S M"], errors="coerce")

# Plot
plt.figure()
plt.plot(df["time (Myr)"], df["d34S M"])
plt.xlabel("Time (Myr)")
plt.ylabel("δ34S Mantle (‰)")
plt.title("Mantle δ34S vs Time")
plt.show()
