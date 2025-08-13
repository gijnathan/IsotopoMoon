import time
import isotope_evolution as ie

start_time = time.perf_counter()

def fmt(val):
    """Smart formatter: 6 dp if stable, else up to 17 dp, strip trailing zeros."""
    rounded = round(val, 6)
    if abs(val - rounded) < 1e-12:
        s = f"{val:.6f}"
    else:
        s = f"{val:.17f}"
    return s.rstrip("0").rstrip(".")  # strip trailing zeros and lone decimal point


a34 = 0.9998
X = 33

print("Testing isotope evolution functions...")
print("Testing a_3X():")
result = ie.a_3X(a34, X)
print(f"\ta_{X} = {fmt(result)}")

print("Testing alphas():")
process = "gr"
nofrac = "no"
MAF = "no"
result = ie.alphas(process, X, nofrac, MAF)
print(f"\talpha_{X} = {fmt(result)}")

a34_ec = ie.alphas("ec", 34, nofrac, MAF)
a33_pi = ie.alphas("pi", 33, nofrac, MAF)
a36_fr = ie.alphas("fr", 36, nofrac, MAF)
print("Testing selected alphas()...")
print(f"\ta34_ec = {fmt(a34_ec)}")
print(f"\ta33_pi = {fmt(a33_pi)}")
print(f"\ta36_fr = {fmt(a36_fr)}")

print("=== Test rate_kgs ===")
proc = "pi"
if f"{proc}_min" in ie.rate_kgs and f"{proc}_max" in ie.rate_kgs:
    rate_min = ie.rate_kgs[f"{proc}_min"]
    rate_max = ie.rate_kgs[f"{proc}_max"]
    print(f"\t{proc}: min = {fmt(rate_min)}, max = {fmt(rate_max)}")

print("\n=== Test VCDT ===")
for iso, ratio in ie.VCDT.items():
    print(f"\tdelta {iso} ratio = {fmt(ratio)}")

# Example isotope ratio
R = 0.045
n = 34
delta = ie.R2d(R, n, ie.VCDT)
print("\n=== Test R2d() ===")
print(f"\tInput ratio R = {fmt(R)}")
print(f"\tIsotope n = {n}")
print(f"\tDelta_{n} = {fmt(delta)}")

print("=== Test d2R ===")
d = 10.0  # +10‰
R = ie.d2R(d, n, ie.VCDT)
print(f"\td = {fmt(d)} ‰, n = {n} -> R = {fmt(R)}")

print("\n=== Test ROT ===")
a, c, b, d_ = 0.0079, 1.0, 0.0442, 0.000153
X = ie.ROT(a, c, b, d_)
print(f"\tX (33S/ST) = {fmt(X)}")

print("\n=== Test kgs2mSs ===")
kgs = 205.0
molSps = ie.kgs2mSs(kgs)
print(f"\t{fmt(kgs)} kg/s SO2 -> {fmt(molSps)} mol S/s")

def fmt(val: float) -> str:
    rounded = round(val, 6)
    s = f"{val:.6f}" if abs(val - rounded) < 1e-12 else f"{val:.17f}"
    s = s.rstrip("0").rstrip(".")
    if "." not in s:
        s += ".0"
    return s

print("\n=== Test MIF ===")
MIF_eq = "FW2003_pg3"

# 1) Standards -> ~0,0
R33s, R34s, R36s = ie.VCDT["33"], ie.VCDT["34"], ie.VCDT["36"]
D33s, D36s = ie.MIF(R33s, R34s, R36s)
print(f"\tAt standards: D33 = {fmt(D33s)}, D36 = {fmt(D36s)}")

# 2) Build ratios from deltas
d34, d33, d36 = 10.0, 5.0, 20.0
R34 = ie.d2R(d34, 34, ie.VCDT)
R33 = ie.d2R(d33, 33, ie.VCDT)
R36 = ie.d2R(d36, 36, ie.VCDT)
D33x, D36x = ie.MIF(R33, R34, R36)
print(f"\tFor d34={fmt(d34)}‰, d33={fmt(d33)}‰, d36={fmt(d36)}‰  -> "
      f"D33 = {fmt(D33x)}, D36 = {fmt(D36x)}")

print("\n=== Test calc_rate ===")
pd_min, pd_max = ie.rate_kgs["pd_min"], ie.rate_kgs["pd_max"]
print(f"\tf='N' -> {fmt(ie.calc_rate('N', pd_min, pd_max))} kg/s")
print(f"\tf=0.5 -> {fmt(ie.calc_rate(0.5, pd_min, pd_max))} kg/s")
print(f"\tf=2.0 -> {fmt(ie.calc_rate(2.0, pd_min, pd_max))} kg/s")

# === Test initial_mass ===
min_mass = 2.95e21
max_mass = 2.24e22
print("=== Test initial_mass ===")
print(f"fraction=0.0 -> {ie.initial_mass(min_mass, max_mass, 0.0):.6f}")
print(f"fraction=1.0 -> {ie.initial_mass(min_mass, max_mass, 1.0):.6f}")
print(f"fraction=0.5 -> {ie.initial_mass(min_mass, max_mass, 0.5):.6f}")

# === Test reservoir_totals ===
M, DM, F, SS, S = 1.0, 2.0, 3.0, 4.0, 5.0
print("\n=== Test reservoir_totals ===")
print(f"Sum = {ie.reservoir_totals(M, DM, F, SS, S):.6f}")

# === Test mass_balance ===
init, now = 100.0, 80.0
print("\n=== Test mass_balance ===")
print(f"Mass balance fraction = {ie.mass_balance(init, now):.6f}")

# === Test S_conc ===
ST = 1.0e6   # moles S
mass = 5.0e9 # kg
print("\n=== Test S_conc ===")
print(f"Concentration = {ie.S_conc(ST, mass):.6f}")

# === Test molSs_mm ===
mag_kgs = 200.0
mag_S = 0.01
f = 0.5
print("\n=== Test molSs_mm ===")
print(f"mol S/s = {ie.molSs_mm(mag_kgs, mag_S, f):.6f}")

print("\n=== Test reservoir_R_S ===")
S32, S33, S34, S36 = 100.0, 1.0, 4.4, 0.15
R33, R34, R36, T32, T33, T34, T36, ST = ie.reservoir_R_S(S32, S33, S34, S36)
print(f"Inputs: S32={S32}, S33={S33}, S34={S34}, S36={S36}")
print(f"R33 = {R33:.6f}, R34 = {R34:.6f}, R36 = {R36:.6f}")
print(f"T32 = {T32:.6f}, T33 = {T33:.6f}, T34 = {T34:.6f}, T36 = {T36:.6f}")
print(f"ST  = {ST:.6f}")

# edge case: S32 == 0
R33z, R34z, R36z, T32z, T33z, T34z, T36z, STz = ie.reservoir_R_S(0.0, 1.0, 1.0, 1.0)
print("\nS32 == 0 edge case:")
print(f"R33={R33z}, R34={R34z}, R36={R36z}, T32={T32z}, T33={T33z}, T34={T34z}, T36={T36z}, ST={STz}")

print("\n=== Test reservoir_isotope ===")
S32, S33, S34, S36 = 100.0, 1.0, 4.4, 0.15
R33, R34, R36, d33, d34, d36, D33, D36, TS, log_ST = ie.reservoir_isotope(S32, S33, S34, S36)
print(f"Inputs: S32={S32}, S33={S33}, S34={S34}, S36={S36}")
print(f"R33={R33:.6f}, R34={R34:.6f}, R36={R36:.6f}")
print(f"d33={d33:.6f}‰, d34={d34:.6f}‰, d36={d36:.6f}‰")
print(f"D33={D33:.6f}‰, D36={D36:.6f}‰")
print(f"ST={TS:.6f}, log_ST={log_ST:.6f}")

# edge case S32 == 0
print("\nS32 == 0 edge case:")
vals = ie.reservoir_isotope(0.0, 1.0, 1.0, 1.0)
print(vals)  # Python returns strings and zeros here by original design

print("\n=== Test initial_reservoir ===")
ST_in, d33_in, d34_in, d36_in = 1.0e9, 0.0, 10.0, 20.0
vals = ie.initial_reservoir(ST_in, d33_in, d34_in, d36_in)
(log_ST_ir, R33_ir, R34_ir, R36_ir,
 T32_ir, T33_ir, T34_ir, T36_ir,
 D33_ir, D36_ir,
 S32_ir, S33_ir, S34_ir, S36_ir) = vals

try:
    fmt  # if you've defined the smart formatter earlier
    f = fmt
except NameError:
    f = lambda x: f"{x:.6f}"

print(f"Inputs: ST={f(ST_in)}, d33={f(d33_in)}‰, d34={f(d34_in)}‰, d36={f(d36_in)}‰")
print(f"R33={f(R33_ir)}, R34={f(R34_ir)}, R36={f(R36_ir)}")
print(f"T32={f(T32_ir)}, T33={f(T33_ir)}, T34={f(T34_ir)}, T36={f(T36_ir)}")
print(f"D33={f(D33_ir)}‰, D36={f(D36_ir)}‰")
print(f"S32={f(S32_ir)}, S33={f(S33_ir)}, S34={f(S34_ir)}, S36={f(S36_ir)}")
print(f"log_ST={f(log_ST_ir)}")

print("\nST <= 0 edge case:")
print(ie.initial_reservoir(0.0, 0.0, 0.0, 0.0))

print("\n=== Test process_R_S ===")

# Example inputs
F = 1.5
res_33R = 0.0079
res_34R = 0.0442
res_36R = 0.000153
a33 = 0.999
a34 = 0.998
a36 = 1.001

result = ie.process_R_S(F, res_33R, res_34R, res_36R, a33, a34, a36)
print("Inputs:")
print(f"  F={F}, res_33R={res_33R}, res_34R={res_34R}, res_36R={res_36R}")
print(f"  a33={a33}, a34={a34}, a36={a36}")
print("Outputs:")
print(f"  R33={result[0]}, R34={result[1]}, R36={result[2]}")
print(f"  S32={result[3]}, S33={result[4]}, S34={result[5]}, S36={result[6]}")

# Edge case: F = 0.0
F_zero = 0.0
res_zero = ie.process_R_S(F_zero, res_33R, res_34R, res_36R, a33, a34, a36)
print("\nEdge case: F=0.0")
print(f"  Output: {res_zero}")


print("\n=== Test mantle, deep_mantle, frost, silsulf, space ===")

b_M  = 1.0e22  # kg
b_mo = 1.0e21  # kg
b_pl = 5.0e20  # kg
b_rs = 2.0e20  # kg
b_dm = 3.0e20  # kg
m_val = ie.mantle(b_M, b_mo, b_pl, b_rs, b_dm)
print("mantle:", m_val)

b_DM = 4.0e20  # kg
dmix_val = ie.deep_mantle(b_DM, b_dm)
print(" deep mantle:", dmix_val)

b_F  = 1.0e20  # kg
b_fr = 2.0e19  # kg
b_bu = 1.0e19  # kg
b_hg = 5.0e18  # kg
frost_val = ie.frost(b_F, b_fr, b_bu, b_hg)
print(" frost:", frost_val)

b_SS = 6.0e19  # kg
b_sq = 3.0e19  # kg
silsulf_val = ie.silsulf(b_SS, b_rs, b_sq, b_pl)
print(" silsulf:", silsulf_val)

b_S  = 7.0e19  # kg
b_pu = 2.0e19  # kg
space_val = ie.space(b_S, b_pu)
print(" space:", space_val)



print("\n=== Test a_bu and a_pr ===")

sn_F = 2.0e19
pd_F = 3.0e19
a_vp = 0.995
a_pd = 1.002
a_bu_val = ie.a_bu(sn_F, pd_F, a_vp, a_pd)
print("a_bu:", a_bu_val)

alphas = {
    "pu":0.98, "gr":1.01, "th":0.97, "ed":1.003,
    "ac":0.999, "ei":1.002, "rc":1.001, "ec":0.998,
    "pi":1.004
}
rates = {"pr":10.0, "ed":2.0, "ac":1.5, "ei":3.0, "rc":1.0, "ec":0.5, "pi":2.0}
a_pr_val = ie.a_pr(alphas, rates)
print("a_pr:", a_pr_val)

print("\n=== Test newton_raphson (module version) ===")

# Solve x^2 - A = 0 with A=2  (root = sqrt(2))
constants = {"A": 2.0}

def eqs(x, c):
    return x*x - c["A"]

def deriv(x, c):
    return 2.0*x

x0   = 1.0    # initial guess
e1   = 1e-12  # tolerance
step = 1.0    # NR step scaling

t0 = time.perf_counter()
root = ie.newton_raphson(x0, constants, e1, step, eqs, deriv)
t1 = time.perf_counter()

print("root =", fmt(root))
print("elapsed_seconds =", fmt(t1 - t0))



import time
import isotope_evolution as ie

print("\n=== Test fractionate (Python) ===")

# constants tuple order must match your function:
# (STT, STE, STC, S32T, S33T, S34T, S36T, a33, a34, a36)
constants = (
    0.0,    # STT (unused in f/df here)
    0.0,    # STE (unused in f/df here)
    10.0,   # STC
    100.0,  # S32T
    1.0,    # S33T
    4.4,    # S34T
    0.15,   # S36T
    0.999,  # a33
    0.998,  # a34
    1.001   # a36
)

nr_step = 1.0
nr_tol  = 1e-12
guessx  = 50.0

t0 = time.perf_counter()
result = ie.fractionate(constants, nr_step, nr_tol, guessx)
t1 = time.perf_counter()

(S32C, S33C, S34C, S36C, S32E, S33E, S34E, S36E) = result

def fmt(val: float) -> str:
    rounded = round(val, 6)
    s = f"{val:.6f}" if abs(val - rounded) < 1e-12 else f"{val:.17f}"
    s = s.rstrip("0").rstrip(".")
    if "." not in s:
        s += ".0"
    return s

print("Crust:")
print(f"  S32C={fmt(S32C)}, S33C={fmt(S33C)}, S34C={fmt(S34C)}, S36C={fmt(S36C)}")
print("Ejecta:")
print(f"  S32E={fmt(S32E)}, S33E={fmt(S33E)}, S34E={fmt(S34E)}, S36E={fmt(S36E)}")
print("elapsed_seconds =", fmt(t1 - t0))





#########################
end_time = time.perf_counter()
print(f"\nExecution time: {fmt(end_time - start_time)} seconds")

