import numpy as np
import matplotlib.pyplot as plt
from math import pi, log

# -----------------------------
#  Monte-Carlo – fuga difusiva
# -----------------------------
N = 10_000                      # número de ensayos
rng = np.random.default_rng(42) # semilla reproducible

# --- Geometría de la tubería ---
ri = 0.05          # radio interior (m)
ro = 0.06          # radio exterior (m)
L  = 100.0         # longitud analizada (m)
ln_ro_ri = log(ro / ri)

# --- Constantes termodinámicas ---
R_gas = 8.314      # J·mol⁻¹·K⁻¹
T     = 293.0      # K
P     = 101_325.0  # Pa
LFL   = 0.04       # 4 % vol → fracción 0.04
L_prec = 0.4 * LFL # 40 % del LFL (1.6 % vol)

# --- Incertidumbres (distribuciones) ---
# 1) Difusividad efectiva (HDPE) – log-normal
D_med  = 1.5e-11         # mediana (m²/s)
g_sd   = 1.8             # desviación geométrica
ln_sig = np.log(g_sd)
D_eff  = rng.lognormal(np.log(D_med), ln_sig, size=N)

# 2) Concentración interna de H2 – normal (COV 5 %)
Ci_mu, Ci_cov = 40.0, 0.05
Ci = rng.normal(Ci_mu, Ci_cov * Ci_mu, size=N)
Ci = np.clip(Ci, 0, None)              # sin negativos

# 3) Volumen del recinto – normal (COV 10 %)
V_mu, V_cov  = 0.20, 0.10              # m³
V = rng.normal(V_mu, V_cov * V_mu, size=N)
V = np.clip(V, 1e-3, None)

# 4) Intervalo de inspección τ – normal (μ = 24 h, COV 50 %)
tau_mu, tau_cov = 24.0, 0.50           # horas
tau_h = rng.normal(tau_mu, tau_cov * tau_mu, size=N)
tau_h = np.clip(tau_h, 0.1, None)
tau_s = tau_h * 3600.0                 # s

# --- Cálculos de fuga ---
J      = D_eff * Ci / ln_ro_ri         # mol·m⁻²·s⁻¹  (Co = 0)
A      = 2 * pi * L * ri               # área lateral (m²)
dot_n  = J * A                         # mol·s⁻¹
n_acc  = dot_n * tau_s                 # mol acumulados

phi = (n_acc * R_gas * T) / (P * V)    # fracción volumétrica

# --- Resultados estadísticos ---
P_prec = np.mean((phi >= L_prec) & (phi < LFL))
P_fl   = np.mean(phi >= LFL)

print("\n***  Monte-Carlo: fuga difusiva de H₂  ***")
print(f"Ensayos simulados:           {N:,}")
print(f"Concentración media:         {phi.mean()*100:.3f} % vol")
print(f"Percentil 95 %:              {np.percentile(phi,95)*100:.3f} % vol")
print(f"P(φ ≥ 40 % LFL = 1.6 %) :   {P_prec:.4f}")
print(f"P(φ ≥ 100 % LFL = 4 %)  :   {P_fl  :.4f}")

# --- Histograma ---
plt.figure(figsize=(7,4))
plt.hist(phi*100, bins=60, edgecolor='black')
plt.axvline(L_prec*100, linestyle='--', label='40 % LFL')
plt.axvline(LFL*100,  linestyle='--', label='LFL (4 %)')
plt.xlabel('Concentración H₂ (% vol)')
plt.ylabel('Frecuencia')
plt.title('Distribución Monte-Carlo de la concentración de H₂')
plt.legend()
plt.tight_layout()
plt.show()
