import json
import math

with open("data/corrosion_database.json") as f:
    db = json.load(f)
p = db["calibrated_parameters"]

R_GAS = 8.31446261815324   # J/(mol*K), matches CorrosionChemistry.h
T_REF = 650.0 + 273.15     # = 923.15 K, matches CorrosionChemistry.h

# Kinetic branch
def r_kin(T, material, salt, redox, position, flow, dT):
    theta0 = p["log_rate0_um_y"]
    Ea = p["Ea_corr_kJ_mol"] * 1000.0  # kJ/mol -> J/mol
    a_M = p.get(f"mat_{material}", 0.0)
    a_S = p.get(f"salt_{salt}", 0.0)
    a_R = p.get(f"redox_{redox}", 0.0)
    a_P = p["hot_leg_bonus"] if position == "hot_leg" else (
          p["cold_leg_corrosion_penalty"] if position == "cold_leg" else 0.0)
    gamma_corr = p["gamma_flow_corr"]
    theta_dT_corr = p["theta_dT_corr"]

    ln_rkin = (theta0
               + (Ea / R_GAS) * (1.0 / T_REF - 1.0 / T)
               + a_M + a_S + a_R + a_P
               + gamma_corr * math.log(flow)
               + theta_dT_corr * math.log1p(dT / 100.0))
    return math.exp(ln_rkin)

# Mass transfer cap
def r_mt(T, salt, flow, dT):
    theta_mt = p["log_mt_cap_um_y"]
    Ea_mt = p["Ea_mt_kJ_mol"] * 1000.0  # kJ/mol -> J/mol
    a_S = p.get(f"salt_{salt}", 0.0)
    gamma_mt = p["gamma_flow_mt"]
    theta_dT_mt = p["theta_dT_mt"]

    ln_rmt = (theta_mt
              + (Ea_mt / R_GAS) * (1.0 / T_REF - 1.0 / T)
              + 0.35 * a_S
              + gamma_mt * math.log(flow)
              + theta_dT_mt * math.log1p(dT / 100.0))
    return math.exp(ln_rmt)

# Harmonic mean
def r_corr(T, material, salt, redox, position, flow, dT):
    rk = r_kin(T, material, salt, redox, position, flow, dT)
    rm = r_mt(T, salt, flow, dT)
    return 1.0 / (1.0 / rk + 1.0 / rm)
    # Quick test: MOD-F-07, corrected temperature (973.15 K)

# Get dimensionless multiplier
def damage_multiplier(salt, redox):
    value = p["damage_base_log"]
    if salt == "chloride" or redox in ("impure_moisture", "chloride_unspecified"):
        value += p["damage_chloride_log"]
    if redox == "tellurium":
        value += p["damage_tellurium_log"]
    if redox == "stressed":
        value += p["damage_stress_log"]
    if redox == "multi_alloy":
        value += p["damage_multi_alloy_log"]
    # clip to [-10, 10] before exponentiating, matching expClip in the C++ code
    value = max(-10.0, min(10.0, value))
    return math.exp(value)

# Uniform attack depth
def x_corr(r_corr_val, time_years):
    return r_corr_val * time_years

# IGC crack depth
def x_igc(r_corr_val, time_years, salt, redox):
    xc = x_corr(r_corr_val, time_years)
    dmg = damage_multiplier(salt, redox)
    return dmg * xc + dmg * math.sqrt(max(xc, 0.0))

# Deposition rate
def r_dep(T, salt, redox, position, flow, dT, surface="metal"):
    beta0 = p["log_dep0_um_y"]
    Ea_dep = p["Ea_dep_kJ_mol"] * 1000.0  # kJ/mol -> J/mol
    a_S = p.get(f"salt_{salt}", 0.0)
    a_dep_S = p.get(f"dep_salt_{salt}", 0.0)
    a_R = p.get(f"redox_{redox}", 0.0)
    redox_cathodic = -0.30 * a_R
    gamma_dep = p["gamma_flow_dep"]
    theta_dT_dep = p["theta_dT_dep"]

    cold_bonus = p["dep_cold_bonus"] if position == "cold_leg" else 0.0

    surface_offset = 0.0
    if surface == "graphite":
        surface_offset += p["dep_graphite_offset"]
    elif surface == "turbulent_metal":
        surface_offset += p["dep_turbulent_bonus"]
    elif surface == "laminar_metal":
        surface_offset += p["dep_laminar_bonus"]
    if redox == "multi_alloy":
        surface_offset += p["dep_multi_alloy_bonus"]

    ln_rdep = (beta0
               + (Ea_dep / R_GAS) * (1.0 / T_REF - 1.0 / T)
               + 0.20 * a_S + a_dep_S
               + redox_cathodic 
               + gamma_dep * math.log(flow)
               + theta_dT_dep * math.log1p(dT / 100.0)
               + cold_bonus + surface_offset)
    return math.exp(ln_rdep)


# Corrosion depth
def corrosion_depth_um(r_corr_val, time_years):
    return r_corr_val * time_years  # same as x_corr


# Deposition depth
def deposition_depth_um(r_dep_val, time_years):
    return r_dep_val * max(time_years, 0.0)

def exp_clip(x, lo=-60.0, hi=60.0):
    return math.exp(max(lo, min(hi, x)))

def salt_cr_ppm_base(r_corr_val, time_years, material, experiment_family, source_id):
    depth = max(corrosion_depth_um(r_corr_val, time_years), 1.0e-9)
    cr_wt_frac = db["alloy_cr_wt_frac"][material]
    cr_factor = cr_wt_frac / 0.07
    is_msre = ("MSRE" in experiment_family) or source_id.startswith("ORNL-TM-3")
    scale_key = "log_ppm_scale_msre" if is_msre else "log_ppm_scale_loop"
    scale = exp_clip(p[scale_key])
    return scale * depth * cr_factor

def salt_cr_ppm(r_corr_val, time_years, material, salt, redox, experiment_family, source_id):
    ppm = salt_cr_ppm_base(r_corr_val, time_years, material, experiment_family, source_id)
    if source_id == "ORNL-TM-4188" and salt == "fluoride_fuel" and redox == "purified_baseline":
        ppm *= exp_clip(p["log_ncl16_cr_inventory_bonus"])
    return ppm


FARADAY = 96485.33212

def bv_overpotential_equivalent_v(redox, T, alpha_n=1.0):
    a_R = p.get(f"redox_{redox}", 0.0)
    return (a_R * R_GAS * T) / (alpha_n * FARADAY)



# Main -- Case MOD-F-07
time_years = 0.02737850787132101  # taken directly from the CSV column, not recomputed
#time_years = 240/8760 

test_rcorr = r_corr(
    T=973.15, material="hastelloy_n", salt="flinak",
    redox="purified_baseline", position="nominal",
    flow=0.35, dT=0.0
)
print("r_corr (corrected, 973.15 K):", test_rcorr)

test_rcorr_orig = r_corr(
    T=923.15, material="hastelloy_n", salt="flinak",
    redox="purified_baseline", position="nominal",
    flow=0.35, dT=0.0
)
print("r_corr (original, 923.15 K):", test_rcorr_orig)

xc_corrected = x_corr(test_rcorr, time_years)
xigc_corrected = x_igc(test_rcorr, time_years, "flinak", "purified_baseline")
print("x_corr (corrected):", xc_corrected)
print("x_igc (corrected):", xigc_corrected)

xc_orig = x_corr(test_rcorr_orig, time_years)
xigc_orig = x_igc(test_rcorr_orig, time_years, "flinak", "purified_baseline")
print("x_corr (original):", xc_orig)
print("x_igc (original):", xigc_orig)

test_rdep_orig = r_dep(923.15, "flinak", "purified_baseline", "nominal", 0.35, 0.0, surface="metal")
print("r_dep (original):", test_rdep_orig)

test_rdep_corrected = r_dep(973.15, "flinak", "purified_baseline", "nominal", 0.35, 0.0, surface="metal")
print("r_dep (corrected):", test_rdep_corrected)

test_ppm_orig = salt_cr_ppm(test_rcorr_orig, time_years, "hastelloy_n", "flinak", "purified_baseline", "Modern fluoride tests", "NPJ-HN-STRESS-FLINAK-2022")
print("salt_cr_ppm (original):", test_ppm_orig)

test_depdepth_orig = deposition_depth_um(test_rdep_orig, time_years)
print("deposition_depth_um (original):", test_depdepth_orig)

test_eta_orig = bv_overpotential_equivalent_v("purified_baseline", 923.15)
print("effective_eta_V (original):", test_eta_orig)

test_eta_corrected = bv_overpotential_equivalent_v("purified_baseline", 973.15)
print("effective_eta_V (corrected):", test_eta_corrected)

test_depdepth_corrected = deposition_depth_um(test_rdep_corrected, time_years)
print("deposition_depth_um (corrected):", test_depdepth_corrected)

test_ppm_corrected = salt_cr_ppm(test_rcorr, time_years, "hastelloy_n", "flinak", "purified_baseline", "Modern fluoride tests", "NPJ-HN-STRESS-FLINAK-2022")
print("salt_cr_ppm (corrected):", test_ppm_corrected)


r_corr_m035 = r_corr(973.15, "hastelloy_n", "flinak", "stressed", "nominal", 0.35, 0.0)
r_dep_m035 = r_dep(973.15, "flinak", "stressed", "nominal", 0.35, 0.0, surface="metal")
eta_m035 = bv_overpotential_equivalent_v("stressed", 973.15)
print("M-035 r_corr:", r_corr_m035)
print("M-035 r_dep:", r_dep_m035)
print("M-035 eta_eff:", eta_m035)

r_corr_m036 = r_corr(973.15, "gh3535", "flinak", "purified_baseline", "nominal", 0.35, 0.0)
r_dep_m036 = r_dep(973.15, "flinak", "purified_baseline", "nominal", 0.35, 0.0, surface="metal")
eta_m036 = bv_overpotential_equivalent_v("purified_baseline", 973.15)
print("M-036 r_corr:", r_corr_m036)
print("M-036 r_dep:", r_dep_m036)
print("M-036 eta_eff:", eta_m036)


# for flow in [0.1, 0.2, 0.35, 0.5, 0.75, 1.0]:
#     rc = r_corr(973.15, "hastelloy_n", "flinak", "purified_baseline", "nominal", flow, 0.0)
#     xigc = x_igc(rc, time_years, "flinak", "purified_baseline")
#     print(f"flow={flow:.2f}  r_corr={rc:.4f}  x_igc={xigc:.4f}")


# for dT in [0, 5, 10, 20]:
#     rc = r_corr(973.15, "hastelloy_n", "flinak", "purified_baseline", "nominal", 0.35, dT)
#     xigc = x_igc(rc, time_years, "flinak", "purified_baseline")
#     print(f"dT={dT:3d}  r_corr={rc:.4f}  x_igc={xigc:.4f}")

# dmg = damage_multiplier("flinak", "purified_baseline")
# print("damage_multiplier for MOD-F-07 (M-037):", dmg)

# dmg_values = [0.79, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0]  # sweeping the log-value directly


# for log_dmg in dmg_values:
#     dmg = math.exp(log_dmg)
#     xigc = dmg * xc_corrected + dmg * math.sqrt(max(xc_corrected, 0.0))
#     factor_error = 100.0 / xigc  # target is 100 (upper bound)
#     print(f"log_dmg={log_dmg:.2f}  dmg={dmg:8.2f}  x_igc={xigc:10.3f}  factor_error(vs 100)={factor_error:.3f}")

