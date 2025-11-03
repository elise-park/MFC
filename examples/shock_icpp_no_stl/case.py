#!/usr/bin/env python3
import json
import math

# --------------------------------------------------------------------
# Primitive variables
# --------------------------------------------------------------------
patmos = 101325.0  # Pa

# --------------------------------------------------------------------
# Material parameters
# --------------------------------------------------------------------
# Gas
gammag = 1.4
rhog = 1.18
Bg = 0.0

# Lung / liquid
gammal = 5.5
Bl = 492.0e6
rhol = 996.0
G_l = 5.0e3

# Shock amplitude
delta_P_s = 2.0e6 #2.0e6
P_amp = patmos + delta_P_s

# Compute shock velocity for lung-like fluid
c_l = math.sqrt(gammal * (patmos + Bl) / rhol)
g_r = gammal
p_r = patmos
p_l = P_amp
B_r = Bl
r_r = rhol

rho_shock_l = (
    (((g_r + 1) / (g_r - 1)) * ((p_l + B_r) / (p_r + B_r)) + 1)
    / ((g_r + 1) / (g_r - 1) + (p_l + B_r) / (p_r + B_r))
    * r_r
)
ushock_l = (c_l / g_r) * ((p_l / p_r) - 1) * (p_r / (p_r + B_r)) / math.sqrt(
    ((g_r + 1) / (2 * g_r)) * ((p_l / p_r) - 1) * (p_r / (p_r + B_r)) + 1
)

# --------------------------------------------------------------------
# Domain settings
# --------------------------------------------------------------------
D = 0.1
Nx, Ny = 199, 99
cfl = 0.25

# --------------------------------------------------------------------
# JSON Case
# --------------------------------------------------------------------
case = {
    "run_time_info": "T",
    "x_domain%beg": -3 * D,
    "x_domain%end": 3 * D,
    "y_domain%beg": -1.5 * D,
    "y_domain%end": 1.5 * D,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "dt": 1.0e-8,
    "t_step_start": 0,
    "t_step_stop": 10000,
    "t_step_save": 100,
    #"cfl_adap_dt": "T", #try adaptive
    "cfl_target": cfl,
    "num_patches": 3,
    "model_eqns": 2,
    "num_fluids": 2,
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "viscous": "F",
    "ib": "F",
    "num_ibs": 0,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",

    # Patch 1: Background medium
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": 1000 * D,
    "patch_icpp(1)%length_y": 1000 * D,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": patmos,
    "patch_icpp(1)%alpha_rho(1)": 0.0, #making this liquid for now
    "patch_icpp(1)%alpha_rho(2)": rhol, #
    "patch_icpp(1)%alpha(1)": 0.0,
    "patch_icpp(1)%alpha(2)": 1.0,

    # Patch 2: Shocked region
    "patch_icpp(2)%geometry": 3,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": -0.2745, #0.27
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%length_x": 0.051, #0.006
    "patch_icpp(2)%length_y": 1000 * D,
    "patch_icpp(2)%vel(1)": ushock_l, #made positive
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": P_amp,
    "patch_icpp(2)%alpha_rho(1)": 0.0,
    "patch_icpp(2)%alpha_rho(2)": rho_shock_l, #rho_shock_l
    "patch_icpp(2)%alpha(1)": 0.0,
    "patch_icpp(2)%alpha(2)": 1.0,

    # Patch 3: Gas (no stl)
    "patch_icpp(3)%geometry": 3,
    "patch_icpp(3)%alter_patch(1)": "T",
    "patch_icpp(3)%x_centroid": 0.0,
    "patch_icpp(3)%y_centroid": 0.0,
    "patch_icpp(3)%length_x": 0.5, #0.1
    "patch_icpp(3)%length_y": 0.1,
    "patch_icpp(3)%vel(1)": 0.0,
    "patch_icpp(3)%vel(2)": 0.0,
    "patch_icpp(3)%pres": patmos,
    "patch_icpp(3)%alpha_rho(1)": rhog,
    "patch_icpp(3)%alpha_rho(2)": 0.0,
    "patch_icpp(3)%alpha(1)": 1.0,
    "patch_icpp(3)%alpha(2)": 0.0,

    # Fluid properties
    "fluid_pp(2)%gamma": 1.0 / (gammal - 1.0),
    "fluid_pp(2)%pi_inf": gammal * Bl / (gammal - 1.0),
    "fluid_pp(2)%G": G_l,
    "fluid_pp(1)%gamma": 1.0 / (gammag - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%G": 0.0,
}

print(json.dumps(case, indent=2))
