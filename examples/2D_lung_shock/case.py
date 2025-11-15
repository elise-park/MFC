#!/usr/bin/env python3
import math
import json

pi = 3.141592653589
patmos = 101325.0

#gas parameters
gammag = 1.4    # gas
Bg = 0.0        # Pa
rhog = 1.18     # kg/m^3
c_g = 347.2     # m/s
G_g = 0.0       # Pa

#liquid parameters
gammal = 5.5
Bl = 492.0e6    # Pa
rhol = 996.0    # kg/m^3
G_l = 5.0e3     # Pa
Cv_l = 1816.0

#shock parameters
lambda_wave = 200.0e-6
delta_P_s = 200000.0
P_amp = delta_P_s + patmos

c_l = math.sqrt(gammal * (patmos + Bl) / rhol)

g_r = gammal
g_l = gammal
p_r = patmos
p_l = P_amp
B_r = Bl
B_l = Bl
r_r = rhol
c_r = c_l

rho_shock_l = (
    (((g_r + 1) / (g_r - 1)) * ((p_l + B_r) / (p_r + B_r)) + 1)
    / (((g_r + 1) / (g_r - 1)) + (p_l + B_r) / (p_r + B_r))
) * r_r

ushock_l = (c_r / g_r) * ((p_l / p_r) - 1) * (p_r / (p_r + B_r)) / math.sqrt(
    ((g_r + 1) / (2 * g_r)) * ((p_l / p_r) - 1) * (p_r / (p_r + B_r)) + 1
)

ps = 248758.567
gam = 1.4
rho = 1.0
c_l_for_dt = math.sqrt(1.4 * ps / rho)
vel = 500.0

leng = 1.0
Ny = 200.0
Nx = Ny * 3.0
dx = leng / Nx

time_end = 5.0 * leng / vel
cfl = 0.02 #0.1
dt = cfl * dx / c_l_for_dt
Nt = int(time_end / dt)

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            "x_domain%beg": -leng / 2.0,
            "x_domain%end": leng / 2 + 2.0 * leng,
            "y_domain%beg": -leng / 2.0,
            "y_domain%end": leng / 2.0,
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": 10000, #10000
            "t_step_save": 100,#100
            # simulation algorithm parameters
            "num_patches": 4,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 1,
            "wave_speeds": 1,
            "avg_state": 2,
            "fd_order": 2,
            "bc_x%beg": -7,
            "bc_x%end": -8,
            "bc_y%beg": -2,
            "bc_y%end": -2,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Turning on hypoelasticity
            #"hypoelasticity": "T",
            # Patch 1: Background Patch (Liquid)
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 10 * leng,
            "patch_icpp(1)%length_y": leng,
            "patch_icpp(1)%vel(1)": 0,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": patmos,
            "patch_icpp(1)%alpha_rho(1)": rhol,
            "patch_icpp(1)%alpha_rho(2)": 0.0e00,
            "patch_icpp(1)%alpha(1)": 1.0e00,
            "patch_icpp(1)%alpha(2)": 0.0e00,
            # Patch 2: Shock (Liquid)
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -3 * leng / 8.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%length_x": leng / 4.0,
            "patch_icpp(2)%length_y": leng,
            "patch_icpp(2)%vel(1)": ushock_l,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%pres": P_amp,
            "patch_icpp(2)%alpha_rho(1)": rho_shock_l,
            "patch_icpp(2)%alpha_rho(2)": 0.0e00,
            "patch_icpp(2)%alpha(1)": 1.0e00,
            "patch_icpp(2)%alpha(2)": 0.0e00,
            # Patch 3: Gas Background (Behind STL)
            "patch_icpp(3)%geometry": 3,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%x_centroid": 1.25,
            "patch_icpp(3)%y_centroid": 0.0,
            "patch_icpp(3)%length_x": 2.5 * leng,
            "patch_icpp(3)%length_y": leng,
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%pres": patmos,
            "patch_icpp(3)%alpha_rho(1)": 0.0,
            "patch_icpp(3)%alpha_rho(2)": rhog, #gas
            "patch_icpp(3)%alpha(1)": 0.0e00,
            "patch_icpp(3)%alpha(2)": 1.0e00,
            # Patch 4: STL (Liquid)
            "patch_icpp(4)%geometry": 21,
            "patch_icpp(4)%model_filepath": "voro_flat.stl",
            "patch_icpp(4)%model_spc": 10,
            "patch_icpp(4)%model_scale(1)": 12, #8
            "patch_icpp(4)%model_scale(2)": 12, #8
            "patch_icpp(4)%model_translate(1)":2.5, #1.25
            "patch_icpp(4)%model_translate(2)":-0.25,   
            "patch_icpp(4)%model_threshold": 0.99,
            "patch_icpp(4)%x_centroid": 0.0e00,
            "patch_icpp(4)%y_centroid": 0.0e00,
            "patch_icpp(4)%alter_patch(1)": "T",
            "patch_icpp(4)%vel(1)": 0.0,
            "patch_icpp(4)%vel(2)": 0.0e00,
            "patch_icpp(4)%pres": patmos,
            "patch_icpp(4)%alpha_rho(1)": rhol,
            "patch_icpp(4)%alpha_rho(2)": 0.0,
            "patch_icpp(4)%alpha(1)": 1.0e00,
            "patch_icpp(4)%alpha(2)": 0.0e00,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0 / (gammal - 1.0),
            "fluid_pp(1)%pi_inf": gammal * Bl / (gammal - 1.0),
            "fluid_pp(1)%G": G_l,
            "fluid_pp(2)%gamma": 1.0 / (gammag - 1.0),
            "fluid_pp(2)%pi_inf": gammag * Bg / (gammag - 1.0),
            "fluid_pp(2)%G": G_g,
        }
    )
)
