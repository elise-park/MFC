#!/usr/bin/env python3
import math
import json

pi = 3.141592653589
patmos = 101325.0  # Pa

# gas properties
gammag = 1.4
Bg = 0.0
rhog = 1.18
c_g = 347.2
G_g = 0.0

# lung properties
gammal = 5.5
Bl = 492.0e6
rhol = 996.0
G_l = 10.0e3
c_l = math.sqrt(gammal * (patmos + Bl) / rhol)

P_amp = 1.0e6  # Shock pressure

leng = 1.0
Ny = 200
Nx = Ny * 3
dx = leng / Nx
vel = 500.0

time_end = 5 * leng / vel
cfl = 0.1
dt = cfl * dx / c_l
Nt = int(time_end / dt)
eps = 1e-8

alphal_back = 1.0 - eps
alphag_back = 1.0 - alphal_back
alphag_lung = 1.0 - eps
alphal_lung = 1.0 - alphag_lung

g_r = gammal
g_l = gammal
p_r = patmos
p_l = P_amp
B_r = Bl
r_r = rhol
c_r = c_l

rho_shock_l = (((g_r + 1)/(g_r - 1))*((p_l + B_r)/(p_r + B_r)) + 1)/(
    (g_r + 1)/(g_r - 1) + (p_l + B_r)/(p_r + B_r))*r_r
ushock_l = (c_r / g_r) * ((p_l / p_r) - 1) * (p_r / (p_r + B_r)) / math.sqrt(
    ((g_r + 1)/(2 * g_r)) * ((p_l / p_r) - 1) * (p_r / (p_r + B_r)) + 1
)


# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -leng / 2.0,
            "x_domain%end": leng / 2 + 2 * leng,
            "y_domain%beg": -leng / 2.0,
            "y_domain%end": leng / 2.0,
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": 10000,
            "t_step_save": 200,
            # Simulation Algorithm Parameters
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
            "riemann_solver": 1, #2
            "fd_order": 4,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -6,
            "bc_x%end": -6,
            "bc_y%beg": -6,
            "bc_y%end": -6,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Turning on hypoelasticity
            "hypoelasticity": "T",
            # Patch 1: Background
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 10 * leng,
            "patch_icpp(1)%length_y": leng,
            "patch_icpp(1)%vel(1)": 0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": patmos,
            "patch_icpp(1)%alpha_rho(1)": rhol * alphal_back,
            "patch_icpp(1)%alpha_rho(2)": rhog * alphag_back,
            "patch_icpp(1)%alpha(1)": alphal_back,
            "patch_icpp(1)%alpha(2)": alphag_back,
            # Patch 2: Shocked state
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -3 * leng / 8.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%length_x": leng / 4.0,
            "patch_icpp(2)%length_y": leng,
            "patch_icpp(2)%vel(1)": 0.12179783573189823, #ushock_l,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%pres": 301325.0, #P_amp,
            "patch_icpp(2)%alpha_rho(1)": 996.0735768378366, #rho_shock_l * alphal_back,
            "patch_icpp(2)%alpha_rho(2)": 1.1800000059292159e-08, #rhog * alphag_back,
            "patch_icpp(2)%alpha(1)": 0.99999999,  #alphal_back,
            "patch_icpp(2)%alpha(2)": 1.0000000050247593e-08, #alphag_back,
            # Patch 3: Lung
            "patch_icpp(3)%geometry": 21,
            "patch_icpp(3)%model_filepath": "voro_flat.stl",
            "patch_icpp(3)%model_spc": 10,
            "patch_icpp(3)%model_scale(1)": 8,
            "patch_icpp(3)%model_scale(2)": 8,
            "patch_icpp(3)%model_translate(1)":1.50,
            "patch_icpp(3)%model_translate(2)":-0.25,   
            "patch_icpp(3)%model_threshold": 0.99,
            "patch_icpp(3)%x_centroid": 0.0e00,
            "patch_icpp(3)%y_centroid": 0.0e00,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0e00,
            "patch_icpp(3)%pres": patmos,
            "patch_icpp(3)%alpha_rho(1)": rhol * alphal_back,
            "patch_icpp(3)%alpha_rho(2)": rhog * alphag_back,
            "patch_icpp(3)%alpha(1)": alphal_back,
            "patch_icpp(3)%alpha(2)": alphag_back,
            #Gas Background
            "patch_icpp(4)%geometry": 3,
            "patch_icpp(4)%alter_patch(1)": "T",
            "patch_icpp(4)%x_centroid": 1.25,
            "patch_icpp(4)%y_centroid": 0.0,
            "patch_icpp(4)%length_x":2.5* leng,
            "patch_icpp(4)%length_y": leng,
            "patch_icpp(4)%vel(1)": 0.0,
            "patch_icpp(4)%vel(2)": 0.0,
            "patch_icpp(4)%pres": patmos,
            "patch_icpp(4)%alpha_rho(1)": rhol * alphal_lung,
            "patch_icpp(4)%alpha_rho(2)": rhog * alphag_lung, #gas
            "patch_icpp(4)%alpha(1)": alphal_lung,
            "patch_icpp(4)%alpha(2)": alphag_lung,
            # Fluid physical parameters
            "fluid_pp(1)%gamma": 1.0 / (gammal - 1.0),
            "fluid_pp(1)%pi_inf": gammal * Bl / (gammal - 1.0),
            "fluid_pp(1)%G": G_l,
            "fluid_pp(2)%gamma": 1.0 / (gammag - 1.0),
            "fluid_pp(2)%pi_inf": gammag * Bg / (gammag - 1.0),
            "fluid_pp(2)%G": G_g,
        }
    )
)
