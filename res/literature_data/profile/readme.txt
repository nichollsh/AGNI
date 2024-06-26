These files correspond to Figure 1 in Selsis+23.
The paper may be found here: https://www.nature.com/articles/s41586-023-06258-3
The data may be found here: https://zenodo.org/records/6877001

They represent the temperature profiles at RCE with two values for the OLR.
For the solar cases:
 - OTR274 was configured with Seff = 0.98 and Alb = 0.18.
 - OTR10000 was configured with Seff = 35.6 and Alb = 0.18.
Their Seff = 1 corresponds to a value of 341.5 Wm−2, which comes from:
    F_inst = S0 / 4, where S0=1366 Wm-2 is the solar constant.

This can be achieved by setting the following parameters in AGNI:
 - s0_fact = 0.6652   
 - zenith_angle = 60.0
 - albedo_b = 0.18
 - instellation = 1366 * S_eff

To solve for these profiles, we need to set the net flux boundary condition such
that the same OLR is achieved. This may be done via:
    F_int = sigma (T_int)^4 = OTR + F_UP_SW - F_DN_SW 

