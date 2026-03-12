import matplotlib.pyplot as plt
import numpy as np
from GetUEFC import UEFC

# SM.5(b) Code - Ricardo Ochoa
# Effective Elevator Range
elevator_range_eff = np.linspace(-10*(np.pi/180), 10*(np.pi/180), 200) # Radians

# Plane Vanilla Constants from Page 4
S = 0.354 # Meters
AR = 9
S_h = UEFC().Sh

# Plane Vanilla Constants from Page 10
c_bar = 0.2 # Meters
b = 1.77 # Meters
C_LW_nom = 0.65
C_MW_nom = -0.15
c_bar_h = 0.0775 # Meters
b_h = 0.525 # Meters
f_e = 0.6
l_h = 0.65 # Meters
m_pay = 250 # Grams

# Derived Quantities from Constants
a_w = (2*np.pi) / (1 + 2/AR)  # Equation 15

AR_h = b_h**2 / S_h           # From definition of AR (AR = b**2/S)

a_h = (2*np.pi) / (1 + 2/AR_h)  # Equation 20
zeta_e = np.arccos(1 - 2*f_e)   # Equation 22
a_e = 2*(np.pi - zeta_e + np.sin(zeta_e)) / (1 + 2/AR_h) # Equation 21
V_h = (S_h*l_h) / (S*c_bar)         # Equation 25

def get_alpha(alpha_e):
    alpha = -(V_h*(c_bar/l_h)*alpha_e) / ((a_w/a_e) + V_h*(c_bar/l_h)*(a_h/a_e))
    return alpha # Radians

plt.figure(1)
plt.plot((180/(np.pi))*elevator_range_eff, (180/(np.pi))*get_alpha(elevator_range_eff))
plt.title(r"$\alpha$ vs $\alpha_{e}$")
plt.xlabel(r"$\alpha_e$ (degrees)")
plt.ylabel(r"$\alpha$ (degrees)")
plt.grid()

# SM.5(c) Code - Ricardo Ochoa

def get_C_LW(alpha_e):
    C_LW = C_LW_nom + a_w*get_alpha(alpha_e) # Equation 14
    return C_LW

def get_C_LH(alpha_e):
    C_LH = a_h*get_alpha(alpha_e) + a_e*alpha_e # Equation 19
    return C_LH

plt.figure(2)
plt.plot((180/(np.pi))*elevator_range_eff, get_C_LW(elevator_range_eff),
         label="$C_{LW}$", color="Blue")
plt.plot((180/(np.pi))*elevator_range_eff, get_C_LH(elevator_range_eff),
         label="$C_{LH}$", color="Orange")

plt.title(r"$C_{LW}$ and $C_{LH}$ vs $\alpha_e$")
plt.xlabel(r"$\alpha_e$ (degrees)")
plt.ylabel("$C_{LW}$ and $C_{LH}$")
plt.legend()
plt.grid()

# SM.5(d) Code - Ricardo Ochoa

C_MW = C_MW_nom # Equation 17

def get_x_cg_over_c_bar(alpha_e):
    x_cg_over_c_bar = 1/4 - C_MW/get_C_LW(alpha_e) # SM.5(a) Result
    return x_cg_over_c_bar

plt.figure(3)
plt.plot((180/(np.pi))*elevator_range_eff, get_x_cg_over_c_bar(elevator_range_eff),
        color="Green")

plt.title(r"$x_{cg}/\bar{c}$ vs $\alpha_e^{trim}$")
plt.xlabel(r"$\alpha_e^{trim}$ (degrees)")
plt.ylabel(r"$x_{cg}/\bar{c}$")
plt.grid()
plt.show()
