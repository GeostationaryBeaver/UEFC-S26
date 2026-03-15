import matplotlib.pyplot as plt
import numpy as np
from GetUEFC import UEFC

# SM.5(b) Code - Ricardo Ochoa
# Effective Elevator Range
elevator_range_eff = np.linspace(-10*(np.pi/180), 10*(np.pi/180), 1000) # Radians

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
a_w = (2*np.pi) / (1 + 2/AR)    # Equation 15

AR_h = b_h**2 / S_h             # From definition of AR (AR = b**2/S)

a_h = (2*np.pi) / (1 + 2/AR_h)  # Equation 20
zeta_e = np.arccos(1 - 2*f_e)   # Equation 22
a_e = 2*(np.pi - zeta_e + np.sin(zeta_e)) / (1 + 2/AR_h) # Equation 21
V_h = (S_h*l_h) / (S*c_bar)     # Equation 25

def get_alpha(alpha_e):
    # SM.5(a) Result
    alpha = -(V_h*(c_bar/l_h)*alpha_e) / ((a_w/a_e) + V_h*(c_bar/l_h)*(a_h/a_e))
    return alpha # Radians

plt.figure(1)
plt.plot((180/(np.pi))*elevator_range_eff, (180/(np.pi))*get_alpha(elevator_range_eff))
plt.title(r"$\alpha$ vs. $\alpha_{e}$")
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
plt.title(r"$C_{LW}$ and $C_{LH}$ vs. $\alpha_e$")
plt.xlabel(r"$\alpha_e$ (degrees)")
plt.ylabel("$C_{LW}$ and $C_{LH}$")
plt.legend()
plt.grid()


# SM.5(d) Code - Ricardo Ochoa

C_MW = C_MW_nom # Equation 17

def get_x_cg_frac(alpha_e):
    C_LH = get_C_LH(alpha_e)
    C_LW = get_C_LW(alpha_e)
    x_cg_frac = 1/4 + (V_h*C_LH - C_MW) / (C_LW + c_bar*V_h*C_LH/l_h) # The SM.1
    # Result when not Nominal
    return x_cg_frac

plt.figure(3)
plt.plot((180/(np.pi))*elevator_range_eff, get_x_cg_frac(elevator_range_eff),
        color="Green")
plt.title(r"$x_{cg}/\bar{c}$ vs. $\alpha_e^{trim}$")
plt.xlabel(r"$\alpha_e^{trim}$ (degrees)")
plt.ylabel(r"$x_{cg}/\bar{c}$")
plt.grid()


# SM.6 Code - Ricardo Ochoa

component_locations = [-0.16, 0.05, 0.20, 0.00, -0.13,
             0.28, 0.45, 0.03, 0.07, 0.75] # Meters
component_masses = [75, 75, 15, 8, 12,
          36, 24, 90, 25, 20] # Grams

m_empty = np.sum(component_masses)

x_cg_empty = np.dot(component_locations, component_masses)/m_empty # Meters
x_cg_empty_frac = x_cg_empty/c_bar
print(f"SM.6 | x_cg_empty = {x_cg_empty} Meters")
print(f"SM.6 | x_cg_empty_frac = {x_cg_empty_frac}")
# x_{cg}^{c_bar} is 0.44105, based on Figure SM.5d,
# the plane WILL FLY TRIMMED WHEN EMPTY around a_e^{trim} = -1.225 - Ricardo


# SM.7 Code - Ricardo Ochoa

m_total = m_pay + m_empty

def get_x_pay(alpha_e):
    x_cg = c_bar*get_x_cg_frac(alpha_e)
    x_pay = (x_cg*(m_total) - x_cg_empty*(m_empty)) / m_pay # Weighted Average
    return x_pay

x_nom_pay = get_x_pay(alpha_e=0) # alpha_e = 0 at Nominal Condition:
x_nom_pay_frac = x_nom_pay / c_bar

print(f"SM.7 | x_nom_pay: {x_nom_pay} Meters")
print(f"SM.7 | x_nom_frac: {x_nom_pay_frac}")


# SM.8 Code - Ricardo Ochoa

def get_delta_x_cg_frac(alpha_e):
    return get_x_cg_frac(alpha_e) - get_x_cg_frac(alpha_e=0) # alpha_e = 0 at
    # the Nominal Condition:

def get_delta_x_pay_frac(alpha_e):
    x_pay_frac = get_x_pay(alpha_e) / c_bar
    return x_pay_frac - x_nom_pay_frac

plt.figure(4)
plt.plot((180/(np.pi))*elevator_range_eff, get_delta_x_cg_frac(elevator_range_eff),
        color="Violet", label=r"$\Delta x_{cg}/\bar{c}$")
plt.plot((180/(np.pi))*elevator_range_eff, get_delta_x_pay_frac(elevator_range_eff),
        color="Green", label=r"$\Delta x_{pay}/\bar{c}$")
plt.title(r"$\Delta x_{pay}/\bar{c}$ vs. $\alpha_e^{trim}$")
plt.xlabel(r"$\alpha_e^{trim}$ (degrees)")
plt.ylabel(r"$\Delta x_{cg}/\bar{c}$ and $\Delta x_{pay}/\bar{c}$")
plt.grid()
plt.legend()


# SM.9 Code - Ricardo Ochoa
def get_SM(alpha_e):
    x_np_frac = (1/4 + V_h/(a_w/a_h + V_h*c_bar/l_h)) # From SM.3 Result
    x_cg_frac = get_x_cg_frac(alpha_e)
    SM = x_np_frac - x_cg_frac
    return SM

SM = get_SM(elevator_range_eff)

plt.figure(5)
plt.plot((180/(np.pi))*elevator_range_eff, SM,
        color="Brown")
plt.title(r"$SM$ vs. $\alpha_e^{trim}$")
plt.xlabel(r"$\alpha_e^{trim}$ (degrees)")
plt.ylabel(r"$SM$")
plt.grid()

# Identify alpha_a_neutral
alpha_e_neutral = elevator_range_eff[np.argmin(abs(SM))] # Radians

print(f"SM.9 | SM at Nominal Conditions: {get_SM(alpha_e=0)}")
print(f"SM.9 | alpha_e_neutral: {(180/np.pi)*alpha_e_neutral} Degrees")
delta_pay_neutral = c_bar*get_delta_x_pay_frac(alpha_e=alpha_e_neutral)
print(f"SM.9 | delta_pay: {delta_pay_neutral} Meters") # Meters

plt.show() # Displays all the 5 Figures
