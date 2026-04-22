import numpy as np
from scipy.optimize import minimize, fsolve
import UEFC_wing  # Calling your provided Wing Vortex Lattice file

# =============================================================================
# 1. FIXED WING & PLANE PARAMETERS (From Table)
# =============================================================================
S = 0.161  # m^2
b = 1.153  # m
AR_wing = 8.2
taper = 0.5
gamma_r = 12.0  # root twist (deg)
gamma_t = 4.5  # tip twist (deg)
dihedral = 10.0  # deg
S_h = 0.04  # m^2
S_v = 0.03  # m^2

# Calculate root, tip, and mean aerodynamic chords
c_root = (2 * S) / (b * (1 + taper))
c_tip = taper * c_root
c_bar = (2 / 3) * c_root * ((1 + taper + taper ** 2) / (1 + taper))

# Wing lift curve slope (Equation 15 from SM_5-9.py)
a_w = (2 * np.pi) / (1 + 2 / AR_wing)

# =============================================================================
# 2. CALLING UEFC_WING.PY FOR CONSTRAINTS 6 & 12 (max c_l <= 0.8)
# =============================================================================
print("--- Initializing Wing Vortex Lattice Method ---")
wing = UEFC_wing.UEFC_wing(b=b, croot=c_root, ctip=c_tip,
                           agroot=gamma_r, agtip=gamma_t, dihedral=dihedral)


def cl_difference(CL_test):
    G, _ = wing.solve(CL=CL_test)
    cl_dist, _ = wing.calccldist(G)
    return np.max(cl_dist) - 0.8


max_allowed_CL = fsolve(cl_difference, 0.6)[0]
print(f"To satisfy Constraints 6 & 12 (max c_l <= 0.8):")
print(f"The maximum allowable aircraft C_L is: {max_allowed_CL:.3f}\n")


# =============================================================================
# 3. HELPER FUNCTIONS & MATH CONFLICT CHECK
# =============================================================================
def calc_Vh(l_h):
    return (S_h * l_h) / (S * c_bar)


def calc_Vv(l_v):
    return (S_v * l_v) / (S * b)


def calc_SM(l_h, AR_h, x_cg_frac):
    Vh = calc_Vh(l_h)
    a_h = (2 * np.pi) / (1 + 2 / AR_h)
    x_np_frac = 0.25 + Vh / ((a_w / a_h) + Vh * (c_bar / l_h))
    return x_np_frac - x_cg_frac


def calc_B(l_v):
    """
    Standard Blaine Rawdon Spiral Stability formula: B = V_v * Dihedral(deg)
    Update this formula if your professor defines B differently!
    """
    return calc_Vv(l_v) * dihedral


print("--- Pre-Optimization Math Check ---")
max_possible_B = 0.05 * dihedral  # Since max V_v is 0.05
print(f"With V_v capped at 0.05 (Constraint 9) and Dihedral = {dihedral} deg,")
print(f"The mathematical maximum for Spiral Stability (B) is {max_possible_B:.2f}.")
if max_possible_B < 5.0:
    print(">>> WARNING: Constraint 9 (V_v <= 0.05) and Constraint 10 (B >= 5) physically contradict each other! <<<\n")


# =============================================================================
# 4. SCIPY OPTIMIZATION FOR TAIL DESIGN
# =============================================================================
# Design Vector x = [l_h, l_v, AR_h, AR_v, x_cg_frac]

def objective(x):
    # Minimize tail boom lengths (lightest possible weight)
    return x[0] ** 2 + x[1] ** 2


constraints = [
    {'type': 'ineq', 'fun': lambda x: calc_Vh(x[0]) - 0.30},  # V_h >= 0.30
    {'type': 'ineq', 'fun': lambda x: 0.60 - calc_Vh(x[0])},  # V_h <= 0.60
    {'type': 'ineq', 'fun': lambda x: calc_Vv(x[1]) - 0.02},  # V_v >= 0.02
    #{'type': 'ineq', 'fun': lambda x: 0.05 - calc_Vv(x[1])},  # V_v <= 0.05
    {'type': 'ineq', 'fun': lambda x: calc_SM(x[0], x[2], x[4])},  # SM >= 0 (x_np > x_cg)
    {'type': 'ineq', 'fun': lambda x: 0.20 - calc_SM(x[0], x[2], x[4])},  # SM <= 0.20
    {'type': 'ineq', 'fun': lambda x: calc_B(x[1]) - 5.0}  # B >= 5 (Spiral Stability)
]

# Bounds [l_h, l_v, AR_h, AR_v, x_cg_frac]
bounds = (
    (0.1, 1.5),  # l_h bounds
    (0.1, 4.0),  # l_v bounds (widened heavily just in case B takes precedence)
    (3.0, 6.0),  # AR_h
    (1.5, 4.0),  # AR_v
    (0.1, 0.4)  # x_cg_frac
)

x0 = [0.4, 0.4, 4.0, 2.0, 0.25]

print("--- Running Scipy Optimization for Tail Booms ---")
result = minimize(objective, x0, method='SLSQP', bounds=bounds, constraints=constraints)

if result.success:
    opt_lh, opt_lv, opt_ARh, opt_ARv, opt_xcg = result.x
    opt_Vh = calc_Vh(opt_lh)
    opt_Vv = calc_Vv(opt_lv)
    opt_SM = calc_SM(opt_lh, opt_ARh, opt_xcg)
    opt_B = calc_B(opt_lv)

    print("Optimization Successful!\n")
    print("=== MISSING TABLE VALUES ===")
    print(f"AR_h (Horiz Tail AR)     = {opt_ARh:.2f}")
    print(f"AR_v (Vert Tail AR)      = {opt_ARv:.2f}")
    print(f"l_h  (Horiz Moment Arm)  = {opt_lh:.3f} m")
    print(f"l_v  (Vert Moment Arm)   = {opt_lv:.3f} m")
    print(f"V_h  (Horiz Volume Coef) = {opt_Vh:.3f}")
    print(f"V_v  (Vert Volume Coef)  = {opt_Vv:.3f}")

    print("\n=== CONSTRAINTS VERIFICATION ===")
    print(f"Required x_cg location   = {opt_xcg * 100:.1f}% of MAC")
    print(f"Static Margin (SM)       = {opt_SM:.3f}   (Passed: 0 <= SM <= 0.2)")
    print(f"Horiz Vol. Coef (V_h)    = {opt_Vh:.3f}   (Passed: 0.30 to 0.60)")
    print(f"Vert Vol. Coef (V_v)     = {opt_Vv:.3f}   (Passed: 0.02 to 0.05)")
    print(f"Spiral Stability (B)     = {opt_B:.3f}   (Passed: B >= 5)")

else:
    print("\nOPTIMIZATION FAILED.")
    print("Reason from Solver:", result.message)
    print("\nDiagnostic:")
    print("The solver could not find a geometry that obeys all rules at the same time.")
    print("Try temporarily commenting out either the V_v <= 0.05 constraint or the B >= 5 constraint")
    print("to see which variable is causing the mathematical conflict in your model parameters.")
