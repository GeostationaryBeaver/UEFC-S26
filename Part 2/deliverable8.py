"""
Table Calculation Script - Stability and Control Analysis
Calculates aerodynamic and stability parameters for different flight conditions
"""

import numpy as np
import matplotlib.pyplot as plt
from SM_5_9 import *

# ============================================================================
# AIRCRAFT PARAMETERS
# ============================================================================

l_v = 0.485  # Vertical tail moment arm (meters)
l_h = 0.527  # Horizontal tail moment arm (meters)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def calculate_x_np():
    """
    Calculate neutral point position in meters
    x_np/c_bar = 1/4 + V_h / (a_w/a_h + V_h*c_bar/l_h)
    """
    x_np_frac = 1/4 + V_h / (a_w/a_h + V_h*c_bar/l_h)
    return x_np_frac * c_bar  # Convert to meters

def calculate_spiral_stability(CL_trim):
    """
    Calculate spiral stability parameter B
    B = (l_v / b) * (dihedral / CL_trim)

    Parameters:
    - CL_trim: trimmed lift coefficient

    Returns:
    - B: spiral stability parameter (positive = stable)
    """
    b_wing = b  # wingspan (1.34 m from SM_5_9.py)
    dihedral_deg = 12.2  # dihedral angle in degrees

    B = (l_v / b_wing) * (dihedral_deg / CL_trim)

    return B

def calculate_flight_condition(alpha_e_rad, condition_name, payload_mass=None):
    """
    Calculate all stability and control parameters for a flight condition at a given α_e

    Parameters:
    - alpha_e_rad: elevator angle in radians
    - condition_name: string describing condition
    - payload_mass: payload mass in grams (optional, uses global m_pay if None)

    Returns:
    - Dictionary with all parameters
    """

    print(f"\n{'='*60}")
    print(f"Calculating: {condition_name}")
    print(f"{'='*60}")
    print(f"Elevator angle: {np.degrees(alpha_e_rad):.2f}°")

    # Calculate aerodynamic coefficients at this α_e
    CL_trim = get_C_LW(alpha_e_rad)

    print(f"CL (wing): {CL_trim:.4f}")

    # Calculate CG, payload, and neutral point positions
    x_cg = c_bar * get_x_cg_frac(alpha_e_rad)
    x_pay = get_x_pay(alpha_e_rad)
    x_np = calculate_x_np()
    SM = get_SM(alpha_e_rad)
    B = calculate_spiral_stability(CL_trim)

    print(f"x_cg: {x_cg:.4f} m")
    print(f"x_pay: {x_pay:.4f} m")
    print(f"x_np: {x_np:.4f} m")
    print(f"SM: {SM:.4f}")
    print(f"B: {B:.4f}")

    return {
        'Condition': condition_name,
        'α_e (deg)': np.degrees(alpha_e_rad),
        'B': B,
        'x_cg (m)': x_cg,
        'x_pay (m)': x_pay,
        'x_np (m)': x_np,
        'SM': SM,
        'CL_trim': CL_trim
    }

# ============================================================================
# MAIN CALCULATION
# ============================================================================

if __name__ == "__main__":

    print("\n" + "="*80)
    print("FLIGHT CONDITION ANALYSIS - PLANE VANILLA")
    print("="*80)

    # Calculate results for each flight condition
    results = []

    # 1. Nominal Condition (α_e = 0)
    print("\n>>> Nominal Condition: Using α_e = 0°")
    results.append(calculate_flight_condition(
        alpha_e_rad=0.0,
        condition_name='Nominal'
    ))

    # 2. Speed Flight (no payload) - Find trim α_e
    print("\n>>> Speed Flight: Finding trim α_e (no payload)...")
    # For speed flight with no payload, find α_e that gives reasonable trim
    SM_values_speed = np.array([get_SM(ae) for ae in elevator_range_eff])
    idx_speed = np.argmin(np.abs(SM_values_speed))
    alpha_e_speed = elevator_range_eff[idx_speed]
    results.append(calculate_flight_condition(
        alpha_e_rad=alpha_e_speed,
        condition_name='Speed Flight'
    ))

    # 3. Aft Pay Flight (payload moved back) - Find trim α_e
    print("\n>>> Aft Pay Flight: Finding trim α_e (payload aft)...")
    # For aft pay flight, find α_e that gives reasonable trim
    SM_values_aft = np.array([get_SM(ae) for ae in elevator_range_eff])
    idx_aft = np.argmin(np.abs(SM_values_aft))
    alpha_e_aft = elevator_range_eff[idx_aft]
    results.append(calculate_flight_condition(
        alpha_e_rad=alpha_e_aft,
        condition_name='Aft Pay Flight'
    ))

    # ========================================================================
    # CREATE OUTPUT TABLE
    # ========================================================================

    print("\n" + "="*80)
    print("RESULTS TABLE")
    print("="*80 + "\n")

    # Print formatted table matching your original format
    print(f"{'Condition':<20} {'α_e (°)':<12} {'B':<12} {'x_cg (m)':<12} "
          f"{'x_pay (m)':<12} {'x_np (m)':<12} {'SM':<12}")
    print("-" * 100)

    for result in results:
        print(f"{result['Condition']:<20} {result['α_e (deg)']:<12.2f} "
              f"{result['B']:<12.4f} {result['x_cg (m)']:<12.4f} "
              f"{result['x_pay (m)']:<12.4f} {result['x_np (m)']:<12.4f} "
              f"{result['SM']:<12.4f}")

    # Verify Nominal matches your known values
    print("\n" + "="*80)
    print("VERIFICATION - Nominal Condition")
    print("="*80)
    nominal = results[0]
    print(f"Your known values:")
    print(f"  x_cg = 0.0733 m  (calculated: {nominal['x_cg (m)']:.4f} m)")
    print(f"  x_pay = 0.035 m  (calculated: {nominal['x_pay (m)']:.4f} m)")
    print(f"  x_np = 0.0351 m  (calculated: {nominal['x_np (m)']:.4f} m)")
    print(f"  SM = 0.118       (calculated: {nominal['SM']:.4f})")

    # Save to CSV for easy table formatting
    try:
        import pandas as pd
        df = pd.DataFrame(results)
        df.to_csv('stability_table.csv', index=False)
        print("\n✓ Results saved to 'stability_table.csv'")
    except ImportError:
        print("\n(pandas not available - skipping CSV export)")

    # Summary comparison
    print("\n" + "="*80)
    print("KEY OBSERVATIONS")
    print("="*80)

    nominal = results[0]
    speed = results[1]
    aft = results[2]

    print(f"\nSpeed Flight vs Nominal:")
    print(f"  Δα_e = {speed['α_e (deg)'] - nominal['α_e (deg)']:+.2f}°")
    print(f"  ΔSM = {speed['SM'] - nominal['SM']:+.4f}")
    print(f"  ΔB = {speed['B'] - nominal['B']:+.4f}")

    print(f"\nAft Pay Flight vs Nominal:")
    print(f"  Δα_e = {aft['α_e (deg)'] - nominal['α_e (deg)']:+.2f}°")
    print(f"  Δx_pay = {(aft['x_pay (m)'] - nominal['x_pay (m)'])*1000:+.2f} mm")
    print(f"  ΔSM = {aft['SM'] - nominal['SM']:+.4f}")
    print(f"  ΔB = {aft['B'] - nominal['B']:+.4f}")

    print("\n" + "="*80)
