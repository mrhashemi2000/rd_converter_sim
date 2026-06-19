import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import csv

# =============================================================================
# PROJECT: Quantum-Mechanistic and Stochastic Modeling for Ribose-to-2-Deoxyribose
# VERSION: 2.1 (Fixed Eyring + Improved Gillespie + Better Output)
# AUTHOR: Reza Hashemi (Updated by Grok)
# DESCRIPTION: Fixed version - realistic rates and stochastic behavior
# =============================================================================

# --- Physical Constants ---
R_GAS = 0.001987          # kcal/(mol·K)
H_PLANCK = 6.62607015e-34 # J·s
K_BOLTZMANN = 1.380649e-23 # J/K

# --- Simulation Parameters ---
T_LIST = [298.15, 310.15, 320.15]  # Temperatures to analyze (K)
DELTA_G_ACT = 28.0                 # Main activation free energy (kcal/mol)
F_PROT = 5.0                       # Protection factor of clay interlayer
K_DEG_RAW = 1e-8                   # Unprotected degradation rate (s^-1)
K_DEG = K_DEG_RAW / F_PROT         # Effective protected degradation rate
T_MAX_H = 72                       # Simulation time (hours)
T_MAX_S = T_MAX_H * 3600           # Simulation time (seconds)
N_TRAJECTORIES = 2000              # Trajectories per N_R0 (reduce to 500 for faster testing)
N_R0_RANGE = np.arange(1, 21)      # Range of initial molecules
SUCCESS_THRESHOLD = 0.20           # f_D > 0.20 for success


def calculate_k_main(barrier_energy, temp):
    """Calculates reaction rate constant using Eyring equation (fixed units)."""
    pre_exponential = (K_BOLTZMANN * temp) / H_PLANCK
    exponent = -barrier_energy / (R_GAS * temp)
    return pre_exponential * np.exp(exponent)


def run_gillespie_simulation(N_R0, k_main, k_deg):
    """Runs a single stochastic trajectory. Returns (success, final_fD)"""
    t = 0.0
    NR = N_R0
    ND = 0
    while t < T_MAX_S and (NR + ND) > 0:
        a1 = k_main * NR          # R → D
        a2 = k_deg * NR           # R degradation
        a3 = k_deg * ND           # D degradation
        a_total = a1 + a2 + a3
        
        if a_total <= 0:
            break
        
        # Stochastic time step
        dt = -np.log(np.random.random()) / a_total
        t += dt
        
        # Choose reaction
        r = np.random.random() * a_total
        if r < a1:
            NR -= 1
            ND += 1
        elif r < a1 + a2:
            NR -= 1
        else:
            ND -= 1
    
    final_total = NR + ND
    final_fD = ND / final_total if final_total > 0 else 0.0
    success = final_fD > SUCCESS_THRESHOLD
    return success, final_fD


def logistic_model(x, L, k, x0):
    """Logistic function for curve fitting."""
    return L / (1 + np.exp(-k * (x - x0)))


if __name__ == "__main__":
    print("Starting Multi-Temperature Stochastic Simulation (v2.1 Fixed)...\n")
    all_probs = []
    all_mean_fD = []
    
    plt.figure(figsize=(11, 7))
    colors = ['black', 'blue', 'red']
    
    for idx, T in enumerate(T_LIST):
        print(f"Simulating for T = {T:.2f} K...")
        k_main = calculate_k_main(DELTA_G_ACT, T)
        print(f"   k_main = {k_main:.2e} s⁻¹")
        
        successes = []
        mean_fDs = []
        
        for N_R0 in N_R0_RANGE:
            results = [run_gillespie_simulation(N_R0, k_main, K_DEG) for _ in range(N_TRAJECTORIES)]
            succ_count = sum(1 for s, _ in results if s)
            mean_fD = np.mean([f for _, f in results])
            
            successes.append(succ_count / N_TRAJECTORIES)
            mean_fDs.append(mean_fD)
        
        probs = np.array(successes)
        all_probs.append(probs)
        all_mean_fD.append(mean_fDs)
        
        # Logistic fit
        try:
            popt, _ = curve_fit(logistic_model, N_R0_RANGE, probs, p0=[1.0, 1.0, 8.0], maxfev=5000)
            threshold_50 = popt[2]
        except:
            threshold_50 = np.nan
        
        label_str = f'T = {T:.0f} K (L50 ≈ {threshold_50:.1f})'
        plt.plot(N_R0_RANGE, probs, 'o', color=colors[idx], label=label_str)
        plt.plot(N_R0_RANGE, logistic_model(N_R0_RANGE, *popt), '-', color=colors[idx], alpha=0.6)

    # --- Export Data to CSV ---
    csv_filename = 'simulation_results_T_analysis_v2.1.csv'
    with open(csv_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        header = ['N_R0'] + [f'P_success_T{t:.0f}K' for t in T_LIST] + [f'Mean_fD_T{t:.0f}K' for t in T_LIST]
        writer.writerow(header)
        for i, nr in enumerate(N_R0_RANGE):
            row = [nr]
            for j in range(len(T_LIST)):
                row.append(all_probs[j][i])
            for j in range(len(T_LIST)):
                row.append(all_mean_fD[j][i])
            writer.writerow(row)

    # --- Plot Formatting ---
    plt.axhline(y=0.5, color='gray', linestyle=':', alpha=0.6, label='50% Threshold')
    plt.xlabel('Initial Ribose Count (N_R0)', fontsize=13)
    plt.ylabel('Probability of Success (f_D > 0.20)', fontsize=13)
    plt.title('Effect of Temperature on Deoxyribose Formation Probability\n(72h Simulation on Fe-Montmorillonite)', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig('temp_analysis_300dpi_v2.1.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("\n" + "="*60)
    print(f"SUCCESS: Data exported to {csv_filename}")
    print(f"SUCCESS: Plot saved as temp_analysis_300dpi_v2.1.png")
    print("Note: You can tune DELTA_G_ACT (28–30) or F_PROT to better match the paper table.")
    print("="*60)
