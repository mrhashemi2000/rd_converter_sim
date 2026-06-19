Quantum-Mechanistic and Stochastic Modeling for Ribose-to-2-Deoxyribose Conversion on Fe-Doped Montmorillonite: Implications for Prebiotic DNA Precursor Emergence under Matter World Hypothesis (Theoretical)

The transition from RNA to DNA genetics requires the selective formation of 2-deoxyribose from ribose, a mechanistically challenging step in prebiotic chemistry. We propose a multi-scale computational model integrating periodic density functional theory (DFT) and stochastic Gillespie simulations to explore C2'-deoxygenation of ribose adsorbed on Fe²⁺/Fe³⁺-doped montmorillonite under UV irradiation. Periodic DFT calculations predict a proton-coupled electron transfer (PCET) pathway with an activation free energy ΔG‡ ≈ 28 kcal/mol for deoxyribose formation. Competing pathways (epimerization, ring-opening, and oxidation) exhibit barriers 7–12 kcal/mol higher, implying strong kinetic selectivity (rate ratios ~10⁵–10⁹ at 298 K). Time-dependent DFT identifies ligand-to-metal charge-transfer (LMCT) excited states in the UV-C range (~250–280 nm) suitable for photochemical initiation. Stochastic simulations in 1 fL clay interlayer compartments reveal non-linear threshold behavior: initial ribose molecule counts N_R0 ≈ 8 yield ~52% probability of deoxyribose fraction f_D > 0.20 after 72 hours (protection factor 5×, with effective degradation rate k_deg ≈ 10^{-8} s^{-1}). This corresponds to a local concentration of ~1.3 × 10^{-8} M. Sensitivity to ±2 kcal/mol barrier uncertainty shifts success probability by ~±22% at N_R0=10. This model, grounded in known clay-sugar interactions and UV photochemistry, generates testable predictions for Fe-montmorillonite catalysis and compartment effects. It supports a "Chemical Darwinism" framework of inevitable molecular selection under the Matter World Hypothesis.
Stochastic Simulation of Ribose-to-2-Deoxyribose Conversion

Project Overview
This repository contains a stochastic simulation based on the Gillespie algorithm to model the conversion of Ribose to 2-Deoxyribose on Fe-doped montmorillonite clay. The model accounts for the competitive rates of conversion and degradation under different temperature regimes.

Key Features
- Eyring-Polanyi Equation: Precise calculation of reaction rate constants based on activation free energy (Delta G^ddagger).
- Stochastic Modeling: Implementation of the Gillespie Direct Method to capture molecular fluctuations in small-volume compartments (femtoliiters).
- Temperature Analysis: Multi-temperature simulation to determine the L_50 threshold.
- Data Export:Automatic generation of CSV results for publication-quality plotting.

Simulation Parameters
- Delta G^\ddagger = 28.0  kcal/mol
- T = [298.15, 310.15, 320.15]  K
- Simulation Time: 72  hours
- Success Threshold: f_D > 0.20

How to Run
1. Install dependencies: `pip install -r requirements.txt`
2. Run the simulator: `python rd_converter_sim_v2_1.py`
3. The results will be saved as `simulation_results_T_analysis_v2.1.csv` and a high-res plot `temp_analysis_300dpi_v2.1.png`.

Citation
If you use this code in your research, please cite it using the Zenodo DOI provided in the release section.
