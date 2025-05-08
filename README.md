# 📉 Robust Kalman-Bucy Filter: Continuous-Time State Estimation under Uncertainty

This project implements a Robust Kalman-Bucy Filter (RKBF) using MATLAB to replicate and validate the theoretical estimator proposed by Jemin George (IEEE Trans. Auto. Control, 2013). The robust formulation ensures optimal performance under persistent excitation, parametric uncertainty, and noise, outperforming the classical Kalman-Bucy Filter in uncertain conditions.

The project includes matrix setup, simulation of stochastic disturbances, Riccati equation solving, Linear Matrix Inequality (LMI) formulation, and detailed graph-based analysis.

## 🎯 Project Objective
To simulate and validate the robust estimator proposed in the referenced paper, demonstrating:

* Robust state estimation under parametric uncertainty

* LMI-based computation of observer gains

* Estimation error convergence with and without disturbances

* Theoretical consistency with standard Kalman-Bucy filter in nominal conditions

## ⚙️ Key Features
✅ Kalman Gain Derivation
* Classical steady-state Kalman gain via CARE (using icare)

* LMI-based validation of stability and robustness

✅ Stochastic Signal Simulation
* Persistent disturbances (W1, W2) generated from Dirichlet and square waveforms

* Non-zero mean excitation mimicking real sensor/system uncertainties

✅ Estimator Dynamics
* Implements error dynamics for both Kalman-Bucy and robust estimator

* Uses Riccati ODE integration to compute error covariance over time

✅ Graphical Outputs
* Covariance plots for states x₁–x₄

* Estimation error norms (Kalman vs. Robust)

* Time-series evolution of robustness under noisy conditions

## 🧱 System Overview
```
Stochastic System: dx = (A x + Wt)dt + dBt
               dy = C x dt + dVt

Add Uncertainties → Define Am, Cm
           ↓
Derive Kalman Gain (CARE)
           ↓
Form Augmented State: Z = [X; NoiseStates]
           ↓
Apply Riccati ODE Solver
           ↓
Solve LMI for Robust Filter Conditions
           ↓
Plot Estimation Covariance and Norms
```
## 📋 Key Code Features
* Modular block structure with clear labels

* Vectorized matrix operations

* Simulates continuous-time stochastic systems

* Uses MATLAB built-in care, ode45, and optionally lmi solvers

* Saves plots and overlays RKBF vs Kalman filters visually

## 🧪 Results Summary
 Kalman filter brings estimation error to 0 in nominal case

* Robust estimator shows comparable accuracy under noise-free conditions

* Under disturbance, RKBF performs significantly better than standard Kalman

* Some divergence from original paper's figures due to undocumented matrix paddings and ambiguous notation
