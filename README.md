# G.I.D.E.O.N. 🛰️
**GU**idance & **I**nertial **D**etermination **E**stimator for **O**rbital **N**avigation

A comprehensive and robust 3D simulation framework for spacecraft translational navigation written in MATLAB. This project implements a 9-state Position and Velocity Extended Kalman Filter (PV-EKF) designed to fuse high-rate inertial data (IMU) with low-rate GNSS observations. It serves as a state-of-the-art baseline for Guidance, Navigation, and Control (GNC) and Attitude and Orbit Control Systems (AOCS) analysis.

## 🚀 Features
* **Multi-Rate Sensor Fusion:** Seamlessly integrates 120 Hz IMU predictions with 1 Hz GNSS measurement updates.
* **High-Fidelity Orbital Dynamics:** The prediction step propagates the state using Runge-Kutta 4th order (RK4) integration, accounting for Central Body gravity and J2 (oblateness) perturbations.
* **Realistic IMU Error Modeling:** Simulates space-grade accelerometer and gyroscope behaviors, including stochastic noise (Velocity/Angle Random Walk, Bias Instability) and deterministic errors (Scale Factor, Non-orthogonality, and Mounting Misalignments).
* **Advanced GNSS Noise & Lever-Arm:** GNSS measurements are strictly modeled with high-frequency thermal noise and slowly-varying correlated errors (Gauss-Markov processes) to simulate ionospheric delay and clock wander. It also fully compensates for the kinematic Lever-Arm effect between the spacecraft's Center of Gravity and the antenna.
* **Dynamic Bias Estimation:** Real-time continuous estimation of the accelerometer's dynamic random walk bias.
* **Rigorous Statistical Validation:** Features a comprehensive visualization suite that automatically plots estimation errors against $\pm 3\sigma$ covariance bounds, computes the Normalized Estimation Error Squared (NEES) for filter consistency, and visualizes 3D orbital trajectories.

## 🛠️ Technology Stack
This project is built entirely in **MATLAB**, utilizing its core matrix manipulation capabilities. It does not strictly require any external or paid toolboxes (like the Navigation or Aerospace toolboxes), as all quaternion mathematics, Jacobians, and orbital mechanics are explicitly coded from scratch.

## 📊 Filter Architecture
The PV-EKF estimates a 9-state vector:
* **Position:** 3 states (ECI frame)
* **Velocity:** 3 states (ECI frame)
* **Accelerometer Bias:** 3 states (Body frame)

The discrete process noise matrix ($Q$) utilizes cross-correlated kinematic formulas to properly map the injection of specific force uncertainty into the velocity and position domains.

## ⚙️ Prerequisites & Execution

### Prerequisites
* **MATLAB:** Any modern release (R2020a or newer recommended). No specialized toolboxes are mandatory.

### How to Run
1. Clone the repository to your local machine.
2. Open MATLAB and navigate to the project root directory.
3. Open and run `mainEKF_Navigation.m`.
4. The simulation will generate the ground truth, synthesize the sensor data, execute the EKF loop, and automatically spawn 7 comprehensive diagnostic figures.

### ⚠️ Note on Tuning
If you modify the simulation duration or the GNSS Gauss-Markov parameters in `initializeGNSS.m`, you may need to adjust the noise floor and process noise tuning factors inside `initializeEKF_Nav.m` to prevent covariance collapse and maintain NEES consistency.

## 🎨 Acknowledgements
This architecture is inspired by standard aerospace industry practices for LEO (Low Earth Orbit) spacecraft navigation filtering.