# G.I.D.E.O.N. 🛰️
**G**uidance & **I**nertial **D**etermination **E**stimator for **O**rbital **N**avigation

A comprehensive, high-fidelity 3D simulation framework for spacecraft Guidance, Navigation, and Control (GNC) written entirely in MATLAB. This repository provides a complete suite of algorithms for translational navigation, attitude estimation, static attitude determination, and rigorous sensor error modeling.

## 🚀 Core Features
### 1. Multi-Rate Sensor Fusion (MEKF & PV-EKF)
* **Attitude Estimation (MEKF):** A Multiplicative Extended Kalman Filter estimating a 6-state error vector (Angle Error and Gyro Bias). It seamlessly fuses high-rate 120 Hz Gyro kinematics (using Hamilton passive quaternions) with asynchronous updates from a 10 Hz Magnetometer and 4 Hz Star Trackers.
* **Translational Navigation (PV-EKF):** A 9-state Position and Velocity Extended Kalman Filter fusing 120 Hz IMU specific forces with 1 Hz GNSS updates. It features Runge-Kutta 4 (RK4) orbital integration including Central Body and J2 oblateness perturbations.

### 2. Static Attitude Determination (Wahba's Problem)
* **QUEST Algorithm:** Implements the optimal QUaternion ESTimator eigenvalue solver for multi-star observations using real stellar data from the Hipparcos catalog. Validates performance against the Cramér-Rao Lower Bound (CRLB).
* **TRIAD Algorithm:** Provides a deterministic, two-vector algebraic solution. Includes a performance comparison against QUEST, demonstrating TRIAD's optimality when the primary vector weight heavily dominates the secondary.

### 3. High-Fidelity Sensor Modeling & Validation
* **Inertial Measurement Unit (IMU):** Simulates space-grade MEMS/FOG specifications including Angle/Velocity Random Walk (ARW/VRW), Rate Random Walk (RRW), and Bias Instability. Includes an empirical **Allan Deviation** analysis script to characterize noise and recommend EKF tuning parameters.
* **Global Navigation Satellite System (GNSS):** Models the physical Lever Arm kinematic coupling, high-frequency white noise, and temporally correlated errors (Gauss-Markov processes representing ionospheric and clock drifts).
* **Magnetometer (MAG):** Simulates 3-axis magnetic field measurements with hard/soft iron deterministic errors, non-orthogonality, and bias instability.

### 4. Rigorous Statistical Analysis
* Automatically computes and plots estimation errors against ±3σ covariance bounds.
* Evaluates filter consistency using the Normalized Estimation Error Squared (NEES) for chi-squared confidence intervals.

## 🛠️ Technology Stack
This project uses **MATLAB** to perform high-fidelity 3D simulations, matrix manipulations, and complex filtering algorithms. While the core GNC algorithms (MEKF, PV-EKF, QUEST, TRIAD) and quaternion mathematics are explicitly coded from scratch without external dependencies, the **Aerospace Toolbox** is utilized exclusively to generate the rigorous environmental ground truth (such as the IGRF-13 magnetic field model and NRLMSISE-00 atmospheric density).

## 📂 Repository Structure
* `/Data` - Stores the cached Hipparcos Stellar Catalog used for the Star Tracker simulations.
* `/Documentation` - Mathematical derivations, algorithm formulations, and academic references.
* `/Figures` - Destination folder for the automatically generated diagnostic plots and performance metrics (saved in `.fig`, `.png`, and `.svg` formats).
* `/Functions` - Helper functions containing all core mathematical operations, quaternion kinematics, state transition matrices, Jacobians, and environmental models.
* `mainEKF_*.m` - Main executable scripts for the Multiplicative EKF (Attitude) and Position-Velocity EKF (Navigation).
* `mainQUEST.m` & `mainTRIAD.m` - Executable scripts to validate the static attitude determination algorithms.
* `mainGNSS.m`, `mainIMU.m`, `mainMAG.m` - Executable scripts to synthesize sensor measurements, characterize noise (e.g., Allan Variance), and validate error models.

## ⚙️ Installation & Usage
Since this project is developed entirely in MATLAB, no external compilation or complex dependency management is required. 

1. **Clone the repository:**
   ```bash
   git clone https://github.com/BFG508/GIDEON.git
2. **Open MATLAB** and navigate to the cloned `GIDEON` directory.
3. **Add to Path**: Ensure that all subdirectories (`/Functions` and `/Data`) are added to your MATLAB path so the main scripts can access the helper functions and the stellar catalog. You can do this by right-clicking the `GIDEON` folder in the Current Folder browser and selecting Add to Path > Selected Folders and Subfolders.
4. **Run the Analysis**: Open and run any of the `main*.m` scripts located in the root directory. The scripts will automatically generate the orbital and attitude ground truth, synthesize the respective sensor data with noise, execute the estimation algorithms, and spawn comprehensive diagnostic figures.