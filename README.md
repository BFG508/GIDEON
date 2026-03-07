# G.I.D.E.O.N. 🛰️
**GU**idance & **I**nertial **D**etermination **E**stimator for **O**rbital **N**avigation

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
* Automatically computes and plots estimation errors against $\pm3\sigma$ covariance bounds.
* Evaluates filter consistency using the Normalized Estimation Error Squared (NEES) for chi-squared confidence intervals.

## 📂 Repository Structure
The suite is modularly divided into several executable `main*.m` scripts. Run any of these to execute a specific simulation:
* `mainEKF_Navigation.m`: Runs the PV-EKF for orbital position, velocity, and accelerometer bias estimation.
* `mainEKF_Attitude.m`: Runs the MEKF for 3D orientation and gyroscope bias estimation.
* `mainQUEST.m`: Validates Star Tracker attitude determination using the QUEST algorithm.
* `mainTRIAD.m`: Validates Star Tracker attitude determination using the TRIAD algorithm.
* `mainIMU.m`: Generates IMU measurements and performs Allan Variance analysis for drift characterization.
* `mainGNSS.m`: Simulates GNSS measurements and analyzes the Lever Arm effect and Gauss-Markov noise.
* `mainMAG.m`: Simulates environmental magnetic fields and magnetometer errors.

## 🛠️ Technology Stack
This project is built entirely in **MATLAB**, utilizing its core matrix manipulation and plotting capabilities. It does not strictly require any external or paid toolboxes, as all quaternion mathematics, Jacobian matrices, sensor kinematic models, and eigenvalue decompositions are explicitly coded from scratch.

## ⚙️ How to Run
1. Clone the repository to your local machine.
2. Open MATLAB and navigate to the project root directory.
3. Open any of the `main*.m` scripts (e.g., `mainEKF_Navigation.m`) and click **Run**.
4. The simulation will generate the ground truth, synthesize the respective sensor data, execute the algorithms, and automatically spawn comprehensive diagnostic figures.