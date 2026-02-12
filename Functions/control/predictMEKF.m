function ekf = predictMEKF(ekf, omega_meas, dt)
    % 1. Propagar cuaternión nominal con giroscopio corregido
    omega_corrected = omega_meas - ekf.x(4:6);  % Compensar bias estimado
    ekf.q_nom = integrateQuaternion(ekf.q_nom, omega_corrected, dt);
    
    % 2. Propagar covarianza: P = Φ*P*Φ' + Q
    Phi = computeSTM(omega_corrected, dt);  % State Transition Matrix
    ekf.P = Phi * ekf.P * Phi' + ekf.Q * dt;
    
    % 3. Reset estado de error (MEKF property)
    ekf.x = zeros(6,1);
end
