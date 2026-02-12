function H = computeMeasJacobian_Att(vPred)
%==========================================================================
% computeMeasJacobian: Computes measurement Jacobian for MEKF vector-based
%                      measurements (magnetometer, star tracker, sun sensor).
%                      Linearizes vector measurement model around current
%                      attitude estimate.
%
% INPUTS:
%    vPred - Predicted measurement vector in body frame (3x1)
%            Examples: magnetic field [nT], star direction [unit], sun direction [unit]
%
% OUTPUTS:
%    H     - Measurement Jacobian (3x6)
%            H = [∂v_body/∂deltaTheta, ∂v_body/∂biasGyro]
%
% METHOD:
%    Measurement model: vBody = C(q) * vECI
%    
%    For small attitude error δθ:
%        C(δq ⊗ qNom) ≈ C(qNom) * (I - [δθ]×)
%    
%    Therefore: vBody ≈ vPred - [vPred]× * δθ
%    
%    Jacobian: ∂v_body/∂δθ    = -[vPred]×
%              ∂v_body/∂bGyro = 0 (bias doesn't affect vector measurements)
%
% REFERENCE:
%    Markley & Crassidis, "Fundamentals of Spacecraft Attitude 
%    Determination and Control", Chapter 7
%==========================================================================

    % --- 1. Skew-Symmetric Matrix of Predicted Vector ---
    vSkew = [        0, -vPred(3),  vPred(2);
              vPred(3),         0, -vPred(1);
             -vPred(2),  vPred(1),        0];
    
    % --- 2. Assemble Jacobian ---
    H = [-vSkew, zeros(3,3)];

end
