function q_prod = quatmultiply_custom(q1, q2)
    % Quaternion multiplication (Hamilton product)
    % Convention: [qw; qx; qy; qz]
    qw1 = q1(1); qv1 = q1(2:4);
    qw2 = q2(1); qv2 = q2(2:4);
    
    qw = qw1*qw2 - dot(qv1, qv2);
    qv = qw1*qv2 + qw2*qv1 + cross(qv1, qv2);
    
    q_prod = [qw; qv];
end