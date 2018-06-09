function R = EAngles321(angles)
phi = angles(1,1);
theta = angles(2,1);
psi = angles(3,1);
R11 = cos(theta)*cos(psi);
R12 = -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
R13 = sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
R21 = cos(theta)*sin(psi);
R22 = cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi);
R23 = -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);
R31 = -sin(theta);
R32 = sin(psi)*cos(theta);
R33 = cos(phi)*cos(theta);
R = [R11 R12 R13;R21 R22 R23;R31 R32 R33];
end