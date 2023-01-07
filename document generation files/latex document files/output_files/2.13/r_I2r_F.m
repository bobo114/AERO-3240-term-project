function r_F_vector = r_I2r_F(rotation_angle,r_I_vector)
%r_I2r_F Rotates r_I_vector about z axis by specified rotation_angle in degrees
r_F_vector = C_3(rotation_angle)*(r_I_vector');
end