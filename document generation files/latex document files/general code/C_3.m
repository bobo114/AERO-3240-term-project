function rotation_matrix_z = C_3(rotation_angle)
%C_3 roatates by given angle in degrees around z axis

rotation_matrix_z = [cosd(rotation_angle) sind(rotation_angle) 0; 
                    -sind(rotation_angle) cosd(rotation_angle) 0;
                        0                       0              1
    ];
end