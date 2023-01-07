function rotation_matrix_x = C_1(rotation_angle)
%C_1 roatates by given angle in degrees around x axis

rotation_matrix_x = [1           0                       0; 
                    0 cosd(rotation_angle) sind(rotation_angle);
                    0 -sind(rotation_angle) cosd(rotation_angle)
    ];
end