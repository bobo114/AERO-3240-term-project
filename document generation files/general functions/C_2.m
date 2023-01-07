function rotation_matrix_y = C_2(rotation_angle)
%C_2 roatates by given angle in degrees around y axis

rotation_matrix_y = [cosd(rotation_angle) 0 -sind(rotation_angle); 
                    0                     1                 0   ;
                    sind(rotation_angle) 0 cosd(rotation_angle)
    ];
end