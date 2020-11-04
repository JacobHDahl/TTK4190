function course_d = guidance(pos,path_info)
    
    x = pos(1); y = pos(2);
    
    [x_p, y_p, y_e] = crosstrack(x_t, y_t, x_ref, y_ref, x, y);
    
    delta_idx = 0;
    delta = path_info(delta_idx);
    Kp = 1/delta;
    cord_angle = 0; %Pi_p
    
    course_d = cord_angle - atan(Kp*y_e);
    


end

