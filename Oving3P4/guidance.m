function [course_d,y_e] = guidance(pos,path_info)
    
    x = pos(1); y = pos(2);
    x1 = path_info(1); y1 = path_info(2);
    x2 = path_info(3); y2 = path_info(4);
    
    y_e = crosstrackWpt(x2, y2, x1, y1, x, y);
    delta = path_info(5);
    Kp = 1/delta;
    Pi_p = path_info(6);
    
    course_first = Pi_p - atan(Kp*y_e);
    course_d = wrapTo2Pi(course_first);


end

