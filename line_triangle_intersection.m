function intersection_point = line_triangle_intersection(P, Q, triangle_vertices)
    eps = 0.000005;
    % P: Coordinates of point P [Px ; Py] 2x1
    % Q: Coordinates of point Q [Qx ; Qy] 2x1
    % triangle_vertices: Matrix containing the vertices of the triangle [xA xB xC;yA yB yC] 2x3;
    
    intersection_point=[];
    % Rotate and shift system of coordiantes
    L = norm(Q-P);
    transformed = [Q(1)-P(1),Q(2)-P(2);P(2)-Q(2),Q(1)-P(1)]*(triangle_vertices-[P(1);P(2)])/L;
    
    Purturbation=[1 2;1 3;2 3];
    %Find the intersecting side
    for i=1:3
        Cur_Line = [transformed(:,Purturbation(i,1)),transformed(:,Purturbation(i,2))];
        if Cur_Line(2,1) * Cur_Line(2,2) > 0 
            continue
        elseif Cur_Line(2,1) == 0 && Cur_Line(2,2) == 0
            X=0.5*(Cur_Line(1,1)+Cur_Line(1,2));
        elseif Cur_Line(2,1) == 0
            X=Cur_Line(1,1);
        elseif Cur_Line(2,2) == 0
            X=Cur_Line(1,2);
        else
            X=det(Cur_Line)/(Cur_Line(2,2)-Cur_Line(2,1));
        end
        if X >= -1*eps && X <=L+eps
            intersection_point=(Q-P)*X/L+P;
            break
        end
    end
end


