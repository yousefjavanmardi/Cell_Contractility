function K=Boundary_membrane_stiffness(coordinates,E,A)
%function K_condensed=Boundary_membrane_stiffness(coordinates,E,A)
x1=coordinates(1,1);x2=coordinates(2,1);y1=coordinates(1,2);y2=coordinates(2,2);
t=atan2((y2-y1),(x2-x1));
L=((x2-x1)^2+(y2-y1)^2)^0.5;
C=cos(t);
S=sin(t);
K=[C^2    C*S    -1*C^2 -1*C*S;...
   C*S    S^2    -1*C*S -1*S^2;...
   -1*C^2 -1*C*S  C^2   C*S;...
   -1*C*S -1*S^2  C*S S^2]*A*E/L;
end