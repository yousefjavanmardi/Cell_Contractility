function Ki = Element_K_Fil(Initial_Ki, Len_Fiber,Angle_Fiber,A_Fiber,E_Fiber,Points,Nodes)
    %This function calculates the stiffness matrix of element due to
    %filaments
    %Points: Coordinates of Fiber ends, [x1,x2;y1,y2]
    %Nodes: Coordinates of the triangluar element:[x1,y1;x2,y2;x3,y3]

    k=E_Fiber*A_Fiber/Len_Fiber*[cos(Angle_Fiber)^2 sin(Angle_Fiber)*cos(Angle_Fiber) -cos(Angle_Fiber)^2 -sin(Angle_Fiber)*cos(Angle_Fiber);
         sin(Angle_Fiber)*cos(Angle_Fiber) sin(Angle_Fiber)^2 -sin(Angle_Fiber)*cos(Angle_Fiber) -sin(Angle_Fiber)^2;   
         -cos(Angle_Fiber)^2 -sin(Angle_Fiber)*cos(Angle_Fiber) cos(Angle_Fiber)^2 sin(Angle_Fiber)*cos(Angle_Fiber);
         -sin(Angle_Fiber)*cos(Angle_Fiber) -sin(Angle_Fiber)^2 sin(Angle_Fiber)*cos(Angle_Fiber) sin(Angle_Fiber)^2];

    %Cytoplasm Stiffness
    xi=Nodes(1,1);yi=Nodes(1,2);
    xj=Nodes(2,1);yj=Nodes(2,2);
    xm=Nodes(3,1);ym=Nodes(3,2);

    ai=xj*ym-yj*xm;aj=yi*xm-xi*ym;am=xi*yj-yj*xi;
    bi=yj-ym;bj=ym-yi;bm=yi-yj;
    gi=xm-xj;gj=xi-xm;gm=xj-xi;
    Double_A=[xi,xj,xm]*[bi;bj;bm];

    Points=[1 1;Points];
    N=[ai bi gi;aj bj gj;am bm gm]*Points;
    Shape_Functions=[N(1,1),0,N(2,1),0,N(3,1),0;0,N(1,1),0,N(2,1),0,N(3,1);N(1,2),0,N(2,2),0,N(3,2),0;0,N(1,2),0,N(2,2),0,N(3,2)]/Double_A;

    Ki=Initial_Ki+Shape_Functions'*k*Shape_Functions;
end




