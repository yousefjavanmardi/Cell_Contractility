% clc
% clear
% close all
% 
% E=30e6;
% Nu=0.25;
% 
% Nodes=[2,0;0,-1;0,1];

function K_Elem=Element_K(E,Nu,Nodes)
    %Cell Stiffness
    xi=Nodes(1,1);yi=Nodes(1,2);
    xj=Nodes(2,1);yj=Nodes(2,2);
    xm=Nodes(3,1);ym=Nodes(3,2);

    %ai=xj*ym-yj*xm;aj=yi*xm-xi*ym;am=xi*yj-yj*xi;
    bi=yj-ym;bj=ym-yi;bm=yi-yj;
    gi=xm-xj;gj=xi-xm;gm=xj-xi;
    Double_A=[xi,xj,xm]*[bi;bj;bm];

    B=[bi,0,bj,0,bm,0;0,gi,0,gj,0,gm;gi,bi,gj,bj,gm,bm]/Double_A;
    %%Plane Stress
    %D=[1,Nu,0;Nu,1,0;0,0,(1-Nu)/2]*E/(1-Nu^2);
    %Plane Strain
    D=[1-Nu,Nu,0;Nu,1-Nu,0;0,0,0.5-Nu]*E/((1+Nu)*(1-2*Nu));
    K_Elem=B'*D*B*Double_A/2;
end