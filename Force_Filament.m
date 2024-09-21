function Nodal_Force = Force_Filament(node1,node2,Coord_Elem,Force_ax)
    mat0=[[1 1 1];Coord_Elem];
    mat1=mat0;
    mat1(2:3,1)=node1;
    mat2=mat0;
    mat2(2:3,2)=node1;
    mat3=mat0;
    mat3(2:3,3)=node1;
    Double_A=det(mat0);
    N1=det(mat1)/Double_A;
    N2=det(mat2)/Double_A;
    N3=det(mat3)/Double_A;
    Lambda=(node2-node1)/norm(node2-node1);
    Nodal_Force=[N1*Lambda,N2*Lambda,N3*Lambda]*Force_ax;
end