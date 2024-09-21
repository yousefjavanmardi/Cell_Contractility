function [New_Eps_dot,New_X,New_Y]=Update_Filaments(Connectivity,Nodes,U,Loc,Points_X,Points_Y,Time_Step)
    N_Major=size(Loc,1);
    New_Length=zeros(N_Major,1);
    New_Eps_dot=zeros(N_Major,1);
    New_X=Points_X;
    New_Y=Points_Y;
    for i_maj=1:N_Major
        Coord_Old=zeros(2);
        Coord_New=zeros(2);
        for i_node=1:2
            %find the coordinates of the filaments end
            node1=[Points_X(i_maj,i_node);Points_Y(i_maj,i_node)];
            %find the associated element
            i_elem=Loc(i_maj,i_node);
            Nodes_of_Elem = Connectivity(:,i_elem);
            X_Elem = Nodes(1,Nodes_of_Elem);
            Y_Elem = Nodes(2,Nodes_of_Elem);
            %find shape functions
            mat0=[[1 1 1];X_Elem;Y_Elem];
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
            %find nodal displacements of the element
            U_Elem = [U(2*Nodes_of_Elem-1),U(2*Nodes_of_Elem)];
            Coord_Old(i_node,:)=node1';
            Coord_New(i_node,:)=node1'+[N1,N2,N3]*U_Elem;
            New_X(i_maj,i_node)=Coord_New(i_node,1);
            New_Y(i_maj,i_node)=Coord_New(i_node,2);
        end
        %update L and Eps dot(rate of strain)
        L_Old=norm(Coord_Old(2,:)-Coord_Old(1,:));
        L_New=norm(Coord_New(2,:)-Coord_New(1,:));
        New_Length(i_maj)=L_New;
        New_Eps=(L_New-L_Old)/L_Old;
        New_Eps_dot(i_maj)=New_Eps/Time_Step;
    end
end