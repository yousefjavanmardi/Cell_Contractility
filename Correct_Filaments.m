function [New_X,New_Y]=Correct_Filaments(Connectivity,Nodes,Points_X,Points_Y,N_Cyto_Elem)

    %find updated boundary (membrane)
    TR = triangulation(Connectivity',Nodes');
    [~,Boundary_Coord] = freeBoundary(TR);

    %find updated boundary (nucleus)
    Connectivity_Nuc = Connectivity(:,N_Cyto_Elem+1:end);
    TR = triangulation(Connectivity_Nuc',Nodes');
    [~,Nucleus_Coord] = freeBoundary(TR);
    
    %find nodes in cytoplasm
    Connectivity_Cyto=Connectivity(:,1:N_Cyto_Elem);
    Nodes_Cyto=Nodes(:,unique(Connectivity_Cyto));
    for i=1:size(Connectivity_Cyto,2)
        Center_Elem_Cyto(:,i)=mean(Nodes_Cyto(:,Connectivity_Cyto(:,i)),2);
    end
    eps=0.0001;

    N_Fil=size(Points_X,1);
    Is_IN_Cyto=zeros(N_Fil,2);
    New_X=Points_X;
    New_Y=Points_Y;
    for i_node=1:2
        %find coordinates of ends of the filament 2x1
        X_Coord_End =Points_X(:,i_node);
        Y_Coord_End =Points_Y(:,i_node);
        %Check if the nodes are inside the cell and outside the nucleus
        %This is based on polyshape not mesh
        Is_In_Cell = inpolygon(X_Coord_End,Y_Coord_End,Boundary_Coord(:,1),Boundary_Coord(:,2));
        Is_Out_Nucleus = ~inpolygon(X_Coord_End,Y_Coord_End,Nucleus_Coord(:,1),Nucleus_Coord(:,2));
        Is_IN_Cyto(:,i_node) = Is_In_Cell & Is_Out_Nucleus;
    end

    %filaments with only one end out
    Is_Fil_Out=xor(~Is_IN_Cyto(:,1),~Is_IN_Cyto(:,2));
    ID=find(Is_Fil_Out == true);
    S_Out=size(ID,1);
    for i_Out=1:S_Out
        Cur_Fil=ID(i_Out);
        %node1 is outside, node2 is inside 
        node1=[Points_X(Cur_Fil,1);Points_Y(Cur_Fil,1)];
        node2=[Points_X(Cur_Fil,2);Points_Y(Cur_Fil,2)];
        The_node=1;
        if Is_IN_Cyto(Cur_Fil,1)
            node1=[Points_X(Cur_Fil,2);Points_Y(Cur_Fil,2)];
            node2=[Points_X(Cur_Fil,1);Points_Y(Cur_Fil,1)];
            The_node=2;
        end
        k=1;
        dT=0.05;
        while true
            Dir=node1-node2;
            L=norm(Dir);
            T=atan2(Dir(2),Dir(1));
            T=T+k*dT;
            k=k+1;
            New_Coord=node2+[L*cos(T);L*sin(T)];
            I1= inpolygon(New_Coord(1),New_Coord(2),Boundary_Coord(:,1),Boundary_Coord(:,2));
            I2 = ~inpolygon(New_Coord(1),New_Coord(2),Nucleus_Coord(:,1),Nucleus_Coord(:,2));
            I=I1 & I2;
            if I
                break
            end
        end

        New_X(Cur_Fil,The_node)=New_Coord(1);
        New_Y(Cur_Fil,The_node)=New_Coord(2);
    end
    
    %filaments with both ends out
    Is_Fil_Out= (~Is_IN_Cyto(:,1)) & (~Is_IN_Cyto(:,2));
    ID=find(Is_Fil_Out == true);
    S_Out=size(ID,1);
    for i_Out=1:S_Out
        Cur_Fil=ID(i_Out);
         
        node1=[Points_X(Cur_Fil,1);Points_Y(Cur_Fil,1)];
        node2=[Points_X(Cur_Fil,2);Points_Y(Cur_Fil,2)];
        %we first bring node2 inside
        Near_Node2 = dsearchn(Center_Elem_Cyto',node2');
        node2=[Center_Elem_Cyto(1,Near_Node2);Center_Elem_Cyto(2,Near_Node2)];
        New_X(Cur_Fil,2)=node2(1);
        New_Y(Cur_Fil,2)=node2(2);
        %Now rotate node1
        k=1;
        dT=0.05;
        while true
            Dir=node1-node2;
            L=norm(Dir);
            T=atan2(Dir(2),Dir(1));
            T=T+k*dT;
            k=k+1;
            New_Coord=node2+[L*cos(T);L*sin(T)];
            I1= inpolygon(New_Coord(1),New_Coord(2),Boundary_Coord(:,1),Boundary_Coord(:,2));
            I2 = ~inpolygon(New_Coord(1),New_Coord(2),Nucleus_Coord(:,1),Nucleus_Coord(:,2));
            I=I1 & I2;
            if I
                break
            end
        end
        New_X(Cur_Fil,1)=New_Coord(1);
        New_Y(Cur_Fil,1)=New_Coord(2);
    end
end


