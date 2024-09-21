function model=Mesh_Generation (P1,P2,AMS)
    %P1 => Cytoplasmic Boundary
    %P2 => Nuclear Boudary
    
    P=[P2;P1];
    %===========Mesh================================ 
    DT1 = triangulation(polyshape(P));
    model = fegeometry(DT1.Points,DT1.ConnectivityList);
    %model only for cytoplasm
    model0 = generateMesh(model,"Hmax",AMS,"Hmin",AMS/10,"GeometricOrder","linear");
    save("Cyto_Mesh.mat","model0");
    %
    Nucl_Edge = nearestEdge(model,P2);
    model=addFace(model,Nucl_Edge);
    model = generateMesh(model,"Hmax",AMS,"Hmin",AMS/10,"GeometricOrder","linear");
    
    %===========Find Nodes on Membrane and Nucleus border
    %nucleus Nodes   
    Edge_Nodes={};
    P2=unique(P2,'rows','stable');
    Vect_Nuc_Edg=[];
    Counter=0;
    for i=1:size(P2,1)
       Nuc_Edge = nearestEdge(model,P2(i,:));
       if ~ismember(Nuc_Edge,Vect_Nuc_Edg)
            Vect_Nuc_Edg=[Vect_Nuc_Edg,Nuc_Edge];
            Counter=Counter+1;
            Edge_Nodes{Counter}=findNodes(model.Mesh,"region","Edge",Nuc_Edge);
       end
    end
    
    %Correct the order of Nodes
    Nucleus_Nodes=Edge_Nodes{1};
    Nodes_on_Nuc=zeros(1,size(Edge_Nodes,2));
    for i=1:size(Edge_Nodes,2)-1
        N2=Edge_Nodes{i}(end);
        Coord_end=model.Mesh.Nodes(:,N2);
        N1=Edge_Nodes{i+1}(1);
        N2=Edge_Nodes{i+1}(end);
        Coord_ini2=model.Mesh.Nodes(:,N1);
        Coord_end2=model.Mesh.Nodes(:,N2);
        norm1=norm(Coord_ini2-Coord_end);
        norm2=norm(Coord_end2-Coord_end);
        %Correct the direction of edge 1
        if i==1
            N1=Edge_Nodes{i}(1);
            Coord_ini=model.Mesh.Nodes(:,N1);
            norm11=norm(Coord_ini2-Coord_ini);
            norm21=norm(Coord_end2-Coord_ini);
            if min(norm11,norm21)<min(norm1,norm2)
                Nucleus_Nodes=fliplr(Edge_Nodes{1});
                norm2=norm21;
                norm1=norm11;
            end
            Nodes_on_Nuc(1)=size(unique(Edge_Nodes{1},'stable'),2);
        end
        if norm2 < norm1
            Edge_Nodes{i+1}=fliplr(Edge_Nodes{i+1});
        end
        Nucleus_Nodes=[Nucleus_Nodes,Edge_Nodes{i+1}];
        Nodes_on_Nuc(i+1)=size(unique(Edge_Nodes{i+1},'stable'),2);
    end 
    %Nucleus_Nodes=unique(Nucleus_Nodes,'stable');

    %Membrane Nodes
    Edge_Nodes={};
    P1=unique(P1,'rows','stable');
    Vect_Mem_Edg=[];
    Counter=0;
    for i=1:size(P1,1)
       Mem_Edge = nearestEdge(model,P1(i,:));
       if ~ismember(Mem_Edge,Vect_Mem_Edg)
            Vect_Mem_Edg=[Vect_Mem_Edg,Mem_Edge];
            Counter=Counter+1;
            Edge_Nodes{Counter}=findNodes(model.Mesh,"region","Edge",Mem_Edge);
       end
    end
    %Correct the order of Nodes
    Membrane_Nodes=Edge_Nodes{1};
    Nodes_on_Edge=zeros(1,size(Edge_Nodes,2));
    for i=1:size(Edge_Nodes,2)-1
        N2=Edge_Nodes{i}(end);
        Coord_end=model.Mesh.Nodes(:,N2);
        N1=Edge_Nodes{i+1}(1);
        N2=Edge_Nodes{i+1}(end);
        Coord_ini2=model.Mesh.Nodes(:,N1);
        Coord_end2=model.Mesh.Nodes(:,N2);
        norm1=norm(Coord_ini2-Coord_end);
        norm2=norm(Coord_end2-Coord_end);
        %Correct the direction of edge 1
        if i==1
            N1=Edge_Nodes{i}(1);
            Coord_ini=model.Mesh.Nodes(:,N1);
            norm11=norm(Coord_ini2-Coord_ini);
            norm21=norm(Coord_end2-Coord_ini);
            if min(norm11,norm21)<min(norm1,norm2)
                Membrane_Nodes=fliplr(Edge_Nodes{1});
                norm2=norm21;
                norm1=norm11;
            end
            Nodes_on_Edge(1)=size(unique(Edge_Nodes{1},'stable'),2);
        end
        if norm2 < norm1
            Edge_Nodes{i+1}=fliplr(Edge_Nodes{i+1});
        end
        Membrane_Nodes=[Membrane_Nodes,Edge_Nodes{i+1}];
        Nodes_on_Edge(i+1)=size(unique(Edge_Nodes{i+1},'stable'),2);
    end 
    %Membrane_Nodes=unique(Membrane_Nodes,'stable');

    save('Boundary_Nodes.mat','Membrane_Nodes',"Nucleus_Nodes","Nodes_on_Edge","Nodes_on_Nuc");   
end