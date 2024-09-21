function [New_X,New_Y]=Update_Location_Orientation(Connectivity,Nodes,Loc,Points_X,Points_Y,Rad,Strain,T_ang,Orientation,Red_Fac)
    %Orientation =1 for Actin fliaments
    %Orientation =-1 for Intermediate fliaments
    
    %Adjust location
    load("Cyto_Mesh.mat","model0");
    load('Parameters.mat','Disp_Factor_Maj','Rotat_Factor_Maj','Disp_Factor_Min','Rotat_Factor_Min');
    Range=max(Strain)-min(Strain);
    
    % Orientation =1 >> Major Filaments
    % Orientation =-1 >> Minor Filaments
    if Orientation==1
        Disp_Factor=Disp_Factor_Maj;
        Rotat_Factor=Rotat_Factor_Maj*(1-Red_Fac^4);
    else
        Disp_Factor=Disp_Factor_Min;
        Rotat_Factor=Rotat_Factor_Min*(1-Red_Fac^4);
    end


    N_Filament=size(Loc,1);
    Dir=zeros(N_Filament,2);
    for i_filament=1:N_Filament
        for i_node=1:2
            %find the associated element
            i_elem=Loc(i_filament,i_node);
            Neigh_Elem = i_elem;
            %find nodes of the element
            Nodes_of_Elem = Connectivity(:,i_elem);
            %find centre of the element
            Center=mean(Nodes(:,Nodes_of_Elem),2);
            %find neighnor element (including the element itself)
            Neigh_Elem = [Neigh_Elem,findElements(model0.Mesh ,"radius",Center,Rad)];
            %max major / min minor principal strain in the neighbourhood:E1_max
            if Orientation == 1
                [E1_Max,I]=min(Strain(Neigh_Elem));
                %major principal strains in the current element:E1_Cur
                E1_Cur=Strain(i_elem);
            else
                [E1_Max,I]=max(Strain(Neigh_Elem));
                %major principal strains in the current element:E1_Cur
                E1_Cur=Strain(i_elem);
            end
            %center of element with Max strain
            Center_Max_Elem = mean(Nodes(:,Connectivity(:,Neigh_Elem(I))),2);
            %direction of dislocation at each node
            Dir(i_filament,:)=Dir(i_filament,:)+0.5*((Center_Max_Elem-Center)*Disp_Factor*Orientation*(E1_Max-E1_Cur)/Range)';
            %Dir(i_filament,:)=Dir(i_filament,:)+0.5*((Center_Max_Elem-Center)*Disp_Factor*(E1_Max-E1_Cur)/Range)';
        end
    end

    %Adjust orientation
    TR=triangulation(Connectivity',Nodes(1,:)',Nodes(2,:)');
    XC = mean(Points_X,2);
    YC = mean(Points_Y,2);
    %Element containing centre
    ID = pointLocation(TR,XC,YC);
    A=find(isnan(ID));
    if ~isnan(A)
        ID(A)=Loc(A,1);
    end

    %Length
    Len=sqrt((Points_X(:,2)-Points_X(:,1)).^2 +(Points_Y(:,2)-Points_Y(:,1)).^2);
    %Angle
    Phi=pi()/2 - atan2(Points_Y(:,2)-Points_Y(:,1),Points_X(:,2)-Points_X(:,1));
    %Adjust Phi
    for i=1:size(Phi,1)
        if Phi(i)>=-pi()/2 && Phi(i)<0
            Phi(i)=Phi(i)+pi();
        elseif Phi(i)>pi() && Phi(i)<=3*pi()/2
            Phi(i)=Phi(i)-pi();
        end
    end


    if Orientation == -1
        T_ang2=zeros(size(T_ang,1),1);
        T_ang2(T_ang < pi /2)=T_ang(T_ang < pi /2) + pi/2;
        T_ang2(T_ang >= pi /2)=T_ang(T_ang >= pi /2) - pi/2;
        T_ang=T_ang2;
    end

    Range=max(T_ang(ID)-Phi);
    Phi=Phi+Rotat_Factor*(T_ang(ID)-Phi)/Range;

    XC=XC+Dir(:,1);
    YC=YC+Dir(:,2);
    New_X=[XC+0.5*Len.*sin(Phi),XC-0.5*Len.*sin(Phi)];
    New_Y=[YC+0.5*Len.*cos(Phi),YC-0.5*Len.*cos(Phi)];
end
