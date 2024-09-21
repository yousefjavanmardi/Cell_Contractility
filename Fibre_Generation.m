function [Major_Points_X,Major_Points_Y,Minor_Points_X,Minor_Points_Y]=Fibre_Generation(Boundary_Coord,Nucleus_Coord)
    load('Parameters.mat',"L_Filament","Area_Ratio_Major","Area_Ratio_Minor");
    load('FE_Data.mat','model');
    Nodes=model.Mesh.Nodes;
    Connectivity=model.Mesh.Elements;
    N_Cyto_Elem=size(findElements(model.Mesh,"region","Face",1),2);
    
    % Length of filaments
    L_Min = 0.2*L_Filament;
    L_Max = 5*L_Filament;
    
    % Cell area = Area of Boundary - Area of Nucleus
    Cell_Area = polyarea(Boundary_Coord(:,1),Boundary_Coord(:,2))-polyarea(Nucleus_Coord(:,1),Nucleus_Coord(:,2));
    
    Length_Total_Major = Area_Ratio_Major * Cell_Area;
    Length_Total_Minor = Area_Ratio_Minor * Cell_Area;
    
    %Find domain's limits
    min_Boundary_Coord = min(Boundary_Coord);
    max_Boundary_Coord = max(Boundary_Coord);
    
    Major_Points_X=zeros(1,2);
    Major_Points_Y=zeros(1,2);
    Minor_Points_X=zeros(1,2);
    Minor_Points_Y=zeros(1,2);
    
    for i_Fillament = 1:2
        % Major >> i_Fillament = 1
        % Minor >> i_Fillament = 2
        Current_Length = 0;
        Counter = 0;
        if i_Fillament == 1
            Length_Total=Length_Total_Major;
        else
            Length_Total=Length_Total_Minor;
        end
        
        while Current_Length < Length_Total
            %Find a random line inside cell and outside nucleus
            Is_Line_In = 0;
            while ~Is_Line_In
                Is_Point_In = zeros(1,2);
            
                %Select a random point (as the first point)
                Rand_Point = min_Boundary_Coord + (max_Boundary_Coord - min_Boundary_Coord).*rand(1,2);
                %Check if it's inside CELL and outside NUCLEUS
                Is_Point_In (1) = ((inpolygon(Rand_Point(1,1),Rand_Point(1,2),Boundary_Coord(:,1),Boundary_Coord(:,2)) == true) && ...
                            (inpolygon(Rand_Point(1,1),Rand_Point(1,2),Nucleus_Coord(:,1),Nucleus_Coord(:,2)) == false));
                
                if Is_Point_In(1)
                    Length=L_Min + (L_Max - L_Min)*rand;
                    Angle = 2*pi()*rand;
                    End_Point = Rand_Point +Length *[cos(Angle),sin(Angle)];
                    Is_Point_In (2) = ((inpolygon(End_Point(1,1),End_Point(1,2),Boundary_Coord(:,1),Boundary_Coord(:,2)) == true) && ...
                                    (inpolygon(End_Point(1,1),End_Point(1,2),Nucleus_Coord(:,1),Nucleus_Coord(:,2)) == false));
                else
                    Is_Line_In = false;
                    continue
                end
                if Is_Point_In (2)
                    Is_Line_In =true;
                    Counter = Counter + 1;
                    if i_Fillament == 1
                        Major_Points_X (Counter,:)=[Rand_Point(1), End_Point(1)];
                        Major_Points_Y (Counter,:)=[Rand_Point(2), End_Point(2)];
                    else
                       Minor_Points_X (Counter,:)=[Rand_Point(1), End_Point(1)];
                       Minor_Points_Y (Counter,:)=[Rand_Point(2), End_Point(2)];
                    end
                end
            end
            if i_Fillament == 1
                Current_Length = Current_Length + sqrt((Major_Points_X (Counter,2)-Major_Points_X (Counter,1))^2 +...
                             (Major_Points_Y (Counter,2)-Major_Points_Y (Counter,1))^2);
            else
                Current_Length = Current_Length + sqrt((Minor_Points_X (Counter,2)-Minor_Points_X (Counter,1))^2 +...
                             (Minor_Points_Y (Counter,2)-Minor_Points_Y (Counter,1))^2);
            end
        end
    end
    %==Final control
    %Major
    [Major_Points_X,Major_Points_Y]=Correct_Filaments(Connectivity,Nodes,Major_Points_X,Major_Points_Y,N_Cyto_Elem);
    %MMinor
    [Minor_Points_X,Minor_Points_Y]=Correct_Filaments(Connectivity,Nodes,Minor_Points_X,Minor_Points_Y,N_Cyto_Elem);
end
