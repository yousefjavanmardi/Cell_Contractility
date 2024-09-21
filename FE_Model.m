function FE_Model()

    load('FE_Data.mat','model','Boundary_Coord','Support_Coord');
    load('Fillement_Data.mat','Major_Points_X','Major_Points_Y','Minor_Points_X','Minor_Points_Y');
    load('Parameters.mat',"Search_Rad_Maj","E_Cell","Nu_Cell","E_Nucl","Nu_Nucl","E_Fiber",...
        "A_Fiber","E_Membrane","A_Membrane","E_Nuc_Mem","A_Nuc_Mem","K_Support",...
        "Sigma_Max","Theta","kf","kb","Eps0_dot","kv","c_signal","Time_Step","N_Step","Search_Rad_Min");
    
    load('Boundary_Nodes.mat','Membrane_Nodes',"Nucleus_Nodes","Nodes_on_Edge","Nodes_on_Nuc");
    
    Connectivity=model.Mesh.Elements;
    Nodes=model.Mesh.Nodes;
    N_Cyto_Elem=size(findElements(model.Mesh,"region","Face",1),2);
    
    E_Cell=1000*E_Cell; % kPa=> Pa
    E_Nucl=1000*E_Nucl; % kPa=> Pa
    E_Fiber=1000*E_Fiber; % kPa=> Pa
    E_Membrane=1000*E_Membrane; % kPa=> Pa
    E_Nuc_Mem=1000*E_Nuc_Mem; % kPa=> Pa
    Sigma_Max=1000*Sigma_Max; % kPa=> Pa
    
    N_Elem=size(Connectivity,2);
    N_Nodes=size(Nodes,2);
    
    %Find the location of filaments' ends
    %Major_Loc >> size: N_Major*2 ; contains: Elem No. begin, Elem No. end
    %Minor_Loc >> size: N_Major*2 ; contains: Elem No. begin, Elem No. end
    N_Major=size(Major_Points_X,1);
    N_Minor=size(Minor_Points_X,1);
    Major_Loc=zeros(N_Major,2);
    Minor_Loc=zeros(N_Minor,2);
    
    %Initializing Filament model parameters
    Eps_dot=zeros(N_Major,1);
    Eta=ones(N_Major,1);
    Sigma_by_Sigma0=zeros(N_Major,1);
    
    %Plot Initial State
    figure('Name',num2str(0));
    triplot(Connectivity',Nodes(1,:)',Nodes(2,:)','color',[0.7 0.7 0.7]);
    hold on
    axis square
    axis off
    %set(gcf,'color','w','Units','normalized','OuterPosition',[0 0 1 1]);
    set(gcf,'color','w','Units','normalized');
    x_min = floor(min(Boundary_Coord(:,1))/10)*10;
    x_max = ceil(max(Boundary_Coord(:,1))/10)*10;
    y_min = floor(min(Boundary_Coord(:,2))/10)*10;
    y_max = ceil(max(Boundary_Coord(:,2))/10)*10;
    xlim([x_min x_max]) 
    ylim([y_min y_max])

    %plot supports
    Supports = dsearchn(Nodes',Support_Coord);
    Supports=unique(Supports,'stable');
    plot(Nodes(1,Supports),Nodes(2,Supports),'o','MarkerEdgeColor','k','MarkerFaceColor','g')
    
    %plot membrane
    Mem_Elem=[Membrane_Nodes;circshift(Membrane_Nodes,-1)];
    plot(Nodes(1,Mem_Elem)',Nodes(2,Mem_Elem)','-k')
    
    %plot nucleus
    Nuc_Elem=[Nucleus_Nodes;circshift(Nucleus_Nodes,-1)];
    plot(Nodes(1,Nuc_Elem)',Nodes(2,Nuc_Elem)','-k','LineWidth',0.75)
    print(gcf,'-vector','-dsvg','Geometry.svg')
    
    
    
    %plot filaments
    plot (Major_Points_X',Major_Points_Y','b-','LineWidth',0.75)
    plot (Minor_Points_X',Minor_Points_Y','r-')
    
    print(gcf,'-vector','-dsvg','Initial.svg')
    
    %title('t= 0:0:0 hr:min:sec')
    title('Iteration = 0')
    U3_M(1)=getframe(gcf);
    drawnow
    
    Total_Strain_MEAN=zeros(N_Cyto_Elem,1);

    Sigma=zeros(N_Major,1);

    Supports = dsearchn(Nodes',Support_Coord);
    Supports=unique(Supports,'stable');

    Total_Disp_Maj=zeros(N_Step,1);
    Total_Disp_Min=zeros(N_Step,1);
    Total_Rot_Maj=zeros(N_Step,1);
    Total_Rot_Min=zeros(N_Step,1);
    
    WB=waitbar(0,'Starting');

    for i_Time=1:N_Step
        TR=triangulation(Connectivity',Nodes');
            
        %Find elements associated with Major filaments
        for i_Major=1:N_Major
            Elem_ID = pointLocation(TR,Major_Points_X(i_Major,:)',Major_Points_Y(i_Major,:)');
            Major_Loc(i_Major,:)=Elem_ID';
        end
        
        %Find elements associated with Minor filaments
        for i_Minor=1:N_Minor
            Elem_ID = pointLocation(TR,Minor_Points_X(i_Minor,:)',Minor_Points_Y(i_Minor,:)');
            Minor_Loc(i_Minor,:)=Elem_ID';
        end
        
        %=======================Global Stiffness Matrix============================
        K=zeros(2*N_Nodes);   
        for i_Elem =1:N_Elem
            %===========Stiffness Due To Cell/Nucelus Stiffness
            Elem_Nodes = Connectivity(:,i_Elem);
            Coord_Nodes=Nodes(:,Elem_Nodes)';
            %Cytoplasm or Nucleus
            if i_Elem <= N_Cyto_Elem
                E= E_Cell;
                Nu=Nu_Cell;
            else
                E= E_Nucl;
                Nu=Nu_Nucl;
            end
        
            %Element Stiffness Matrix
            Ki=Element_K(E,Nu,Coord_Nodes); % t=1um => unit of Ki: nN/m
        
            %================Stiffness due to Filaments
            %==Major
            %Is element associated with filaments?
            Elem_Filam = find (Major_Loc == i_Elem);
            %Find all filaments
            if ~isnan(Elem_Filam)
                Elem_Filam = Elem_Filam - (Elem_Filam > N_Major)*N_Major;
                Elem_Filam=unique(Elem_Filam,'stable');
            end
            flag =0;
            for i_Film=1:size(Elem_Filam,1)
                curr_Film=Elem_Filam(i_Film);
                x1=Major_Points_X(curr_Film,1);
                x2=Major_Points_X(curr_Film,2);
                y1=Major_Points_Y(curr_Film,1);
                y2=Major_Points_Y(curr_Film,2);
                Len_Fiber=sqrt((x2-x1)^2+(y2-y1)^2);
                Angle_Fiber=atan2(y2-y1,x2-x1);
                if Major_Loc(curr_Film,1)==Major_Loc(curr_Film,2)
                    Points=[x1,x2;y1,y2];
                else
                    if(inpolygon(x2,y2,Coord_Nodes(:,1),Coord_Nodes(:,2)))
                        x3=x1;
                        y3=y1;
                        x1=x2;
                        y1=y2;
                        x2=x3;
                        y2=y3;
                    end
                    intersection_point = line_triangle_intersection([x1;y1], [x2;y2], Coord_Nodes');
                    if isempty(intersection_point)
                        flag=1;
                        break
                    end
                    Points=[x1,intersection_point(1);y1,intersection_point(2)];
    
                end
                Ki = Element_K_Fil(Ki, Len_Fiber,Angle_Fiber,A_Fiber,E_Fiber,Points,Coord_Nodes);
            end
        
            %==Minor
            %Is element associated with filaments?
            Elem_Filam = find (Minor_Loc == i_Elem);
            %Find all filaments
            if ~isnan(Elem_Filam)
                Elem_Filam = Elem_Filam - (Elem_Filam > N_Minor)*N_Minor;
                Elem_Filam=unique(Elem_Filam,'stable');
            end
            for i_Film=1:size(Elem_Filam,1)
                curr_Film=Elem_Filam(i_Film);
                x1=Minor_Points_X(curr_Film,1);
                x2=Minor_Points_X(curr_Film,2);
                y1=Minor_Points_Y(curr_Film,1);
                y2=Minor_Points_Y(curr_Film,2);
                Len_Fiber=sqrt((x2-x1)^2+(y2-y1)^2);
                Angle_Fiber=atan2(y2-y1,x2-x1);
                if Minor_Loc(curr_Film,1)==Minor_Loc(curr_Film,2)
                    Points=[x1,x2;y1,y2];
                else
                    if(inpolygon(x2,y2,Coord_Nodes(:,1),Coord_Nodes(:,2)))
                        x3=x1;
                        y3=y1;
                        x1=x2;
                        y1=y2;
                        x2=x3;
                        y2=y3;
                    end
                    intersection_point = line_triangle_intersection([x1;y1], [x2;y2], Coord_Nodes');
                    if isempty(intersection_point)
                        flag=1;
                        break
                    end
                    Points=[x1,intersection_point(1);y1,intersection_point(2)];
    
                end
                Ki = Element_K_Fil(Ki, Len_Fiber,Angle_Fiber,A_Fiber,E_Fiber,Points,Coord_Nodes);
            end
            if flag ==1
                break
            end
            
            %=================Assembly=============================================
            R=2*Elem_Nodes-1;
            for i_Row =1:3
                for i_Col =1:3
                    K(R(i_Row):R(i_Row)+1,R(i_Col):R(i_Col)+1) = K(R(i_Row):R(i_Row)+1,R(i_Col):R(i_Col)+1) + Ki(2*i_Row-1:2*i_Row,2*i_Col-1:2*i_Col);
                end
            end
        end
        if flag ==1
            break
        end
        %==================End of Global Stiffness Matrix==========================
        
        %================Displacement Boundary Conditions==========================
        
        %==1-Plasma Membrane Stiffness
        N_Edge=size(Nodes_on_Edge,2);
        for i=1:N_Edge
            %Find nodes on each edge
            if i==1
                begin_edge=1;
            else
                begin_edge=begin_edge+Nodes_on_Edge(i-1);
            end
            end_edge=begin_edge+Nodes_on_Edge(i)-1;
            Edge_Nodes = Membrane_Nodes(begin_edge:end_edge);
            NNode_Edge = size(Edge_Nodes,2);
            nod=zeros(1,2);
            for i_Case=1:NNode_Edge-1
                % for the nodes on each edge
                for j=1:NNode_Edge-i_Case
                    %find each element
                    nod(1)=Edge_Nodes(j);
                    nod(2)=Edge_Nodes(j+i_Case);
                    Coodinate=[Nodes(:,nod(1))';Nodes(:,nod(2))'];
                    %find K of each truss element
                    Elem_Stiff=Boundary_membrane_stiffness(Coodinate,E_Membrane,A_Membrane);
                    %asseble global K
                    for i_row=1:2
                        for j_row=1:2
                            K(2*nod(i_row)-1:2*nod(i_row),2*nod(j_row)-1:2*nod(j_row))=...
                            K(2*nod(i_row)-1:2*nod(i_row),2*nod(j_row)-1:2*nod(j_row)) + ...
                            Elem_Stiff(2*i_row-1:2*i_row,2*j_row-1:2*j_row);
                        end
                    end
                end
            end
        end
        
        %==2-Nucleus Stiffness
        N_Edge=size(Nodes_on_Nuc,2);
        for i=1:N_Edge
            %Find nodes on each edge
            if i==1
                begin_edge=1;
            else
                begin_edge=begin_edge+Nodes_on_Nuc(i-1);
            end
            end_edge=begin_edge+Nodes_on_Nuc(i)-1;
            Edge_Nodes = Nucleus_Nodes(begin_edge:end_edge);
            NNode_Edge = size(Edge_Nodes,2);
            nod=zeros(1,2);
            % for the nodes on each edge
            for j=1:NNode_Edge-1
                %find each element
                nod(1)=Edge_Nodes(j);
                nod(2)=Edge_Nodes(j+1);
                Coodinate=[Nodes(:,nod(1))';Nodes(:,nod(2))'];
                %find K of each truss element
                Elem_Stiff=Boundary_membrane_stiffness(Coodinate,E_Nuc_Mem,A_Nuc_Mem);
                %asseble global K
                for i_row=1:2
                    for j_row=1:2
                        K(2*nod(i_row)-1:2*nod(i_row),2*nod(j_row)-1:2*nod(j_row))=...
                        K(2*nod(i_row)-1:2*nod(i_row),2*nod(j_row)-1:2*nod(j_row)) + ...
                        Elem_Stiff(2*i_row-1:2*i_row,2*j_row-1:2*j_row);
                    end
                end
            end
        end

        %==3-Spring Supports        
        for i_Supp=1:size(Supports,1)
            j_Supp=2*Supports(i_Supp,1)-1;
            K(j_Supp,j_Supp)=K(j_Supp,j_Supp)+K_Support;
            K(j_Supp+1,j_Supp+1)=K(j_Supp+1,j_Supp+1)+K_Support;
        end
        %==========================================================================
        
        %%========Force boundary conditions========================================
        %===Forces in Major filaments
        Total_Force=zeros(2*N_Nodes,1);
        
        
        Time=(i_Time-1)*Time_Step;

        % if i_Time == 1
        %     Sigma=(Sigma_Max/N_Step)*ones(N_Major,1);
        % end
        %Find C
        C=c_signal*exp(-1*Time/Theta);
        %Find Sigma / Sigma0
        Eps_dot=-Eps_dot; %To convert extensive strain to compressive strain
        for i_Major = 1:N_Major
            if Eps_dot(i_Major)/Eps0_dot <-Eta(i_Major)/kv
                Sigma_by_Sigma0(i_Major)=0;
            elseif Eps_dot(i_Major)/Eps0_dot <= 0
                if Eta(i_Major)==0
                    Eta(i_Major)=0.00001;
                end
                Sigma_by_Sigma0(i_Major)=1+(kv/Eta(i_Major))*Eps_dot(i_Major)/Eps0_dot;
            else
                Sigma_by_Sigma0(i_Major)=1;
            end
        end
        %Find Eta dot
        Eta_dot=(1-Eta)*C*kf/Theta-(1-Sigma_by_Sigma0).*Eta*kb/Theta;
        %Update Eta
        Eta=Eta + Eta_dot*Time_Step;
        Eta(Eta<0)=0;
        Eta(Eta>1)=1;
        Sigma=Eta.*Sigma_by_Sigma0*Sigma_Max/N_Step;
        %Find d_Sigma
        %Sigma=Sigma1-Sigma;
        Sigma(Sigma>Sigma_Max)=Sigma_Max;
        Sigma(Sigma<-Sigma_Max)=-Sigma_Max;
        %Find Force in the filament
        Force_Fil=Sigma*A_Fiber;
    
        for cur_Major=1:N_Major
            %Apply force to nodes of first element
            Elem1=Major_Loc(cur_Major,1);
            Node_of_Elem =Connectivity(:,Elem1);
            Coord_Elem = Nodes(:,Node_of_Elem);
            node1=[Major_Points_X(cur_Major,1);Major_Points_Y(cur_Major,1)];
            node2=[Major_Points_X(cur_Major,2);Major_Points_Y(cur_Major,2)];
            Nodal_Force = Force_Filament(node1,node2,Coord_Elem,Force_Fil(cur_Major));
            for i_NE=1:3
                Total_Force(2*Node_of_Elem(i_NE)-1,1)=Total_Force(2*Node_of_Elem(i_NE)-1,1)+Nodal_Force(1,i_NE);
                Total_Force(2*Node_of_Elem(i_NE),1)=Total_Force(2*Node_of_Elem(i_NE),1)+Nodal_Force(2,i_NE);
            end
            %Apply force to nodes of second element
            Elem2=Major_Loc(cur_Major,2);
            Node_of_Elem =Connectivity(:,Elem2);
            Coord_Elem = Nodes(:,Node_of_Elem);
            node1=[Major_Points_X(cur_Major,2);Major_Points_Y(cur_Major,2)];
            node2=[Major_Points_X(cur_Major,1);Major_Points_Y(cur_Major,1)];
            Nodal_Force = Force_Filament(node1,node2,Coord_Elem,Force_Fil(cur_Major));
            for i_NE=1:3
                Total_Force(2*Node_of_Elem(i_NE)-1,1)=Total_Force(2*Node_of_Elem(i_NE)-1,1)+Nodal_Force(1,i_NE);
                Total_Force(2*Node_of_Elem(i_NE),1)=Total_Force(2*Node_of_Elem(i_NE),1)+Nodal_Force(2,i_NE);
            end
        end
        %==========================================================================
        
        %======================Solving Equations===================================
        U=K\Total_Force;
        %==========================================================================
        
        %==========================Update Parameters===============================
        
        %===Find Mean Principal Strain in Elements
        [Strain_MEAN,T_ang]=Strain_Finder(Connectivity(:,1:N_Cyto_Elem),Nodes,U);
        Total_Strain_MEAN=Total_Strain_MEAN+Strain_MEAN;
    
        %===Update Filaments Coordinates, Update Eps dot
        [Eps_dot,Major_Points_X,Major_Points_Y]=Update_Filaments(Connectivity(:,1:N_Cyto_Elem),Nodes,U,...
            Major_Loc,Major_Points_X,Major_Points_Y,Time_Step);
        
        [~,Minor_Points_X,Minor_Points_Y]=Update_Filaments(Connectivity(:,1:N_Cyto_Elem),Nodes,U,...
            Minor_Loc,Minor_Points_X,Minor_Points_Y,Time_Step);
    
        %Update Mesh Coordinate
        Nodes(1,:)=Nodes(1,:)+U(1:2:end,1)';
        Nodes(2,:)=Nodes(2,:)+U(2:2:end,1)';
    
        %===Update Location and Orientation of filaments
        %Major
        Ini_X=Major_Points_X;
        Ini_Y=Major_Points_Y;
        %Rad=Search_Rad*(1-((i_Time-1)/(N_Step-1))^2)+1;
        n_Fac=1.5;
        k_Fac=log(Search_Rad_Maj)*n_Fac^2/N_Step^2;
        Rad=Search_Rad_Maj*exp(-1*k_Fac*i_Time^2);
        [Major_Points_X,Major_Points_Y]=Update_Location_Orientation(Connectivity(:,1:N_Cyto_Elem),Nodes,Major_Loc,Major_Points_X,...
            Major_Points_Y,Rad,Total_Strain_MEAN,T_ang,1,i_Time/N_Step);
        
        if i_Time > 1
            Total_Disp_Maj(i_Time)=Total_Disp_Maj(i_Time-1)+sum(sqrt(mean(Ini_X-Major_Points_X,2).^2+mean(Ini_Y-Major_Points_Y,2).^2));
            Vec_Before=[Ini_X(:,2)-Ini_X(:,1),Ini_Y(:,2)-Ini_Y(:,1)];
            Vec_After=[Major_Points_X(:,2)-Major_Points_X(:,1),Major_Points_Y(:,2)-Major_Points_Y(:,1)];
            Total_Rot_Maj(i_Time)=Total_Rot_Maj(i_Time-1)+sum(acosd(dot(Vec_Before',Vec_After')./(vecnorm(Vec_Before').*vecnorm(Vec_After'))));
        end
    
        %Minor
        Ini_X=Minor_Points_X;
        Ini_Y=Minor_Points_Y;
        n_Fac=1.5;
        k_Fac=log(Search_Rad_Min)*n_Fac^2/N_Step^2;
        Rad=Search_Rad_Min*exp(-1*k_Fac*i_Time^2);
        [Minor_Points_X,Minor_Points_Y]=Update_Location_Orientation(Connectivity(:,1:N_Cyto_Elem),Nodes,Minor_Loc,Minor_Points_X,...
            Minor_Points_Y,Rad,Total_Strain_MEAN,T_ang,-1,i_Time/N_Step);
        if i_Time >1
            Total_Disp_Min(i_Time)=Total_Disp_Min(i_Time-1)+sum(sqrt(mean(Ini_X-Minor_Points_X,2).^2+mean(Ini_Y-Minor_Points_Y,2).^2));
            Vec_Before=[Ini_X(:,2)-Ini_X(:,1),Ini_Y(:,2)-Ini_Y(:,1)];
            Vec_After=[Minor_Points_X(:,2)-Minor_Points_X(:,1),Minor_Points_Y(:,2)-Minor_Points_Y(:,1)];
            Total_Rot_Min(i_Time)=Total_Rot_Maj(i_Time-1)+sum(acosd(dot(Vec_Before',Vec_After')./(vecnorm(Vec_Before').*vecnorm(Vec_After'))));
        end
    
        %===Correct location of filaments to make sure they remain in cytoplasm
        [Major_Points_X,Major_Points_Y]=Correct_Filaments(Connectivity,Nodes,Major_Points_X,Major_Points_Y,N_Cyto_Elem);
        [Minor_Points_X,Minor_Points_Y]=Correct_Filaments(Connectivity,Nodes,Minor_Points_X,Minor_Points_Y,N_Cyto_Elem);
        
        %===Output presentation
        %plot mesh
        close all
        figure('Name',num2str(floor(i_Time*Time_Step/60)),'Visible','off');
        triplot(Connectivity',Nodes(1,:)',Nodes(2,:)','color',[0.7 0.7 0.7]);
        %t=Time+Time_Step;
        %hr=int2str(floor(t/3600));
        %t_rem=rem(t,3600);
        %minut=int2str(floor(t_rem/60));
        %secon=int2str(round(rem(t_rem,60)));
        %title(['t= ',hr,':',minut,':',secon,' hr:min:sec'])
        title (['Iteration= ', int2str(i_Time)])
        hold on
        axis square
        axis off
        %set(gcf,'color','w','Units','normalized','OuterPosition',[0 0 1 1]);
        set(gcf,'color','w','Units','normalized');
        xlim([x_min x_max]) 
        ylim([y_min y_max])
        
        %plot filaments
        plot (Major_Points_X',Major_Points_Y','b-','LineWidth',0.5)
        plot (Minor_Points_X',Minor_Points_Y','r-','LineWidth',0.5)
        
        %plot supports
        plot(Nodes(1,Supports),Nodes(2,Supports),'o','MarkerEdgeColor','k','MarkerFaceColor','g')
    
        %plot membrane
        Mem_Elem=[Membrane_Nodes;circshift(Membrane_Nodes,-1)];
        plot(Nodes(1,Mem_Elem)',Nodes(2,Mem_Elem)','-k')
    
        %plot nucleus
        Nuc_Elem=[Nucleus_Nodes;circshift(Nucleus_Nodes,-1)];
        plot(Nodes(1,Nuc_Elem)',Nodes(2,Nuc_Elem)','-k','LineWidth',0.75)
        
        U3_M(i_Time+1)=getframe(gcf);
        drawnow 
        waitbar(i_Time/N_Step,WB,sprintf('Progress: %.2f %%',round(i_Time/N_Step*100,2)));
    end
    close (WB);
    print(gcf,'-vector','-dsvg','Major_Minor.svg')
    video=VideoWriter('Filament Evolution.avi', 'Uncompressed AVI');
    video.FrameRate=10;
    open(video)
    writeVideo(video,U3_M);
    close(video);

    save('Displacemet_Rot.mat',"Total_Disp_Maj","Total_Disp_Min","Total_Rot_Maj","Total_Rot_Min");
    %==============================Plot figures for paper==================
    %==Only Major filaments
    figure ('Visible','off')
    triplot(Connectivity',Nodes(1,:)',Nodes(2,:)','color',[0.7 0.7 0.7]);
    hold on
    axis square
    axis off
    %set(gcf,'color','w','Units','normalized','OuterPosition',[0 0 1 1]);
    set(gcf,'color','w','Units','normalized');
    xlim([x_min x_max]) 
    ylim([y_min y_max])
    
    %plot filaments
    plot (Major_Points_X',Major_Points_Y','b-','LineWidth',0.75)
    
    %plot supports
    plot(Nodes(1,Supports),Nodes(2,Supports),'o','MarkerEdgeColor','k','MarkerFaceColor','g')
    
    %plot membrane
    Mem_Elem=[Membrane_Nodes;circshift(Membrane_Nodes,-1)];
    plot(Nodes(1,Mem_Elem)',Nodes(2,Mem_Elem)','-k')
    
    %plot nucleus
    Nuc_Elem=[Nucleus_Nodes;circshift(Nucleus_Nodes,-1)];
    plot(Nodes(1,Nuc_Elem)',Nodes(2,Nuc_Elem)','-k','LineWidth',0.75)
    print(gcf,'-vector','-dsvg','Major.svg')

    %==Only Minor filaments
    figure ('Visible','off')
    triplot(Connectivity',Nodes(1,:)',Nodes(2,:)','color',[0.7 0.7 0.7]);
    hold on
    axis square
    axis off
    %set(gcf,'color','w','Units','normalized','OuterPosition',[0 0 1 1]);
    set(gcf,'color','w','Units','normalized');
    xlim([x_min x_max]) 
    ylim([y_min y_max])
    
    %plot filaments
    plot (Minor_Points_X',Minor_Points_Y','r-','LineWidth',0.75)
    
    %plot supports
    plot(Nodes(1,Supports),Nodes(2,Supports),'o','MarkerEdgeColor','k','MarkerFaceColor','g')
    
    %plot membrane
    Mem_Elem=[Membrane_Nodes;circshift(Membrane_Nodes,-1)];
    plot(Nodes(1,Mem_Elem)',Nodes(2,Mem_Elem)','-k')
    
    %plot nucleus
    Nuc_Elem=[Nucleus_Nodes;circshift(Nucleus_Nodes,-1)];
    plot(Nodes(1,Nuc_Elem)',Nodes(2,Nuc_Elem)','-k','LineWidth',0.75)
    print(gcf,'-vector','-dsvg','Minor.svg')

    %Plot Mean Strain and Direction of Major principal strain
    Xn=[Nodes(1,Connectivity(1,:))',Nodes(1,Connectivity(2,:))',Nodes(1,Connectivity(3,:))'];
    Yn=[Nodes(2,Connectivity(1,:))',Nodes(2,Connectivity(2,:))',Nodes(2,Connectivity(3,:))'];
    figure ('Visible','off')
    Cn=[Total_Strain_MEAN;zeros(45,1)];
    T=[T_ang;zeros(45,1)];
    patch(Xn',Yn',Cn,'edgecolor',[0.5,0.5,0.5])
    colorbar
    colormap(jet)
    axis square
    axis off
    Xm=mean(Xn,2);
    Ym=mean(Yn,2);
    L=0.5;
    XL=[Xm+L*sin(T),Xm-L*sin(T)];
    YL=[Ym+L*cos(T),Ym-L*cos(T)];
    hold on
    plot(XL',YL','k')
    print(gcf,'-vector','-dsvg','Strain.svg')

    %Plot Average Accumulative Displacement and Rotation
    figure ("Visible","off")
    Iteration=0:N_Step-1;
    plot(Iteration,Total_Disp_Maj/N_Major,'b-','LineWidth',1.5)
    hold on
    plot(Iteration,Total_Disp_Min/N_Minor,'r-','LineWidth',1.5)
    legend({"Major","Minor"},'Location','Northwest')
    xlabel("Iteration")
    ylabel("Accumulative Displacement (\mum)")
    fontsize(14,"Points")
    print(gcf,'-vector','-dsvg','Disp.svg')

    figure ("Visible","off")
    plot(Iteration,Total_Rot_Maj/N_Major,'b-','LineWidth',1.5)
    hold on
    plot(Iteration,Total_Rot_Min/N_Minor,'r-','LineWidth',1.5)
    legend({"Major","Minor"},'Location','Northwest')
    xlabel("Iteration")
    ylabel("Accumulative Rotation (\circ)")
    fontsize(14,"Points")
    print(gcf,'-vector','-dsvg','Rotate.svg')

end