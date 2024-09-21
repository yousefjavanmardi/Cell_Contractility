function [Stress_MEAN,T]=Strain_Finder(Connectivity,Nodes,U)%E,Nu
    N_Elem=size(Connectivity,2);
    Stress_MEAN=zeros(N_Elem,1);
    T=zeros(N_Elem,1);
    for i_elem=1:N_Elem
        Nodes_of_Elem = Connectivity(:,i_elem);
        X_Elem = Nodes(1,Nodes_of_Elem);
        Y_Elem = Nodes(2,Nodes_of_Elem);
        b1=Y_Elem(2)-Y_Elem(3);
        b2=Y_Elem(3)-Y_Elem(1);
        b3=Y_Elem(1)-Y_Elem(2);
        g1=X_Elem(3)-X_Elem(2);
        g2=X_Elem(1)-X_Elem(3);
        g3=X_Elem(2)-X_Elem(1);
        Double_A=X_Elem(1)*b1+X_Elem(2)*b2+X_Elem(3)*b3;
        B=[b1,0,b2,0,b3,0;0,g1,0,g2,0,g3;g1,b1,g2,b2,g3,b3]/Double_A;
        U_Elem = [U(2*Nodes_of_Elem(1)-1);U(2*Nodes_of_Elem(1));...
                 U(2*Nodes_of_Elem(2)-1);U(2*Nodes_of_Elem(2));...
                 U(2*Nodes_of_Elem(3)-1);U(2*Nodes_of_Elem(3))];
        S=B*U_Elem;
        Stress_MEAN(i_elem)=0.5*(S(1)+S(2));

        T(i_elem)=pi()/2-0.5*atan2(0.5*S(3),0.5*(S(1)-S(2)));
    end 
end