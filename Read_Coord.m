function [Boundary_Coord,Nucleus_Coord,Support_Coord]= Read_Coord(Boundary_File,Nucleus_File,Support_File)
    fileID = fopen(Boundary_File, 'r');
    Boundary_Coord = fscanf(fileID,'%f',[2 Inf])';
    fclose(fileID);
    
    fileID = fopen(Nucleus_File, 'r');
    Nucleus_Coord = fscanf(fileID,'%f',[2 Inf])';
    fclose(fileID);

    fileID = fopen(Support_File, 'r');
    Support_Coord = fscanf(fileID,'%f',[2 Inf])';
    fclose(fileID);
end