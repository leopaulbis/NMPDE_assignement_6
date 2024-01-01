function contourPlot(u,X,T)

p=X';
switch size(T,2)
    case 3 %linear triangles
        t=[T';ones(1,size(T,1))];
    case 4 %linear quadrilateral
        t=[[T(:,[1,2,3]);T(:,[1,3,4])]';ones(1,2*size(T,1))]; 
    case 9; %quadratic quadrilateral
        t=[[T(:,[1,5,9]);...
            T(:,[1,9,8]);...
            T(:,[5,2,6]);...
            T(:,[5,6,9]);...
            T(:,[8,9,7]);...
            T(:,[8,7,4]);...
            T(:,[9,6,3]);...
            T(:,[9,3,7]);...
            ]';ones(1,8*size(T,1))];         
    otherwise
        error('This function work only for linear and Q2 FEM :-(')
end
  
pdecont(p,t,u,40); %to see vortices, augmenter le nombre de surfaces de niveau 
%pdesurf(p,t,u);