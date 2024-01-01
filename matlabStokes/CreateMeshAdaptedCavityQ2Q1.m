function [X,T,XP,TP] = CreateMeshAdaptedCavityQ2Q1(nelem1d)

[X,T] = CreateMesh(0,9,0,1,0,1,2*nelem1d+1,2*nelem1d+1);
[XP,TP] = CreatePresNod(0,4,X,T,2*nelem1d+1,2*nelem1d+1);

end
%_____________________________________________________________
function [X,T] = CreateMesh(elem,nen,x1,x2,y1,y2,npx,npy)
% [X,T] = CreateMesh(elem,nen,x1,x2,y1,y2,nx,ny)
% Creates the topology of an structured mesh over a
% rectangular domain [x1,x2]x[y1,y2]
%
% Input:    
%   elem:           element type (0 for quadrilaterals, 1 for triangles)
%   nen:            number of element nodes
%   x1,x2,y1,y2:    vertices' coordinates
%   npx,npy:        number of nodes on each direction 
% Output:   
%   X:              nodal coordinates
%   T:              connectivities


% Number of elements in each direction
nx = npx-1; ny = npy-1;

% Allocate of storage for the nodal coordinates matrix
X = zeros((npx)*(npy),2);

hx = (x2-x1)/nx;
hy = (y2-y1)/ny;
xs = linspace(x1,x2,npx)'; 
xs(2:nx)  = (x2-x1)*(tanh(5/2)+tanh(5*(xs(2:nx) -(x1+x2)/2)))/(2*tanh(5/2))+x1;
unos = ones(npx,1);
% Nodes' coordinates
yys = linspace(y1,y2,npy);
yys(2:ny) = (y2-y1)*(tanh(5/2)+tanh(5*(yys(2:ny)-(y1+y2)/2)))/(2*tanh(5/2))+y1;
for i=1:npy
     ys = yys(i)*unos; 
    posi = [(i-1)*(npx)+1:i*(npx)]; 
    X(posi,:)=[xs,ys];
end

% Connectivities
if elem == 0            % Quadrilaterals
    if nen == 4             % Q1
        T = zeros(nx*ny,4);
        for i=1:ny
            for j=1:nx
                ielem = (i-1)*nx+j;
                inode = (i-1)*(npx)+j;
                T(ielem,:) = [inode inode+1 inode+(nx+2) inode+(npx)];
            end   
        end
    elseif nen == 9         % Q2
        if (nx-2*floor(nx/2) ~=0) & (ny-2*floor(ny/2) ~=0)
            error('Number of nodes in X or Y direction is not odd')
        end
        for i=1:ny/2
            for j=1:nx/2
                ielem=(i-1)*nx/2+j;
                inode=(i-1)*2*(npx)+2*(j-1)+1;
                T(ielem,:)=[inode inode+2 inode+2*(npx)+2 inode+2*(npx) inode+1 inode+2+(npx) inode+2*(npx)+1 inode+(npx) inode+(npx)+1 ];
            end
        end
    end
elseif elem == 1        % Triangles
    if nen == 3             % P1
        T = zeros(nx*ny,3);
        for i=1:ny
            for j=1:nx
                ielem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*(npx)+j;
                T(ielem,:) = [inode   inode+1   inode+(npx)];
                T(ielem+1,:) = [inode+1   inode+1+npx   inode+npx];
            end   
        end
        % Modification of left lower and right upper corner elements to avoid them 
        % having all their nodes on the boundary
        T(1,:) = [1  npx+2   npx+1];
        T(2,:) = [1    2     npx+2];
        aux = size(T,1);
        T(aux,:) = [npx*ny-1    npx*npy   npx*npy-1];
        T(aux-1,:)   = [npx*ny-1    npx*ny    npx*npy];
    elseif nen == 4         % P1+  (Note that "bubble" coordinates are not considered)
        T = zeros(nx*ny,4);
        for i=1:ny
            for j=1:nx
                ielem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*(npx)+j;
                n_ad = npx*npy + 2*((i-1)*nx+j)-1;
                T(ielem,:) = [inode   inode+1   inode+(npx)  n_ad];
                T(ielem+1,:) = [inode+1   inode+1+npx   inode+npx n_ad+1];
            end   
        end
        % Modification of left lower and right upper corner elements to avoid them 
        % having all their nodes on the boundary
        aux = size(T,1);
        T(1,:) = [1  npx+2   npx+1  npx*npy+1];
        T(2,:) = [1    2     npx+2  npx*npy+2];
        T(aux-1,:) = [npx*ny-1    npx*npy   npx*npy-1   npx*npy+2*nx*ny-1];
        T(aux,:)   = [npx*ny-1    npx*ny    npx*npy   npx*npy+2*nx*ny];
    elseif nen == 6         % P2
        for i=1:ny/2
            for j=1:nx/2
                ielem=2*((i-1)*nx/2+j)-1;
                inode=(i-1)*2*(npx)+2*(j-1)+1;
                T(ielem,:) = [inode   inode+2   inode+2*npx   inode+1    inode+1+npx   inode+npx];
                T(ielem+1,:) = [inode+2    inode+2+2*npx   inode+2*npx   inode+2+npx   inode+1+2*npx   inode+1+npx];
            end    
        end
        % Modification of left lower and right upper corner elements to avoid them 
        % having all their nodes on the boundary
        T(1,:) = [1   2*npx+3    2*npx+1   npx+2   2*npx+2   npx+1];
        T(2,:) = [1      3       2*npx+3     2      npx+3    npx+2];
        aux = size(T,1);
        T(aux-1,:) = [npx*(ny-1)-2   npx*(ny-1)    npx*npy    npx*(ny-1)-1    npx*ny     npx*ny-1];
        T(aux,:)   = [npx*(ny-1)-2    npx*npy     npx*npy-2     npx*ny-1     npx*npy-1   npx*ny-2 ];
    end
end   

end
%_________________________________________________________________
function [XP,TP] = CreatePresNod(elem,nenP,X,T,npx,npy)
% [XP,TP] = CreatePresNod(elem,nenP,X,T,npx,npy)
% Creates a mesh for the pressure field with nenP nodes per element,
% based on a velocity mesh (X,T)
%
% Input:    
%   elem:       element type (0 for quadrilaterals, 1 for triangles)
%   nenP:       number of element nodes
%   X,T:        nodal coordinates and connectivities of the velocity mesh
%   npx,npy:    number of nodes on each direction 
%
% Output:   
%   XP:         nodal coordinates forthe pressure mesh
%   TP:         connectivities of the pressure mesh
% 


[numel,nen] = size(T); 
numnp = size(X,1);

nx = npx-1;  ny = npy-1;

if elem == 0  % Quadrilateral
    if nenP == 1  % P0
        XP = zeros(numel,2); TP = zeros(numel,nenP);
        for ielem=1:numel
            XP(ielem,:) = mean(X(T(ielem,:),:));
            TP(ielem)   = ielem;
        end        
    elseif nenP == 4 
        if nen == 4  % Q1Q1
            XP = X; TP = T;
        elseif  nen==9 % Q2Q1
            % Only for structured meshes!!!
            npx=nx/2+1; npy=ny/2+1;
            XP = zeros(npx*npy,2); TP = zeros(numel,nenP);
            for irow=1:npy
                irowodd = 2*(irow-1)+1;
                XP((irow-1)*npx+1:irow*npx,:) = X((irowodd-1)*(nx+1)+1:2:irowodd*(nx+1),:);
            end
            for i=1:ny/2
                for j=1:nx/2
                    ielem = (i-1)*nx/2+j;
                    inode = (i-1)*(npx)+j;
                    TP(ielem,:) = [inode inode+1 inode+(npx+1) inode+(npx)];
                end   
            end
        end
    elseif nenP==9 & nen==9 % Q2Q2
        XP = X; TP = T;
    end
elseif elem == 1   % Triangles
    if nenP == 1 %P0
        XP = zeros(numel,2); TP = zeros(numel,nenP);
        for ielem=1:numel
            XP(ielem,:) = mean(X(T(ielem,:),:));
            TP(ielem)   = ielem;
        end
    elseif nenP == 3
        if nen == 3  % P1P1
            XP = X; 
            TP = T;
        elseif nen == 4  % MINI (P1+P1)
            XP = X(1:(nx+1)*(ny+1),:);
            TP = zeros(nx*ny,3); TP = T(:,1:3);
        elseif nen == 6  % P2P1
            npx=nx/2+1; npy=ny/2+1;
            XP = zeros(npx*npy,2);
            TP = zeros(nx*ny/4,3);
            for irow=1:npy
                irowodd = 2*(irow-1)+1;
                XP((irow-1)*npx+1:irow*npx,:) = X((irowodd-1)*(nx+1)+1:2:irowodd*(nx+1),:);
            end
            for i=1:ny/2
                for j=1:nx/2
                    ielem = 2*((i-1)*(nx/2)+j)-1;
                    inode = (i-1)*(npx)+j;
                    TP(ielem,:) = [inode   inode+1   inode+(npx)];
                    TP(ielem+1,:) = [inode+1   inode+1+npx   inode+npx];
                end
            end
            TP(1,:) = [1  npx+2   npx+1];
            TP(2,:) = [1     2    npx+2];
            aux = size(TP,1);
            TP(aux,:) =     [npx*(npy-1)-1   npx*(npy-1)    npx*npy];
            TP(aux-1,:)   = [npx*(npy-1)-1     npx*npy     npx*npy-1];
        end
    elseif nenP == 6 & nen == 6  % P2P2
        XP = X;
        TP = T;
    end
end

end