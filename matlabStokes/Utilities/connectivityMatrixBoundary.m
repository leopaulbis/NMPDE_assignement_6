function Tboundary=connectivityMatrixBoundary(T,referenceElement)

if size(referenceElement.faceNodes,1)==3 %triangle
    aux=[T(:,1) T(:,2)
        T(:,2) T(:,3)
        T(:,3) T(:,1)];
else %quadrilateral
    aux=[T(:,1) T(:,2)
        T(:,2) T(:,3)
        T(:,3) T(:,4)
        T(:,4) T(:,1)];
end
nOfNodes=max(max(T));

n1=aux(:,1); n2=aux(:,2); %1st and 2nd node of each face (internal faces are repeated with 2 orientations)

A=sparse(n1,n2,ones(size(n1)),nOfNodes,nOfNodes); %A(i,j)=1 if face xi-xj exists with proper orientation
B=A+A'; %A(i,j)=2 if face xi-xj if interior face, =1 if exterior face
[n1,n2]=find(B==1 & not(A==0));
[n1,permutation] = sort(n1); n2=n2(permutation);

% (n1,n2) are the vertexes of the faces on the boundary
% Tboundary=[n1,n2]; if linear
nOfElementFaces=size(referenceElement.faceNodes,1);
Tboundary=[];
for i=1:length(n1)
   [e1,kk]=find(T==n1(i)); [e2,kk]=find(T==n2(i));
   elem = intersect(e1,e2);
   Te=T(elem,:); 
   for f=1:nOfElementFaces
      Tf=Te(referenceElement.faceNodes(f,:));
      if any(Tf==n1(i)) & any(Tf==n2(i)) 
          Tboundary=[Tboundary;Tf]; break;
      end
   end
end



