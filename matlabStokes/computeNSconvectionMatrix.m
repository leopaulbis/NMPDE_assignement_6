function C=computeNSconvectionMatrix(u,X,T,referenceElement)

GaussPoints=referenceElement.IPcoordinates;
GaussWeights=referenceElement.IPweights;
N=referenceElement.N;
Nxi=referenceElement.Nxi;
Neta=referenceElement.Neta;
Np = referenceElement.Npressure;

nOfNodes = size(X,1);
nOfElements = size(T,1);

C=spalloc(2*nOfNodes,2*nOfNodes,8*nOfNodes);%global matrix initiallization

%Loop in elements
for i=1:nOfElements
    Te=T(i,:); %nodes in the element
    Xe=X(Te,:);
    ass = [Te,nOfNodes+Te];
    ue=u(ass);
    Ce=computeElementalMatrix(ue,Xe,GaussPoints,GaussWeights,N,Nxi,Neta);
    C(ass,ass)=C(ass,ass)+Ce; %assembly
end

%_______________________________________
%Computation of elemental matrix and vector
function Ce=computeElementalMatrix(ue,Xe,GaussPoints,GaussWeights,N,Nxi,Neta)


nOfNodes = size(Xe,1);
Ce=zeros(2*nOfNodes);
xe = Xe(:,1); ye = Xe(:,2); %x and y coordinates of the element nodes
%Loop in integration points
for k=1:length(GaussWeights)
    Nk=N(k,:);
    dNkdxi=Nxi(k,:);
    dNkdeta=Neta(k,:); 
    xk = Nk*Xe;
    %Jacobian of the isoparametric transformation 
    J = [dNkdxi*xe dNkdxi*ye;dNkdeta*xe dNkdeta*ye];
    %Computation of the derivatives of the shape functions with respecto to x and y
    % B*ue: approximation of (u1 u2)' at Gauss point
    % S*ue: approximation of (du1/dx du1/dy du2/dx du2/dy)' at point
    Grad1k = J\[dNkdxi;dNkdeta];
    S = zeros(4,2*nOfNodes);
    S(1:2,1:nOfNodes) = Grad1k;
    S(3:4,nOfNodes+1:2*nOfNodes) = Grad1k;
    B = zeros(2,2*nOfNodes);
    B(1,1:nOfNodes) = Nk;
    B(2,1+nOfNodes:2*nOfNodes) = Nk;
    %integration weight at physical element
    dxy=GaussWeights(k)*det(J);
    %contribution of the integration point to the elemental matrix
    aux=S*ue;
    gradu=[aux(1) aux(2); aux(3) aux(4)];
    Ce =  Ce + B'*gradu*B*dxy;
end
