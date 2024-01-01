function [K,G,f]=computeSystemStokes(X,T,Xp,Tp,referenceElement,source,nu)

GaussPoints=referenceElement.IPcoordinates;
GaussWeights=referenceElement.IPweights;
N=referenceElement.N;
Nxi=referenceElement.Nxi;
Neta=referenceElement.Neta;
Np = referenceElement.Npressure;

nOfNodes = size(X,1);
nOfElements = size(T,1);
nOfNodesp = size(Xp,1);

K=spalloc(2*nOfNodes,2*nOfNodes,8*nOfNodes);%global matrix initiallization
G=spalloc(2*nOfNodes,nOfNodesp,2*nOfNodes);
f=zeros(2*nOfNodes,1);

%Loop in elements
for i=1:nOfElements
    Te=T(i,:); %nodes in the element
    Tep = Tp(i,:);
    Xe=X(Te,:);
    Xep = Xp(Tep,:);%coordinates of the element nodes
    [Ke,Ge,fe]=computeElementalMatrices(Xe,Xep,GaussPoints,GaussWeights,N,Nxi,Neta,Np,source);
    ass = [Te,nOfNodes+Te];
    K(ass,ass)=K(ass,ass)+nu*Ke; %assembly
    G(ass,Tep) = G(ass,Tep) + Ge;
    f(ass) = f(ass) + fe;
end

%_______________________________________
%Computation of elemental matrix and vector
function [Ke,Ge,fe]=computeElementalMatrices(Xe,Xep,GaussPoints,GaussWeights,N,Nxi,Neta,Np,source)

nOfNodes = size(Xe,1);
nOfNodesp = size(Xep,1);
Ke=zeros(2*nOfNodes);
Ge = zeros(2*nOfNodes,nOfNodesp);
fe=zeros(2*nOfNodes,1);
xe = Xe(:,1); ye = Xe(:,2); %x and y coordinates of the element nodes
%Loop in integration points
for k=1:length(GaussWeights)
    Nk=N(k,:);
    dNkdxi=Nxi(k,:);
    dNkdeta=Neta(k,:); 
    Nkp = Np(k,:);
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
    Ke =  Ke + S'*S*dxy;%info de la forma debil
    Ge = Ge - ([1,0,0,1]*S)'*Nkp*dxy;
    fe = fe + B'*source(xk)*dxy;             
end
