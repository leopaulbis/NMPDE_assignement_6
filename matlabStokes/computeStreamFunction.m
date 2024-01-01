function phi=computeStreamFunction(ux,uy,X,T,Tboundary,referenceElement)

% Computation of the matrix (discretization of Laplacian)
IPweights = referenceElement.IPweights;
IPcoord = referenceElement.IPcoordinates;
N=referenceElement.N; Nxi=referenceElement.Nxi; Neta=referenceElement.Neta;

nOfNodes = size(X,1); nOfElements = size(T,1);

K=spalloc(nOfNodes,nOfNodes,9*nOfNodes);
f=zeros(nOfNodes,1);
%Loop in elements
for i=1:nOfElements
    Te=T(i,:); %index of the nodes in the element
    Xe=X(Te,:); %coordinates of the nodes in the element
    [Ke,fe]=elementalMatrix(ux(Te),uy(Te),Xe,IPcoord,IPweights,N,Nxi,Neta);
    K(Te,Te)=K(Te,Te)+Ke; %assembly of elemental matrix
    f(Te) = f(Te) + fe;
end

%Computation of r.h.s.
nOfElements = size(Tboundary,1);

IPcoord = referenceElement.IPcoordinates1d; 
IPweights = referenceElement.IPweights1d; 
N = referenceElement.N1d; 

for i=1:nOfElements
    Te=Tboundary(i,:); %index of the nodes in the element
    Xe=X(Te,:); %coordinates of the nodes in the element
    fe=elementalVector(ux(Te),uy(Te),Xe,IPcoord,IPweights,N);
    f(Te) = f(Te) + fe;
end
 
% %phi(1)=0
% K=K(2:end,2:end); f=f(2:end);
% phi=[0;K\f];
%phi(end)=0
K=K(1:end-1,1:end-1); f=f(1:end-1);
phi=[K\f;0];


%____________________________________
function fe=elementalVector(uxe,uye,Xe,IPcoord,IPweights,N)

nnodes = size(Xe,1);
fe=zeros(nnodes,1);

x1 = Xe(1,:); x2 = Xe(end,:);
t = x2 - x1; h = norm(t);  
n = [t(2), -t(1)]; n = n/h; 
%Loop in integration points
for k=1:length(IPweights)
    Nk=N(k,:); %basis functions at integration point
    xk = Nk*Xe; %xy-coordinates of the integration point
    dl=IPweights(k)*h/2; %diferential of line
    uxk=Nk*uxe; uyk=Nk*uye; %velocity at integration point
    fe = fe + Nk'*(n(1)*uyk-n(2)*uxk)*dl;
end
  

%_______________________________________
%Calcul de la matriu elemental
function [Ke,fe]=elementalMatrix(uxe,uye,Xe,IPcoord,IPweights,N,Nxi,Neta)

nnodes = size(Xe,1);
Ke=zeros(nnodes); fe=zeros(nnodes,1);
xe = Xe(:,1); ye = Xe(:,2);
%Bucle en punts d'integraci?
for k=1:length(IPweights)
    Nk=N(k,:);
    Nkxi=Nxi(k,:);
    Nketa=Neta(k,:); 
    xk = Nk*Xe; 
    %Jacobia 
    J = [Nkxi*xe Nkxi*ye;Nketa*xe Nketa*ye];
    % Derivadas de las funciones de forma respecto a (x,y)
    Nkxy = J\[Nkxi;Nketa];
    Nkx = Nkxy(1,:); Nky = Nkxy(2,:);
    %diferencial de volum
    dxy=IPweights(k)*det(J);
    Gk=Nkxy;
    Ke = Ke + Gk'*Gk*dxy;
    fe = fe + Nk'*(Nky*uxe-Nkx*uye)*dxy;
end
