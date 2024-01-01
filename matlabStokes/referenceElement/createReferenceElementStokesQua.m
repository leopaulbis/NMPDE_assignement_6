function theReferenceElement=createReferenceElementStokesQua(degree)

% theReferenceElement=createReferenceElementStokes(degree)
% Input:
%  degree: k for degree k for velocity and k-1 for pressure
% Output:
%  theReferenceElement: struct containing
%     .IPcoordinates: coordinates of the integration points for 2D elemens
%     .IPweights: weights of the integration points for 2D elements
%     .N: shape functions at the IP
%     .Nxi,.Neta: derivatives of the shape functions at the IP
%     .NP: shape functions for pressure
%     .IPcoordinates1d: coordinates of the integration points for 1D boundary elemens
%     .IPweights1d: weights of the integration points for 1D boundary elements
%     .N1d: 1D shape functions at the IP
%     .N1dxi: derivatives of the 1D shape functions at the IP
%     .faceNodes: matrix [nOfFaces nOfNodesPerFace] with the edge nodes numbering
%     .innerNodes: vector [1 nOfInnerNodes] with the inner nodes numbering
%     .faceNodes1d: vector [1 nOfNodesPerElement] with the 1D nodes numbering
%     .NodesCoord: spatial coordinates of the element nodes
%     .NodesCoord1d: spatial coordinates of the 1D element nodes

theReferenceElement=createReferenceElementQua(degree);
theReferenceElementP=createReferenceElementQua(degree-1);

theReferenceElement.NodesCoordP=theReferenceElementP.NodesCoord;
NP=evaluateNodalBasisQuawithoutDerivatives(theReferenceElement.IPcoordinates,theReferenceElementP.NodesCoord,degree-1);
theReferenceElement.Npressure=NP; 
theReferenceElement.IPweights1d=(theReferenceElement.IPweights1d)';


