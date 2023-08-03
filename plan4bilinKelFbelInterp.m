function [ Kel, fe ] = plan4bilinKelFbelInterp( ex, ey, ep,eq)
%
% Syntax:
% [ Ke, fe ] = plan4bilinKelFbelInterp( ex, ey, ep, eq )
%-------------------------------------------------------------
% PURPOSE
%  Compute the stiffness matrix and element external force vector
%  for a bilinear plane element
%
% INPUT:  ex = [ x1 x2 x3 x4 ]         element nodal x-coordinates
%         ey = [ y1 y2 y3 y4 ]         element nodal y-coordinates
%
%         ep = [ngp,lambda-x,lambda-y]
%
%              ngp:  number of Gauss points 
%                    in each direction ( ksi 
%                    and eta) (ngp = 1 or 2 
%                    or 3)
%
%         eq = [ bz]                    bz: force z-dir
%--------------------------------------------------------------
% OUTPUT: Kel : element stiffness matrix (4 x 4)
%         fel : equivalent nodal forces (4 x 1)
%--------------------------------------------------------------
% 
%
% MODIFIED for MEF by Luis Verduzco 2022-11-01
%--------------------------------------------------------------

ngp   = ep(1); % Total gauss points = ( NoGaussPoits per direction )^2

%  Tolerance for the Jacobian determinant
minDetJ = 1.e-16;

ngp   = ngp^2; % Total gauss points = ( NoGaussPoits per direction )^2

if ngp == 1
 
    intWeight=[2,2];
    GaussPoints=[0,0];

elseif ngp == 4
 
     intWeight=[1,1;
                1,1];


     GaussPoints=[-0.57735026918962,-0.57735026918962;
                  0.57735026918962,0.57735026918962];
                
         
elseif ngp == 9
 
     intWeight=[0.55555555555,0.55555555555
                0.888888888888,0.888888888888;
                0.55555555555,0.55555555555];

     GaussPoints=[-0.774596666,-0.774596666;
                 0,0;
                 0.774596666,0.774596666];
               
else
 
 error('Only 1,2 or 3 Gauss Points in each direction apply')

end

gauss_point=0;

%  Initialize Ke and fe with zeros for all of their elements
lambdax=ep(2); lambday=ep(3);

fe=zeros(4,1);
Kex=zeros(4); Key=zeros(4);
for punto_gaus_xsi=1:ngp^0.5
    for punto_gaus_eta=1:ngp^0.5
        
        gauss_point=gauss_point+1;
        
        % Compute derivatives (with respect to xsi and eta) of the
        % shape functions at coordinate (xsi,eta). Since the element is
        % isoparametic, these are also the derivatives of the basis functions.
        xsi=GaussPoints(punto_gaus_xsi,1);
        weightXsi=intWeight(punto_gaus_xsi,1);
        
        eta=GaussPoints(punto_gaus_eta,2);
        weightEta=intWeight(punto_gaus_eta,2); 
    
        Ne=[(xsi-1).*(eta-1)/4  -(1+xsi).*(eta-1)/4 ... 
            (xsi+1).*(eta+1)/4 -(xsi-1).*(1+eta)/4]; % 1x4

        dNr=[-(1-eta)/4 (1-eta)/4 (1+eta)/4 -(1+eta)/4;
             -(1-xsi)/4 -(1+xsi)/4 (1+xsi)/4 (1-xsi)/4];
     
        %  Use shape function derivatives and element vertex coordinates 
        %  to establish the Jacobian matrix.
    
        JT=dNr*[ex;ey]';

        detJ=det(JT);

        % Determinant seems OK - invert the transpose of the Jacobian
        JTinv=inv(JT);
        % Compute derivatives with respect to x and y, of all basis 
        % functions Ninv=inv(Ne);
        dNxy=JTinv*dNr;

        Bxsieta=[dNxy(1,1) dNxy(1,2) dNxy(1,3) dNxy(1,4);
                 dNxy(2,1) dNxy(2,2) dNxy(2,3) dNxy(2,4)];
             
        Kex=Kex+lambdax*Bxsieta(1,:)'*Bxsieta(1,:)*detJ*weightXsi*weightEta;
        Key=Key+lambday*Bxsieta(2,:)'*Bxsieta(2,:)*detJ*weightXsi*weightEta;
        
        % Compute the contribution to element stiffness and load matrix 
        % from current Gauss point
        Kel=Kex+Key;

        % Forces
        fe=fe+Ne'*detJ*weightXsi*weightEta*eq;
    end

end
