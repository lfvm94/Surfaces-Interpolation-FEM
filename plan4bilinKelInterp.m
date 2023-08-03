function [Kel] = plan4bilinKelInterp(ex,ey, dx, dy, xk, yk)
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
%--------------------------------------------------------------
% OUTPUT: Kel : element stiffness matrix (4 x 4)
%--------------------------------------------------------------
%
% MODIFIED for MEF by Luis Verduzco 2022-11-01
%--------------------------------------------------------------

% Compute derivatives of shape functions:
Ae=dx*dy;
Ne=1/Ae*[(xk-ex(2))*(yk-ey(4)), (xk-ex(1))*(yk-ey(3)),...
         (xk-ex(4))*(yk-ey(2)), (xk-ex(3))*(yk-ey(1))];
Kel=Ne'*Ne;
        