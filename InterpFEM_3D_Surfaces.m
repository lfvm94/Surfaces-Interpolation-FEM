% ---------------------------------------------------------------------
% HW-12 Interpolation with the Finite Element Method for 3D surfaces
% ---------------------------------------------------------------------
% By Luis Fernando Verduzco Martínez
% All Rights Reserved (c)
% ---------------------------------------------------------------------
clear all
clc

% ---------------------------------------------------------------------
% Plotting original surface:
% ---------------------------------------------------------------------

rangexy=[-pi,pi]; % squared range

% ---------------------------------------------------------------------
% Generating sample points:
% ---------------------------------------------------------------------
samplePoints=100;
typeSampleData=1; % 1: random, 2: uniform
if typeSampleData==1
    % To generate a random sample of data points (more realistic)

    x=rangexy(1)+rand(samplePoints,1)*(rangexy(2)-rangexy(1));
    y=rangexy(1)+rand(samplePoints,1)*(rangexy(2)-rangexy(1));

elseif typeSampleData==2
    % To generate a uniform sample of data points
    npsxy=fix(sqrt(samplePoints));
    nps=npsxy^2;
    samplePoints=nps;
    dxx=(rangexy(2)-rangexy(1))/(npsxy-1);
    dyy=dxx;
    yi=rangexy(1);
    for i=1:npsxy
        xi=rangexy(1);
        for j=1:npsxy
            x((i-1)*npsxy+j)=xi+(j-1)*dxx;
            y((i-1)*npsxy+j)=yi;
        end
        yi=yi+dyy;
    end
end
z=fxy(x,y);

% ---------------------------------------------------------------------
% Generate discretization mesh
% ---------------------------------------------------------------------

domxy=rangexy(2)-rangexy(1); % this is the range of the sample
npxy=12; % the n nodes on each direction x,y to generate the FE mesh

% Generate discretization mesh

nelxy=npxy-1; % this is the number of finite elements
dx=domxy/(nelxy); % this is the length of the elements (uniform for all 
			  % elements)
dy=domxy/(nelxy); % this is the length of the elements (uniform for all 
			  % elements)
              
mallax=rangexy(1):dx:rangexy(2);
mallay=rangexy(1):dy:rangexy(2);

% ---------------------------------------------------------------------
% Coordinates of mesh
% ---------------------------------------------------------------------

npoints_mesh=npxy^2;
Coord=zeros(npoints_mesh,2);
for i=1:npxy % to traverse mesh in y
    coordy=mallay(i);
    for j=1:npxy % to traverse mesh in x
        Coord((i-1)*npxy+j,1)=mallax(j);
        Coord((i-1)*npxy+j,2)=coordy;
        
    end
end

% --------------------------------------------------------------------
% Topology
% --------------------------------------------------------------------
nelem=nelxy^2;
ex=zeros(nelem,4);
ey=zeros(nelem,4);
Edof=zeros(nelem,5);
for i=1:nelxy
    for j=1:nelxy
        % ------------------------------------------------------------
        % Elemental node coordinates:
        % ------------------------------------------------------------
        ex((i-1)*nelxy+j,1)=Coord((i-1)*nelxy+j,1);
        ex((i-1)*nelxy+j,2)=Coord((i-1)*nelxy+j+1,1);
        ex((i-1)*nelxy+j,3)=Coord(i*nelxy+j+1,1);
        ex((i-1)*nelxy+j,4)=Coord(i*nelxy+j,1);
        
        ey((i-1)*nelxy+j,1)=Coord((i-1)*nelxy+j,2);
        ey((i-1)*nelxy+j,2)=Coord((i-1)*nelxy+j+1,2);
        ey((i-1)*nelxy+j,3)=Coord(i*nelxy+j+1,2);
        ey((i-1)*nelxy+j,4)=Coord(i*nelxy+j,2);
        
        % ------------------------------------------------------------
        % Topology matrix:
        % ------------------------------------------------------------
        Edof((i-1)*nelxy+j,1)=(i-1)*nelxy+j;
        Edof((i-1)*nelxy+j,2)=((i-1)*npxy+j);
        Edof((i-1)*nelxy+j,3)=((i-1)*npxy+j+1);
        Edof((i-1)*nelxy+j,4)=((i-1)*npxy+j+1+npxy);
        Edof((i-1)*nelxy+j,5)=((i-1)*npxy+j+npxy);
    end
end
   
% --------------------------------------------------------------------
% FE model
% --------------------------------------------------------------------

% ---------------------------------------------------------------------
% This section is to recollect the values of the sample data points
% for each range of each finite element in the discrete model:
% ---------------------------------------------------------------------
ndof=length(Coord(:,1));

ya=rangexy(1);
sumZk=zeros(nelem,1);
npointsk=zeros(nelem,1);
for i=1:nelxy
    xa=rangexy(1);
    for j=1:nelxy
        exel=ex((i-1)*nelxy+j,:);
        eyel=ey((i-1)*nelxy+j,:);
        
        npoints_elem=0;
        for k=1:samplePoints
            if (x(k)>=xa && x(k)<xa+dx) && ...
                    (y(k)>=ya && y(k)<ya+dy)
                npoints_elem=npoints_elem+1;
                sumZk((i-1)*nelxy+j)=sumZk((i-1)*nelxy+j)+z(k);
            end
        end
        npointsk((i-1)*nelxy+j)=npoints_elem;
        if npoints_elem~=0
            sumZk((i-1)*nelxy+j)=sumZk((i-1)*nelxy+j)/...
                npointsk((i-1)*nelxy+j);
        end
        xa=xa+dx;
    end
    ya=ya+dy;
end
ngp=2;

lambdax=2;
lambday=2; 
ep=[ngp,lambdax,lambday];
eq=sumZk;

nElements=length(Edof(:,1));

% ---------------------------------------------------------------------
% Generating system of equations
% ---------------------------------------------------------------------
Kxy=zeros(ndof);
f=zeros(ndof,1);
for i=1:nElements
    [ Ke, fe ] = plan4bilinKelFbelInterp(ex(i,:),ey(i,:),ep,eq(i));
    [Kxy,f]=assem(Edof,Kxy,Ke,f,fe);
end
K=Kxy;

% ----------------------------------------------------------------------
% Original surface
% ----------------------------------------------------------------------
limit=npxy;
delta=(rangexy(2)-rangexy(1))/(limit-1);
[xValues,yValues]=meshgrid(rangexy(1):delta:rangexy(2),...
    rangexy(1):delta:rangexy(2));

zValues=zeros(limit,limit);
for j=1:limit
    for k=1:limit
        zValues(j,k)=fxy(xValues(j,k), yValues(j,k));
    end
end

% ---------------------------------------------------------------------
% Solving the system
% ---------------------------------------------------------------------

% Boundary conditions -------------------------------------------------
bounds=[1 2 3 4];
if length(bounds)==4
    bc=zeros(npxy*2+2*(npxy-2),2);
    % Bottom boundary:
    for i=1:npxy
        bc(i,1)=i;
        bc(i,2)=zValues(1,i);
    end
    % Left boundary:
    for i=1:npxy-2
        bc(npxy+i,1)=(i)*npxy+1;
        bc(npxy+i,2)=zValues(i+1,1);
    end

    % Right boundary:
    for i=1:npxy-2
        bc(npxy+npxy-2+i,1)=(i+1)*npxy;
        bc(npxy+npxy-2+i,2)=zValues(i+1,npxy);
    end

    % Upper boundary:
    for i=1:npxy
        bc(npxy*2+2*(npxy-2)-i+1,1)=npxy^2-i+1;
        bc(npxy*2+2*(npxy-2)-i+1,2)=zValues(npxy,npxy-i+1);
    end
    [u,r]=solveq(K,f,bc);
end
% ----------------------------------------------------------------------
% Plotting results - Interpolating surface
% ----------------------------------------------------------------------
figure(1)
surf(xValues,yValues,zValues);
colormap default;
shading interp;
view([-7 -9 10]);
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')
title('Original surface')
colorbar;
uValues=zeros(limit,limit);
for i=1:limit
    for j=1:limit
        uValues(i,j)=u((i-1)*limit+j,1);
    end
end

figure(2);
mesh(xValues,yValues,uValues);
colormap default;
shading interp;
view([-7 -9 10]);
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')
title('Mesh for the interpolating surface')
hold on

figure(3);
surf(xValues,yValues,uValues);
colormap default;
shading interp;
colorbar;
view([-7 -9 10]);
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')
title('Interpolating surface')
hold on


% ----------------------------- FIN ----------------------------------
%