clear all
close all
clc
B = 5;
n=20;                   %% number of solution asked 
Mass = 0.067;           %% effective mass, constant over all the structure...
Nx=2;                  %% Meshing point in x-direction
Ny=3;                  %% Meshing point in y-direction
Nz=4;                  %% Meshing point in x-direction
Mx=1;              
My=1;     
Mz=1;   
x=linspace(-Mx/2,Mx/2,Nx);
y=linspace(-My/2,My/2,Ny);
z=linspace(-Mz/2,Mz/2,Nz);
Nx=length(x);
Ny=length(y);
Nz=length(z);
dx=x(2)-x(1);
dy=y(2)-y(1);
dz=z(2)-z(1);
Y=kron(kron(diag(ones(1,Nx)), diag(y)), diag(ones(1,Nz)));
A1 = ones(1,Nx*Ny*Nz-Ny*Nz);
DX = (-1)*diag(A1,-Ny*Nz) + (1)*diag(A1,Ny*Nz);
DX2 = (-2)*diag(ones(1,Nx*Ny*Nz))+(1)*diag(A1,-Ny*Nz) + (1)*diag(A1,Ny*Nz);
B1 = ones(1,Ny*Nz-Nz);
B2 = (-2)*diag(ones(1,Nz*Ny))+(1)*diag(B1,-Nz) + (1)*diag(B1,Nz);
B3 = diag(ones(1,Nx));
DY2 = kron(B3, B2);
C1 = ones(1,Nz-1);
C2 = (-2)*diag(ones(1,Nz))+(1)*diag(C1,-1) + (1)*diag(C1,1);
C3 = diag(ones(1,Nx*Ny));
DZ2 = kron(C3, C2);
% H=(-1/(2*Mass)) * ( DX2/dx^2 + B.^2*Y.^2 - 2*B*Y*DX/(2*dx) + DY2/dy^2 + DZ2/dz^2 );
H=(-1/(2*Mass)) * ( DX2/dx^2 + B.^2*Y.^2 - 2*B*Y*DX/(2*dx) + DY2/dy^2 );
% H=sparse(H);
[PSI,Energy] = eigs(H,n,'SM');
E1 = diag(Energy);
