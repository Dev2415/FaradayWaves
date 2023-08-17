clear; close all; clc;

% diff = 90/2;
% theta = 0 : diff : 360 - diff;
% 
% kx = cosd(theta);
% ky = sind(theta);
% 
% K = [kx', ky'];
% 
% Nx = 100;
% Ny = 100;
% 
% p = zeros(Ny, Nx);
% for i = 1 : Ny
%     for j = 1 : Nx
%         for kk = 1 : length(K)
%             p(i,j) = p(i,j) + cos(j*K(kk,1) + i*K(kk,2));
%         end
%     end
% end
% 
% contourf(p);

% diff = 90/3;
% theta = 0 : diff : 360 - diff;
% kc = 7.76e2;
% Lx = 60/kc;
% Ly = 60/kc;
% kx = kc*cosd(theta);
% ky = kc*sind(theta);
% K = [kx', ky'];
 Nx = 101;
 Ny = 101;
% lx = linspace(0,Lx, Nx);
% ly = linspace(0,Ly,Ny);
% 
% 
% p = zeros(Ny, Nx);
% for i = 1 : Ny
%     for j = 1 : Nx
%         for k = 1 : length(K)
%             x = lx(j);
%             y = ly(i);
%             for kk = 1 : length(K)
%                 p(i,j) = p(i,j) + cos(x*K(kk,1) + y*K(kk,2));
%             end
%         end
%     end
% end
% contourf(p);
Lx = 1.5*4.27e-2;
Ly = 1.5*4.27e-2;
lx = linspace(0,Lx,Nx);
ly = linspace(0,Ly,Ny);
kc = 3/4.27e-2;
K = kc*[1, 0; -1/2, sqrt(3)/2; -1/2, -sqrt(3)/2];
K1 = [K(:,1)*cosd(22) - K(:,2)*sind(22), K(:,1)*sind(22) + K(:,2)*cosd(22)];
p = zeros(Ny, Nx);
for i = 1 : Ny
    for j = 1 : Nx
        x = lx(j); y = ly(i);
        for kk = 1 : 3
            p(i,j) = p(i,j) + 0.1*cos(K(kk,1)*x/Lx + K(kk,2)*y/Ly) + 0.1*cos(K1(kk,1)*x/Lx + K1(kk,2)*y/Ly);
        end
    end
end
Lz = 1e-2;
Nz = 21;
lz = linspace(0,Lz, Nz);
PHI = zeros(Ny, Nx, Nz);
ep = 1e-3;
h = Lz/5;
for i = 1 : Ny
    for j = 1 : Nx
        for k = 1 : Nz
           x = lx(j); y = ly(i); z = lz(k);
           for kk = 1 : 3
            PHI(i,j,k) = PHI(i,j,k) + h*(1 + ep*cos(K(kk,1)*x/Lx + K(kk,2)*y/Ly) + ep*cos(K1(kk,1)*x/Lx + K1(kk,2)*y/Ly));
           end
           PHI(i,j,k) = PHI(i,j,k) - z;
        end
    end
end
isosurface(PHI);

