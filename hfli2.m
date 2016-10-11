function zhfli2 = hfli2(sx,sy,x,y,z)
%HFLI2 Higher-order Finite-difference-based Least-squares Integration.
%   D * Z = G.
%
%   Reference: 
%   Guanghui Li, Yanqiu Li, Ke Liu, Xu Ma, and Hai Wang, "Improving wavefront
%   reconstruction accuracy by using integration equations with higher-order 
%   truncation errors in the Southwell geometry," J. Opt. Soc. Am. A 30, 
%   1448-1459 (2013) 

%   Copyright since 2013 Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2013-12-01 Original Version
%   2014-01-31 First revision on matrix
%   2014-11-13 revision on speed
%   2016-01-29 Solve "Warning: Rank deficient" for complete dataset.

% Validate number of input arguments.
narginchk(4,5);
% Validate number of output arguments.
nargoutchk(1,1); 

% Calculate size and mask.
[Ny, Nx] = size(sx);
mask = isnan(sx) | isnan(sy);

% Expand sx, sy, x, and y.
sx = [NaN(Ny,1),sx,NaN(Ny,2)];
x  = [NaN(Ny,1),x ,NaN(Ny,2)];
sy = [NaN(1,Nx);sy;NaN(2,Nx)];
y  = [NaN(1,Nx);y ;NaN(2,Nx)];

ExpandMaskx = isnan(sx);
se = [1 1 0 1 0];
DilatedExpandMaskx = imdilate(ExpandMaskx,se);
Maskx = DilatedExpandMaskx(:,2:end-2) & ~ExpandMaskx(:,2:end-2);

ExpandMasky = isnan(sy);
se = [1;1;0;1;0];
DilatedExpandMasky = imdilate(ExpandMasky,se);
Masky = DilatedExpandMasky(2:end-2,:) & ~ExpandMasky(2:end-2,:);

% Generate matrix.
Num = Ny*Nx;
ee = ones(Num,1);
Dx = spdiags([-ee,ee],[0,Ny],Num,Num);
Dy = spdiags([-ee,ee],[0,1],Num,Num);

% O(h^5)
Gx = (-1/13*sx(:,1:end-3)+sx(:,2:end-2)+sx(:,3:end-1)-1/13*sx(:,4:end)) ...
    .*(x(:,3:end-1)-x(:,2:end-2))*13/24;
Gy = (-1/13*sy(1:end-3,:)+sy(2:end-2,:)+sy(3:end-1,:)-1/13*sy(4:end,:)) ...
    .*(y(3:end-1,:)-y(2:end-2,:))*13/24;

% O(h^3)
Gx3 = (sx(:,2:end-2)+sx(:,3:end-1)).*(x(:,3:end-1)-x(:,2:end-2))/2;
Gy3 = (sy(2:end-2,:)+sy(3:end-1,:)).*(y(3:end-1,:)-y(2:end-2,:))/2;
Gx(Maskx) = Gx3(Maskx);
Gy(Masky) = Gy3(Masky);

if nargin==5
    Dz = spdiags(ee,0,Ny*Nx,Ny*Nx);
end
clear sx sy x y Gx3 Gy3;

% Delete NaN pixels.
if nargin==4
    % D
    D = [Dx(isfinite(Gx(:)),:);Dy(isfinite(Gy(:)),:)];
    clear Dx Dy;
    % G
    G = [Gx(isfinite(Gx(:)));Gy(isfinite(Gy(:)))];
%     clear Gx Gy;

    % Calculate Z with least squares method.
    % Z=(D'*D)\D'*G;

    % Solve "Warning: Rank deficient" for complete dataset by assuming Z(Ind) is 0.  
    Ind = find(D(1,:)==-1,1);
    D(:,Ind) = [];   
    Z = D\G; 
    Z = [Z(1:Ind-1);0;Z(Ind:end)];      
    
elseif nargin==5
    % D
    D = [Dx(isfinite(Gx(:)),:);Dy(isfinite(Gy(:)),:);Dz(isfinite(z (:)),:)];
    clear Dx Dy Dz;
    % G
    G = [Gx(isfinite(Gx(:)));Gy(isfinite(Gy(:)));z(isfinite(z))];
    clear Dx Dy z;

    % Calculate Z with least squares method.
    % Z=(D'*D)\D'*G;
    Z = D\G;
    
end

clear D G;

% Compose 2D matrix of z.
zhfli2 = reshape(Z,Ny,Nx);
zhfli2(mask)= NaN;

end