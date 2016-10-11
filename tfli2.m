function ztfli2  = tfli2(sx,sy,x,y,z)
%TFLI2 Traditional Finite-difference-based Least-squares Integration.
%   D * Z = G.
% 
%   Reference:
%   W.H. Southwell, "Wave-front estimation from wave-front slope 
%   measurements," J. Opt. Soc. Am. 70, 998-1006 (1980) 


%   Copyright since 2010, Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2010-10-01 Original Version
%   2014-02-01 Revision on data structure


% Validate number of input arguments.
narginchk(4,5);
% Validate number of output arguments.
nargoutchk(1,1); 

% Calculate size and mask.
[ny, nx] = size(sx);
mask = isnan(sx) | isnan(sy);

% Expand sy and y.
sy = [sy;NaN(1,nx)];
y  = [y ;NaN(1,nx)];

% Make matrix.
ee = ones(ny*nx,1);
Dx = spdiags([-ee,ee],[0,ny],ny*(nx-1),ny*nx); 
Dy = spdiags([-ee,ee],[0,1],ny*nx,ny*nx);
Gx = (sx(:,1:end-1)+sx(:,2:end)).*(x(:,2:end)-x(:,1:end-1))/2;  
Gy = (sy(1:end-1,:)+sy(2:end,:)).*(y(2:end,:)-y(1:end-1,:))/2; 
if nargin==5
    Dz = spdiags(ee,0,ny*nx,ny*nx);
end
clear sx sy x y;

% delete NaN points
if nargin==4
    % D
    Dx(isnan(Gx(:)),:)=[];
    Dy(isnan(Gy(:)),:)=[];
    D = [Dx;Dy];
    clear Dx Dy;
    % G
    Gx(isnan(Gx(:)))=[];
    Gy(isnan(Gy(:)))=[];
    G = [Gx(:);Gy(:)];
    clear Gx Gy;
    
    % Solve "Warning: Rank deficient" for complete dataset by assuming Z(Ind) is 0.  
    Ind = find(D(1,:)==-1,1);
    D(:,Ind) = [];   
    Z = D\G; 
    Z = [Z(1:Ind-1);0;Z(Ind:end)];      

elseif nargin==5
    % D
    Dx(isnan(Gx(:)),:)=[];
    Dy(isnan(Gy(:)),:)=[];
    Dz(isnan(z (:)),:)=[];
    D = [Dx;Dy;Dz];
    clear Dx Dy Dz;
    % G
    Gx(isnan(Gx(:)))=[];
    Gy(isnan(Gy(:)))=[];
    z(isnan(z))=[];
    G = [Gx(:);Gy(:);z(:)];
    clear Gx Gy z;
    
    % Calculate Z with least squares method.
    Z = D\G; 

end

clear D G;

% zls
ztfli2 = reshape(Z,ny,nx);
clear Z;
ztfli2(mask)= NaN;

% End of Function.
end

