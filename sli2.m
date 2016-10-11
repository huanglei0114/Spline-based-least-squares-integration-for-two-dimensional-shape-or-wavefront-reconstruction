function zsli2 = sli2(sx, sy, x, y)
%SLI2 Spline-based Least-squares integration.
%   D * Z = G (G is mainly composed by spline estimated values).

%   Copyright since 2016, Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%
%   2016-09-29 Original Version


% Validate number of input arguments.
narginchk(4,4);
% Validate number of output arguments.
nargoutchk(1,1);


% Calculate size and mask..................................................
[Ny, Nx] = size(sx);
ValidMask = isfinite(sx) & isfinite(sy);

% Compose matrix D.........................................................
% Expand sy and y.
sy = [sy;NaN(1,Nx)];
y  = [y ;NaN(1,Nx)];

% Compose matrix.
ee = ones(Ny*Nx,1);
Dx = spdiags([-ee,ee],[0,Ny],Ny*(Nx-1),Ny*Nx); 
Dy = spdiags([-ee,ee],[0,1],Ny*Nx,Ny*Nx);
gx = (sx(:,1:end-1)+sx(:,2:end)).*(x(:,2:end)-x(:,1:end-1))/2;  
gy = (sy(1:end-1,:)+sy(2:end,:)).*(y(2:end,:)-y(1:end-1,:))/2; 

% Delete NaN points.
% D
Dx(isnan(gx(:)),:)=[];
Dy(isnan(gy(:)),:)=[];
D = [Dx;Dy];
clear Dx Dy;


% Compose matrix g.........................................................

% Compose matrix SpGx.
spGx = ComposeSpGx(x,sx,ValidMask,Nx,Ny);
% Compose matrix SpGy.
spGy = ComposeSpGy(y,sy,ValidMask,Nx,Ny);
clear sx sy x y;

% Replace with spline values, if available.
gy(end,:)=[];
gx(isfinite(spGx)) = spGx(isfinite(spGx));
gy(isfinite(spGy)) = spGy(isfinite(spGy));

% g
g = [gx(:);gy(:)];
g = g(isfinite(g));
clear gx gy;

% Solve "Warning: Rank deficient" for complete data by assuming Z(Ind)=0.  
Ind = find(D(1,:)==-1,1);
D(:,Ind) = [];
Z = D\g;
Z = [Z(1:Ind-1);0;Z(Ind:end)];

% Compose 2D matrix of z.
zsli2 = reshape(Z,Ny,Nx);
clear Z;
zsli2(~ValidMask)= NaN;

end




%% Subfunctions.

% Compose matrix spGx
function SpGx = ComposeSpGx(x,sx,ValidMask,Nx,Ny)

SpGx = NaN(Ny,Nx-1);
for ny = 1:Ny
    xl = x(ny,:)';
    vl = sx(ny,:)';
    
    % Check the number of sections.
    ml = ValidMask(ny,:)';   
    if all(ml)==true      
        Ns = 1;
    else
        ss = regionprops(ml,'PixelIdxList');
        Ns = length(ss);
    end
    
    % Spline fitting section by section.
    gs = cell(Ns,1);
    for ns = 1:Ns
        if all(ml)==true
            idx = 1:Nx;
        else
            idx = ceil(ss(ns).PixelIdxList);
        end
        xx = xl(idx);
        vv = vl(idx);
        if length(xx)>1
            pp = spline(xx,vv); % "not-a-knot end condition"
            c = pp.coefs;   
            switch(size(c,2))
                case 4  % 4 points for piecewise cubic spline fitting.
                    dx = diff(xx);
                    gs{ns} = dx.*(c(:,4) + dx.*(c(:,3)./2 + dx.*(c(:,2)./3 + dx.*c(:,1)./4)));
                
                case 3  % 3 points for 2nd order polynominal fitting.
                    % Here we do not use the polynomials. 
                    % We are going to use the Southwell expression instead
                    gs{ns} = diff(xx)*NaN; 
                
                case 2  % 2 points for 1st order polynominal fitting.
                    % Here we do not use the polynomials. 
                    % We are going to use the Southwell expression instead
                    gs{ns} = diff(xx)*NaN;
                
                case 1
                    % Logically impossible.
                    error('Only one point for fitting in x direction!');
                
                otherwise
                    % Logically impossible.
                    error('Unexpected number of points for fitting in x direction!');
            end
        end
    end
    sg = cat(1,gs{:});
    Valid = ml(1:end-1) & ml(1+1:end);
    pt = 1;
    for nx = 1 : Nx-1
        if Valid(nx) == 1
            SpGx(ny, nx) = sg(pt);
            pt = pt + 1;
        end
    end
end
end


% Compose matrix spGy
function SpGy = ComposeSpGy(y,sy,ValidMask,Nx,Ny)
SpGy = NaN(Ny-1,Nx);
for nx = 1:Nx
    yl = y(:,nx);
    vl = sy(:,nx);
    
    % Check the number of sections.
    ml = ValidMask(:,nx);
    if all(ml)==true      
        Ns = 1;
    else
        ss = regionprops(ml,'PixelIdxList');
        Ns = length(ss);
    end
    
    % Spline fitting section by section.
    gs = cell(Ns,1);
    for ns = 1:Ns
        if all(ml)==true
            idx = 1:Nx;
        else
            idx = ceil(ss(ns).PixelIdxList);
        end
        yy = yl(idx);
        vv = vl(idx);
        if length(yy)>1
            pp = spline(yy,vv); % "not-a-knot end condition"
            c = pp.coefs;   
            switch(size(c,2))
                case 4  % 4 points for piecewise cubic spline fitting.
                    dy = diff(yy);
                    gs{ns} = dy.*(c(:,4) + dy.*(c(:,3)./2 + dy.*(c(:,2)./3 + dy.*c(:,1)./4)));
                    
                case 3  % 3 points
                    % Here we do not use the polynomials. 
                    % We are going to use the Southwell expression instead
                    gs{ns} = diff(yy)*NaN;

                case 2  % 2 points
                    % Here we do not use the polynomials. 
                    % We are going to use the Southwell expression instead
                    gs{ns} = diff(yy)*NaN;
                
                case 1
                    % Logically impossible.
                    error('Only one point for fitting in y direction!');
                
                otherwise
                    % Logically impossible.
                    error('Unexpected number of points for fitting in y direction!');
            end
        end
    end
    sg = cat(1,gs{:});
    Valid = ml(1:end-1) & ml(1+1:end);
    pt = 1;
    for ny = 1 : Ny-1
        if Valid(ny) == 1
            SpGy(ny, nx) = sg(pt);
            pt = pt + 1;
        end
    end
end
end