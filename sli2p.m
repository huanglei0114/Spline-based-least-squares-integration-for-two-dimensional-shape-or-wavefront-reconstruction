function zsli2p = sli2p(sx, sy, x, y)
%SLI2P Spline-based Least-squares integration. (Polynomial Version)
%   D * Z = G (G is composed by spline or polynomial estimated values).

%   Copyright since 2016, Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2016-09-29 Original Version

% Check the number of arguments............................................
% Validate number of input arguments.
narginchk(4,4);
% Validate number of output arguments.
nargoutchk(1,1); 

% Generate Matrix D and G..................................................
% Calculate size and mask.
[Ny, Nx] = size(sx);
ValidMask = isfinite(sx) & isfinite(sy);

% Expand sy.
sy = [sy; NaN(1,Nx)];

% Compose matrix.
ee = ones(Ny*Nx,1);
Dx = spdiags([-ee,ee],[0,Ny],Ny*(Nx-1),Ny*Nx); 
Dy = spdiags([-ee,ee],[0,1],Ny*Nx,Ny*Nx);
mx = isnan(sx(:,1:end-1)+sx(:,2:end));
my = isnan(sy(1:end-1,:)+sy(2:end,:));

% Compose D.
Dx(mx(:),:)=[];
Dy(my(:),:)=[];
D = [Dx;Dy];
clear Dx Dy mx my;

% Compose matrix G.
SpGx = ComposeSpGx(x,sx,ValidMask,Nx,Ny);
SpGy = ComposeSpGy(y,sy,ValidMask,Nx,Ny);
clear sx sy x y;

% G
G = [SpGx(:);SpGy(:)];
G = G(isfinite(G));
clear Gx Gy;

% Solve "Warning: Rank deficient" for complete data by assuming Z(Ind)=0.  
Ind = find(D(1,:)==-1,1);
D(:,Ind) = [];   
Z = D\G; 
Z = [Z(1:Ind-1);0;Z(Ind:end)];

% Compose 2D matrix of z.
zsli2p = reshape(Z,Ny,Nx);
zsli2p(~ValidMask)= NaN;

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
    [Ns, Indices] = CheckSection(ml, Nx);
    
    % Spline fitting section by section.
    gs = cell(Ns,1);
    for ns = 1:Ns
        idx = Indices{ns};
        xx = xl(idx);
        vv = vl(idx);
        if length(xx)>1
            pp = spline(xx,vv); % "not-a-knot end condition"
            c = pp.coefs;   
            switch(size(c,2))
                case 4  % 4 points for piecewise cubic spline fitting.
                    dx = diff(xx);
                    if sign(mean(dx))==1
                        gs{ns} = dx.*(c(:,4) + dx.*(c(:,3)./2 + dx.*(c(:,2)./3 + dx.*c(:,1)./4)));
                    else
                        dx = -flipud(dx);
                        gs{ns} = dx.*(c(:,4) + dx.*(c(:,3)./2 + dx.*(c(:,2)./3 + dx.*c(:,1)./4)));
                        gs{ns} = -flipud(gs{ns});
                    end
                    
                case 3  % 3 points for 2nd order polynominal fitting.
                    xf = xx-xx(1);
                    ff = c(:,3) + xf.*(c(:,2)./2 + xf.*c(:,1)./3);
                    g = ff(2:end)-ff(1:end-1)+vv(1:end-1);
                    gs{ns} = g(:).*diff(xx);
                
                case 2  % 2 points for 1st order polynominal fitting.
                    dx = diff(xx);
                    gs{ns} = dx.*(c(:,2) + dx.*c(:,1)./2);
                
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
    [Ns, Indices] = CheckSection(ml, Ny);
    
    % Spline fitting section by section.
    gs = cell(Ns,1);
    for ns = 1:Ns
        idx = Indices{ns};
        yy = yl(idx);
        vv = vl(idx);
        if length(yy)>1
            pp = spline(yy,vv); % "not-a-knot end condition"
            c = pp.coefs;   
            switch(size(c,2))
                case 4  % 4 points for piecewise cubic spline fitting.
                    dy = diff(yy);
                    if sign(mean(dy))==1
                        gs{ns} = dy.*(c(:,4) + dy.*(c(:,3)./2 + dy.*(c(:,2)./3 + dy.*c(:,1)./4)));
                    else
                        dy = -flipud(diff(yy));
                        gs{ns} = dy.*(c(:,4) + dy.*(c(:,3)./2 + dy.*(c(:,2)./3 + dy.*c(:,1)./4)));
                        gs{ns} = -flipud(gs{ns});
                    end
                    
                case 3  % 3 points for 2nd order polynominal fitting.
                    yf = yy-yy(1);
                    ff = c(:,3) + yf.*(c(:,2)./2 + yf.*c(:,1)./3);
                    g = ff(2:end)-ff(1:end-1)+vv(1:end-1);
                    gs{ns} = g(:).*diff(yy);

                case 2  % 2 points for 1st order polynominal fitting.
                    dy = diff(yy);
                    gs{ns} = dy.*(c(:,2) + dy.*c(:,1)./2);
                
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


% Check Sections.
function [Ns, Indices] = CheckSection(ml, N)
if all(ml)==true      
    Ns = 1;
    Indices{Ns} = 1:N;
else
    Indices = cell(N,1);
    first = nan;
    last = nan;
    Ns = 0;
    for n = 1:N
        % Find the first.
        if n==1
            if ml(n)==true
                first = n;
            end
        else
            if ml(n)==true && ml(n-1)==false
                first = n;
            end
        end

        % Find the last.
        if n==N
            if ml(n)==true
                last = n;
            end
        else
            if ml(n)==true && ml(n+1)==false
                last = n;
            end
        end

        % Sum up the total number of sections and compose the Indices.
        if isfinite(first) && isfinite(last)
            Ns = Ns + 1;
            Indices{Ns} = first:last;
            first = nan;
            last = nan;
        end
    end
end
end
