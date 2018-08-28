function [u] = dMapEnergy(u,r,lambdaAdjuster,visualize)
%[u] = dMapEnergy(u,respMap,lambdaAdjuster,visualize)
% This function assumes ISOTROPIC inputs 
%   r: some kind of filter response map, e.g., difference map
%   lambda: weight for smoothness term

intFlag_mu = 3;

% pad matrices to alliviate the boundary effect
u = padarray(u, [1 1], 'replicate');
r = padarray(r, [1 0 1], 'replicate');

[row,col,slc] = size(r);

[slc_sub,row_sub] = meshgrid(1:slc,1:row);

dr = gradient(r); % r for response

lambda = 10^lambdaAdjuster;
% % Initial external and internal energy and their magnitudes
% ind = sub2ind([row,col,slc],row_sub(:),round(u(:)),slc_sub(:));
% E_ext0 = sum(r(ind));
% mag_ext = log10(E_ext0);
% [ux,uy]=gradient(u);
% du_square_norm = ux.^2 + uy.^2;
% E_int0 = sum(du_square_norm(:))/2;
% mag_int = log10(E_int0);
% % make the internal energy term N magnitude(s) higher than the external 
% % energy term in the begining
% lambda = 10^ceil(lambdaAdjuster+mag_ext-mag_int); % weight for smoothing term

maxIter = 10000*intFlag_mu;
E_ext = zeros(1,maxIter);
E_int = zeros(1,maxIter);
E_seg = zeros(1,maxIter);
dE = zeros(1,maxIter-1);
par_der = zeros(size(u));
mu = 1;
flag_converge = false;
for iterCnt=1:maxIter   
    % external energy
    ind = sub2ind([row,col,slc],row_sub(:),round(u(:)),slc_sub(:));
    E_ext(iterCnt) = sum(r(ind));
    par_der(:) = dr(ind);
    
    % internal energy
    u=NeumannBoundCond(u);
    [ux,uy]=gradient(u); 
    du_square_norm = ux.^2 + uy.^2;
    E_int(iterCnt) = sum(du_square_norm(:))/2;
%     normDu=sqrt(ux.^2 + uy.^2 + 1e-10);
%     Nx=ux./normDu;
%     Ny=uy./normDu;
%     div=curvature_central(Nx,Ny);
    lap=4*del2(u);
    
    E_seg(iterCnt) = E_ext(iterCnt) + lambda*E_int(iterCnt);
    % adaptive update step mu
    if iterCnt>1
        dE(iterCnt-1)=E_int(iterCnt)-E_int(iterCnt-1);
    end
    if iterCnt>intFlag_mu
        divideMu = true;
        for j = 1:intFlag_mu
            if dE(iterCnt-j)<=0
                divideMu = false;
                break;
            end
        end
    else
        divideMu = false;
    end
    if divideMu
        fprintf('%d consecutive energy increases at Iter %d. Dividing mu...\n', ...
            intFlag_mu, iterCnt);
        mu = mu/2;
    end
    
    delta_u = mu*(lambda*lap-par_der);
    u = u + delta_u;
    % ad-hoc fix for out-of-bound error
    u(u<1) = 1;
    u(u>col) = col;
    
    if max(abs(delta_u(:)))<1e-3
        disp(['CONVERGED AT ITER ' num2str(iterCnt)]);
        flag_converge = true;
        break;
    end
end
if ~flag_converge
    warning('FAILED TO CONVERGE!!!');
end

% visual debug codes
if exist('visualize', 'var') && visualize
    figure; imshow(mat2gray(col-u));
    figure;
    subplot(3,1,1); plot(E_seg(1:iterCnt)); title('E_{seg}');
    subplot(3,1,2); plot(E_ext(1:iterCnt)); title('E_{ext}');
    subplot(3,1,3); plot(E_int(1:iterCnt)); title('E_{int}');
end

% strip the padded frame
u = u(2:end-1,2:end-1);

end

function K = curvature_central(nx,ny)
    [nxx,~]=gradient(nx);
    [~,nyy]=gradient(ny);
    K=nxx+nyy;
end

