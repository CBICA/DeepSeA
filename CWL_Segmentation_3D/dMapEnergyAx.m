function [u] = dMapEnergyAx(u,diffMap,lambdaAdjuster,visualize)
%DMAPTV 此处显示有关此函数的摘要
%   CAUTION: This function assumes ISOTROPIC input.

intFlag_mu = 3;

% pad matrices to alliviate the boundary effect
u = padarray(u, [1 1], 'replicate');
diffMap = padarray(diffMap, [0 1 1], 'replicate');

[row,col,slc] = size(diffMap);

[~,dM,~] = gradient(diffMap);

[slc_sub,col_sub] = meshgrid(1:slc,1:col);

lambda = 10^lambdaAdjuster;
% % determine lambda value adaptively based on ratio between initial external
% % and internal energy
% ind = sub2ind([row,col,slc],round(u(:)),col_sub(:),slc_sub(:));
% E_ext0 = sum(diffMap(ind));
% [ux,uy]=gradient(u);
% du_square_norm = ux.^2 + uy.^2;
% E_int0 = sum(du_square_norm(:))/2;
% % make the internal energy term N magnitude(s) higher than the external 
% % energy term in the begining
% % mag_ext = floor(log10(E_ext0));
% % mag_int = floor(log10(E_int0));
% % lambda = 10^(lambdaAdjuster+mag_ext-mag_int); % weight for smoothing term
% mag_ext = log10(E_ext0);
% mag_int = log10(E_int0);
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
    ind = sub2ind([row,col,slc],round(u(:)),col_sub(:),slc_sub(:));
    E_ext(iterCnt) = sum(diffMap(ind));
    par_der(:) = dM(ind);
    
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
    u(u>row) = row; 
    
    if max(abs(delta_u(:)))<1e-3
        disp(['CONVERGED AT ITER ' num2str(iterCnt)]);
        flag_converge = true;
        break;
    end
end
if ~flag_converge
    warning('FAILED TO CONVERGE!!!');
end

% strip the padded frame
u = u(2:end-1,2:end-1);

% visual debug codes
if exist('visualize','var') && visualize
    figure;
    imshow(mat2gray(row-u));
    figure;
    subplot(3,1,1); plot(E_seg(1:iterCnt)); title('E(D)');
    subplot(3,1,2); plot(E_ext(1:iterCnt)); title('E_{ext}');
    subplot(3,1,3); plot(E_int(1:iterCnt)); title('E_{int}');
end

end

function K = curvature_central(nx,ny)
    [nxx,~]=gradient(nx);
    [~,nyy]=gradient(ny);
    K=nxx+nyy;
end

