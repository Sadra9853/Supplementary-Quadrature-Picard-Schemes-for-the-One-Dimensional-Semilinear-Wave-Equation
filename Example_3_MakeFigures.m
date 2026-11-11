%% Example_3_MakeFigures.m
% Example 3: sine-Gordon BREATHER, both supplementary families.
%   Newton-Cotes (order kNC)  and  Gauss-Legendre (kGL points per block),
%   on the same list of mesh sizes m, h = T/m (m = number of time panels
%   for NC, number of time blocks for GL).
%   Output: console table, ex3_table.tex (rows for tab:ex3) and
%   ex3_convergence_error.pdf (left: convergence of both families,
%   right: pointwise error at t=T for the finest mesh).
%
% Exact solution (omega in (0,1), s = sqrt(1-omega^2)):
%   u(t,x) = 4*atan( s*cos(omega*t) / (omega*cosh(s*x)) ),  u_t(0,x)=0,
% so that the d'Alembert linear part is L(t,x) = (f(x+t)+f(x-t))/2.
%
% DOMAIN NOTE. The breather decays only like exp(-s|x|) (s=0.866 for
% omega=0.5), much more slowly than the kink's source term would suggest
% in a small window, so the computational domain must extend far beyond
% the reported window: padExtra is large on purpose (Xmax = Xcore+T+padExtra).
% If you reduce padExtra, check that e_m does not change.

clear; clc; close all;

omega = 0.5;  s = sqrt(1-omega^2);
Ffun      = @(u) sin(u);
uExactFun = @(t,x) 4*atan( s*cos(omega*t)./(omega*cosh(s*x)) );
f0   = @(x) uExactFun(0,x);
Lfun = @(t,x) 0.5*( f0(x+t) + f0(x-t) );

T = 1;  Xcore = 6;  padExtra = 23;  tol = 1e-12;  maxIter = 100;

kNC = 4;                  % Newton-Cotes order (Boole)
kGL = 5;                  % Gauss points per time block (k=5: matches Table 9 and Figure 9)
mList = [4 8 16 32];      % mesh sizes (multiples of kNC)
hxGL  = 0.5;              % spatial Gauss block width for the GL runs
kX    = 10;               % Gauss points per spatial block (GL)
nBlocksX = round(2*(Xcore+T+padExtra)/hxGL);

E = zeros(numel(mList),2);   % column 1: NC, column 2: GL
for i = 1:numel(mList)
    o1 = solveWaveNC(Ffun,Lfun,uExactFun,T,kNC,mList(i),Xcore,padExtra,tol,maxIter);
    E(i,1) = o1.e;
    o2 = solveWaveGL(Ffun,Lfun,uExactFun,T,mList(i),kGL,nBlocksX,kX,Xcore,padExtra,tol,maxIter,false);
    E(i,2) = o2.e;
    fprintf('m=%3d  NC(k=%d) e=%.3e   GL(k=%d) e=%.3e\n',mList(i),kNC,E(i,1),kGL,E(i,2));
end

Q = NaN(size(E));
for j = 1:2
    for i = 2:numel(mList)
        Q(i,j) = log2(E(i-1,j)/E(i,j));
    end
end

fprintf('\n%6s   %12s %8s   %12s %8s\n','m',sprintf('NC k=%d',kNC),'order',sprintf('GL k=%d',kGL),'order');
for i = 1:numel(mList)
    if i==1, q1='   --   '; q2='   --   '; else, q1=sprintf('%8.2f',Q(i,1)); q2=sprintf('%8.2f',Q(i,2)); end
    fprintf('%6d   %12.3e %s   %12.3e %s\n',mList(i),E(i,1),q1,E(i,2),q2);
end

fid = fopen('ex3_table.tex','w');
for fo = [1 fid]
    fprintf(fo,'%% ---- paste into tab:ex3 in wave_paper_1plus1.tex ----\n');
    for i = 1:numel(mList)
        if i==1, q1='--'; q2='--'; else, q1=sprintf('%.2f',Q(i,1)); q2=sprintf('%.2f',Q(i,2)); end
        fprintf(fo,'$m_%d=%d$ & %s & %s & %s & %s \\\\\n',i,mList(i),texSci(E(i,1)),q1,texSci(E(i,2)),q2);
    end
end
fclose(fid);

%% Finest meshes: pointwise error at t = T
mF = max(mList);
oN = solveWaveNC(Ffun,Lfun,uExactFun,T,kNC,mF,Xcore,padExtra,tol,maxIter);
xN = oN.x(oN.core);  errN = abs(oN.U(end,oN.core) - uExactFun(T,xN));  errN(errN==0) = eps;
oG = solveWaveGL(Ffun,Lfun,uExactFun,T,mF,kGL,nBlocksX,kX,Xcore,padExtra,tol,maxIter,true);
errG = abs(oG.uT - uExactFun(T,oG.xT));  errG(errG==0) = eps;

h = T./mList;
fig = figure('Color','w','Units','centimeters','Position',[2 2 14 6]);
subplot(1,2,1);
loglog(h,E(:,1),'o-','LineWidth',1.6,'MarkerSize',6); hold on;
loglog(h,E(:,2),'s-','LineWidth',1.6,'MarkerSize',6);
loglog(h,E(1,1)*(h/h(1)).^(kNC+2),'k--','LineWidth',1);
loglog(h,E(1,2)*(h/h(1)).^(kGL+1),'k:','LineWidth',1.2);
grid on; xlabel('h = T/m'); ylabel('e_m');
legend(sprintf('NC, k=%d',kNC),sprintf('GL, k=%d',kGL),sprintf('h^{%d}',kNC+2),sprintf('h^{%d}',kGL+1),'Location','southeast');
subplot(1,2,2);
semilogy(xN,errN,'b.','LineWidth',1.2); hold on;
semilogy(oG.xT,errG,'r.','MarkerSize',5);
grid on; xlabel('x'); ylabel('|u_h(T,x) - u(T,x)|');
legend('NC','GL','Location','south');
set(findall(fig,'-property','FontSize'),'FontSize',10);
exportgraphics(fig,'ex3_convergence_error.pdf','ContentType','vector');
exportgraphics(fig,'ex3_convergence_error.eps','ContentType','vector');
fprintf('Saved ex3_table.tex and ex3_convergence_error.pdf\n');

%% ------------------------------------------------------------------
function str = texSci(x)
    e = floor(log10(abs(x)));  mant = x/10^e;
    str = sprintf('$%.3f\\times10^{%d}$',mant,e);
end

function out = solveWaveNC(Ffun,Lfun,uExactFun,T,k,nt,Xcore,padExtra,tol,maxIter)
    tGrid = linspace(0,T,nt+1);
    h     = tGrid(2)-tGrid(1);
    pad   = T + padExtra;
    NxHalf = round((Xcore+pad)/h);
    xEval  = (-NxHalf:NxHalf)*h;
    Nx = numel(xEval);
    coreIdx = abs(xEval) <= Xcore + 1e-9;

    Wt = supplementaryNCVolterraWeights(tGrid,k);
    spaceKernel = cell(nt,1);
    for m = 1:nt
        spaceKernel{m} = fixedDomainNCKernel(m,k,h);
    end

    L = Lfun(tGrid.',xEval);
    U = L;
    for iter = 1:maxIter
        Unew = L;
        for r = 1:nt
            idx = find(Wt(r+1,:) ~= 0);
            corr = zeros(1,Nx);
            for ii = idx
                m = r - ii + 1;
                if m == 0, continue; end
                mm = abs(m);  sgn = sign(m);
                kernel = spaceKernel{mm};
                Fvals = Ffun(U(ii,:));
                Fpad  = [zeros(1,mm), Fvals, zeros(1,mm)];
                Sint  = zeros(1,Nx);
                for p = 0:2*mm
                    Sint = Sint + kernel(p+1)*Fpad(p+1:p+Nx);
                end
                corr = corr + Wt(r+1,ii)*sgn*Sint;
            end
            Unew(r+1,:) = L(r+1,:) - 0.5*corr;
        end
        change = max(abs(Unew(:)-U(:)));
        U = Unew;
        if change < tol, break; end
    end

    [Tg,Xg] = ndgrid(tGrid,xEval);
    Err = abs(U-uExactFun(Tg,Xg));
    out.e = max(max(Err(:,coreIdx)));
    out.U = U;  out.x = xEval;  out.core = coreIdx;
end

function W = supplementaryNCVolterraWeights(z,k)
    n = numel(z)-1;  hloc = z(2)-z(1);
    full = lagrangeIntegralWeights(k,k);
    rem = zeros(k,k+1);
    for q = 1:k-1
        rem(q,:) = lagrangeIntegralWeights(k,q);
    end
    W = zeros(n+1,n+1);
    for r = 0:n
        ell = floor(r/k);  q = mod(r,k);
        for block = 0:ell-1
            idxB = block*k + (0:k) + 1;
            W(r+1,idxB) = W(r+1,idxB) + hloc*full;
        end
        if q > 0
            idxB = ell*k + (0:k) + 1;
            W(r+1,idxB) = W(r+1,idxB) + hloc*rem(q,:);
        end
    end
end

function kernel = fixedDomainNCKernel(m,k,h)
    N = 2*m;  nFull = floor(N/k);  remPanels = N - nFull*k;
    kernel = zeros(1,N+1);
    if nFull > 0
        fullW = lagrangeIntegralWeights(k,k);
        for b = 0:nFull-1
            idxB = b*k + (0:k) + 1;
            kernel(idxB) = kernel(idxB) + fullW;
        end
    end
    if remPanels > 0
        remW = lagrangeIntegralWeights(remPanels,remPanels);
        idxB = nFull*k + (0:remPanels) + 1;
        kernel(idxB) = kernel(idxB) + remW;
    end
    kernel = h*kernel;
end

function weights = lagrangeIntegralWeights(k,q)
    weights = zeros(1,k+1);  nodes = 0:k;
    for j = 1:k+1
        others = nodes([1:j-1,j+1:end]);
        coeff = poly(others) / prod(nodes(j)-others);
        anti  = polyint(coeff);
        weights(j) = polyval(anti,q)-polyval(anti,0);
    end
end
function out = solveWaveGL(Ffun,Lfun,uExactFun,T,nBlocksT,mT,nBlocksX,mX,Xcore,padExtra,tol,maxIter,wantNodes)
    Xmax = Xcore + T + padExtra;  aX = -Xmax;  bX = Xmax;
    [Wt,tEval,tauSource] = supplementaryGaussLegendreVolterraWeights(0,T,nBlocksT,mT);
    nTau = numel(tauSource);
    tauRowIdx = zeros(1,nTau);
    for cc = 1:nTau
        [~,tauRowIdx(cc)] = min(abs(tEval-tauSource(cc)));
    end
    WtSS = Wt(tauRowIdx,:);

    [xiX,wgX] = gaussLegendreRule(mX);
    deltaX = (bX-aX)/nBlocksX;
    xSource = reshape((aX + (0:nBlocksX-1)*deltaX) + deltaX*(xiX(:)+1)/2,1,[]);
    nXs = numel(xSource);
    [AntiCoeffs,antiAtMinus1] = buildAntiderivativeCache(xiX);

    L = Lfun(tauSource(:),xSource);  U = L;
    for iter = 1:maxIter
        [FallBlocks,CumBAll] = blockData(Ffun(U),mX,nBlocksX,nTau,wgX,deltaX);
        Unew = L;
        for cRow = 1:nTau
            corr = zeros(1,nXs);
            for cp = find(WtSS(cRow,:) ~= 0)
                wd = tauSource(cRow) - tauSource(cp);
                if wd == 0, continue; end
                Gp = evalGvec(xSource+wd,aX,deltaX,nBlocksX,AntiCoeffs,antiAtMinus1,CumBAll(:,cp),FallBlocks(:,:,cp));
                Gm = evalGvec(xSource-wd,aX,deltaX,nBlocksX,AntiCoeffs,antiAtMinus1,CumBAll(:,cp),FallBlocks(:,:,cp));
                corr = corr + WtSS(cRow,cp)*(Gp-Gm);
            end
            Unew(cRow,:) = L(cRow,:) - 0.5*corr;
        end
        change = max(abs(Unew(:)-U(:)));  U = Unew;
        if change < tol, break; end
    end
    out.iters = iter;

    xOut = xSource(abs(xSource) <= Xcore);
    [FallBlocks,CumBAll] = blockData(Ffun(U),mX,nBlocksX,nTau,wgX,deltaX);
    Uout = zeros(numel(tEval),numel(xOut));
    for r = 1:numel(tEval)
        corr = zeros(1,numel(xOut));
        for cp = find(Wt(r,:) ~= 0)
            wd = tEval(r) - tauSource(cp);
            if wd == 0, continue; end
            Gp = evalGvec(xOut+wd,aX,deltaX,nBlocksX,AntiCoeffs,antiAtMinus1,CumBAll(:,cp),FallBlocks(:,:,cp));
            Gm = evalGvec(xOut-wd,aX,deltaX,nBlocksX,AntiCoeffs,antiAtMinus1,CumBAll(:,cp),FallBlocks(:,:,cp));
            corr = corr + Wt(r,cp)*(Gp-Gm);
        end
        Uout(r,:) = Lfun(tEval(r),xOut) - 0.5*corr;
    end
    [Tg,Xg] = ndgrid(tEval,xOut);
    out.e = max(max(abs(Uout - uExactFun(Tg,Xg))));
    if wantNodes
        out.uT = Uout(end,:);  out.xT = xOut;
    end
end

function [FallBlocks,CumBAll] = blockData(Fall,mX,nBlocksX,nTau,wgX,deltaX)
    FallBlocks = reshape(Fall.',mX,nBlocksX,nTau);
    BlockSums  = reshape(sum(FallBlocks.*wgX(:),1),nBlocksX,nTau);
    CumBAll    = [zeros(1,nTau); deltaX/2*cumsum(BlockSums,1)];
end


function [W,evalNodes,sourceNodes] = supplementaryGaussLegendreVolterraWeights(a,b,nBlocks,m)
% (from Example_1_Supplementary_Gauss_Legendre_1D.m, unmodified)
    [xiL,w] = gaussLegendreRule(m);
    delta = (b-a)/nBlocks;
    blockLeft = a + (0:nBlocks-1)*delta;
    sourceNodes = reshape(blockLeft + delta*(xiL(:)+1)/2,1,[]);
    targetNodes = linspace(a,b,nBlocks+1);
    evalNodes = unique(sort([targetNodes,sourceNodes]));
    W = zeros(numel(evalNodes),numel(sourceNodes));
    for r = 1:numel(evalNodes)
        t = evalNodes(r);
        if abs(t-b) <= 100*eps(max(1,abs(b)))
            complete = nBlocks; partial = 0;
        else
            complete = floor((t-a)/delta);
            left = a + complete*delta;
            partial = (t-left)/delta;
        end
        for block = 0:complete-1
            idx = block*m + (1:m);
            W(r,idx) = delta*w/2;
        end
        if partial > 100*eps
            eta = 2*partial-1;
            idx = complete*m + (1:m);
            W(r,idx) = delta/2 * lagrangePartialWeights(xiL,eta);
        end
    end
end

function [AntiCoeffs,antiAtMinus1] = buildAntiderivativeCache(xi)
    m = numel(xi);
    AntiCoeffs = zeros(m,m+1);  antiAtMinus1 = zeros(m,1);
    for j = 1:m
        other = xi([1:j-1,j+1:end]);
        coeff = poly(other) / prod(xi(j)-other);
        anti = polyint(coeff);
        AntiCoeffs(j,:) = anti;  antiAtMinus1(j) = polyval(anti,-1);
    end
end

function PW = evalPartialWeightsVec(AntiCoeffs,antiAtMinus1,etaVec)
    m = size(AntiCoeffs,1);
    PW = zeros(numel(etaVec),m);
    for j = 1:m
        PW(:,j) = polyval(AntiCoeffs(j,:),etaVec) - antiAtMinus1(j);
    end
end

function Gvals = evalGvec(Xvec,aX,deltaX,nBlocksX,AntiCoeffs,antiAtMinus1,CumB,FvecBlocks)
% G(X) = int_{aX}^{X} phi ds for all X in Xvec at once (O(mX) per point);
% points outside [aX,bX] are clamped.
    bX = aX + nBlocksX*deltaX;
    Xc = min(max(Xvec,aX),bX);
    rel = (Xc-aX)/deltaX;
    complete = min(floor(rel),nBlocksX);
    partial = rel - complete;
    atRight = complete >= nBlocksX;
    complete(atRight) = nBlocksX;  partial(atRight) = 0;
    Gvals = CumB(complete(:)+1).';
    need = partial > 100*eps & ~atRight;
    if any(need)
        eta = 2*partial(need)-1;
        PW = evalPartialWeightsVec(AntiCoeffs,antiAtMinus1,eta(:));
        Fsel = FvecBlocks(:,complete(need)+1);
        Gvals(need) = Gvals(need) + (deltaX/2)*sum(PW.'.*Fsel,1);
    end
end

function [x,w] = gaussLegendreRule(m)
    j = (1:m-1)';
    off = j ./ sqrt(4*j.^2-1);
    [V,D] = eig(diag(off,1)+diag(off,-1));
    [x,order] = sort(diag(D));
    V = V(:,order);
    w = 2*(V(1,:).^2).';
    x = x.';  w = w.';
end

function weights = lagrangePartialWeights(nodes,eta)
    m = numel(nodes);
    weights = zeros(1,m);
    for j = 1:m
        other = nodes([1:j-1,j+1:end]);
        coeff = poly(other) / prod(nodes(j)-other);
        anti = polyint(coeff);
        weights(j) = polyval(anti,eta)-polyval(anti,-1);
    end
end
