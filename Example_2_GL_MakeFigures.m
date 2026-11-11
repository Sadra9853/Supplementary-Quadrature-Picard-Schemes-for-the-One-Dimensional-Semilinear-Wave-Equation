%% Example_2_GL_MakeFigures.m
% Example 2 (sine-Gordon kink) with the SUPPLEMENTARY GAUSS-LEGENDRE rule.
%   m = number of TIME blocks (h = T/m), k = Gauss points per time block.
%   The spatial discretisation is kept fixed and fine (nBlocksX x kX) so
%   that the TIME rule is what is measured.
%   Output: console table, ex2gl_table.tex (rows for tab:ex2gl) and
%   ex2gl_convergence_error.pdf (left: convergence, right: pointwise error
%   at t=T, plotted at the spatial Gauss nodes).
%
% EXPECTED ORDER: a k-point Gauss rule is exact to degree 2k-1, but on the
% unfinished (supplementary) block the degree-(k-1) interpolant is
% integrated only up to t, with local error O(h^(k+1)). So at general
% evaluation times the observed rate may be close to k+1, not 2k-1.
% Read the printed orders and then set refOrder below; do not assume it.

clear; clc; close all;

v = 0.5;  gamma = sqrt(1 - v^2);
Ffun      = @(u) sin(u);
uExactFun = @(t,x) 4*atan( exp((x - v*t)/gamma) );
f0 = @(x) 4*atan( exp(x/gamma) );
G0 = @(x) -2*v*atan( sinh(x/gamma) );
Lfun = @(t,x) 0.5*( f0(x-t) + f0(x+t) ) + 0.5*( G0(x+t) - G0(x-t) );

T = 1;  Xcore = 10;  tol = 1e-12;  maxIter = 100;

kList = [2 3];            % Gauss points per TIME block (two table columns)
mList = [2 4 8 16];       % number of TIME blocks, h = T/m
nBlocksX = 52;  kX = 10;  % fixed fine spatial discretisation (width 26)
refOrder = @(k) k+1;      % reference slope h^refOrder(k): adjust after seeing the orders

nk = numel(kList);
E = zeros(numel(mList),nk);
for j = 1:nk
    for i = 1:numel(mList)
        out = solveWaveGL(Ffun,Lfun,uExactFun,T,mList(i),kList(j),nBlocksX,kX,Xcore,tol,maxIter,false);
        E(i,j) = out.e;
        fprintf('k=%d m=%d  e=%.3e  (iters=%d)\n',kList(j),mList(i),out.e,out.iters);
    end
end

Q = NaN(size(E));
for j = 1:nk
    for i = 2:numel(mList)
        Q(i,j) = log2(E(i-1,j)/E(i,j));
    end
end

fprintf('\n%6s','m');
for j = 1:nk, fprintf('   %12s %8s',sprintf('e_m (k=%d)',kList(j)),'order'); end
fprintf('\n');
for i = 1:numel(mList)
    fprintf('%6d',mList(i));
    for j = 1:nk
        if i==1, qs='   --   '; else, qs=sprintf('%8.2f',Q(i,j)); end
        fprintf('   %12.3e %s',E(i,j),qs);
    end
    fprintf('\n');
end

fid = fopen('ex2gl_table.tex','w');
for fo = [1 fid]
    fprintf(fo,'%% ---- paste into tab:ex2gl in wave_paper_1plus1.tex ----\n');
    fprintf(fo,'$m$ & $e_m$, $k=%d$ & order $q_m$ & $e_m$, $k=%d$ & order $q_m$ \\\\\n',kList(1),kList(2));
    fprintf(fo,'\\midrule\n');
    for i = 1:numel(mList)
        if i==1, q1='--'; q2='--'; else, q1=sprintf('%.2f',Q(i,1)); q2=sprintf('%.2f',Q(i,2)); end
        fprintf(fo,'$m_%d=%d$ & %s & %s & %s & %s \\\\\n',i,mList(i),texSci(E(i,1)),q1,texSci(E(i,2)),q2);
    end
end
fclose(fid);

%% Finest run: pointwise error at t = T at the spatial Gauss nodes
kF = max(kList);  mF = max(mList);
out = solveWaveGL(Ffun,Lfun,uExactFun,T,mF,kF,nBlocksX,kX,Xcore,tol,maxIter,true);
err = abs(out.uT - uExactFun(T,out.xT));  err(err==0) = eps;

h = T./mList;  colors = lines(nk);  markers = {'o','s'};
fig = figure('Color','w','Units','centimeters','Position',[2 2 14 6]);
subplot(1,2,1);
for j = 1:nk
    loglog(h,E(:,j),[markers{j} '-'],'Color',colors(j,:),'LineWidth',1.6,'MarkerSize',6); hold on;
end
for j = 1:nk
    loglog(h,E(1,j)*(h/h(1)).^refOrder(kList(j)),'--','Color',colors(j,:)*0.55,'LineWidth',1);
end
grid on; xlabel('h = T/m'); ylabel('e_m');
legend([arrayfun(@(k) sprintf('k=%d',k),kList,'UniformOutput',false), ...
        arrayfun(@(k) sprintf('h^{%d}',refOrder(k)),kList,'UniformOutput',false)], 'Location','southeast');
subplot(1,2,2);
semilogy(out.xT,err,'.-','Color',[0 0 0.8],'LineWidth',1,'MarkerSize',6);
grid on; xlabel('x'); ylabel('|u_h(T,x) - u(T,x)|');
set(findall(fig,'-property','FontSize'),'FontSize',10);
exportgraphics(fig,'ex2gl_convergence_error.pdf','ContentType','vector');
fprintf('Saved ex2gl_table.tex and ex2gl_convergence_error.pdf\n');

%% ------------------------------------------------------------------
function out = solveWaveGL(Ffun,Lfun,uExactFun,T,nBlocksT,mT,nBlocksX,mX,Xcore,tol,maxIter,wantNodes)
    Xmax = Xcore + 2*T + 1;  aX = -Xmax;  bX = Xmax;
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

function str = texSci(x)
    e = floor(log10(abs(x)));  mant = x/10^e;
    str = sprintf('$%.3f\\times10^{%d}$',mant,e);
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
