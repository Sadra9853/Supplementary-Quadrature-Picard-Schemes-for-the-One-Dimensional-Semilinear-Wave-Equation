%% Example_1_Wave_Equation_Supplementary_Gauss_Legendre_Figures.m

%
%  Panel 1 : log-log convergence e_m vs h = T/m with reference slopes
%  Panel 2 : pointwise error |u_h(T,x) - u(T,x)| at the finest run
%  Panel 3 : 3D surface of log10|u_h(t,x) - u(t,x)| on (t,x)

clear; clc; close all;

%% ---------- Problem setup ----------
v   = 0.4;
mu  = 1 - v^2;

Ffun      = @(u) mu*(u - 2*u.^3);
uExactFun = @(t,x) 1 ./ cosh(x - v*t);

f0 = @(x) 1 ./ cosh(x);
G0 = @(x) v*(1 - 1 ./ cosh(x));

Lfun = @(t,x) 0.5*( f0(x-t) + f0(x+t) ) + ...
              0.5*( G0(x+t) - G0(x-t) );

T = 1;  Xcore = 5;  tol = 1e-12;  maxIter = 60;

% k labels identical to NC paper table:
%   k=2  <->  mT=2 (time order 4)
%   k=4  <->  mT=3 (time order 6)
kList  = [2 10];
mTList = [2 3];
mList  = [8 40 60];
nk     = numel(kList);

%% ---------- Convergence study ----------
E = zeros(numel(mList), nk);
for j = 1:nk
    for i = 1:numel(mList)
        out = solveWaveGL(Ffun,Lfun,uExactFun,T, ...
                          mList(i),mTList(j),Xcore,tol,maxIter,false);
        E(i,j) = out.e;
        fprintf('k=%d (mT=%d)  m=%2d  iter=%2d  e=%.3e\n', ...
                kList(j),mTList(j),mList(i),out.iter,out.e);
    end
end

%% ---------- Empirical orders ----------
Q = NaN(numel(mList), nk);
for j = 1:nk
    for i = 2:numel(mList)
        Q(i,j) = log2(E(i-1,j)/E(i,j));
    end
end

%% ---------- Plain-text + LaTeX table ----------
fprintf('\n%6s','m');
for j = 1:nk
    fprintf('   %12s %8s', sprintf('e_m (k=%d)',kList(j)), 'order');
end
fprintf('\n');
for i = 1:numel(mList)
    fprintf('%6d', mList(i));
    for j = 1:nk
        if i == 1, qs = '   --   '; else, qs = sprintf('%8.2f',Q(i,j)); end
        fprintf('   %12.3e %s', E(i,j), qs);
    end
    fprintf('\n');
end

fid = fopen('ex1_table_GL.tex','w');
for fo = [1 fid]
    fprintf(fo,'$m$ & $e_m$, $k=%d$ & order $q_m$ & $e_m$, $k=%d$ & order $q_m$ \\\\\n', ...
            kList(1), kList(2));
    fprintf(fo,'\\midrule\n');
    for i = 1:numel(mList)
        if i == 1, q1 = '--'; q2 = '--';
        else,      q1 = sprintf('%.2f',Q(i,1)); q2 = sprintf('%.2f',Q(i,2)); end
        fprintf(fo,'$m_%d=%d$ & %s & %s & %s & %s \\\\\n', ...
            i, mList(i), texSci(E(i,1)), q1, texSci(E(i,2)), q2);
    end
end
fclose(fid);
fprintf('LaTeX rows written to ex1_table_GL.tex\n\n');

%% ---------- Finest run for panels 2 & 3 ----------
kFinest  = kList(end);
mTFinest = mTList(end);
mFinest  = mList(end);

out = solveWaveGL(Ffun,Lfun,uExactFun,T, ...
                  mFinest,mTFinest,Xcore,tol,maxIter,true);

x   = out.x;
uh  = out.U(end,:);
ue  = uExactFun(T,x);
err = abs(uh-ue);  err(err==0) = eps;

tGrid = out.t;
Uc    = out.U;
[Tg,Xg] = ndgrid(tGrid,x);
ErrST = abs(Uc - uExactFun(Tg,Xg));
ErrST(ErrST==0) = eps;

h = T ./ mList;

%% ============================================================
%%  Figure  (identical layout to Newton–Cotes file)
%% ============================================================
markers = {'o','s','d','^'};
colors  = lines(nk);

f = figure('Color','w','Units','centimeters','Position',[2 2 22 6.5]);

% ---- Panel 1 ----
ax1 = subplot(1,3,1);
hLines = gobjects(2*nk,1);
labels = cell(2*nk,1);
for j = 1:nk
    mk = markers{mod(j-1,numel(markers))+1};
    hLines(j) = loglog(ax1, h, E(:,j), [mk '-'], ...
        'Color', colors(j,:), 'LineWidth',2, 'MarkerSize',10, ...
        'MarkerFaceColor','w');
    hold(ax1,'on');
    labels{j} = sprintf('k=%d', kList(j));

    refSlope = E(1,j) * (h/h(1)).^(2*mTList(j));
    hLines(nk+j) = loglog(ax1, h, refSlope, '--', ...
        'Color', colors(j,:)*0.55, 'LineWidth',2);
    labels{nk+j} = sprintf('h^{%d}', 2*mTList(j));
end
set(ax1,'XScale','log','YScale','log'); grid(ax1,'on');
xlabel(ax1,'h = T/m'); ylabel(ax1,'e_m');
legend(ax1,hLines,labels,'Location','southeast');
ax1.FontSize = 12;

infoLeft = { ...
    sprintf('Cubic soliton,  v = %.2f', v), ...
    sprintf('T = %g,  X_{core} = %g', T, Xcore), ...
    sprintf('Gauss–Legendre supplementary'), ...
    sprintf('m \\in \\{%s\\}', strjoin(string(mList), ', '))};
for j = 1:nk
    infoLeft{end+1} = sprintf('e(k=%d,m=%d) = %.2e', ...
        kList(j), mList(end), E(end,j));
end
text(ax1,0.03,0.97,infoLeft,'Units','normalized', ...
    'VerticalAlignment','top','HorizontalAlignment','left', ...
    'FontSize',9,'BackgroundColor',[1 1 1 0.85], ...
    'EdgeColor',[0.5 0.5 0.5],'Margin',4,'Interpreter','tex');

% ---- Panel 2 ----
ax2 = subplot(1,3,2);
semilogy(ax2, x, err, 'b-', 'LineWidth',2);
set(ax2,'YScale','log'); grid(ax2,'on');
xlabel(ax2,'x'); ylabel(ax2,'|u_h(T,x) - u(T,x)|');
ax2.FontSize = 12;

[emax,imax] = max(err);
infoRight = { ...
    sprintf('Finest: k=%d, m=%d', kFinest, mFinest), ...
    sprintf('||u_h - u||_\\infty = %.2e', emax), ...
    sprintf('at x = %.3f', x(imax))};
text(ax2,0.03,0.97,infoRight,'Units','normalized', ...
    'VerticalAlignment','top','HorizontalAlignment','left', ...
    'FontSize',9,'BackgroundColor',[1 1 1 0.85], ...
    'EdgeColor',[0.5 0.5 0.5],'Margin',4,'Interpreter','tex');

% ---- Panel 3 ----
ax3 = subplot(1,3,3);
Z = log10(ErrST);
surf(ax3, x, tGrid, Z, 'EdgeColor','none');
hold(ax3,'on');
plot3(ax3, x, T*ones(size(x)), Z(end,:), 'r', 'LineWidth',1.5);
shading(ax3,'interp'); colormap(ax3, parula);
cb = colorbar(ax3);
cb.Label.String = 'log_{10}|u_h - u|';
cb.Label.Interpreter = 'tex';
view(ax3,40,30); grid(ax3,'on');
xlabel(ax3,'x'); ylabel(ax3,'t'); zlabel(ax3,'log_{10}|u_h - u|');
ax3.FontSize = 12;

infoST = { ...
    sprintf('Finest: k=%d, m=%d', kFinest, mFinest), ...
    sprintf('max|u_h-u| = %.2e', emax)};
text(ax3,-0.05,0.97,infoST,'Units','normalized', ...
    'VerticalAlignment','top','HorizontalAlignment','left', ...
    'FontSize',9,'BackgroundColor',[1 1 1 0.85], ...
    'EdgeColor',[0.5 0.5 0.5],'Margin',4,'Interpreter','tex');


%% ============================================================
%%  Gauss–Legendre solver  (SPACE GRID TIED TO TIME GRID)
%% ============================================================
function out = solveWaveGL(Ffun,Lfun,uExactFun,T,nBlocksT,mT,Xcore,tol,maxIter,computeFull)

    % ---------------- Time ----------------
    [Wt,tEval,tauSource] = ...
        supplementaryGaussLegendreVolterraWeights(0,T,nBlocksT,mT);
    nTau = numel(tauSource);

    tauRowIdx = zeros(1,nTau);
    for cc = 1:nTau
        [~,tauRowIdx(cc)] = min(abs(tEval - tauSource(cc)));
    end
    WtSS = Wt(tauRowIdx,:);

    % ---------------- Space (FIXED: refine with time) ----------------
    %   delta_x = 2*delta_t  (space step proportional to time step)
    %   mX      = 8          (low degree => no polynomial cancellation)
    delta_t  = T / nBlocksT;
    delta_x  = 2 * delta_t;
    Xmax     = Xcore + 2*T + 1;

    nBlocksX = max(2, ceil(2*Xmax / delta_x));
    mX       = 8;

    aX = -Xmax;  bX = Xmax;

    [xiX,wgX] = gaussLegendreRule(mX);
    deltaX    = (bX-aX)/nBlocksX;
    blockLeftX = aX + (0:nBlocksX-1)*deltaX;
    xSource   = reshape(blockLeftX + deltaX*(xiX(:)+1)/2, 1, []);
    nXs       = numel(xSource);

    [AntiCoeffs, antiAtMinus1] = buildAntiderivativeCache(xiX);

    if computeFull
        fprintf('  space: %d blocks x %d pts = %d nodes (delta_x=%.4g)\n', ...
                nBlocksX, mX, nXs, deltaX);
    end

    % ---------------- Picard iteration ----------------
    L = Lfun(tauSource(:), xSource);
    U = L;
    iter = 0;
    for iter = 1:maxIter
        Fall = Ffun(U);
        FallBlocks = reshape(Fall.',mX,nBlocksX,nTau);
        BlockSums  = squeeze(sum(FallBlocks.*wgX(:),1));
        if nTau==1, BlockSums = BlockSums(:); end
        CumBAll = [zeros(1,nTau); deltaX/2*cumsum(BlockSums,1)];

        Unew = L;
        for cRow = 1:nTau
            idx  = find(WtSS(cRow,:) ~= 0);
            corr = zeros(1,nXs);
            for cp = idx
                wd = tauSource(cRow) - tauSource(cp);
                if wd==0, continue; end
                FvecBlocks = FallBlocks(:,:,cp);
                CumB       = CumBAll(:,cp);
                Gp = evalGvec(xSource+wd,aX,deltaX,nBlocksX, ...
                     AntiCoeffs,antiAtMinus1,CumB,FvecBlocks);
                Gm = evalGvec(xSource-wd,aX,deltaX,nBlocksX, ...
                     AntiCoeffs,antiAtMinus1,CumB,FvecBlocks);
                corr = corr + WtSS(cRow,cp)*(Gp-Gm);
            end
            Unew(cRow,:) = L(cRow,:) - 0.5*corr;
        end
        change = max(abs(Unew(:)-U(:)));
        U = Unew;
        if change < tol, break; end
    end

    % ---------------- Error at t = T ----------------
    Fall = Ffun(U);
    FallBlocks = reshape(Fall.',mX,nBlocksX,nTau);
    BlockSums  = squeeze(sum(FallBlocks.*wgX(:),1));
    if nTau==1, BlockSums = BlockSums(:); end
    CumBAll = [zeros(1,nTau); deltaX/2*cumsum(BlockSums,1)];

    coreIdx = abs(xSource) <= Xcore + 1e-9;

    [~,iT] = min(abs(tEval - T));
    UT = computeAtTime(T,Lfun,Wt(iT,:),tauSource,xSource, ...
                       aX,deltaX,nBlocksX,AntiCoeffs,antiAtMinus1, ...
                       CumBAll,FallBlocks);
    out.e    = max(abs(UT(coreIdx) - uExactFun(T,xSource(coreIdx))));
    out.iter = iter;

    if computeFull
        nEval = numel(tEval);
        Uout  = zeros(nEval,nXs);
        for r = 1:nEval
            Uout(r,:) = computeAtTime(tEval(r),Lfun,Wt(r,:), ...
                tauSource,xSource,aX,deltaX,nBlocksX, ...
                AntiCoeffs,antiAtMinus1,CumBAll,FallBlocks);
        end
        out.U = Uout(:,coreIdx);
        out.x = xSource(coreIdx);
        out.t = tEval;
    end
end

function Urow = computeAtTime(t,Lfun,wrow,tauSource,xSource, ...
    aX,deltaX,nBlocksX,AntiCoeffs,antiAtMinus1,CumBAll,FallBlocks)
    idx  = find(wrow ~= 0);
    Lrow = Lfun(t,xSource);
    corr = zeros(1,numel(xSource));
    for cp = idx
        wd = t - tauSource(cp);
        if wd==0, continue; end
        FvecBlocks = FallBlocks(:,:,cp);
        CumB       = CumBAll(:,cp);
        Gp = evalGvec(xSource+wd,aX,deltaX,nBlocksX, ...
             AntiCoeffs,antiAtMinus1,CumB,FvecBlocks);
        Gm = evalGvec(xSource-wd,aX,deltaX,nBlocksX, ...
             AntiCoeffs,antiAtMinus1,CumB,FvecBlocks);
        corr = corr + wrow(cp)*(Gp-Gm);
    end
    Urow = Lrow - 0.5*corr;
end


%% ============================================================
%%  Helper functions (unchanged from original GL file)
%% ============================================================
function [W,evalNodes,sourceNodes] = ...
    supplementaryGaussLegendreVolterraWeights(a,b,nBlocks,m)

    [xiL,w] = gaussLegendreRule(m);
    delta   = (b-a)/nBlocks;
    blockLeft = a + (0:nBlocks-1)*delta;
    sourceNodes = reshape(blockLeft + delta*(xiL(:)+1)/2,1,[]);
    targetNodes = linspace(a,b,nBlocks+1);
    evalNodes   = unique(sort([targetNodes,sourceNodes]));

    W = zeros(numel(evalNodes),numel(sourceNodes));
    for r = 1:numel(evalNodes)
        t = evalNodes(r);
        if abs(t-b) <= 100*eps(max(1,abs(b)))
            complete = nBlocks; partial = 0;
        else
            complete = floor((t-a)/delta);
            left     = a + complete*delta;
            partial  = (t-left)/delta;
        end
        for block = 0:complete-1
            idx = block*m + (1:m);
            W(r,idx) = delta*w/2;
        end
        if partial > 100*eps
            eta = 2*partial-1;
            idx = complete*m + (1:m);
            W(r,idx) = delta/2*lagrangePartialWeights(xiL,eta);
        end
    end
end

function [AntiCoeffs,antiAtMinus1] = buildAntiderivativeCache(xi)
    m = numel(xi);
    AntiCoeffs = zeros(m,m+1);
    antiAtMinus1 = zeros(m,1);
    for j = 1:m
        other = xi([1:j-1,j+1:end]);
        coeff = poly(other)/prod(xi(j)-other);
        anti  = polyint(coeff);
        AntiCoeffs(j,:)   = anti;
        antiAtMinus1(j)   = polyval(anti,-1);
    end
end

function PW = evalPartialWeightsVec(AntiCoeffs,antiAtMinus1,etaVec)
    m  = size(AntiCoeffs,1);
    PW = zeros(numel(etaVec),m);
    for j = 1:m
        PW(:,j) = polyval(AntiCoeffs(j,:),etaVec) - antiAtMinus1(j);
    end
end

function Gvals = evalGvec(Xvec,aX,deltaX,nBlocksX, ...
    AntiCoeffs,antiAtMinus1,CumB,FvecBlocks)

    bX = aX + nBlocksX*deltaX;
    Xc = min(max(Xvec,aX),bX);
    rel = (Xc-aX)/deltaX;
    complete = min(floor(rel),nBlocksX);
    partial  = rel - complete;

    atRightEdge = complete >= nBlocksX;
    complete(atRightEdge) = nBlocksX;
    partial(atRightEdge)  = 0;

    Gvals = CumB(complete(:)+1).';

    needPartial = partial > 100*eps & ~atRightEdge;
    if any(needPartial)
        eta  = 2*partial(needPartial)-1;
        PW   = evalPartialWeightsVec(AntiCoeffs,antiAtMinus1,eta(:));
        blks = complete(needPartial)+1;
        Fsel = FvecBlocks(:,blks);
        contrib = sum(PW.'.*Fsel,1);
        Gvals(needPartial) = Gvals(needPartial) + (deltaX/2)*contrib;
    end
end

function [x,w] = gaussLegendreRule(m)
    j   = (1:m-1)';
    off = j./sqrt(4*j.^2-1);
    [V,D] = eig(diag(off,1)+diag(off,-1));
    [x,ord] = sort(diag(D));
    V = V(:,ord);
    w = 2*(V(1,:).^2).';
    x = x.';  w = w.';
end

function weights = lagrangePartialWeights(nodes,eta)
    m = numel(nodes);
    weights = zeros(1,m);
    for j = 1:m
        other = nodes([1:j-1,j+1:end]);
        coeff = poly(other)/prod(nodes(j)-other);
        anti  = polyint(coeff);
        weights(j) = polyval(anti,eta) - polyval(anti,-1);
    end
end

function str = texSci(x)
    e    = floor(log10(abs(x)));
    mant = x/10^e;
    str  = sprintf('$%.3f\\times10^{%d}$',mant,e);
end
