%% Examples_Unified_Sweep.m
% Runs, for the three examples of the wave paper (cubic soliton, sine-Gordon
% kink, sine-Gordon breather), BOTH supplementary families with matched
% orders:
%     Newton-Cotes   k = 2, 4   (predicted order p0 = k+2 = 4, 6)
%     Gauss-Legendre k = 3, 5   (predicted order p0 = k+1 = 4, 6)
% and writes everything needed to complete Section 4:
%     ex<n>_nc_table.tex, ex<n>_gl_table.tex   rows for the NC / GL tables
%     summary_rows.tex                         rows for tab:summary
%     iters_time_rows.tex                      Picard iterations and CPU time
%     ex<n>gl_convergence_error.pdf (n=1,2)    GL figure: convergence + pointwise error
%     efficiency_Nt.pdf                        error versus number of TIME nodes, 3 panels
% Set runEx to run only some examples (Example 3 is the slowest).
%
% Number of time nodes N_t: Newton-Cotes m+1, Gauss-Legendre k*m.
% Errors are max errors on |x| <= Xcore at all time nodes.  The spatial
% discretisation differs between the families (NC: h_x = h; GL: fixed fine
% Gauss mesh), so N_t compares TIME nodes only.

clear; clc; close all;

runEx = [1 2 3];
kNC = [2 4];   kGL = [3 5];   kX = 10;   hxGL = 0.5;
tol = 1e-12;   maxIter = 100;

EX = setupExamples();
colors = lines(4);
effFig = figure('Color','w','Units','centimeters','Position',[2 2 22 6.5]);
summaryFid = fopen('summary_rows.tex','w');
itFid      = fopen('iters_time_rows.tex','w');

for ie = runEx
    P = EX(ie);
    nBlocksX = round(2*(P.Xcore+P.T+P.padExtra)/hxGL);
    nm = numel(P.mNC);  ng = numel(P.mGL);
    ENC = zeros(nm,2); ITNC = zeros(nm,2); TNC = zeros(nm,2);
    EGL = zeros(ng,2); ITGL = zeros(ng,2); TGL = zeros(ng,2);

    for j = 1:2
        for i = 1:nm
            tic; o = solveWaveNC(P.F,P.L,P.uex,P.T,kNC(j),P.mNC(i),P.Xcore,P.padExtra,tol,maxIter);
            TNC(i,j) = toc; ENC(i,j) = o.e; ITNC(i,j) = o.iters;
            fprintf('Ex%d NC k=%d m=%3d  e=%.3e  iters=%d  t=%.1fs\n',ie,kNC(j),P.mNC(i),o.e,o.iters,TNC(i,j));
        end
        for i = 1:ng
            tic; o = solveWaveGL(P.F,P.L,P.uex,P.T,P.mGL(i),kGL(j),nBlocksX,kX,P.Xcore,P.padExtra,tol,maxIter,false);
            TGL(i,j) = toc; EGL(i,j) = o.e; ITGL(i,j) = o.iters;
            fprintf('Ex%d GL k=%d m=%3d  e=%.3e  iters=%d  t=%.1fs\n',ie,kGL(j),P.mGL(i),o.e,o.iters,TGL(i,j));
        end
    end
    QNC = orders(ENC);  QGL = orders(EGL);

    writeTable(sprintf('ex%d_nc_table.tex',ie),P.mNC,ENC,QNC,true);
    writeTable(sprintf('ex%d_gl_table.tex',ie),P.mGL,EGL,QGL,false);

    for j = 1:2
        fprintf(summaryFid,'%d & NC & %d & %d & %.2f, %.2f & %s \\\\\n',ie,kNC(j),kNC(j)+2,QNC(end-1,j),QNC(end,j),texSci(ENC(end,j)));
    end
    for j = 1:2
        fprintf(summaryFid,'%d & GL & %d & %d & %.2f, %.2f & %s \\\\\n',ie,kGL(j),kGL(j)+1,QGL(end-1,j),QGL(end,j),texSci(EGL(end,j)));
    end
    for j = 1:2
        fprintf(itFid,'%d & NC & %d & %d--%d & %.1f \\\\\n',ie,kNC(j),min(ITNC(:,j)),max(ITNC(:,j)),TNC(end,j));
    end
    for j = 1:2
        fprintf(itFid,'%d & GL & %d & %d--%d & %.1f \\\\\n',ie,kGL(j),min(ITGL(:,j)),max(ITGL(:,j)),TGL(end,j));
    end

    % efficiency panel: error versus number of TIME nodes
    figure(effFig); subplot(1,3,ie);
    loglog(P.mNC+1,ENC(:,1),'o--','Color',colors(1,:),'LineWidth',1.3,'MarkerSize',5); hold on;
    loglog(P.mNC+1,ENC(:,2),'o-','Color',colors(2,:),'LineWidth',1.3,'MarkerSize',5);
    loglog(kGL(1)*P.mGL,EGL(:,1),'s--','Color',colors(3,:),'LineWidth',1.3,'MarkerSize',5);
    loglog(kGL(2)*P.mGL,EGL(:,2),'s-','Color',colors(4,:),'LineWidth',1.3,'MarkerSize',5);
    grid on; xlabel('time nodes N_t'); title(sprintf('Example %d',ie));
    if ie == runEx(1), ylabel('e'); legend('NC k=2','NC k=4','GL k=3','GL k=5','Location','southwest'); end

    % GL figure for Examples 1 and 2 (convergence + pointwise error at t=T)
    if ie <= 2
        oG = solveWaveGL(P.F,P.L,P.uex,P.T,P.mGL(end-1),kGL(2),nBlocksX,kX,P.Xcore,P.padExtra,tol,maxIter,true);
        err = abs(oG.uT - P.uex(P.T,oG.xT)); err(err==0) = eps;
        h = P.T./P.mGL;
        f = figure('Color','w','Units','centimeters','Position',[2 2 14 6]);
        subplot(1,2,1);
        loglog(h,EGL(:,1),'o-','LineWidth',1.6,'MarkerSize',6); hold on;
        loglog(h,EGL(:,2),'s-','LineWidth',1.6,'MarkerSize',6);
        loglog(h,EGL(1,1)*(h/h(1)).^(kGL(1)+1),'k--','LineWidth',1);
        loglog(h,EGL(1,2)*(h/h(1)).^(kGL(2)+1),'k:','LineWidth',1.2);
        grid on; xlabel('h = T/m'); ylabel('e_m');
        legend(sprintf('k=%d',kGL(1)),sprintf('k=%d',kGL(2)),sprintf('h^{%d}',kGL(1)+1),sprintf('h^{%d}',kGL(2)+1),'Location','southeast');
        subplot(1,2,2);
        semilogy(oG.xT,err,'.-','Color',[0 0 0.8],'LineWidth',1,'MarkerSize',6);
        grid on; xlabel('x'); ylabel('|u_h(T,x) - u(T,x)|');
        set(findall(f,'-property','FontSize'),'FontSize',10);
        exportgraphics(f,sprintf('ex%dgl_convergence_error.pdf',ie),'ContentType','vector');
    end
end
fclose(summaryFid); fclose(itFid);
set(findall(effFig,'-property','FontSize'),'FontSize',10);
exportgraphics(effFig,'efficiency_Nt.pdf','ContentType','vector');
fprintf('Done. Files written: ex*_table.tex, summary_rows.tex, iters_time_rows.tex, *.pdf\n');

%% ------------------------------------------------------------------
function EX = setupExamples()
    % Example 1: cubic soliton
    v = 0.4;  mu = 1 - v^2;
    EX(1).F   = @(u) mu*(u - 2*u.^3);
    EX(1).uex = @(t,x) sech(x - v*t);
    EX(1).L   = @(t,x) 0.5*(1+v)*sech(x-t) + 0.5*(1-v)*sech(x+t);
    EX(1).T = 1;  EX(1).Xcore = 5;  EX(1).padExtra = 2;
    EX(1).mNC = [8 16 32 64];  EX(1).mGL = [2 4 8 16];
    % Example 2: sine-Gordon kink
    v2 = 0.5;  g2 = sqrt(1-v2^2);
    f0 = @(x) 4*atan(exp(x/g2));  G0 = @(x) -2*v2*atan(sinh(x/g2));
    EX(2).F   = @(u) sin(u);
    EX(2).uex = @(t,x) 4*atan(exp((x - v2*t)/g2));
    EX(2).L   = @(t,x) 0.5*(f0(x-t)+f0(x+t)) + 0.5*(G0(x+t)-G0(x-t));
    EX(2).T = 1;  EX(2).Xcore = 10;  EX(2).padExtra = 2;
    EX(2).mNC = [8 16 32 64];  EX(2).mGL = [2 4 8 16];
    % Example 3: sine-Gordon breather
    w = 0.5;  s = sqrt(1-w^2);
    EX(3).F   = @(u) sin(u);
    EX(3).uex = @(t,x) 4*atan(s*cos(w*t)./(w*cosh(s*x)));
    fb = @(x) 4*atan(s./(w*cosh(s*x)));
    EX(3).L   = @(t,x) 0.5*(fb(x+t)+fb(x-t));
    EX(3).T = 1;  EX(3).Xcore = 6;  EX(3).padExtra = 23;
    EX(3).mNC = [4 8 16 32];  EX(3).mGL = [4 8 16 32];
end

function Q = orders(E)
    Q = NaN(size(E));
    for j = 1:size(E,2)
        for i = 2:size(E,1)
            Q(i,j) = log2(E(i-1,j)/E(i,j));
        end
    end
end

function writeTable(fname,mList,E,Q,withNt)
    fid = fopen(fname,'w');
    for i = 1:numel(mList)
        if i == 1, q1 = '--'; q2 = '--'; else, q1 = sprintf('%.2f',Q(i,1)); q2 = sprintf('%.2f',Q(i,2)); end
        if withNt
            fprintf(fid,'$m_%d=%d$ & %d & %s & %s & %s & %s \\\\\n',i,mList(i),mList(i)+1,texSci(E(i,1)),q1,texSci(E(i,2)),q2);
        else
            fprintf(fid,'$m_%d=%d$ & %s & %s & %s & %s \\\\\n',i,mList(i),texSci(E(i,1)),q1,texSci(E(i,2)),q2);
        end
    end
    fclose(fid);
end

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
    out.iters = iter;
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
