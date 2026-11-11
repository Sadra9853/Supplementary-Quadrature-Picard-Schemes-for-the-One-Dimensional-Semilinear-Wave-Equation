%% Example_2_MakeFigures.m
%   - Panel 1: log-log convergence e_m vs h with reference slopes h^(k+2)
%   - Panel 2: pointwise absolute error |u_h(T,x) - u(T,x)|
%   - Panel 3: 3D surface of log10|u_h(t,x) - u(t,x)| over the space-time grid
% (Case 2: sine-Gordon kink, supplementary Newton-Cotes).

clear; clc; close all;

%% ---------- Problem setup: Sine-Gordon Kink ----------
% Equation : u_tt - u_xx + sin(u) = 0
% Exact    : u(t,x) = 4*arctan( exp((x - v*t)/gamma) ),  gamma = sqrt(1 - v^2)
% Initial  : f(x) = 4*arctan(exp(x/gamma))
%            g(x) = -(2v/gamma) * sech(x/gamma)

v     = 0.5;
gamma = sqrt(1 - v^2);

Ffun      = @(u) sin(u);
uExactFun = @(t,x) 4*atan( exp((x - v*t)/gamma) );

% d'Alembert data for the linear wave operator
f0 = @(x) 4*atan( exp(x/gamma) );               % u(0,x)
G0 = @(x) -2*v*atan( sinh(x/gamma) );           % ∫_0^x g(s) ds

Lfun = @(t,x) 0.5*( f0(x-t) + f0(x+t) ) + ...
              0.5*( G0(x+t) - G0(x-t) );

T = 1;  Xcore = 10;  tol = 1e-12;  maxIter = 200;

% ============================================================
kList = [2 4];
mList = [8 16 32 64];
% ============================================================

kFinest = max(kList);
mFinest = max(mList);
nk      = numel(kList);

%% ---------- Convergence study ----------
E = zeros(numel(mList), nk);
for j = 1:nk
    for i = 1:numel(mList)
        out = solveWaveNC(Ffun,Lfun,uExactFun,T,kList(j),mList(i),Xcore,tol,maxIter);
        E(i,j) = out.e;
        fprintf('k=%d m=%d e=%.3e\n', kList(j), mList(i), out.e);
    end
end

%% ---------- Table for the paper (tab:ex2) ----------
% Empirical orders q_m = log2(e_m / e_{2m}), column by column.
Q = NaN(numel(mList), nk);
for j = 1:nk
    for i = 2:numel(mList)
        Q(i,j) = log2(E(i-1,j)/E(i,j));
    end
end

% Plain-text table in the command window
fprintf('\n%6s', 'm');
for j = 1:nk, fprintf('   %12s %8s', sprintf('e_m (k=%d)',kList(j)), 'order'); end
fprintf('\n');
for i = 1:numel(mList)
    fprintf('%6d', mList(i));
    for j = 1:nk
        if i == 1, qs = '   --   '; else, qs = sprintf('%8.2f',Q(i,j)); end
        fprintf('   %12.3e %s', E(i,j), qs);
    end
    fprintf('\n');
end

% Ready-to-paste LaTeX rows (also saved to ex2_table.tex)
fid = fopen('ex2_table.tex','w');
for fo = [1 fid]
    fprintf(fo,'%% ---- paste into tab:ex2 (Newton--Cotes) in wave_paper_1plus1.tex ----\n');
    fprintf(fo,'$m$ & $e_m$, $k=%d$ & order $q_m$ & $e_m$, $k=%d$ & order $q_m$ \\\\\n', kList(1), kList(2));
    fprintf(fo,'\\midrule\n');
    for i = 1:numel(mList)
        if i == 1
            q1 = '--'; q2 = '--';
        else
            q1 = sprintf('%.2f',Q(i,1)); q2 = sprintf('%.2f',Q(i,2));
        end
        fprintf(fo,'$m_%d=%d$ & %s & %s & %s & %s \\\\\n', ...
            i, mList(i), texSci(E(i,1)), q1, texSci(E(i,2)), q2);
    end
end
fclose(fid);
fprintf('LaTeX rows written to ex2_table.tex\n\n');

%% ---------- Finest run (for panels 2 & 3) ----------
out = solveWaveNC(Ffun,Lfun,uExactFun,T,kFinest,mFinest,Xcore,tol,maxIter);
x  = out.x(out.core);
uh = out.U(end,out.core);
ue = uExactFun(T,x);

h   = T./mList;
err = abs(uh - ue);
err(err == 0) = eps;

% --- full space-time error on the core region (for 3D panel) ---
tGrid = linspace(0,T,mFinest+1);
xc    = out.x(out.core);
Uc    = out.U(:, out.core);
[Tg, Xg] = ndgrid(tGrid, xc);
ErrST = abs( Uc - uExactFun(Tg, Xg) );    % (nt+1) x Nx_core
ErrST(ErrST == 0) = eps;

%% ============================================================
%%  Figure
%% ============================================================
markers = {'o','s','d','^','v','>','<','p','h'};
colors  = lines(nk);

f = figure('Color','w','Units','centimeters','Position',[2 2 22 6.5]);

% ------------------------------------------------------------
% Panel 1: log-log convergence
% ------------------------------------------------------------
ax1 = subplot(1,3,1);

hLines = gobjects(2*nk,1);
labels = cell(2*nk,1);
for j = 1:nk
    mk = markers{mod(j-1,numel(markers))+1};

    hLines(j) = loglog(ax1, h, E(:,j), [mk '-'], ...
        'Color', colors(j,:), 'LineWidth', 2, 'MarkerSize', 10, ...
        'MarkerFaceColor','w');
    hold(ax1,'on');
    labels{j} = sprintf('k=%d', kList(j));

    refSlope = E(1,j) * (h/h(1)).^(kList(j)+2);
    hLines(nk+j) = loglog(ax1, h, refSlope, '--', ...
        'Color', colors(j,:)*0.55, 'LineWidth', 2);
    labels{nk+j} = sprintf('h^{%d}', kList(j)+2);
end

set(ax1, 'XScale','log', 'YScale','log');
grid(ax1,'on');
xlabel(ax1,'h = T/m'); ylabel(ax1,'e_m');
legend(ax1, hLines, labels, 'Location','southeast');
ax1.FontSize = 12;

infoLeft = { ...
    sprintf('Sine-Gordon kink,  v = %.2f', v), ...
    sprintf('T = %g,  X_{core} = %g', T, Xcore), ...
    sprintf('k \\in \\{%s\\}', strjoin(string(kList), ', ')), ...
    sprintf('m \\in \\{%s\\}', strjoin(string(mList), ', ')), ...
    sprintf('tol = %.0e', tol)};
for j = 1:nk
    infoLeft{end+1} = sprintf('e(k=%d, m=%d) = %.2e', ...
        kList(j), mList(end), E(end,j));   %#ok<SAGROW>
end

text(ax1, 0.03, 0.97, infoLeft, ...
    'Units','normalized', ...
    'VerticalAlignment','top', 'HorizontalAlignment','left', ...
    'FontSize', 9, ...
    'BackgroundColor', [1 1 1 0.85], ...
    'EdgeColor', [0.5 0.5 0.5], 'Margin', 4, ...
    'Interpreter','tex');

% ------------------------------------------------------------
% Panel 2: pointwise error at t = T
% ------------------------------------------------------------
ax2 = subplot(1,3,2);
semilogy(ax2, x, err, 'b-', 'LineWidth', 2);
set(ax2, 'YScale','log');
grid(ax2,'on');
xlabel(ax2,'x'); ylabel(ax2,'|u_h(T,x) - u(T,x)|');
ax2.FontSize = 12;

[emax, imax] = max(err);
infoRight = { ...
    sprintf('Finest: k = %d, m = %d', kFinest, mFinest), ...
    sprintf('||u_h - u||_\\infty = %.2e', emax), ...
    sprintf('at x = %.3f', x(imax))};

text(ax2, 0.03, 0.97, infoRight, ...
    'Units','normalized', ...
    'VerticalAlignment','top', 'HorizontalAlignment','left', ...
    'FontSize', 9, ...
    'BackgroundColor', [1 1 1 0.85], ...
    'EdgeColor', [0.5 0.5 0.5], 'Margin', 4, ...
    'Interpreter','tex');

% ------------------------------------------------------------
% Panel 3: 3D surface of log10|error| over (t,x)
% ------------------------------------------------------------
ax3 = subplot(1,3,3);

Z = log10(ErrST);

surf(ax3, xc, tGrid, Z, 'EdgeColor','none');
hold(ax3,'on');

% خط t = T روی سطح (اختیاری)
plot3(ax3, xc, T*ones(size(xc)), Z(end,:), 'r', 'LineWidth', 1.5);

shading(ax3, 'interp');
colormap(ax3, parula);
cb = colorbar(ax3);
cb.Label.String = 'log_{10}|u_h - u|';
cb.Label.Interpreter = 'tex';

view(ax3, 40, 30);
grid(ax3,'on');
xlabel(ax3, 'x'); ylabel(ax3, 't'); zlabel(ax3, 'log_{10}|u_h - u|');
ax3.FontSize = 12;

% کادر اطلاعات گوشه
infoST = { ...
    sprintf('Finest: k=%d, m=%d', kFinest, mFinest), ...
    sprintf('max|u_h-u| = %.2e', emax)};
text(ax3, -0.05, 0.97, infoST, ...
    'Units','normalized', ...
    'VerticalAlignment','top', 'HorizontalAlignment','left', ...
    'FontSize', 9, ...
    'BackgroundColor', [1 1 1 0.85], ...
    'EdgeColor', [0.5 0.5 0.5], 'Margin', 4, ...
    'Interpreter','tex');


%% ============================================================
%%  Solver (Newton-Cotes)
%% ============================================================
function out = solveWaveNC(Ffun,Lfun,uExactFun,T,k,nt,Xcore,tol,maxIter)
    tGrid = linspace(0,T,nt+1);
    h     = tGrid(2)-tGrid(1);
    pad   = T + 4*h;
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

function str = texSci(x)
% 1.234e-05  ->  $1.234\times10^{-5}$
    e = floor(log10(abs(x)));
    mant = x/10^e;
    str = sprintf('$%.3f\\times10^{%d}$', mant, e);
end