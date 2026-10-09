function R = run_mri_synthetic(options)
%RUN_MRI_SYNTHETIC Controlled test of the nonlinear/varying-kernel theory.
% The problem is posed in R^n with known solution xStar, nonlinear kernels
% M_k, alternating kernel parameters, and a strongly monotone A+B.
%
% R = run_mri_synthetic();
% R = run_mri_synthetic(mri_options());
%
% Output:
%   <outputRoot>/results_synthetic/synthetic_summary.csv
%   <outputRoot>/results_synthetic/synthetic.mat
%   <outputRoot>/results_synthetic/figures/synthetic_linear_rate.pdf
%   <outputRoot>/results_synthetic/figures/synthetic_linear_rate.png
% Figures use a white background; PDF is vector, PNG is 600 DPI.
% Four distinct line styles and staggered markers distinguish the curves.
%
% This experiment is deterministic and independent of the image data.
% Audit 2026-10-03: enforce final Newton defect, record termination status,
% distinguish graph and formula residuals, and explicitly set a log Y axis.

if nargin < 1
    options = struct();
end
setup_mri();
opts = mri_options(options);
outDir = fullfile(opts.outputRoot,'results_synthetic');
figDir = fullfile(outDir,'figures');
if ~exist(outDir,'dir'), mkdir(outDir); end
if ~exist(figDir,'dir'), mkdir(figDir); end

% Problem data: A(x)=a(x-xStar), B(x)=b(x-xStar).
n = 1000;
a = 1;
b = 1;
lambdaCocoercive = 1 / b;
S_lower = 1;
Lkernel = 0.6;
stepBound = 2 * lambdaCocoercive * S_lower * (1-Lkernel^2);

idx = (1:n)';
xStar = sin(2*pi*idx/n) + 0.5*cos(6*pi*idx/n);
x0 = zeros(n,1);

% Nonlinear varying kernels: M_k x = mu_k^{-1}(x+ell_k*tanh(x)).
% The two parameter pairs alternate forever, so {M_k} does not converge.
muEven = 0.8; ellEven = 0.2;
muOdd  = 1.0; ellOdd  = 0.6;
assert(max(muEven,muOdd) <= stepBound + 100*eps(stepBound), ...
    'MRI:SyntheticStep','Synthetic stepsize exceeds the theorem bound.');

rho = 0.3;
alphaList = [0 0.5 1 2];
assert(all(rho*(1+alphaList) < 1), ...
    'MRI:SyntheticInertia','Synthetic inertia-relaxation condition failed.');

maxIter = 500;
tol = 1e-10;
newtonTol = 1e-13;
newtonMaxIter = 50;

nA = numel(alphaList);
iterations = zeros(nA,1);
finalRelError = zeros(nA,1);
finalRelResidual = zeros(nA,1);
empiricalFactor = NaN(nA,1);
maxInnerResidual = zeros(nA,1);
histories = cell(nA,1);
converged = false(nA,1);
status = repmat({'maxiter'},nA,1);
maxResidualMismatch = zeros(nA,1);
maxBackwardPointErrorBound = zeros(nA,1);
fitFirstIteration = NaN(nA,1);
fitLastIteration = NaN(nA,1);
fitPoints = zeros(nA,1);

xStarNorm = max(1,norm(xStar));

for ia = 1:nA
    alpha = alphaList(ia);
    xPrev = x0;
    x = x0;
    relErr = zeros(maxIter+1,1);
    relRes = NaN(maxIter+1,1);
    muHist = NaN(maxIter,1);
    ellHist = NaN(maxIter,1);
    innerHist = NaN(maxIter,1);
    mismatchHist = NaN(maxIter,1);
    equationCheckHist = NaN(maxIter,1);
    backwardBoundHist = NaN(maxIter,1);
    newtonHist = NaN(maxIter,1);
    relErr(1) = norm(x-xStar)/xStarNorm;
    stopIter = maxIter;

    for k = 0:maxIter-1
        y = x + alpha*(x-xPrev);
        if mod(k,2)==0
            mu = muEven;
            ell = ellEven;
        else
            mu = muOdd;
            ell = ellOdd;
        end
        muHist(k+1) = mu;
        ellHist(k+1) = ell;

        % Backward equation, componentwise:
        % (1+mu*a)t + ell*tanh(t)
        %   = (1-mu*b)y + ell*tanh(y) + mu*(a+b)xStar.
        rhs = (1-mu*b)*y + ell*tanh(y) + mu*(a+b)*xStar;
        t = rhs/(1+mu*a);
        nNewton = 0;
        for it = 0:newtonMaxIter
            F = (1+mu*a)*t + ell*tanh(t) - rhs;
            innerResidual = norm(F,inf);
            if ~isfinite(innerResidual)
                error('MRI:SyntheticNewton','Nonfinite backward defect at k=%d.',k);
            end
            if innerResidual <= newtonTol
                break;
            end
            if it == newtonMaxIter
                break;
            end
            dF = (1+mu*a) + ell*(1-tanh(t).^2);
            t = t-F./dF;
            nNewton = nNewton+1;
        end
        % F is evaluated AT THE RETURNED t, not at the previous Newton point.
        innerHist(k+1) = innerResidual;
        newtonHist(k+1) = nNewton;
        if innerResidual > newtonTol
            error('MRI:SyntheticNewton', ...
                'Backward solve failed at alpha=%g, k=%d: defect=%g > %g.', ...
                alpha,k,innerResidual,newtonTol);
        end
        xTilde = t;

        % Corrected point y - mu(M_k y - M_k xTilde), with S=I.
        xHat = xTilde - ell*(tanh(y)-tanh(xTilde));
        xNext = (1-rho)*x + rho*xHat;

        if any(~isfinite(xNext)) || norm(xNext)>1e12*xStarNorm
            error('MRI:SyntheticNumerics','Nonfinite or unbounded iterate at k=%d.',k);
        end
        % Direct graph element: valid even at a numerical backward point.
        % For exact backward evaluation this equals the theorem residual.
        r = (a+b)*(xTilde-xStar);
        rFormula = (y-xTilde+ell*(tanh(y)-tanh(xTilde)))/mu ...
            + b*(xTilde-y);
        mismatchHist(k+1) = norm(rFormula-r,inf);
        % Algebraic identity: rFormula-r = -F/mu, up to roundoff.
        equationCheckHist(k+1) = norm(rFormula-r+F/mu,inf);
        scale = max([1,norm(rFormula,inf),norm(r,inf),norm(y,inf)]);
        assert(equationCheckHist(k+1)<=1e3*eps(scale), ...
            'MRI:SyntheticResidual','Backward/graph residual identity failed.');
        % Scalar derivative >=1+mu*a gives an infinity-norm inverse bound.
        backwardBoundHist(k+1) = innerResidual/(1+mu*a);
        relErr(k+2) = norm(xNext-xStar)/xStarNorm;
        relRes(k+2) = norm(r)/xStarNorm;

        xPrev = x;
        x = xNext;

        if relErr(k+2) <= tol && relRes(k+2) <= tol
            stopIter = k+1;
            converged(ia) = true;
            status{ia} = 'converged';
            break;
        end
    end

    relErr = relErr(1:stopIter+1);
    relRes = relRes(1:stopIter+1);
    muHist = muHist(1:stopIter);
    ellHist = ellHist(1:stopIter);
    innerHist = innerHist(1:stopIter);
    mismatchHist = mismatchHist(1:stopIter);
    equationCheckHist = equationCheckHist(1:stopIter);
    backwardBoundHist = backwardBoundHist(1:stopIter);
    newtonHist = newtonHist(1:stopIter);

    iterations(ia) = stopIter;
    finalRelError(ia) = relErr(end);
    finalRelResidual(ia) = relRes(end);
    maxInnerResidual(ia) = max(innerHist);
    maxResidualMismatch(ia) = max(mismatchHist);
    maxBackwardPointErrorBound(ia) = max(backwardBoundHist);
    [empiricalFactor(ia),fitFirstIteration(ia),fitLastIteration(ia), ...
        fitPoints(ia)] = local_empirical_factor(relErr);

    histories{ia} = struct('alpha',alpha,'rho',rho, ...
        'relativeError',relErr,'relativeResidual',relRes, ...
        'mu',muHist,'ell',ellHist,'innerResidual',innerHist, ...
        'residualMismatch',mismatchHist,'residualIdentityCheck',equationCheckHist, ...
        'backwardPointErrorBoundInf',backwardBoundHist,'newtonSteps',newtonHist, ...
        'xFinal',x,'xBackwardFinal',xTilde,'converged',converged(ia), ...
        'status',status{ia});
    fprintf('Synthetic: alpha=%g, rho=%g, k=%d, error=%.3e, residual=%.3e, %s\n', ...
        alpha,rho,stopIter,relErr(end),relRes(end),status{ia});
end

summary = table(alphaList(:),repmat(rho,nA,1),iterations,finalRelError, ...
    finalRelResidual,empiricalFactor,maxInnerResidual,converged,status, ...
    maxResidualMismatch,maxBackwardPointErrorBound,fitFirstIteration, ...
    fitLastIteration,fitPoints, ...
    'VariableNames',{'Alpha','Rho','Iterations','FinalRelativeError', ...
    'FinalRelativeResidual','EmpiricalLinearFactor','MaxBackwardResidual', ...
    'Converged','Status','MaxResidualMismatch','MaxBackwardPointErrorBoundInf', ...
    'FitFirstIteration','FitLastIteration','FitPoints'});
writetable(summary,fullfile(outDir,'synthetic_summary.csv'));

metadata = struct();
metadata.codeVersion = 'MRI-SYNTH-2026-10-03-audit1';
metadata.n = n;
metadata.a = a;
metadata.b = b;
metadata.xStarFormula = 'sin(2*pi*i/n)+0.5*cos(6*pi*i/n)';
metadata.muEven = muEven;
metadata.ellEven = ellEven;
metadata.muOdd = muOdd;
metadata.ellOdd = ellOdd;
metadata.Lkernel = Lkernel;
metadata.stepBound = stepBound;
metadata.rho = rho;
metadata.alphaList = alphaList;
metadata.maxIter = maxIter;
metadata.tol = tol;
metadata.newtonTol = newtonTol;
metadata.newtonMaxIter = newtonMaxIter;
metadata.backwardDefectNorm = 'Infinity norm of unscaled scalar equation defects';
metadata.graphResidual = '(a+b)*(xBackward-xStar)';
metadata.empiricalFit = ['Last at most 30 indices in second half of finite errors ' ...
    'above 1e-13; at least 8 fit points; exp(log-error OLS slope)'];
metadata.matlabVersion = version;
metadata.computer = computer;
save(fullfile(outDir,'synthetic.mat'),'summary','histories','metadata','xStar','-v7');

if opts.makePlots
    % Plot the computed histories without smoothing or resampling.
    % The axes, line styles and staggered markers match the reference plot.
    f = figure('Visible','off','Units','inches', ...
        'Position',[1 1 6.2 4.2],'Color','w', ...
        'Name','Synthetic convergence histories','NumberTitle','off');
    figureGuard = onCleanup(@()close(f));

    % Keep this figure light even when the MATLAB desktop uses a dark theme.
    if isprop(f,'Theme')
        theme(f,'light');
    end
    set(f,'Color','w');

    % Set YScale before hold on so that the vertical axis is logarithmic.
    ax = axes('Parent',f,'Units','normalized', ...
        'Position',[0.14 0.15 0.83 0.80], ...
        'Color','w','XColor','k','YColor','k', ...
        'XScale','linear','YScale','log', ...
        'FontName','Arial','FontSize',10,'FontWeight','normal', ...
        'LineWidth',0.75,'TickDir','out', ...
        'XMinorTick','off','YMinorTick','off');
    hold(ax,'on');
    box(ax,'on');

    % Blue, orange, green and red, respectively.
    curveColors = [31 119 180; 255 127 14; 44 160 44; 214 39 40]/255;
    lineStyles = {'-','--','-.',':'};
    markerShapes = {'o','s','^','d'};
    lineHandles = gobjects(nA,1);
    positiveErrors = [];
    lastIteration = 0;

    for ia = 1:nA
        h = histories{ia};
        errorHistory = h.relativeError(:);
        kPlot = (0:numel(errorHistory)-1)';
        lastIteration = max(lastIteration,kPlot(end));

        % Keep zero/nonfinite values out of the logarithmic display only.
        % The stored histories, CSV and MAT data are not modified.
        positive = isfinite(errorHistory) & errorHistory>0;
        positiveErrors = [positiveErrors; errorHistory(positive)]; %#ok<AGROW>
        plotError = errorHistory;
        plotError(~positive) = NaN;

        styleIndex = mod(ia-1,numel(lineStyles))+1;
        % Mark iterations 2,4,6,8 for the first four curves, then every 12.
        % Staggering the markers makes nearly coincident curves readable.
        firstMarker = min(3+2*(styleIndex-1),numel(kPlot));
        markerIndices = firstMarker:12:numel(kPlot);
        lineHandles(ia) = plot(ax,kPlot,plotError, ...
            'Color',curveColors(styleIndex,:), ...
            'LineStyle',lineStyles{styleIndex}, ...
            'Marker',markerShapes{styleIndex}, ...
            'MarkerIndices',markerIndices, ...
            'MarkerSize',4,'MarkerEdgeColor',curveColors(styleIndex,:), ...
            'MarkerFaceColor',curveColors(styleIndex,:), ...
            'LineWidth',1.3, ...
            'DisplayName',sprintf('\\alpha = %g',h.alpha));
    end

    xlabel(ax,'Iteration','FontName','Arial','FontSize',11, ...
        'FontWeight','normal','Color','k','Interpreter','tex');
    ylabel(ax,'Relative error','FontName','Arial','FontSize',11, ...
        'FontWeight','normal','Color','k','Interpreter','tex');

    % With the current data these limits are [0,65] and [3e-11,1.3].
    % Expand automatically if a later run contains a longer/larger history.
    xMax = max(1,lastIteration);
    xTickStep = max(10,10*ceil(xMax/100));
    xlim(ax,[0 xMax]);
    xticks(ax,unique([0:xTickStep:xMax xMax]));
    yLower = 3e-11;
    yUpper = 1.3;
    if ~isempty(positiveErrors)
        yLower = max(realmin,min(yLower,0.4*min(positiveErrors)));
        yUpper = max(yUpper,1.1*max(positiveErrors));
    end
    ylim(ax,[yLower yUpper]);
    firstExponent = 2*ceil(log10(yLower)/2);
    lastExponent = 2*floor(log10(yUpper)/2);
    exponentStep = max(2,2*ceil((lastExponent-firstExponent)/12));
    % MATLAB requires tick values to be in increasing order.
    yticks(ax,10.^(firstExponent:exponentStep:lastExponent));

    grid(ax,'on');
    set(ax,'GridColor','k','GridAlpha',0.2,'GridLineStyle','-', ...
        'XMinorGrid','off','YMinorGrid','off','Layer','top');
    legend(ax,lineHandles,'Location','southwest','Interpreter','tex', ...
        'FontName','Arial','FontSize',10, ...
        'Color','w','TextColor','k','EdgeColor','k', ...
        'Box','on','AutoUpdate','off');
    drawnow;

    % Export this same figure: vector PDF and 600-DPI PNG, both white.
    pdfFile = fullfile(figDir,'synthetic_linear_rate.pdf');
    pngFile = fullfile(figDir,'synthetic_linear_rate.png');
    exportgraphics(f,pdfFile, ...
        'ContentType','vector','BackgroundColor','white');
    exportgraphics(f,pngFile, ...
        'Resolution',600,'BackgroundColor','white');
    fprintf('Figures saved:\n%s\n%s\n',pdfFile,pngFile);
    clear figureGuard;
end

R = struct('summary',summary,'histories',{histories}, ...
    'metadata',metadata,'outDir',outDir);
fprintf('Synthetic experiment completed. Results: %s\n',outDir);
if ~all(converged)
    warning('MRI:SyntheticCap','Some runs reached the cap. Inspect the Status column.');
end
end

function [q,firstIteration,lastIteration,nPoints] = local_empirical_factor(err)
% Descriptive one-step factor; not the q in the theorem's Lyapunov bound.
err = err(:);
k = (0:numel(err)-1)';
valid = find(isfinite(err) & err>1e-13);
q=NaN; firstIteration=NaN; lastIteration=NaN; nPoints=0;
if isempty(valid), return; end
startPos = floor(numel(valid)/2)+1;
tail = valid(startPos:end);
if numel(tail)>30, tail=tail(end-29:end); end
nPoints=numel(tail);
if nPoints<8, return; end
firstIteration=k(tail(1)); lastIteration=k(tail(end));
p = polyfit(k(tail),log(err(tail)),1);
q = exp(p(1));
end
