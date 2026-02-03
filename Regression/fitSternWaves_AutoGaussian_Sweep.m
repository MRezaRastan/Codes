function fitSternWaves_AutoGaussian_Sweep(inputFile, outputBase)
% fitSternWaves_AutoGaussian_Sweep
% Global Gaussian model sweep with Fr-dependent polynomials.
%
% Excel format expected: paired columns:
%   Fr=0.75_x, Fr=0.75_y, Fr=1.00_x, Fr=1.00_y, ...
% or using underscore instead of '=': Fr_0.75_x, Fr_0.75_y, ...
%
% Usage:
%   fitSternWaves_AutoGaussian_Sweep('stern_wave_data_low.xlsx','gaussian_model')
%   -> outputs gaussian_model.tex and gaussian_model.xlsx

%% ---------------- Parameters you may tweak ----------------
N_gauss_range = 9:9;    % number of Gaussians to sweep
polyOrder_range = 1:2;   % polynomial order in Fr to sweep
xGridCount = 200;        % points for smooth profiles (per Fr) in the output file/plots
sigmaFrac = 0.25;        % initial sigma fraction of global x-span
rng(1);                  % deterministic init for reproducibility
%% ----------------------------------------------------------

%% --- Read Excel with preserved headers ---
T = readtable(inputFile, 'VariableNamingRule','preserve');
colNames = T.Properties.VariableNames;

% Detect Fr list from headers like Fr=1.5_x or Fr_1.5_x
Fr_list = [];
xCols = {};
yCols = {};
for i = 1:numel(colNames)
    nm = colNames{i};
    tok = regexp(nm, '^Fr[=_]([\d\.]+)_x$', 'tokens');
    if ~isempty(tok)
        Fr_list(end+1) = str2double(tok{1}{1}); %#ok<AGROW>
        xCols{end+1} = nm; %#ok<AGROW>
        % find matching y
        yname = regexprep(nm, '_x$', '_y');
        if ~ismember(yname, colNames)
            error('Missing y column for "%s". Expected "%s".', nm, yname);
        end
        yCols{end+1} = yname; %#ok<AGROW>
    end
end
nFr = numel(Fr_list);
if nFr==0
    error('No Froude numbers detected. Check Excel headers like "Fr=1.50_x".');
end

% Sort by increasing Fr (keep pairs aligned)
[Fr_list, sortIdx] = sort(Fr_list(:)');
xCols = xCols(sortIdx);
yCols = yCols(sortIdx);

fprintf('Found Froude numbers: %s\n', num2str(Fr_list));

% Build per-Fr data cell and global dataset
dataCell = cell(nFr,1);
X_all = []; y_all = []; Fr_all = [];

for iFr = 1:nFr
    xCol = T.(xCols{iFr});
    yCol = T.(yCols{iFr});
    valid = ~isnan(xCol) & ~isnan(yCol);
    Xi = xCol(valid);
    Yi = yCol(valid);
    dataCell{iFr} = [Xi, Yi];
    X_all = [X_all; Xi]; %#ok<AGROW>
    y_all = [y_all; Yi]; %#ok<AGROW>
    Fr_all = [Fr_all; Fr_list(iFr)*ones(numel(Xi),1)]; %#ok<AGROW>
end

% Global dataset [X, Fr] -> y
XFr_all = [X_all, Fr_all];

% Useful global stats for initial guesses
x_min = min(X_all);
x_max = max(X_all);
x_span = x_max - x_min;
x_mean = mean(X_all);
y_mean = mean(y_all);
y_std  = std(y_all);
if ~isfinite(y_std) || y_std==0, y_std = max(1, std(y_all(~isnan(y_all)))); end

%% --- Sweep models ---
resultsAll = struct('N_gauss',{},'polyOrder',{},'rmse',{},'pFit',{},'fun',{});
rowID = 0;

for N_gauss = N_gauss_range
    for polyOrder = polyOrder_range
        rowID = rowID + 1;

        [modelFun, p0] = makeGaussianModel(N_gauss, polyOrder, x_mean, x_span, y_mean, y_std, sigmaFrac);

        opts = optimoptions('lsqcurvefit', ...
                            'Display','off', ...
                            'MaxIterations', 10000, ...
                            'MaxFunctionEvaluations', 200000);

        % Use robust try/catch; if fails, mark rmse Inf
        try
            pFit = lsqcurvefit(@(p,XFr) modelFun(p, XFr), p0, XFr_all, y_all, [], [], opts);

            y_pred = modelFun(pFit, XFr_all);
            if any(~isfinite(y_pred))
                rmse = Inf; pFit = [];
            else
                rmse = sqrt(mean((y_all - y_pred).^2));
            end
        catch
            rmse = Inf; pFit = [];
        end

        resultsAll(rowID).N_gauss = N_gauss;
        resultsAll(rowID).polyOrder = polyOrder;
        resultsAll(rowID).rmse = rmse;
        resultsAll(rowID).pFit = pFit;
        resultsAll(rowID).fun = @(Z) modelFun(pFit, Z); % freeze params for plotting
    end
end

% Guard: ensure at least one finite RMSE
allRMSE = [resultsAll.rmse];
if all(~isfinite(allRMSE))
    error('All fits failed (RMSE=Inf). Try reducing N_gauss_range/polyOrder_range or adjust sigmaFrac.');
end

%% --- Select best ---
[~, bestIdx] = min(allRMSE);
bestResult = resultsAll(bestIdx);

fprintf('\n===== Best Fit Summary =====\n');
fprintf('N_gauss: %d\n', bestResult.N_gauss);
fprintf('PolyOrder: %d\n', bestResult.polyOrder);
fprintf('RMSE: %.6f\n', bestResult.rmse);

%% --- Equation (LaTeX + Excel) ---
[pStruct, latexEq, excelFormula] = getGaussianEquation(bestResult.pFit, bestResult.N_gauss, bestResult.polyOrder);

% Write LaTeX
texFile = sprintf('%s.tex', outputBase);
fid = fopen(texFile, 'w'); fprintf(fid, '%s\n', latexEq); fclose(fid);
fprintf('LaTeX equation saved to %s\n', texFile);

% Print Excel formula preview
fprintf('\nExcel Formula (put in C2 where A2=Fr, B2=x):\n=%s\n', excelFormula);

%% --- RMSE table (top 15) ---
RMSE_table = struct2table(resultsAll);
RMSE_table = sortrows(RMSE_table, 'rmse');
fprintf('\n===== RMSE Table (top 15) =====\n');
disp(RMSE_table(1:min(15,height(RMSE_table)),:));

%% --- Plot best fit vs. data ---
figure('Name','Best Gaussian Fit vs Data','Color','w');
tiledlayout(nFr,1,'TileSpacing','compact','Padding','compact');
for iFr = 1:nFr
    Xi = dataCell{iFr}(:,1);
    Yi = dataCell{iFr}(:,2);
    % predict on Xi
    Yi_fit = bestResult.fun([Xi, Fr_list(iFr)*ones(size(Xi))]);
    nexttile
    plot(Xi, Yi, 'o', 'DisplayName', sprintf('Fr=%.3g data', Fr_list(iFr)));
    hold on
    plot(Xi, Yi_fit, '-', 'LineWidth', 1.5, 'DisplayName', 'Gaussian Fit');
    legend('Location','best'); xlabel('x'); ylabel('y'); grid on;
    title(sprintf('Fr=%.3g', Fr_list(iFr)));
end

%% --- Save fitted profiles + equation to Excel ---
% Common x grid (global)
x_common = linspace(x_min, x_max, xGridCount)';

% Build LES interpolations and model fits on the common grid
y_les_interp = nan(xGridCount, nFr);
y_fit_interp = nan(xGridCount, nFr);

for iFr = 1:nFr
    Xi = dataCell{iFr}(:,1);
    Yi = dataCell{iFr}(:,2);

    % LES interpolation (linear with extrap)
    y_les_interp(:,iFr) = interp1(Xi, Yi, x_common, 'linear', 'extrap');

    % Model prediction
    y_fit_interp(:,iFr) = bestResult.fun([x_common, Fr_list(iFr)*ones(size(x_common))]);
end

% Assemble table: x | LES_Fr=... | Fit_Fr=...
lesNames = strcat("LES_Fr=", string(Fr_list));
fitNames = strcat("Fit_Fr=", string(Fr_list));
varNames = ["x", lesNames, fitNames];
result_table = array2table([x_common, y_les_interp, y_fit_interp], 'VariableNames', varNames);

xlsxFile = sprintf('%s.xlsx', outputBase);
writetable(result_table, xlsxFile, 'Sheet', 'Profiles', 'Range', 'A1');

% Sheet 2: Equation and Excel formula
writecell( {'Fitted Equation (LaTeX)', latexEq; 'Excel Formula (use in C2; A2=Fr, B2=x)', ['=' excelFormula]}, ...
           xlsxFile, 'Sheet', 'Equation', 'Range', 'A1');

fprintf('Results saved to %s\n', xlsxFile);

end % ====== main ======

%% ==================== Model Builder ====================
function [modelFun, p0] = makeGaussianModel(N_gauss, polyOrder, x_mean, x_span, y_mean, y_std, sigmaFrac)
% Parameter layout:
% p = [ C0_coeffs(nCoeff),
%       A1_coeffs(nCoeff), mu1_coeffs(nCoeff), sigma1_coeffs(nCoeff),
%       ...,
%       AN_coeffs(nCoeff), muN_coeffs(nCoeff), sigmaN_coeffs(nCoeff) ]
nCoeff = polyOrder + 1;
nParams = (1 + 3*N_gauss) * nCoeff; %#ok<NASGU> % for info

% Initial guesses (only constant term non-zero)
p0 = zeros((1 + 3*N_gauss)*nCoeff,1);
idx = 1;

% C0(Fr) ~ mean(y)
c0_const = y_mean;
p0(idx) = c0_const; idx = idx + nCoeff;

% Place Gaussian centers across x-range, amplitudes ~ y_std/N, sigma ~ fraction of span
if x_span <= 0, x_span = 1; end
sig0 = max(1e-3, sigmaFrac * x_span);
for n = 1:N_gauss
    mu0 = x_mean + ( (n - (N_gauss+1)/2) / max(N_gauss-1,1) ) * (0.5*x_span);
    A0  = y_std / max(N_gauss,1);
    % A_n(Fr)
    p0(idx) = A0; idx = idx + nCoeff;
    % mu_n(Fr)
    p0(idx) = mu0; idx = idx + nCoeff;
    % sigma_n(Fr) as exp(poly), so store log(sigma) in constant term ~ log(sig0)
    p0(idx) = log(max(sig0,1e-3)); idx = idx + nCoeff;
end

% Model function: sigma(Fr) = exp(poly(Fr)) ensures positivity
modelFun = @(p, XFr) modelGaussianSum(XFr(:,1), XFr(:,2), p, N_gauss, polyOrder);
end

%% ==================== Gaussian Model ====================
function y = modelGaussianSum(X, Fr, p, N_gauss, polyOrder)
% X, Fr are column vectors
X = X(:)'; Fr = Fr(:)'; nCoeff = polyOrder + 1; idx = 1;

% C0(Fr)
c0 = p(idx:idx+nCoeff-1); idx = idx + nCoeff;
C0 = polyval(flip(c0), Fr);

y = C0;
for n = 1:N_gauss
    a  = p(idx:idx+nCoeff-1); idx = idx + nCoeff;      % A_n(Fr)
    mu = p(idx:idx+nCoeff-1); idx = idx + nCoeff;      % mu_n(Fr)
    sg = p(idx:idx+nCoeff-1); idx = idx + nCoeff;      % log_sigma_n(Fr)

    AFr  = polyval(flip(a),  Fr);
    muFr = polyval(flip(mu), Fr);
    sigFr = exp(polyval(flip(sg), Fr));                % strictly positive

    y = y + AFr .* exp(-((X - muFr)./sigFr).^2);
end

y = y(:);
end

%% ==================== Equation Builder ====================
function [pStruct, latexEq, excelFormula] = getGaussianEquation(p, N_gauss, polyOrder)
nCoeff = polyOrder + 1; idx = 1; pStruct = struct();

% Extract coeff vectors (highest power first inside pStruct for readability)
pStruct.C0    = flip(p(idx:idx+nCoeff-1)); idx = idx + nCoeff;
for n=1:N_gauss
    pStruct.(['A' num2str(n)])     = flip(p(idx:idx+nCoeff-1)); idx = idx + nCoeff;
    pStruct.(['mu' num2str(n)])    = flip(p(idx:idx+nCoeff-1)); idx = idx + nCoeff;
    pStruct.(['logsigma' num2str(n)]) = flip(p(idx:idx+nCoeff-1)); idx = idx + nCoeff;
end

% Build LaTeX: sigma(Fr) = exp(poly(Fr))
latexEq = sprintf('y(Fr,x) = %s', poly2latex(pStruct.C0,'Fr'));
for n=1:N_gauss
    Atex  = poly2latex(pStruct.(['A' num2str(n)]),'Fr');
    mtex  = poly2latex(pStruct.(['mu' num2str(n)]),'Fr');
    lsTex = poly2latex(pStruct.(['logsigma' num2str(n)]),'Fr'); % log sigma
    stex  = sprintf('\\exp\\left(%s\\right)', lsTex);           % sigma(Fr) = exp(poly)
    term = sprintf(' + (%s)\\exp\\left(-\\left(\\frac{x-(%s)}{%s}\\right)^{2}\\right)', Atex, mtex, stex);
    latexEq = [latexEq, term]; %#ok<AGROW>
end

% Build Excel formula: Fr -> A2, x -> B2; sigma = EXP(poly_excel(A2))
excelFormula = poly2excel(pStruct.C0, 'A2');
for n=1:N_gauss
    Axl  = poly2excel(pStruct.(['A' num2str(n)]),'A2');
    mxl  = poly2excel(pStruct.(['mu' num2str(n)]),'A2');
    lsql = poly2excel(pStruct.(['logsigma' num2str(n)]),'A2');  % log sigma
    sxl  = sprintf('EXP(%s)', lsql);                            % sigma(Fr)
    term = sprintf(' + (%s)*EXP(-((B2-(%s))/(%s))^2)', Axl, mxl, sxl);
    excelFormula = [excelFormula, term]; %#ok<AGROW>
end
end

%% ==================== Utility: Polynomial to LaTeX ====================
function s = poly2latex(coeffs, varname)
coeffs = coeffs(:)'; deg = numel(coeffs)-1; parts = {};
for i=1:numel(coeffs)
    c = coeffs(i); pow = deg-(i-1);
    if abs(c) < 1e-12, continue; end
    if pow == 0
        parts{end+1} = sprintf('%.6g', c); %#ok<AGROW>
    elseif pow == 1
        parts{end+1} = sprintf('%.6g\\,%s', c, varname); %#ok<AGROW>
    else
        parts{end+1} = sprintf('%.6g\\,%s^{%d}', c, varname, pow); %#ok<AGROW>
    end
end
if isempty(parts), s = '0'; else, s = strjoin(parts,' + '); end
end

%% ==================== Utility: Polynomial to Excel ====================
function s = poly2excel(coeffs, varname)
coeffs = coeffs(:)'; deg = numel(coeffs)-1; parts = {};
for i=1:numel(coeffs)
    c = coeffs(i); pow = deg-(i-1);
    if abs(c) < 1e-12, continue; end
    if pow == 0
        parts{end+1} = sprintf('%.6g', c); %#ok<AGROW>
    elseif pow == 1
        parts{end+1} = sprintf('%.6g*%s', c, varname); %#ok<AGROW>
    else
        parts{end+1} = sprintf('%.6g*%s^%d', c, varname, pow); %#ok<AGROW>
    end
end
if isempty(parts), s = '0'; else, s = strjoin(parts,'+'); end
end
