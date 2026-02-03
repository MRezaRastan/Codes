function fitSternWaves_AutoFourier_Sweep(inputFile, outputBase, N_range, polyOrder_range, xGridCount)
% fitSternWaves_AutoFourier_Sweep
% Sweeps Fourier-model families (cosine, sine, full) across numbers of harmonics
% and polynomial orders for Fr-dependent coefficients. Exports:
%   - outputBase.tex  (LaTeX equation)
%   - outputBase.xlsx (Profiles sheet + Equation sheet with Excel formula)
%
% Usage:
%   fitSternWaves_AutoFourier_Sweep('stern_wave_data.xlsx','fourier_model')
%
% Optional args:
%   N_range        - vector of N_harmonics to try (default 1:10)
%   polyOrder_range- vector of polynomial orders in Fr (default 0:4)
%   xGridCount     - number of x points in output file (default 200)

if nargin < 3 || isempty(N_range), N_range = 10:10; end
if nargin < 4 || isempty(polyOrder_range), polyOrder_range = 1:2; end
if nargin < 5 || isempty(xGridCount), xGridCount = 200; end

rng(1); % reproducible initial guesses

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

% Useful global stats for plotting/export
x_min = min(X_all);
x_max = max(X_all);
x_span = x_max - x_min;

%% --- Sweep models (cosine, sine, full) ---
modelTypes = {'cosine','sine','full'};
resultsAll = struct('type',{},'N_harmonics',{},'polyOrder',{},'rmse',{},'pFit',{},'fun',{});
rowID = 0;

opts = optimoptions('lsqcurvefit', ...
                    'Display','off', ...
                    'MaxIterations', 10000, ...
                    'MaxFunctionEvaluations', 200000);

for N_harmonics = N_range
    for polyOrder = polyOrder_range
        for mType = modelTypes
            type = mType{1};
            rowID = rowID + 1;

            [modelFun, p0, nParams] = makeModelFunction(type, N_harmonics, polyOrder);

            try
                pFit = lsqcurvefit(@(pp,XFr) modelFun(pp,XFr), p0, XFr_all, y_all, [], [], opts);
                y_pred = modelFun(pFit, XFr_all);
                if any(~isfinite(y_pred))
                    rmse = Inf; pFit = [];
                else
                    rmse = sqrt(mean((y_all - y_pred).^2));
                end
            catch
                rmse = Inf;
                pFit = [];
            end

            resultsAll(rowID).type = type;
            resultsAll(rowID).N_harmonics = N_harmonics;
            resultsAll(rowID).polyOrder = polyOrder;
            resultsAll(rowID).rmse = rmse;
            resultsAll(rowID).pFit = pFit;
            if ~isempty(pFit)
                % store a frozen function for later evaluation
                resultsAll(rowID).fun = @(Z) modelFun(pFit, Z);
            else
                resultsAll(rowID).fun = @(Z) nan(size(Z,1),1);
            end
        end
    end
end

% Guard: ensure at least one finite RMSE
allRMSE = [resultsAll.rmse];
if all(~isfinite(allRMSE))
    error('All fits failed (RMSE=Inf). Try reducing N_range/polyOrder_range or adjust data.');
end

%% --- Select best ---
[~, bestIdx] = min(allRMSE);
bestResult = resultsAll(bestIdx);

fprintf('\n===== Best Fit Summary =====\n');
fprintf('Model: %s\n', bestResult.type);
fprintf('N_harmonics: %d\n', bestResult.N_harmonics);
fprintf('PolyOrder: %d\n', bestResult.polyOrder);
fprintf('RMSE: %.6f\n', bestResult.rmse);

%% --- Equation (LaTeX + Excel) ---
[pStruct, latexEq, excelFormula] = getFourierEquation(bestResult.pFit, bestResult.type, bestResult.N_harmonics, bestResult.polyOrder);

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
figure('Name','Best Fourier Fit vs Data','Color','w');
tiledlayout(nFr,1,'TileSpacing','compact','Padding','compact');
for iFr = 1:nFr
    Xi = dataCell{iFr}(:,1);
    Yi = dataCell{iFr}(:,2);
    % predict on Xi
    Yi_fit = bestResult.fun([Xi, Fr_list(iFr)*ones(size(Xi))]);
    nexttile
    plot(Xi, Yi, 'o', 'DisplayName', sprintf('Fr=%.3g data', Fr_list(iFr)));
    hold on
    plot(Xi, Yi_fit, '-', 'LineWidth', 1.5, 'DisplayName', 'Fourier Fit');
    legend('Location','best'); xlabel('x'); ylabel('y'); grid on;
    title(sprintf('%s model, Fr=%.3g', bestResult.type, Fr_list(iFr)));
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

%% ==================== Model builder =====================================
function [modelFun, p0, nParams] = makeModelFunction(type, N_harmonics, polyOrder)
% Build model function that receives p and XFr and returns y
nCoeff = polyOrder + 1;
switch type
    case 'cosine'
        % A_n (N), phi_n (N), k (1), C0 (1) -> total (2N + 2) * nCoeff
        nParams = (2*N_harmonics + 2) * nCoeff;
        modelFun = @(p,XFr) modelFourier(XFr(:,1), XFr(:,2), p, N_harmonics, polyOrder, 'cosine');
    case 'sine'
        nParams = (2*N_harmonics + 2) * nCoeff;
        modelFun = @(p,XFr) modelFourier(XFr(:,1), XFr(:,2), p, N_harmonics, polyOrder, 'sine');
    case 'full'
        % A_n (N), B_n (N), k (1), C0 (1) -> (2N + 2)*nCoeff (same count)
        nParams = (2*N_harmonics + 2) * nCoeff;
        modelFun = @(p,XFr) modelFourier(XFr(:,1), XFr(:,2), p, N_harmonics, polyOrder, 'full');
    otherwise
        error('Unknown model type: %s', type);
end
p0 = zeros(nParams,1);
% simple initial guess: constant term 1 for first nCoeff entries
p0(1:nCoeff) = 1;
end

%% ==================== General Fourier model ==============================
function y = modelFourier(X, Fr, p, N_harmonics, polyOrder, type)
% X, Fr are column vectors of same length
X = X(:)'; Fr = Fr(:)'; nCoeff = polyOrder + 1; idx = 1;

switch type
    case {'cosine','sine'}
        % A_n
        A = zeros(N_harmonics, length(Fr));
        for n=1:N_harmonics
            coeffs = p(idx:idx+nCoeff-1); idx = idx+nCoeff;
            A(n,:) = polyval(flip(coeffs), Fr);
        end
        % phi_n
        phi = zeros(N_harmonics, length(Fr));
        for n=1:N_harmonics
            coeffs = p(idx:idx+nCoeff-1); idx = idx+nCoeff;
            phi(n,:) = polyval(flip(coeffs), Fr);
        end
    case 'full'
        % A_n
        A = zeros(N_harmonics, length(Fr));
        for n=1:N_harmonics
            coeffs = p(idx:idx+nCoeff-1); idx = idx+nCoeff;
            A(n,:) = polyval(flip(coeffs), Fr);
        end
        % B_n
        B = zeros(N_harmonics, length(Fr));
        for n=1:N_harmonics
            coeffs = p(idx:idx+nCoeff-1); idx = idx+nCoeff;
            B(n,:) = polyval(flip(coeffs), Fr);
        end
end

% k(Fr)
coeffs = p(idx:idx+nCoeff-1); idx = idx+nCoeff;
k = polyval(flip(coeffs), Fr);

% C0(Fr)
coeffs = p(idx:idx+nCoeff-1); idx = idx+nCoeff;
C0 = polyval(flip(coeffs), Fr);

switch type
    case 'cosine'
        % y = C0 + sum_n A_n(Fr)*cos(n * k(Fr) * X + phi_n(Fr))
        y = C0 + sum(A .* cos((1:N_harmonics)'.*(k.*X) + phi), 1);
    case 'sine'
        y = C0 + sum(A .* sin((1:N_harmonics)'.*(k.*X) + phi), 1);
    case 'full'
        y = C0 + sum(A .* cos((1:N_harmonics)'.*(k.*X)) + B .* sin((1:N_harmonics)'.*(k.*X)), 1);
end
y = y(:);
end

%% ==================== Equation Builder ===================================
function [pStruct, latexEq, excelFormula] = getFourierEquation(p, type, N_harmonics, polyOrder)
% Extract polynomial coefficient vectors (highest-power-first inside pStruct)
nCoeff = polyOrder+1; idx = 1; pStruct = struct();

switch type
    case {'cosine','sine'}
        for n=1:N_harmonics
            pStruct.(['A' num2str(n)]) = flip(p(idx:idx+nCoeff-1)); idx = idx + nCoeff;
        end
        for n=1:N_harmonics
            pStruct.(['phi' num2str(n)]) = flip(p(idx:idx+nCoeff-1)); idx = idx + nCoeff;
        end
    case 'full'
        for n=1:N_harmonics
            pStruct.(['A' num2str(n)]) = flip(p(idx:idx+nCoeff-1)); idx = idx + nCoeff;
        end
        for n=1:N_harmonics
            pStruct.(['B' num2str(n)]) = flip(p(idx:idx+nCoeff-1)); idx = idx + nCoeff;
        end
end
pStruct.k = flip(p(idx:idx+nCoeff-1)); idx = idx + nCoeff;
pStruct.C0 = flip(p(idx:idx+nCoeff-1));

% Build LaTeX
latexEq = sprintf('y(Fr,x) = %s', poly2latex(pStruct.C0,'Fr'));
switch type
    case 'cosine'
        for n=1:N_harmonics
            Atex = poly2latex(pStruct.(['A' num2str(n)]),'Fr');
            ktex = poly2latex(pStruct.k,'Fr');
            phtex = poly2latex(pStruct.(['phi' num2str(n)]),'Fr');
            latexEq = sprintf('%s + (%s)\\cos\\left(%d\\cdot(%s) x + (%s)\\right)', ...
                latexEq, Atex, n, ktex, phtex);
        end
    case 'sine'
        for n=1:N_harmonics
            Atex = poly2latex(pStruct.(['A' num2str(n)]),'Fr');
            ktex = poly2latex(pStruct.k,'Fr');
            phtex = poly2latex(pStruct.(['phi' num2str(n)]),'Fr');
            latexEq = sprintf('%s + (%s)\\sin\\left(%d\\cdot(%s) x + (%s)\\right)', ...
                latexEq, Atex, n, ktex, phtex);
        end
    case 'full'
        for n=1:N_harmonics
            Atex = poly2latex(pStruct.(['A' num2str(n)]),'Fr');
            Btex = poly2latex(pStruct.(['B' num2str(n)]),'Fr');
            ktex = poly2latex(pStruct.k,'Fr');
            latexEq = sprintf('%s + (%s)\\cos\\left(%d\\cdot(%s) x\\right) + (%s)\\sin\\left(%d\\cdot(%s) x\\right)', ...
                latexEq, Atex, n, ktex, Btex, n, ktex);
        end
end

% Build Excel formula: Fr -> A2, x -> B2; use Excel COS and SIN
% Start with C0
excelFormula = poly2excel(pStruct.C0, 'A2');
switch type
    case 'cosine'
        for n=1:N_harmonics
            Axl  = poly2excel(pStruct.(['A' num2str(n)]),'A2');     % amplitude
            kxl  = poly2excel(pStruct.k,'A2');                      % k(Fr)
            phxl = poly2excel(pStruct.(['phi' num2str(n)]),'A2');   % phi(Fr)
            % term = (A(Fr))*COS(n*(k(Fr))*B2 + (phi(Fr)))
            term = sprintf(' + (%s)*COS(%d*(%s)*B2 + (%s))', Axl, n, kxl, phxl);
            excelFormula = [excelFormula, term]; %#ok<AGROW>
        end
    case 'sine'
        for n=1:N_harmonics
            Axl  = poly2excel(pStruct.(['A' num2str(n)]),'A2');
            kxl  = poly2excel(pStruct.k,'A2');
            phxl = poly2excel(pStruct.(['phi' num2str(n)]),'A2');
            term = sprintf(' + (%s)*SIN(%d*(%s)*B2 + (%s))', Axl, n, kxl, phxl);
            excelFormula = [excelFormula, term]; %#ok<AGROW>
        end
    case 'full'
        for n=1:N_harmonics
            Axl  = poly2excel(pStruct.(['A' num2str(n)]),'A2');
            Bxl  = poly2excel(pStruct.(['B' num2str(n)]),'A2');
            kxl  = poly2excel(pStruct.k,'A2');
            term = sprintf(' + (%s)*COS(%d*(%s)*B2) + (%s)*SIN(%d*(%s)*B2)', Axl, n, kxl, Bxl, n, kxl);
            excelFormula = [excelFormula, term]; %#ok<AGROW>
        end
end
end

%% ==================== Utilities: poly -> LaTeX / Excel ====================
function s = poly2latex(coeffs, varname)
coeffs = coeffs(:)'; deg = numel(coeffs)-1; terms = {};
for i=1:numel(coeffs)
    c = coeffs(i); pow = deg-(i-1);
    if abs(c) < 1e-12, continue; end
    if pow == 0
        terms{end+1} = sprintf('%.6g', c);
    elseif pow == 1
        terms{end+1} = sprintf('%.6g\\,%s', c, varname);
    else
        terms{end+1} = sprintf('%.6g\\,%s^{%d}', c, varname, pow);
    end
end
if isempty(terms), s = '0'; else, s = strjoin(terms, ' + '); end
end

function s = poly2excel(coeffs, varname)
% coeffs highest-power-first
coeffs = coeffs(:)'; deg = numel(coeffs)-1; parts = {};
for i=1:numel(coeffs)
    c = coeffs(i); pow = deg-(i-1);
    if abs(c) < 1e-12, continue; end
    if pow == 0
        parts{end+1} = sprintf('%.6g', c);
    elseif pow == 1
        parts{end+1} = sprintf('%.6g*%s', c, varname);
    else
        parts{end+1} = sprintf('%.6g*%s^%d', c, varname, pow);
    end
end
if isempty(parts), s = '0'; else, s = strjoin(parts,'+'); end
end
