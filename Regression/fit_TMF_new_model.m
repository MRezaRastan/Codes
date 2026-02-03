% fit_TMF_new_model_constrained.m
% Updated pipeline with coefficient constraints:
% x: (C1 + 0*A + C3*TS - C4*P) * Phi1   (C2=0)
% y: (C1 + C2*A + C3*TS - 0*P) * Phi1 + (C5 + C6*A + C7*TS - 0*P) * Phi3   (C4=0, C8=0)

clear; close all; clc;

%% -------------------- USER SETTINGS --------------------
manualFiles = []; % leave empty to auto-detect TMF_Fr*_nearWake.xlsx

rho_w = 1000;       % kg/m3
U_inf = 1.4;        % m/s
S = rho_w * U_inf;  % physical scaling
Td = 0.05;          % draft for x scaling
xmin = -30; xmax = 20; % x/Td window to keep
binWidth = 0.1;     % 0 => no binning; >0 uses bin averaging denoise
tmf_threshold = 0;  % abs threshold on TMF to drop very small values (0 keep all)
doBootstrap = true;
nboot = 200;        % bootstrap replications (case-level)
doLOFO = true;
savePngDpi = 300;
verbose = true;

opts_lsq = optimoptions('lsqnonlin','Display','iter','MaxIterations',1000,'MaxFunctionEvaluations',20000,'TolX',1e-9,'TolFun',1e-9);

%% -------------------- FILE DETECTION --------------------
if ~isempty(manualFiles)
    files = manualFiles(:);
else
    D = dir('TMF_Fr*_nearWake.xlsx');
    if isempty(D)
        D = dir('TMF_Fr*.xlsx');
    end
    files = sort({D.name});
end
if isempty(files), error('No Excel files found with expected names.'); end
nFiles = numel(files);
if verbose
    fprintf('Using %d files:\n', nFiles);
    for i=1:nFiles, fprintf('  %s\n', files{i}); end
end

%% -------------------- READ & MASK DATA --------------------
cases = struct();
for k=1:nFiles
    T = readtable(files{k});
    col = @(cands) findFirstVar(T, cands);

    ix_x = col({'x','X'});
    ix_xTMF = col({'xTMF','x_TMF','x_TMF_target','x_TMF_FW_Kelli','x_TMF_NW_Kelli'});
    ix_yTMF = col({'yTMF','y_TMF','y_TMF_target','y_TMF_FW_Kelli','y_TMF_NW_Kelli'});
    if isempty(ix_x), error('File %s: cannot find x column.', files{k}); end
    if isempty(ix_xTMF), error('File %s: cannot find xTMF column.', files{k}); end
    if isempty(ix_yTMF), error('File %s: cannot find yTMF column.', files{k}); end

    % Column names - only using Phi1 for x, Phi1 and Phi3 for y
    names_x = {'x_M1_S1','x_M2_S1','x_M3_S1','x_M4_S1'}; % Only Phi1 terms
    names_y = {'y_M1_S1','y_M2_S1','y_M3_S1','y_M4_S1',... % Phi1 terms
               'y_M1_S3','y_M2_S3','y_M3_S3','y_M4_S3'};   % Phi3 terms
    
    for nm = [names_x names_y]
        if isempty(findFirstVar(T, nm{1}))
            error('File %s: missing required column "%s".', files{k}, nm{1});
        end
    end

    xraw = T{:, ix_x};
    xTMFraw = T{:, ix_xTMF};
    yTMFraw = T{:, ix_yTMF};

    cols_x = zeros(numel(xraw), numel(names_x));
    for j=1:numel(names_x), cols_x(:,j) = T{:, findFirstVar(T, names_x{j})}; end
    cols_y = zeros(numel(xraw), numel(names_y));
    for j=1:numel(names_y), cols_y(:,j) = T{:, findFirstVar(T, names_y{j})}; end

    xscaled = xraw ./ Td;
    xTMF_s = xTMFraw ./ S;
    yTMF_s = yTMFraw ./ S;
    cols_x_s = cols_x ./ S;
    cols_y_s = cols_y ./ S;

    mask = (xscaled >= xmin) & (xscaled <= xmax);
    if tmf_threshold > 0
        mask = mask & (abs(xTMF_s) >= tmf_threshold | abs(yTMF_s) >= tmf_threshold);
    end

    cases(k).fname = files{k};
    cases(k).x = xscaled(mask);
    cases(k).xTMF = xTMF_s(mask);
    cases(k).yTMF = yTMF_s(mask);
    cases(k).xcols = cols_x_s(mask,:);
    cases(k).ycols = cols_y_s(mask,:);
    if verbose, fprintf('Loaded %s -> kept %d rows\n', files{k}, sum(mask)); end
end

%% -------------------- POOL & DECODE --------------------
% X-direction: Only using Phi1 terms (4 columns)
Phi1_x_all = []; A_x_all = []; TS_x_all = []; P_x_all = []; Yx_all = []; case_idx_x = []; xcoord_x = [];

% Y-direction: Using both Phi1 and Phi3 terms (8 columns)
Phi1_y_all = []; Phi3_y_all = []; A_y_all = []; TS_y_all = []; P_y_all = []; Yy_all = []; case_idx_y = []; xcoord_y = [];

for k=1:nFiles
    colsx = cases(k).xcols;
    if ~isempty(colsx)
        % x_M1_S1, x_M2_S1, x_M3_S1, x_M4_S1
        Phi1_x = colsx(:,1); % x_M1_S1 = 1 * Phi1x
        APh1 = colsx(:,2);   % x_M2_S1 = A * Phi1x
        TSPh1 = colsx(:,3);  % x_M3_S1 = τnS_n * Phi1x
        PPh1 = colsx(:,4);   % x_M4_S1 = -P/ε * Phi1x
        
        valid = isfinite(Phi1_x) & (abs(Phi1_x) > 1e-16);
        idxs = find(valid);
        if ~isempty(idxs)
            An_x = NaN(size(Phi1_x)); TauSn_x = NaN(size(Phi1_x)); Pn_x = NaN(size(Phi1_x));
            An_x(valid) = APh1(valid) ./ Phi1_x(valid);
            TauSn_x(valid) = TSPh1(valid) ./ Phi1_x(valid);
            Pn_x(valid) = -PPh1(valid) ./ Phi1_x(valid);

            Phi1_x_all = [Phi1_x_all; Phi1_x(idxs)];
            A_x_all    = [A_x_all;    An_x(idxs)];
            TS_x_all   = [TS_x_all;   TauSn_x(idxs)];
            P_x_all    = [P_x_all;    Pn_x(idxs)];
            Yx_all     = [Yx_all;     cases(k).xTMF(idxs)];
            case_idx_x = [case_idx_x; k*ones(numel(idxs),1)];
            xcoord_x   = [xcoord_x; cases(k).x(idxs)];
        end
    end

    colsy = cases(k).ycols;
    if ~isempty(colsy)
        % First 4 columns: Phi1 terms, Next 4 columns: Phi3 terms
        Phi1_y = colsy(:,1); % y_M1_S1 = 1 * Phi1y
        Phi3_y = colsy(:,5); % y_M1_S3 = 1 * Phi3y
        
        % Get modulators from Phi1 terms (same for Phi3)
        APh1y = colsy(:,2);   % y_M2_S1 = A * Phi1y
        TSPh1y = colsy(:,3);  % y_M3_S1 = τnS_n * Phi1y
        PPh1y = colsy(:,4);   % y_M4_S1 = -P/ε * Phi1y
        
        validy = isfinite(Phi1_y) & (abs(Phi1_y) > 1e-16);
        idy = find(validy);
        if ~isempty(idy)
            An_y = NaN(size(Phi1_y)); TauSn_y = NaN(size(Phi1_y)); Pn_y = NaN(size(Phi1_y));
            An_y(validy) = APh1y(validy) ./ Phi1_y(validy);
            TauSn_y(validy) = TSPh1y(validy) ./ Phi1_y(validy);
            Pn_y(validy) = -PPh1y(validy) ./ Phi1_y(validy);

            Phi1_y_all = [Phi1_y_all; Phi1_y(idy)];
            Phi3_y_all = [Phi3_y_all; Phi3_y(idy)];
            A_y_all    = [A_y_all;    An_y(idy)];
            TS_y_all   = [TS_y_all;   TauSn_y(idy)];
            P_y_all    = [P_y_all;    Pn_y(idy)];
            Yy_all     = [Yy_all;     cases(k).yTMF(idy)];
            case_idx_y = [case_idx_y; k*ones(numel(idy),1)];
            xcoord_y   = [xcoord_y; cases(k).x(idy)];
        end
    end
end

% Ensure column vectors
Phi1_x_all = Phi1_x_all(:); A_x_all = A_x_all(:); TS_x_all = TS_x_all(:); P_x_all = P_x_all(:); Yx_all = Yx_all(:);
Phi1_y_all = Phi1_y_all(:); Phi3_y_all = Phi3_y_all(:); A_y_all = A_y_all(:); TS_y_all = TS_y_all(:); P_y_all = P_y_all(:); Yy_all = Yy_all(:);
case_idx_x = case_idx_x(:); case_idx_y = case_idx_y(:); xcoord_x = xcoord_x(:); xcoord_y = xcoord_y(:);

if verbose, fprintf('Pooled rows: x=%d, y=%d\n', numel(Yx_all), numel(Yy_all)); end

%% -------------------- BINNING (optional denoise) --------------------
if binWidth > 0
    edges = xmin:binWidth:xmax; binCenters = (edges(1:end-1)+edges(2:end))/2; nBins = numel(binCenters);
    P1x_b=[]; A_x_b=[]; TS_x_b=[]; P_x_b=[]; Yx_b=[]; case_x_b=[]; xcoord_x_b=[];
    P1y_b=[]; P3y_b=[]; A_y_b=[]; TS_y_b=[]; P_y_b=[]; Yy_b=[]; case_y_b=[]; xcoord_y_b=[];
    for k=1:nFiles
        idx = find(case_idx_x==k);
        if ~isempty(idx)
            [~,~,binIdx] = histcounts(xcoord_x(idx), edges);
            for b=1:nBins
                sel = binIdx==b;
                if any(sel)
                    P1x_b(end+1,1) = mean(Phi1_x_all(idx(sel)));
                    A_x_b(end+1,1) = mean(A_x_all(idx(sel)));
                    TS_x_b(end+1,1) = mean(TS_x_all(idx(sel)));
                    P_x_b(end+1,1) = mean(P_x_all(idx(sel)));
                    Yx_b(end+1,1) = mean(Yx_all(idx(sel)));
                    case_x_b(end+1,1) = k;
                    xcoord_x_b(end+1,1) = binCenters(b);
                end
            end
        end
        idy = find(case_idx_y==k);
        if ~isempty(idy)
            [~,~,binIdy] = histcounts(xcoord_y(idy), edges);
            for b=1:nBins
                sel = binIdy==b;
                if any(sel)
                    P1y_b(end+1,1) = mean(Phi1_y_all(idy(sel)));
                    P3y_b(end+1,1) = mean(Phi3_y_all(idy(sel)));
                    A_y_b(end+1,1) = mean(A_y_all(idy(sel)));
                    TS_y_b(end+1,1) = mean(TS_y_all(idy(sel)));
                    P_y_b(end+1,1) = mean(P_y_all(idy(sel)));
                    Yy_b(end+1,1) = mean(Yy_all(idy(sel)));
                    case_y_b(end+1,1) = k;
                    xcoord_y_b(end+1,1) = binCenters(b);
                end
            end
        end
    end
    if ~isempty(P1x_b)
        Phi1_x_all = P1x_b; A_x_all = A_x_b; TS_x_all = TS_x_b; P_x_all = P_x_b; Yx_all = Yx_b; case_idx_x = case_x_b; xcoord_x = xcoord_x_b;
    end
    if ~isempty(P1y_b)
        Phi1_y_all = P1y_b; Phi3_y_all = P3y_b; A_y_all = A_y_b; TS_y_all = TS_y_b; P_y_all = P_y_b; Yy_all = Yy_b; case_idx_y = case_y_b; xcoord_y = xcoord_y_b;
    end
    if verbose, fprintf('After binning: pooled rows: x=%d, y=%d\n', numel(Yx_all), numel(Yy_all)); end
    edges_used = edges;
else
    edges_used = linspace(xmin, xmax, 51);
    binCenters = (edges_used(1:end-1)+edges_used(2:end))/2;
end

%% -------------------- FIT MODELS (WITH CONSTRAINTS) --------------------
% x-model: y = (C1 + 0*A + C3*TS - C4*P) * Phi1   (C2=0)
% y-model: y = (C1 + C2*A + C3*TS - 0*P) * Phi1 + (C5 + C6*A + C7*TS - 0*P) * Phi3   (C4=0, C8=0)

% X-direction: 3 coefficients (C1, C3, C4) - C2=0
C0x = [0.1; 0; 0];  % C1, C3, C4  (Note: C2 is omitted since it's fixed to 0)

% Y-direction: 6 coefficients (C1, C2, C3, C5, C6, C7) - C4=0, C8=0
C0y = [0.1; 0; 0; 0.1; 0; 0];  % C1, C2, C3, C5, C6, C7

% Residual functions with constraints
% X: C2=0
resid_x = @(C) ((C(1) + 0*A_x_all + C(2).*TS_x_all - C(3).*P_x_all) .* Phi1_x_all - Yx_all);

% Y: C4=0, C8=0
resid_y = @(C) ((C(1) + C(2).*A_y_all + C(3).*TS_y_all - 0*P_y_all) .* Phi1_y_all + ...
                (C(4) + C(5).*A_y_all + C(6).*TS_y_all - 0*P_y_all) .* Phi3_y_all - Yy_all);

if verbose, fprintf('Fitting x-model (3 coefficients, C2=0)...\n'); end
Cx = lsqnonlin(resid_x, C0x, [], [], opts_lsq);
Cx = Cx(:); % ensure column
% Build full x-coefficient vector for reporting: [C1; C2=0; C3; C4]
Cx_full = [Cx(1); 0; Cx(2); Cx(3)];

if verbose, fprintf('Fitting y-model (6 coefficients, C4=0, C8=0)...\n'); end
Cy = lsqnonlin(resid_y, C0y, [], [], opts_lsq);
Cy = Cy(:); % ensure column
% Build full y-coefficient vector for reporting: [C1; C2; C3; C4=0; C5; C6; C7; C8=0]
Cy_full = [Cy(1); Cy(2); Cy(3); 0; Cy(4); Cy(5); Cy(6); 0];

fprintf('\nFitted coefficients:\n');
fprintf('Cx (C1, C2=0, C3, C4) = %.4f, %.4f, %.4f, %.4f\n', Cx_full(1), Cx_full(2), Cx_full(3), Cx_full(4));
fprintf('Cy (C1, C2, C3, C4=0, C5, C6, C7, C8=0) = %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n', ...
        Cy_full(1), Cy_full(2), Cy_full(3), Cy_full(4), Cy_full(5), Cy_full(6), Cy_full(7), Cy_full(8));

%% -------------------- predictions & term contributions --------------------
% X-terms: 3 terms (C1, C3, C4) * modulator * Phi1 (C2=0)
termX = zeros(numel(Yx_all),3);
termX(:,1) = Cx(1) .* Phi1_x_all;               % C1 * Phi1
termX(:,2) = Cx(2) .* (TS_x_all .* Phi1_x_all); % C3 * τS * Phi1
termX(:,3) = -Cx(3) .* (P_x_all .* Phi1_x_all); % -C4 * P * Phi1
Yx_pred = sum(termX,2);

% Y-terms: 6 terms 
termY = zeros(numel(Yy_all),6);
% Phi1 terms (C1, C2, C3) - C4=0
termY(:,1) = Cy(1) .* Phi1_y_all;               % C1 * Phi1
termY(:,2) = Cy(2) .* (A_y_all .* Phi1_y_all);  % C2 * A * Phi1
termY(:,3) = Cy(3) .* (TS_y_all .* Phi1_y_all); % C3 * τS * Phi1
% Phi3 terms (C5, C6, C7) - C8=0
termY(:,4) = Cy(4) .* Phi3_y_all;               % C5 * Phi3
termY(:,5) = Cy(5) .* (A_y_all .* Phi3_y_all);  % C6 * A * Phi3
termY(:,6) = Cy(6) .* (TS_y_all .* Phi3_y_all); % C7 * τS * Phi3
Yy_pred = sum(termY,2);

%% -------------------- Calculate Conditioned Statistics --------------------
cond_corr_x = nan(nFiles,1); cond_corr_y = nan(nFiles,1);
cond_normstd_x = nan(nFiles,1); cond_normstd_y = nan(nFiles,1);

global_corr_x = corr(Yx_all, Yx_pred);
global_corr_y = corr(Yy_all, Yy_pred);
global_normstd_x = std(Yx_pred)/std(Yx_all);
global_normstd_y = std(Yx_pred)/std(Yy_all);

% Calculate per-case statistics
for k=1:nFiles
    % X-direction statistics
    idx = case_idx_x==k;
    if sum(idx) > 1
        Yx_case = Yx_all(idx);
        Yx_pred_case = Yx_pred(idx);
        
        if std(Yx_case) > 0 && std(Yx_pred_case) > 0
            cond_corr_x(k) = corr(Yx_case, Yx_pred_case);
        end
        
        if std(Yx_case) > 0
            cond_normstd_x(k) = std(Yx_pred_case)/std(Yx_case);
        end
    end
    
    % Y-direction statistics
    idy = case_idx_y==k;
    if sum(idy) > 1
        Yy_case = Yy_all(idy);
        Yy_pred_case = Yy_pred(idy);
        
        if std(Yy_case) > 0 && std(Yy_pred_case) > 0
            cond_corr_y(k) = corr(Yy_case, Yy_pred_case);
        end
        
        if std(Yy_case) > 0
            cond_normstd_y(k) = std(Yy_pred_case)/std(Yy_case);
        end
    end
end

% Save conditioned statistics
cond_stats_table = table((1:nFiles)', cond_corr_x, cond_normstd_x, cond_corr_y, cond_normstd_y, ...
    'VariableNames', {'case', 'cond_corr_x', 'cond_normstd_x', 'cond_corr_y', 'cond_normstd_y'});
writetable(cond_stats_table, 'TMF_constrained_model_conditioned_statistics.csv');

% Display statistics
fprintf('\n=== GLOBAL CONDITIONED STATISTICS ===\n');
fprintf('X-direction:\n');
fprintf('  Correlation coefficient (C): %.4f\n', global_corr_x);
fprintf('  Normalized std (σ_model/σ_data): %.4f\n', global_normstd_x);
fprintf('Y-direction:\n');
fprintf('  Correlation coefficient (C): %.4f\n', global_corr_y);
fprintf('  Normalized std (σ_model/σ_data): %.4f\n', global_normstd_y);

fprintf('\n=== PER-CASE CONDITIONED STATISTICS ===\n');
for k=1:nFiles
    fprintf('Case %d:\n', k);
    fprintf('  X: C=%.4f, σ_N=%.4f\n', cond_corr_x(k), cond_normstd_x(k));
    fprintf('  Y: C=%.4f, σ_N=%.4f\n', cond_corr_y(k), cond_normstd_y(k));
end

R2_x = 1 - sum((Yx_all-Yx_pred).^2) / (sum((Yx_all-mean(Yx_all)).^2)+eps);
R2_y = 1 - sum((Yy_all-Yy_pred).^2) / (sum((Yy_all-mean(Yy_all)).^2)+eps);
fprintf('\nGlobal R2: x=%.4f, y=%.4f\n', R2_x, R2_y);

% per-case RMSE
RMSE_x = nan(nFiles,1); RMSE_y = nan(nFiles,1);
for k=1:nFiles
    idx = case_idx_x==k; if any(idx), RMSE_x(k)=sqrt(mean((Yx_all(idx)-Yx_pred(idx)).^2)); end
    idy = case_idx_y==k; if any(idy), RMSE_y(k)=sqrt(mean((Yy_all(idy)-Yy_pred(idy)).^2)); end
end
writetable(table((1:nFiles)',RMSE_x,RMSE_y,'VariableNames',{'case','RMSE_x','RMSE_y'}),'TMF_constrained_model_RMSE_perCase.csv');

%% -------------------- LOFO --------------------
if doLOFO
    LOFO_RMSE_x = nan(nFiles,1); LOFO_RMSE_y = nan(nFiles,1);
    for k=1:nFiles
        % X-direction LOFO (C2=0)
        train = (case_idx_x ~= k); test = (case_idx_x == k);
        if sum(test)>=1 && sum(train)>=10
            Xtrain = [Phi1_x_all(train), TS_x_all(train).*Phi1_x_all(train), -P_x_all(train).*Phi1_x_all(train)];
            b_cv = Xtrain \ Yx_all(train);
            Xtest = [Phi1_x_all(test), TS_x_all(test).*Phi1_x_all(test), -P_x_all(test).*Phi1_x_all(test)];
            ypred = Xtest * b_cv;
            LOFO_RMSE_x(k) = sqrt(mean((Yx_all(test)-ypred).^2));
        end
        
        % Y-direction LOFO (C4=0, C8=0)
        trainY = (case_idx_y ~= k); testY = (case_idx_y == k);
        if sum(testY)>=1 && sum(trainY)>=10
            XtrainY = [Phi1_y_all(trainY), A_y_all(trainY).*Phi1_y_all(trainY), ...
                       TS_y_all(trainY).*Phi1_y_all(trainY), ...
                       Phi3_y_all(trainY), A_y_all(trainY).*Phi3_y_all(trainY), ...
                       TS_y_all(trainY).*Phi3_y_all(trainY)];
            bcvY = XtrainY \ Yy_all(trainY);
            XtestY = [Phi1_y_all(testY), A_y_all(testY).*Phi1_y_all(testY), ...
                      TS_y_all(testY).*Phi1_y_all(testY), ...
                      Phi3_y_all(testY), A_y_all(testY).*Phi3_y_all(testY), ...
                      TS_y_all(testY).*Phi3_y_all(testY)];
            ypredY = XtestY * bcvY;
            LOFO_RMSE_y(k) = sqrt(mean((Yy_all(testY)-ypredY).^2));
        end
    end
    writetable(table((1:nFiles)',LOFO_RMSE_x,LOFO_RMSE_y,'VariableNames',{'case','LOFO_RMSE_x','LOFO_RMSE_y'}),'TMF_constrained_model_LOFO_RMSE.csv');
end

%% -------------------- BOOTSTRAP (case-level) & CI --------------------
if doBootstrap
    rng(1);
    Cx_boot = nan(nboot,3);  % 3 coefficients for x (C1, C3, C4)
    Cy_boot = nan(nboot,6);  % 6 coefficients for y (C1, C2, C3, C5, C6, C7)
    
    for ib=1:nboot
        samp = randsample(nFiles,nFiles,true);
        maskx = false(size(Yx_all)); masky = false(size(Yy_all));
        for s=1:numel(samp)
            maskx = maskx | (case_idx_x==samp(s));
            masky = masky | (case_idx_y==samp(s));
        end
        
        try
            Cb_x = lsqnonlin(@(C) ((C(1) + 0*A_x_all(maskx) + C(2).*TS_x_all(maskx) - C(3).*P_x_all(maskx)) .* Phi1_x_all(maskx) - Yx_all(maskx)), ...
                             C0x, [], [], optimoptions(opts_lsq,'Display','off'));
            Cx_boot(ib,:) = Cb_x(:)';
        catch
            Cx_boot(ib,:) = Cx(:)';
        end
        
        try
            Cb_y = lsqnonlin(@(C) ((C(1) + C(2).*A_y_all(masky) + C(3).*TS_y_all(masky) - 0*P_y_all(masky)) .* Phi1_y_all(masky) + ...
                                    (C(4) + C(5).*A_y_all(masky) + C(6).*TS_y_all(masky) - 0*P_y_all(masky)) .* Phi3_y_all(masky) - Yy_all(masky)), ...
                             C0y, [], [], optimoptions(opts_lsq,'Display','off'));
            Cy_boot(ib,:) = Cb_y(:)';
        catch
            Cy_boot(ib,:) = Cy(:)';
        end
    end
    
    Cx_lo = prctile(Cx_boot,2.5); Cx_hi = prctile(Cx_boot,97.5);
    Cy_lo = prctile(Cy_boot,2.5); Cy_hi = prctile(Cy_boot,97.5);
    
    % Create coefficient tables with constraints
    coefNames_x = {'C1';'C2';'C3';'C4'};
    Cx_est_full = [Cx(1); 0; Cx(2); Cx(3)];
    Cx_CIlo_full = [Cx_lo(1); 0; Cx_lo(2); Cx_lo(3)];
    Cx_CIhi_full = [Cx_hi(1); 0; Cx_hi(2); Cx_hi(3)];
    
    coefNames_y = {'C1';'C2';'C3';'C4';'C5';'C6';'C7';'C8'};
    Cy_est_full = [Cy(1); Cy(2); Cy(3); 0; Cy(4); Cy(5); Cy(6); 0];
    Cy_CIlo_full = [Cy_lo(1); Cy_lo(2); Cy_lo(3); 0; Cy_lo(4); Cy_lo(5); Cy_lo(6); 0];
    Cy_CIhi_full = [Cy_hi(1); Cy_hi(2); Cy_hi(3); 0; Cy_hi(4); Cy_hi(5); Cy_hi(6); 0];
    
    Tcoefs_x = table(coefNames_x, Cx_est_full, Cx_CIlo_full, Cx_CIhi_full, ...
        'VariableNames', {'Coefficient','Estimate','CI_lower','CI_upper'});
    Tcoefs_y = table(coefNames_y, Cy_est_full, Cy_CIlo_full, Cy_CIhi_full, ...
        'VariableNames', {'Coefficient','Estimate','CI_lower','CI_upper'});
    
    writetable(Tcoefs_x,'TMF_constrained_model_coeffs_x_with_CI.csv');
    writetable(Tcoefs_y,'TMF_constrained_model_coeffs_y_with_CI.csv');
    
    save('TMF_constrained_model_fit_results.mat','Cx','Cy','Cx_boot','Cy_boot','Tcoefs_x','Tcoefs_y','-v7.3');
end

%% -------------------- PLOTTING --------------------
edges = edges_used; binCenters = (edges(1:end-1)+edges(2:end))/2; nBins = numel(binCenters);
colors = lines(8);

% Define meaningful term names with constraints
xTermNames = {
    'C1·Φ₁ (Base)';
    'C3·τₛ·Φ₁ (Shear)';
    '-C4·P·Φ₁ (Production)'
    % Note: C2 term is omitted (set to 0)
};

yTermNames = {
    'C1·Φ₁ (Base)';
    'C2·A·Φ₁ (Anisotropy)';
    'C3·τₛ·Φ₁ (Shear)';
    'C5·Φ₃ (Buoyancy)';
    'C6·A·Φ₃ (Buoyancy+Aniso)';
    'C7·τₛ·Φ₃ (Buoyancy+Shear)'
    % Note: C4 and C8 terms are omitted (set to 0)
};

% Define markers and line styles
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*', '+', 'x'};
lineStyles = {'-', '--', ':', '-.', '-', '--', ':', '-.'};

%% ----- POOLED CONDITIONAL MEAN X -----
[~,~,binIdxAll] = histcounts(xcoord_x, edges);
bin_true_all = nan(nBins,1); bin_pred_all = nan(nBins,1); 
bin_terms_all = nan(nBins,3); % 3 terms for x-model (C1, C3, C4)

for b=1:nBins
    sel = binIdxAll==b;
    if any(sel)
        bin_true_all(b) = mean(Yx_all(sel));
        bin_pred_all(b) = mean(Yx_pred(sel));
        bin_terms_all(b,:) = mean(termX(sel,:),1);
    end
end

h = figure('Visible','off', 'Position', [100, 100, 1200, 800]);
plot(binCenters, bin_true_all, 'k-', 'LineWidth', 2.5, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); hold on;
plot(binCenters, bin_pred_all, 'r-', 'LineWidth', 2.5, 'Marker', 's', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% Plot the 3 x-model terms
for tt = 1:3
    plot(binCenters, bin_terms_all(:,tt), ...
        'LineStyle', lineStyles{tt}, ...
        'LineWidth', 1.5, ...
        'Marker', markers{tt}, ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', colors(tt,:), ...
        'Color', colors(tt,:));
end

xlabel('x/T_d', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Scaled xTMF', 'FontSize', 12, 'FontWeight', 'bold');
title('Pooled Conditional Mean: x-direction TMF (C2=0)', 'FontSize', 14, 'FontWeight', 'bold');
legendLabels = {'LES Data', 'Model Prediction', xTermNames{:}};
legend(legendLabels, 'Location', 'bestoutside', 'FontSize', 10);
grid on;
box on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
print('-dpng', 'TMF_constrained_model_pooled_condMean_x.png', sprintf('-r%d', savePngDpi));
close(h);

%% ----- POOLED CONDITIONAL MEAN Y -----
[~,~,binIdxAllY] = histcounts(xcoord_y, edges);
bin_true_all_y = nan(nBins,1); bin_pred_all_y = nan(nBins,1); 
bin_terms_all_y = nan(nBins,6); % 6 terms for y-model

for b=1:nBins
    sel = binIdxAllY==b;
    if any(sel)
        bin_true_all_y(b) = mean(Yy_all(sel));
        bin_pred_all_y(b) = mean(Yy_pred(sel));
        bin_terms_all_y(b,:) = mean(termY(sel,:),1);
    end
end

h = figure('Visible','off', 'Position', [100, 100, 1200, 800]);
plot(binCenters, bin_true_all_y, 'k-', 'LineWidth', 2.5, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); hold on;
plot(binCenters, bin_pred_all_y, 'r-', 'LineWidth', 2.5, 'Marker', 's', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% Plot all 6 y-model terms
for tt = 1:6
    plot(binCenters, bin_terms_all_y(:,tt), ...
        'LineStyle', lineStyles{tt}, ...
        'LineWidth', 1.5, ...
        'Marker', markers{tt}, ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', colors(tt,:), ...
        'Color', colors(tt,:));
end

xlabel('x/T_d', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Scaled yTMF', 'FontSize', 12, 'FontWeight', 'bold');
title('Pooled Conditional Mean: y-direction TMF (C4=0, C8=0)', 'FontSize', 14, 'FontWeight', 'bold');
legendLabels = {'LES Data', 'Model Prediction', yTermNames{:}};
legend(legendLabels, 'Location', 'bestoutside', 'FontSize', 9);
grid on;
box on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
print('-dpng', 'TMF_constrained_model_pooled_condMean_y.png', sprintf('-r%d', savePngDpi));
close(h);

%% ----- PER-CASE PLOTS -----
for k=1:nFiles
    % X-DIRECTION PLOT
    idx = find(case_idx_x==k);
    if ~isempty(idx)
        [~,~,bidx] = histcounts(xcoord_x(idx), edges);
        bin_true = []; bin_pred = []; bcent = [];
        Tcon = [];
        
        for b=1:nBins
            sel = bidx==b;
            if any(sel)
                bin_true(end+1,1) = mean(Yx_all(idx(sel)));
                bin_pred(end+1,1) = mean(Yx_pred(idx(sel)));
                Tcon = [Tcon; mean(termX(idx(sel),:),1)];
                bcent(end+1,1) = binCenters(b);
            end
        end
        
        if ~isempty(bin_true) && ~isempty(Tcon) && size(Tcon,1) == length(bcent)
            fig = figure('Visible','off', 'Position', [100, 100, 1200, 800]);
            plot(bcent, bin_true, 'k-', 'LineWidth', 2.5, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); hold on;
            plot(bcent, bin_pred, 'r-', 'LineWidth', 2.5, 'Marker', 's', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            
            % Plot the 3 x-model terms
            for tt = 1:3
                plot(bcent, Tcon(:,tt), ...
                    'LineStyle', lineStyles{tt}, ...
                    'LineWidth', 1.5, ...
                    'Marker', markers{tt}, ...
                    'MarkerSize', 6, ...
                    'MarkerFaceColor', colors(tt,:), ...
                    'Color', colors(tt,:));
            end
            
            xlabel('x/T_d', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('Scaled xTMF', 'FontSize', 12, 'FontWeight', 'bold');
            title(sprintf('Case %d: x-direction TMF (C2=0)', k), 'FontSize', 14, 'FontWeight', 'bold');
            legendLabels = {'LES Data', 'Model Prediction', xTermNames{:}};
            legend(legendLabels, 'Location', 'bestoutside', 'FontSize', 10);
            grid on;
            box on;
            set(gca, 'FontSize', 11, 'LineWidth', 1.5);
            print('-dpng', sprintf('TMF_constrained_model_case%02d_condMean_x.png', k), sprintf('-r%d', savePngDpi));
            close(fig);
        end
    end
    
    % Y-DIRECTION PLOT
    idy = find(case_idx_y==k);
    if ~isempty(idy)
        [~,~,bidx] = histcounts(xcoord_y(idy), edges);
        bin_true = []; bin_pred = []; bcent = [];
        Tcon = [];
        
        for b=1:nBins
            sel = bidx==b;
            if any(sel)
                bin_true(end+1,1) = mean(Yy_all(idy(sel)));
                bin_pred(end+1,1) = mean(Yy_pred(idy(sel)));
                Tcon = [Tcon; mean(termY(idy(sel),:),1)];
                bcent(end+1,1) = binCenters(b);
            end
        end
        
        if ~isempty(bin_true) && ~isempty(Tcon) && size(Tcon,1) == length(bcent)
            fig = figure('Visible','off', 'Position', [100, 100, 1200, 800]);
            plot(bcent, bin_true, 'k-', 'LineWidth', 2.5, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); hold on;
            plot(bcent, bin_pred, 'r-', 'LineWidth', 2.5, 'Marker', 's', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            
            % Plot all 6 y-model terms
            for tt = 1:6
                plot(bcent, Tcon(:,tt), ...
                    'LineStyle', lineStyles{tt}, ...
                    'LineWidth', 1.5, ...
                    'Marker', markers{tt}, ...
                    'MarkerSize', 6, ...
                    'MarkerFaceColor', colors(tt,:), ...
                    'Color', colors(tt,:));
            end
            
            xlabel('x/T_d', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('Scaled yTMF', 'FontSize', 12, 'FontWeight', 'bold');
            title(sprintf('Case %d: y-direction TMF (C4=0, C8=0)', k), 'FontSize', 14, 'FontWeight', 'bold');
            legendLabels = {'LES Data', 'Model Prediction', yTermNames{:}};
            legend(legendLabels, 'Location', 'bestoutside', 'FontSize', 9);
            grid on;
            box on;
            set(gca, 'FontSize', 11, 'LineWidth', 1.5);
            print('-dpng', sprintf('TMF_constrained_model_case%02d_condMean_y.png', k), sprintf('-r%d', savePngDpi));
            close(fig);
        end
    end
end

%% -------------------- SAVE SUMMARY --------------------
save('TMF_constrained_model_fit_fullpipeline.mat', ...
     'Cx', 'Cy', 'Cx_boot', 'Cy_boot', 'Tcoefs_x', 'Tcoefs_y', ...
     'Yx_all', 'Yx_pred', 'Yy_all', 'Yy_pred', 'termX', 'termY', ...
     'cases', 'files', 'RMSE_x', 'RMSE_y', 'LOFO_RMSE_x', 'LOFO_RMSE_y', ...
     'R2_x', 'R2_y', 'cond_stats_table', 'global_corr_x', 'global_corr_y', ...
     'global_normstd_x', 'global_normstd_y', '-v7.3');

fprintf('\nAll done. Results saved to MAT/CSV and PNG files.\n');
fprintf('Conditioned statistics saved to: TMF_constrained_model_conditioned_statistics.csv\n');

%% -------------------- helper: findFirstVar --------------------
function idx = findFirstVar(Tbl, cands)
    if ischar(cands) || isstring(cands), cands = {char(cands)}; end
    idx = [];
    for i=1:numel(cands)
        nm = cands{i};
        if ismember(nm, Tbl.Properties.VariableNames), idx=find(strcmp(nm, Tbl.Properties.VariableNames),1); return; end
        im = find(strcmpi(nm, Tbl.Properties.VariableNames),1);
        if ~isempty(im), idx=im; return; end
        nv = strrep(nm,'-','_'); nv = strrep(nv,' ','_');
        im2 = find(strcmpi(nv, Tbl.Properties.VariableNames),1);
        if ~isempty(im2), idx=im2; return; end
    end
end