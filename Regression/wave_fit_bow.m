% wave_fit_JFM.m
% Fits eta(x,Fr) as 4th-order polynomial in x with coefficients quadratic in Fr.
% Produces: fitted_wave_profiles.xlsx (coeff table + equation + RMSE),
% wave_fit_summary.tex, validation_plot.png, model_text_version.txt
clear; clc; close all;

%% ---- User settings ----
filename = 'wave_data.xlsx';
Fr_values = [0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5];  % adjust if needed
poly_order_x = 4;   % 4th-order in x
poly_order_Fr = 2;  % quadratic in Fr
num_coeffs = poly_order_x + 1;

%% ---- 1. Load input table ----
if ~isfile(filename)
    error('Input file not found: %s. Place wave_data.xlsx in the current folder.', filename);
end
data = readtable(filename, 'VariableNamingRule', 'preserve');
colNames = data.Properties.VariableNames;

x_data = cell(length(Fr_values),1);
y_data = cell(length(Fr_values),1);
for i = 1:length(Fr_values)
    if abs(Fr_values(i) - round(Fr_values(i))) < 1e-6
        Fr_str = sprintf('%.1f', Fr_values(i));
    else
        Fr_str = sprintf('%.2f', Fr_values(i));
    end
    x_col = ['Fr=' Fr_str '_x'];
    y_col = ['Fr=' Fr_str '_y'];
    if ~ismember(x_col, colNames) || ~ismember(y_col, colNames)
        error('Missing column(s): %s or %s in %s', x_col, y_col, filename);
    end
    x_temp = data{:, x_col};
    y_temp = data{:, y_col};
    valid = ~isnan(x_temp) & ~isnan(y_temp);
    x_temp = x_temp(valid);
    y_temp = y_temp(valid);
    if ~issorted(x_temp)
        [x_temp, idxs] = sort(x_temp);
        y_temp = y_temp(idxs);
    end
    x_data{i} = x_temp(:);
    y_data{i} = y_temp(:);
end

%% ---- 2. Fit polynomial in x for each Fr ----
coeff_matrix = zeros(length(Fr_values), num_coeffs); % each row: descending coeffs
for i = 1:length(Fr_values)
    p = polyfit(x_data{i}, y_data{i}, poly_order_x);
    coeff_matrix(i, :) = p;
end

%% ---- 3. Fit each x-coefficient as quadratic in Fr ----
coeff_func = cell(1, num_coeffs); % each is [a2 a1 a0] for a2*Fr^2 + a1*Fr + a0
for j = 1:num_coeffs
    y_for_Fr = coeff_matrix(:, j);
    coeff_func{j} = polyfit(Fr_values(:)', y_for_Fr(:)', poly_order_Fr);
end

%% ---- 4. Build symbolic expression (for display/latex) ----
useSymbolic = true;
try
    syms x Fr
    G_sym = sym(0);
    powers = poly_order_x:-1:0;
    for j = 1:num_coeffs
        p_vec = coeff_func{j};
        coeff_in_Fr = poly2sym(p_vec, Fr);
        G_sym = G_sym + coeff_in_Fr * x^powers(j);
    end
    G_clean = simplify(vpa(G_sym, 6));
catch
    useSymbolic = false;
    G_clean = []; %#ok<NASGU>
end

%% ---- 5. Prepare coefficient strings for paper-style display ----
% We want A4..A0 where A_n(Fr) = a2 Fr^2 + a1 Fr + a0
powers = poly_order_x:-1:0;
A_strings_tex = cell(num_coeffs,1);  % LaTeX-friendly strings (for tex doc)
A_strings_txt = cell(num_coeffs,1);  % plain-text strings (for model_text_version.txt)
for j = 1:num_coeffs
    p = coeff_func{j}; % [a2 a1 a0]
    A_strings_tex{j} = quadFrTexTex(p(1), p(2), p(3));
    A_strings_txt{j} = quadFrTexPlain(p(1), p(2), p(3));
end

%% ---- 6. Compute RMSE across all original LES points ----
residuals = [];
for i = 1:length(Fr_values)
    Fr_i = Fr_values(i);
    % evaluate coefficient functions at this Fr to get polynomial in x (descending)
    p_local = zeros(1, num_coeffs);
    for j = 1:num_coeffs
        p_local(j) = polyval(coeff_func{j}, Fr_i);
    end
    % predict for each x_data point
    y_pred = polyval(p_local, x_data{i});
    res = y_data{i} - y_pred;
    residuals = [residuals; res(:)]; %#ok<AGROW>
end
RMSE_all = sqrt(mean(residuals.^2));

%% ---- 7. Save coefficient table + RMSE to Excel ----
out_xlsx = 'fitted_wave_profiles.xlsx';
power_labels = arrayfun(@(k) sprintf('x^%d', k), powers, 'UniformOutput', false);
coeff_table_cell = cell(num_coeffs, 3);
for j = 1:num_coeffs
    coeff_table_cell{j,1} = coeff_func{j}(1); % Fr^2
    coeff_table_cell{j,2} = coeff_func{j}(2); % Fr
    coeff_table_cell{j,3} = coeff_func{j}(3); % const
end
T = table(power_labels(:), cell2mat(coeff_table_cell(:,1)), cell2mat(coeff_table_cell(:,2)), cell2mat(coeff_table_cell(:,3)), ...
    'VariableNames', {'Term', 'Fr2', 'Fr', 'Constant'});
writetable(T, out_xlsx, 'Sheet', 'Coefficients', 'Range', 'A1');

% Equation sheet: if symbolic available use its latex; else write built A_n lines
if useSymbolic
    writecell({'Fitted Equation (symbolic)'; char(G_clean)}, out_xlsx, 'Sheet', 'Equation', 'Range', 'A1');
else
    lines = [{'Fitted Equation (A_n listed below):'}; A_strings_txt];
    writecell(lines, out_xlsx, 'Sheet', 'Equation', 'Range', 'A1');
end

% Metrics sheet with RMSE
writecell({'Metric', 'Value'; 'RMSE_all', RMSE_all}, out_xlsx, 'Sheet', 'Metrics', 'Range', 'A1');

%% ---- 8. Write paper-style text version (txt + tex) ----
% Plain text file for quick copy-paste into manuscript
txt_file = 'model_text_version.txt';
fid = fopen(txt_file, 'w');
fprintf(fid, 'Text Version (for paper):\n\n');
fprintf(fid, 'The fitted function is expressed as a fourth-order polynomial in x, with coefficients dependent on the Froude number Fr:\n\n');
fprintf(fid, 'f(x, Fr) = sum_{n=0}^4 A_n(Fr) x^n\n\n');
fprintf(fid, 'where the coefficients are:\n\n');
for j = 1:num_coeffs
    fprintf(fid, 'A_%d(Fr) = %s\n', poly_order_x-(j-1), A_strings_txt{j});
end
fprintf(fid, '\nRoot-mean-square error (RMSE) between LES and fitted model (all points): %.6g\n', RMSE_all);
fclose(fid);

% Build LaTeX doc with the A_n lines and RMSE
tex_filename = 'wave_fit_summary.tex';
fid = fopen(tex_filename, 'w');
if fid == -1, error('Could not open %s for writing.', tex_filename); end
fprintf(fid, '%% Auto-generated LaTeX for JFM-style presentation\n');
fprintf(fid, '\\documentclass[10pt]{article}\n\\usepackage{amsmath,booktabs}\n\\usepackage[margin=1in]{geometry}\n');
fprintf(fid, '\\title{Fitted Bow Wave Profile Equation}\\author{}\\date{}\n\\begin{document}\\maketitle\n\n');
fprintf(fid, 'Text Version (for paper):\n\n');
fprintf(fid, 'The fitted function is expressed as a fourth-order polynomial in $x$, with coefficients dependent on the Froude number $\\mathrm{Fr}$:\n\n');
fprintf(fid, '\\begin{equation*}\\label{eq:fitted}\n f(x,\\mathrm{Fr}) = \\sum_{n=0}^{4} A_n(\\mathrm{Fr}) x^n\n\\end{equation*}\n\n');
fprintf(fid, 'where the coefficients are:\n\n');
for j = 1:num_coeffs
    Aname = sprintf('A_{%d}(\\mathrm{Fr})', poly_order_x-(j-1));
    fprintf(fid, '\\noindent %s = %s \\\\\n', Aname, A_strings_tex{j});
end
fprintf(fid, '\n\\vspace{1em}\nRoot-mean-square error (RMSE) between LES and fitted model (all points): $%.6g$.\n', RMSE_all);
fprintf(fid, '\\end{document}\n');
fclose(fid);
fprintf('Saved: %s and %s\n', txt_file, tex_filename);

%% ---- 9. Plot comparison (LES vs fitted) ----
all_x = vertcat(x_data{:});
x_common = linspace(min(all_x), max(all_x), 300)';
fig = figure('Units','pixels','Position',[200 200 1000 600]); hold on;
for i = 1:length(Fr_values)
    % LES interpolate for smooth display
    y_les_i = interp1(x_data{i}, y_data{i}, x_common, 'linear', 'extrap');
    % evaluate fit
    p_local = zeros(1, num_coeffs);
    for j = 1:num_coeffs
        p_local(j) = polyval(coeff_func{j}, Fr_values(i));
    end
    y_fit_i = polyval(p_local, x_common);
    plot(x_common, y_les_i, '-', 'DisplayName', sprintf('LES Fr=%.2g', Fr_values(i)));
    plot(x_common, y_fit_i, '--', 'LineWidth', 1.2, 'DisplayName', sprintf('Fit Fr=%.2g', Fr_values(i)));
end
xlabel('$x$','Interpreter','latex'); ylabel('$\eta$','Interpreter','latex');
title('LES vs Fitted Wave Profiles'); legend('Location','bestoutside','NumColumns',2); grid on; box off;
print('validation_plot.png','-dpng','-r300');
fprintf('Plot saved: validation_plot.png\n');

%% ---- finished ----
fprintf('All done. Outputs: %s, %s, %s, %s\n', out_xlsx, tex_filename, txt_file, 'validation_plot.png');
fprintf('RMSE (all points): %.6g\n', RMSE_all);

%% ---- Helper functions ----
function s = num2latex(x)
    % scientific or simple formatting into LaTeX-friendly string (4 sig figs)
    if ~isfinite(x)
        s = '0';
        return;
    end
    fs = sprintf('%.4g', x); % compact
    epos = strfind(lower(fs), 'e');
    if isempty(epos)
        s = fs;
    else
        mant = fs(1:epos-1);
        expnt = str2double(fs(epos+1:end));
        s = sprintf('%s \\times 10^{%d}', mant, expnt);
    end
end

function s = quadFrTexTex(a2,a1,a0)
    % Create TeX-ready string for a2*Fr^2 + a1*Fr + a0, with signs & spacing
    parts = {};
    if abs(a2) > 0
        parts{end+1} = [num2latex(a2) ' Fr^2'];
    end
    if abs(a1) > 0
        if isempty(parts)
            parts{end+1} = [num2latex(a1) ' Fr'];
        else
            parts{end+1} = [signStr(a1) ' ' num2latex(abs(a1)) ' Fr'];
        end
    end
    if abs(a0) > 0 || isempty(parts)
        if isempty(parts)
            parts{end+1} = num2latex(a0);
        else
            parts{end+1} = [signStr(a0) ' ' num2latex(abs(a0))];
        end
    end
    % join parts into one expression and replace leading '+ ' if any
    expr = strjoin(parts, ' ');
    expr = strrep(expr, '+ -', '- ');
    % ensure unary plus is removed
    if startsWith(expr, '+ ')
        expr = expr(3:end);
    end
    s = expr;
end

function s = quadFrTexPlain(a2,a1,a0)
    % plain text version with typical decimal formatting (6 significant)
    fmt = @(v) sprintf('%.6g', v);
    parts = {};
    if abs(a2) > 0
        parts{end+1} = [fmt(a2) ' Fr^2'];
    end
    if abs(a1) > 0
        if isempty(parts)
            parts{end+1} = [fmt(a1) ' Fr'];
        else
            parts{end+1} = [signStr(a1) ' ' fmt(abs(a1)) ' Fr'];
        end
    end
    if abs(a0) > 0 || isempty(parts)
        if isempty(parts)
            parts{end+1} = fmt(a0);
        else
            parts{end+1} = [signStr(a0) ' ' fmt(abs(a0))];
        end
    end
    expr = strjoin(parts, ' ');
    expr = strrep(expr, '+ -', '- ');
    if startsWith(expr, '+ ')
        expr = expr(3:end);
    end
    s = expr;
end

function s = signStr(v)
    if v < 0
        s = '-';
    else
        s = '+';
    end
end
