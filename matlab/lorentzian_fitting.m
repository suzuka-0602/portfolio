function results = lorentzian_fitting(dataFolder, configFcn)
%LORENTZIAN_FITTING Fit Lorentzian peaks for configured spectra.
%   RESULTS = LORENTZIAN_FITTING(DATAFOLDER, CONFIGFCN) loads .spe spectra
%   from DATAFOLDER, applies Caliblation() to map indices to cm^-1, and
%   fits Lorentzian peaks defined by CONFIGFCN (default: @lorentzian_regions).
%   The function saves the amplitudes (A_i) for each spectrum into
%   lorentzian_results.mat and returns a struct RESULTS summarizing each fit.
%
%   Each region in CONFIGFCN must include:
%     spectrum          - filename of the .spe file within DATAFOLDER.
%     center_cm1        - central wavenumber for the fitting window.
%     window_cm1        - half-window width around center_cm1 to fit.
%     peak_centers_cm1  - 1xN array with initial center estimates per peak.
%     initial_amplitude - 1xN array with initial amplitudes (optional).
%     initial_gamma_cm1 - 1xN array with initial HWHM values (optional).
%
%   The amplitudes are exposed as RESULTS.<spectrum>.A_k fields (k=1..N)
%   and also stored inside RESULTS.<spectrum>.amplitudes.

if nargin < 1 || isempty(dataFolder)
    dataFolder = '.';
end

if nargin < 2 || isempty(configFcn)
    configFcn = @lorentzian_regions;
elseif ischar(configFcn) || isstring(configFcn)
    configFcn = str2func(configFcn);
end

regions = configFcn();
if isempty(regions)
    error('lorentzian_fitting:NoRegions', 'No fitting regions returned by configuration.');
end

results = struct();

for idx = 1:numel(regions)
    region = regions(idx);
    spectrumPath = fullfile(dataFolder, region.spectrum);
    [wavenumbers, intensity] = load_spectrum_with_calibration(spectrumPath);

    [xFit, yFit] = limit_region(wavenumbers, intensity, region.center_cm1, region.window_cm1);
    params0 = build_initial_guess(region);
    [lb, ub] = build_bounds(region, params0);

    options = optimoptions('lsqcurvefit', 'Display', 'off');
    model = @(p, x) multi_lorentzian(p, x);
    params = lsqcurvefit(model, params0, xFit, yFit, lb, ub, options);

    amplitudes = params(1:3:end);
    fileStem = matlab.lang.makeValidName(strip_extension(region.spectrum));

    for k = 1:numel(amplitudes)
        label = sprintf('lorentz_A_%d', k);
        results.(fileStem).(label) = amplitudes(k);
    end

    results.(fileStem).amplitudes = amplitudes;
    results.(fileStem).parameters = params;
    results.(fileStem).fitted_range_cm1 = [min(xFit), max(xFit)];
    results.(fileStem).center_cm1 = region.center_cm1;
    results.(fileStem).window_cm1 = region.window_cm1;
end

save(fullfile(dataFolder, 'lorentzian_results.mat'), 'results');
end

%% Local helpers
function params0 = build_initial_guess(region)
    numPeaks = numel(region.peak_centers_cm1);

    if isfield(region, 'initial_amplitude') && ~isempty(region.initial_amplitude)
        amps = region.initial_amplitude;
    else
        % Estimate amplitude using the max intensity scaled down by peak count.
        amps = repmat(1, 1, numPeaks);
    end

    if isfield(region, 'initial_gamma_cm1') && ~isempty(region.initial_gamma_cm1)
        gammas = region.initial_gamma_cm1;
    else
        % Roughly set HWHM to 1/10 of the window.
        gammas = repmat(region.window_cm1 ./ 10, 1, numPeaks);
    end

    params0 = zeros(1, 3 * numPeaks);
    for i = 1:numPeaks
        baseIdx = 3 * (i - 1);
        params0(baseIdx + 1) = amps(i);
        params0(baseIdx + 2) = region.peak_centers_cm1(i);
        params0(baseIdx + 3) = gammas(i);
    end

    % Clip negative guesses to small positive numbers to aid the optimizer.
    params0 = max(params0, eps);
end

function [lb, ub] = build_bounds(region, params0)
    numPeaks = numel(region.peak_centers_cm1);
    lb = zeros(size(params0));
    ub = inf(size(params0));

    % Constrain centers within the fitting window.
    for i = 1:numPeaks
        baseIdx = 3 * (i - 1);
        lb(baseIdx + 2) = region.center_cm1 - region.window_cm1;
        ub(baseIdx + 2) = region.center_cm1 + region.window_cm1;

        % Gamma must be positive but allow up to the window width.
        lb(baseIdx + 3) = 0;
        ub(baseIdx + 3) = max(region.window_cm1, params0(baseIdx + 3) * 10);
    end
end

function y = multi_lorentzian(params, x)
    % Sum of Lorentzian profiles: A / (1 + ((x - x0)/gamma).^2)
    numPeaks = numel(params) / 3;
    y = zeros(size(x));
    for i = 1:numPeaks
        baseIdx = 3 * (i - 1);
        A = params(baseIdx + 1);
        x0 = params(baseIdx + 2);
        gamma = params(baseIdx + 3);
        y = y + A ./ (1 + ((x - x0) ./ gamma) .^ 2);
    end
end

function [xRange, yRange] = limit_region(x, y, center, window)
    mask = x >= (center - window) & x <= (center + window);
    xRange = x(mask);
    yRange = y(mask);
    if isempty(xRange)
        error('lorentzian_fitting:EmptyWindow', 'No data points within the specified window.');
    end
end

function [wavenumbers, intensity] = load_spectrum_with_calibration(path)
    raw = readmatrix(path);
    if size(raw, 2) == 1
        intensity = raw(:);
    else
        % Assume intensity is the last column when multiple columns are present.
        intensity = raw(:, end);
    end

    wavenumbers = map_indices_to_cm1(numel(intensity));
end

function cm1 = map_indices_to_cm1(numPoints)
    indices = (1:numPoints)';

    if exist('Caliblation', 'file')
        % Preferred: Caliblation returns the calibrated wavenumber axis when
        % provided pixel indices.
        try
            cm1 = Caliblation(indices);
            if numel(cm1) >= numPoints
                cm1 = cm1(:);
                return;
            end
        catch
            % Some Caliblation scripts may populate a workspace variable instead
            % of returning a value; execute and capture wn_calib when present.
            try
                Caliblation;
                if exist('wn_calib', 'var') && numel(wn_calib) >= numPoints
                    cm1 = wn_calib(:);
                    return;
                end
            catch
                % Continue to fallback if execution fails.
            end
        end
    end

    warning(['Using index-based wavenumber axis; please ensure Caliblation.m ' ...
             'is on the path and returns wn_calib for your spectra.']);
    cm1 = indices;
end

function name = strip_extension(file)
    [~, name, ~] = fileparts(file);
end
