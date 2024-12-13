function result = referenceDeconvoluteFT1(spectra_path, LW_Hertz, LB, TA, magnetStrength, show_plot)

%REFERENCEDECONVOLUTEFT1 Process and correct NMR spectra using reference deconvolution.
%
%   This function loads FT1 NMR data files, synthesizes a reference DSS peak,
%   applies a Lorentzian line broadening correction based on input parameters,
%   and performs exponential decay apodization. It can optionally display a
%   plot of the original, corrected, and apodized spectra for visual inspection.
%
%   Inputs:
%       spectra_path   - A string specifying the path to the spectra data
%                        files. Can be a single file path.
%                        Will accept only one file. Example: ./folder/file.ft1
%       LW_Hertz       - A scalar specifying the Lorentzian width in Hertz
%                        of the target DSS peak
%       LB             - A scalar specifying the line broadening factor in
%                        Hertz applied to the resulting spectra
%       TA             - A scalar specifying the percentage at which the
%                        ramp down of the trapezoidal apodization will be 
%                        applied. A value of 1 will eliminate the 
%                        apodization. The scalar should be between 0 and 1
%       magnetStrength - A scalar representing the magnet strength in MHz.
%       show_plot      - A boolean (optional, default true) to control whether
%                        the function plots the spectra comparison.
%
%   Outputs:
%       None. The function writes the processed spectra to a file.
%
%   Example Usage:
%       referenceDeconvoluteFT1('path/to/data', 1.1, 0.5, 0.75, 600, true);
%
%   The function generates a file with the processed spectrum using the parameters
%   provided, and optionally displays a comparison plot. The output filename is
%   automatically generated based on the input parameters and original file name.
%
%   Note:
%       This function depends on custom helper functions `load_ft1_files` and
%       `Setup1D` to load and prepare the data. Make sure these functions are
%       available in your MATLAB path.
%


    if nargin < 5  % If the number of inputs is less than 4, default show_plot to true
        show_plot = true;
    end

    if TA < 0 || TA > 1
        error('TA must be between 0 and 1');
    end

    % Load NMR data
    [ful_res_spectra, outputPath] = load_ft1_files(spectra_path);
    % Ensure only one file's data is loaded
    if numel(ful_res_spectra) > 1
        error('Multiple files detected. This function is designed to process only one file.');
    end

    % Put spectra into matrix
    [ful_res_X, X_ppm, ~] = Setup1D(ful_res_spectra);

    % Original spectrum and ppm values
    spectrum = ful_res_X;
    ppm = X_ppm;

    % Define and isolate DSS range, create synthetic DSS
    dssRange = [0 - 0.06, 0 + 0.06];
    dssIndices = find(ppm >= dssRange(1) & ppm <= dssRange(2));    
    zeroedSpectrum = zeros(size(spectrum));
    dssRegion = spectrum(dssIndices);
    zeroedSpectrum(dssIndices) = dssRegion;

    lowestNonZero = min(zeroedSpectrum(zeroedSpectrum > 0));
    adjustedDssRegion = max(zeroedSpectrum - lowestNonZero, 0);
    dssSpectrum = adjustedDssRegion;
    % figure();
    % plotr(ppm, dssSpectrum)

    % Synthesize DSS
    integral_value = trapz(ppm,dssSpectrum);
    x0 = 0;  % Peak center
    Gamma =LW_Hertz/magnetStrength;  % Adjust based on your sample spectrum (HWHM)
    A = integral_value;  % Area under the curve
    LorentzianPeak  = (A / pi) * (Gamma ./ ((ppm - x0).^2 + Gamma^2));
    % Apply the Hilbert transform to generate the complete reference signal
    % imag_LorentzianPeak = imag(hilbert(LorentzianPeak));
    % synthDssSpectrum = LorentzianPeak + 1i * imag_LorentzianPeak;
    synthDssSpectrum=LorentzianPeak;
    % plotr(ppm, synthDssSpectrum);

    % Apply corrections
    dssFID_mirrored = ifft(real(dssSpectrum));
    dssFID = ifft(real(dssSpectrum));
    % figure();
    % hold on;
    % plot(real(dssFID_mirrored))
    % disp(length(real(dssFID)))
    dssFID = [dssFID(1:end/2), zeros(1, length(dssFID)/2)];
    % plot(real(dssFID))
    
    % 
    dssSimFID_mirrored = ifft(real(synthDssSpectrum));
    dssSimFID = ifft(real(synthDssSpectrum));
    % figure();
    % hold on;
    % plot(real(dssSimFID_mirrored))
    dssSimFID = [dssSimFID(1:end/2), zeros(1, length(dssSimFID)/2)];
    % plot(real(dssSimFID))

    % % correctionFactor = dssSimFID ./ (dssFID + threshold);
    % 
    % % threshold = 1e-6;  % Threshold below which the denominator is considered too small
    % % correctionFactor = ones(size(dssFID));  % Initialize correctionFactor to zero
    % % validIndices = abs(dssFID) > threshold;  % Indices where dssFID is large enough
    % % correctionFactor(validIndices) = dssSimFID(validIndices) ./ dssFID(validIndices);
    % % % % Smoothing with a moving average filter
    % % windowSize = 100;  % Size of the moving window
    % % correctionFactorSmoothed = movmean(correctionFactor, windowSize);
    % 
    fullFID_mirrored = ifft(real(ful_res_X));
    fullFID = ifft(real(ful_res_X));
    fullFID = [fullFID(1:end/2), zeros(1, length(fullFID)/2)];
    % figure();
    % hold on;
    % plot(real(fullFID_mirrored))
    % plot(real(fullFID))

    
    threshold = 1e-20;  

    correctedFID = fullFID./ (dssFID+threshold) .* dssSimFID;
    correctedFID_mirrored = fullFID_mirrored./ (dssFID_mirrored+threshold) .* dssSimFID_mirrored;
    % 
    % % 
    correctedSpectrum_mirrored = fft(correctedFID_mirrored);
    correctedSpectrum = fft(correctedFID);
    % figure();
    % hold on;
    % plotr(ppm,correctedSpectrum_mirrored/2)
    % plotr(ppm,correctedSpectrum)
    % plotr(ppm,ful_res_X)
    
    % 
    % 
    sw = (max(X_ppm) - min(X_ppm)) * magnetStrength;
    t_half = linspace(0, length(fullFID)/2, length(fullFID)/2);
    expDecay_half = exp( -pi*t_half*LB / sw );
    expDecay_full = [expDecay_half, fliplr(expDecay_half(1:end))];
    % If the number of points is odd, adjust by adding one more element from the mirrored part
    if mod(length(fullFID), 2) == 1
        expDecay_full = [expDecay_half, fliplr(expDecay_half)];
    end
    expDecay = expDecay_full;

    % % Parameters for the trapezoidal apodization
    trapezoidalStart = TA;

    rampDownStart = ceil(trapezoidalStart * length(t_half));  % Start ramp down at 90% of half length
    rampDownLength = length(t_half) - rampDownStart;  % Length of ramp down
    % Create half trapezoidal apodization function
    plateau = ones(1, rampDownStart);
    rampDown = linspace(1, 0, rampDownLength);
    % Combine to form half trapezoid
    trapezoid_half = [plateau, rampDown];
    % Check if fullFID length is odd, adjust trapezoid accordingly
    if mod(length(fullFID), 2) == 1
        % Extend by repeating the middle value for symmetry
        trapezoid_half = [trapezoid_half, trapezoid_half(end)];
    end
    trapezoid_full = [trapezoid_half, fliplr(trapezoid_half)];
    trapezoid_full = trapezoid_full(1:length(fullFID));

    correctedFID_apodized = correctedFID .* expDecay;
    correctedFID_apodized = correctedFID_apodized .* trapezoid_full;


    correctedSpectrum_apodized = fft(correctedFID_apodized);
    % correctedSpectrum = fft(correctedFID);
    correctedSpectrum = correctedSpectrum_mirrored/2;
    % Prepare output file and spectrum plots
    if show_plot
        figure;
        hold on;
        plot(ppm, spectrum, 'b');  % Original spectrum
        % plot(ppm, correctedSpectrum, 'r');  % Corrected spectrum
        plot(ppm,correctedSpectrum,'r')
        title('Comparison of Spectra');
        legend('Original Spectrum', 'Apodized Corrected Spectrum');
        xlabel('ppm');
        ylabel('Intensity');
        set(gca, 'XDir','reverse');
    end

    % Write output file
    if isfolder(spectra_path)
        if contains(outputPath, '/')
            if outputPath(end) ~= '/' 
                input_fname = strcat(outputPath,'/',ful_res_spectra.FileName);
            else 
                input_fname = strcat(outputPath,ful_res_spectra.FileName);
            end
        end
        if contains(outputPath, '\')
            if outputPath(end) ~= '\'
                input_fname = strcat(outputPath,'\',ful_res_spectra.FileName);
            else 
                input_fname = strcat(outputPath,ful_res_spectra.FileName);
            end
        end
    else
        input_fname = spectra_path;
        [folderPath, ~, ~] = fileparts(spectra_path);
        outputPath = folderPath;
    end

    % input_fname 
    spectrum_in = fopen(input_fname, 'r');
    header_in = fread(spectrum_in, 512, 'float32');
    fclose(spectrum_in);
    [~, original_name, ext] = fileparts(ful_res_spectra.FileName);
    LW_Hertz_str = strrep(num2str(LW_Hertz), '.', '_');
    LB_str = strrep(num2str(LB), '.', '_');
    output_filename = fullfile(outputPath, [original_name, '_refDec_', num2str(LW_Hertz_str), 'HzLW_', num2str(LB_str), 'HzLB',ext]);
    fid = fopen(output_filename, 'w');  % 'b' for binary mode
    fwrite(fid, header_in, 'float32');
    wData = flip(correctedSpectrum_apodized');
    fwrite(fid, wData, 'float32');
    fclose(fid);

    result = correctedSpectrum_apodized;
end
