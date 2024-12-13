function [ref_spectra, new_ppm] = reference1D(ful_res_X, X_ppm, ref_max_ppm, ref_min_ppm)
% reference1D - Aligns NMR spectra based on a reference peak within a specified ppm range.
%
% Syntax: [ref_spectra, new_ppm] = reference1D(ful_res_X, X_ppm, ref_max_ppm, ref_min_ppm)
%
% Inputs:
%    ful_res_X - A matrix where each row represents a 1D NMR spectrum.
%    X_ppm - A vector representing the ppm values corresponding to the columns of ful_res_X.
%    ref_max_ppm - The maximum ppm value defining the range to search for the reference peak.
%    ref_min_ppm - The minimum ppm value defining the range to search for the reference peak.
%
% Outputs:
%    ref_spectra - A matrix of the aligned NMR spectra with a common ppm vector.
%    new_ppm - The common ppm vector for the aligned spectra.
%
% Notes:
%    - This function assumes that the input spectra are already processed and ready for alignment.
%    - The specified ppm range should contain a prominent DSS peak that can be used as a reference for alignment.
%    - The function uses linear interpolation to align the spectra to a common ppm vector.
%
% Author: Leandro Ponce (chatgpt helped :P )
% University of Georgia
% CCRC Edison Lab
% Date: June 10, 2024
%
% See also: findpeaks, interp1
    % Initialize variables
    [num_samples, ~] = size(ful_res_X);
    new_ppm_cell = cell(num_samples, 1);

    % Loop through each sample
    for i = 1:num_samples
        % Extract the spectrum for the current sample
        spectrum = ful_res_X(i, :);

        % Identify the reference region indices
        ref_indices = X_ppm >= ref_min_ppm & X_ppm <= ref_max_ppm;
        ref_region = spectrum(ref_indices);
        ref_ppm_region = X_ppm(ref_indices);

        % Detect peaks in the reference region
        [peaks, locs] = findpeaks(ref_region);
        
        % Find the maximum peak in the reference region
        [~, max_index] = max(peaks);
        ref_peak_index = locs(max_index);
        ref_peak_ppm = ref_ppm_region(ref_peak_index);

        % Plot the spectrum and identified peaks
        % figure;
        % plot(X_ppm, spectrum);
        % hold on;
        % plot(ref_ppm_region(locs), peaks, 'ro'); % Mark the peaks with red circles
        % plot(ref_peak_ppm, peaks(max_index), 'go'); % Mark the reference peak with a green circle
        % hold off;
        % title(['Sample ' num2str(i) ': Spectrum with Identified Peaks']);
        % xlabel('PPM');
        % ylabel('Intensity');
        % legend('Spectrum', 'Identified Peaks', 'Reference Peak');
        % set(gca, 'XDir', 'reverse'); % Reverse x-axis for NMR spectra

        % Calculate the shift needed to set the reference peak to 0 ppm
        ppm_shift = ref_peak_ppm;

        % Adjust the ppm vector for the current sample
        new_ppm_cell{i} = X_ppm - ppm_shift;
    end

    % Determine a common ppm vector
    common_ppm_min = max(cellfun(@min, new_ppm_cell));
    common_ppm_max = min(cellfun(@max, new_ppm_cell));
    new_ppm = linspace(common_ppm_min, common_ppm_max, length(X_ppm));

    % Interpolate each spectrum to the common ppm vector
    ref_spectra = zeros(num_samples, length(new_ppm));
    for i = 1:num_samples
        ref_spectra(i, :) = interp1(new_ppm_cell{i}, ful_res_X(i, :), new_ppm, 'linear', 'extrap');
    end
end
