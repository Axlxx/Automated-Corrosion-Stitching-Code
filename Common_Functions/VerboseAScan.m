function VerboseAScan(scan, axes)
%VERBOSEASCAN Summary of this function goes here
%   Detailed explanation goes here
    % Plot signal & peaks
    s1 = subplot(2, 1, 1, axes);
    plot(s1, scan.time_vector(scan.initial_deadzone:length(scan.short_signal) + scan.initial_deadzone -1), scan.short_signal);
    hold on;
    plot(s1, scan.time_vector(scan.short_signal_peaks(2, :) + double(scan.initial_deadzone)), scan.short_signal_peaks(1, :), "*");
    % hb = abs(hilbert(scan.short_signal));
    % plot(s1, scan.time_vector(scan.initial_deadzone:length(scan.short_signal) + scan.initial_deadzone -1), hb);
    hold off;
    title('Relevant part of the original signal')
    xlabel('Time vector (s)');
    ylabel('Amplitude (displacement)');

    s2 = subplot(2, 1, 2);
    % Plot crosscorr
    plot(s2, scan.short_xcorr)
    hold on;
    hb = abs(hilbert(scan.short_xcorr));
    plot(s2, hb(1:end-100))
    plot(s2, scan.short_xcorr_peaks(2, :), scan.short_xcorr_peaks(1, :), "*");
    hold off;
    title('Autocorrelation to get tof')
    xlabel('Displacement (tau)');
    ylabel('Autocorrelation');
end

