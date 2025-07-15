% Generate strips from real scanned data using verasonics array
% stripsFromRealData(steel_db, aperture, 6000 * 10^3, 62.5 * 10^6, false)
function strips = stripsFromRealData(steel_db, aperture, c_material, sampling_freq, verbose)
    miss_counter = 0;
    try 
        % Find all values in dataset
        % strips = A_scan(zeros(1, 1, 1));
        % strips = createArray(1,)
        % strips = 0;
        noStrips = 0;
        for i = 1:(size(steel_db, 1)-1)
            indeces = find(~cellfun(@isempty,steel_db(i, :)));
            no_x = length(indeces);
        
            if no_x
                % A strip 
                noStrips = noStrips + 1;
                for j = 1:no_x
                    % Cscan = steel_db{i, indeces(j)};
                    Cscan = steel_db{i, indeces(j)};

                    for k = 1:(64-aperture)  % 1 to 65 
                        signals = zeros(aperture, size(Cscan, 2));
                        for k_ap = 1:aperture
                            signals(k_ap, :) = Cscan(k + k_ap - 1, :);
                        end
                        signal = mean(signals, 1); 

                        scan = A_scan(signal, sampling_freq);
                        % scan = scan.set_timestamp(Cscan.timestamp);
                        scan = scan.set_coordinates(i, indeces(j), 10);
                        % scan = scan.set_transducer(Cscan.trans_width, Cscan.trans_length, Cscan.trans_longitudinal, Cscan.transducer_operating_freq)
                        scan = scan.set_transducer(15, 0.5, 1, 5*10^6);
                        scan = scan.set_c(c_material);
                        try
                            scan = scan.calculate_depth();

                            if scan.depth > 10.5 || scan.depth < 7.5
                                if verbose
                                    clf(figure(11))
                                    figure(11);
                                    subplot(2, 1, 1);
                                    plot(scan.signal(200:1000));
                    
                                    subplot(2, 1, 2);
                                    plot(scan.short_signal);
                                    hold on;
                                    plot(scan.short_signal_peaks(2, 1), scan.short_signal_peaks(1, 1), '*r');
                                    plot(scan.short_signal_peaks(2, 2), scan.short_signal_peaks(1, 2), '*r'); 
                                    hold off;
                                end
            
                                scan = calculate_tof(scan, 125, 300, 900);
                                time_between_samples = 1 / scan.trans_sampl_freq;
                                scan.depth = (double(scan.tof/2) * time_between_samples) * scan.mat_c;
                
                             if scan.depth > 12 || scan.depth < 7
                                miss_counter = miss_counter + 1;
                                 if verbose
                                    clf(figure(12))
                                    figure(12);
                                    subplot(2, 1, 1);
                                    plot(scan.signal(200:1000));
                    
                                    subplot(2, 1, 2);
                                    plot(scan.short_signal);
                                    hold on;
                                    plot(scan.short_signal_peaks(2, 1), scan.short_signal_peaks(1, 1), '*r');
                                    plot(scan.short_signal_peaks(2, 2), scan.short_signal_peaks(1, 2), '*r'); 
                                    
                                    sgtitle(sprintf("Depth: %0.04f", scan.depth));
                                    hold off;
                                 end
                             end
                            end



                            if verbose
                                pause(0.2)
                                figure(1);
                                % Plot signal & peaks
                                subplot(2, 1, 1);
                                plot(scan.time_vector(scan.initial_deadzone:length(scan.short_signal) + scan.initial_deadzone -1), scan.short_signal);
                                hold on;
                                plot(scan.time_vector(scan.short_signal_peaks(2, :) + double(scan.initial_deadzone)), scan.short_signal_peaks(1, :), "*");
                                % hb = abs(hilbert(scan.short_signal));
                                % plot(scan.time_vector(scan.initial_deadzone:length(scan.short_signal) + scan.initial_deadzone -1), hb);
                                hold off;
    
                                subplot(2, 1, 2);
                                % Plot crosscorr
                                plot(scan.short_xcorr)
                                hold on;
                                hb = abs(hilbert(scan.short_xcorr));
                                plot(hb(1:end-100))
                                plot(scan.short_xcorr_peaks(2, :), scan.short_xcorr_peaks(1, :), "*");
                                hold off;
                            end
                            strips(noStrips, j, k) = scan.depth;   

                        catch ex
                            fprintf("coords: i=%d, j=%d, k=%d", i, indeces(j), k); 
                            strips(noStrips, j, k) = 0;  

                            % if k > 1
                            %     strips(noStrips, j, k) = strips(noStrips, j, k-1);
                            % else
                            %     strips(noStrips, j, k) = strips(noStrips, j-1, 1);
                            % end
                        end
                        
                    end
                end
                % noStrips  
                % i
            end
        end
    catch ex
        fprintf(ex.message);
    end
    fprintf("Miss counter: %d", miss_counter)

end