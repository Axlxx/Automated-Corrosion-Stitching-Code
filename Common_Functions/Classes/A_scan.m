classdef A_scan
    %A_SCAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % The signal received from the A-scan
        signal double                       % array with the signal 
        short_signal double
        short_xcorr double
        initial_deadzone (1,1) int32
        time_vector double                  % array with time at which each individual 
                                            %     sample point from signal array was collected

        % Variables related to the signal
        len (1,1) int32                     % length of the signal
        timestamp (1,1) datetime            % time when the scan was taken (when scan was completed)
        tof (1,1)                           % time of flight between peaks 
        short_signal_peaks (2, 2) double    % row 1 = peak amplitude, row 2 = peak index 
        short_xcorr_peaks (2, 2) double     % row 1 = peak amplitude, row 2 = peak index 
        depth (1,1)                         % depth measured (tof/2 * mat_c)
        no_averages (1,1) int32             % how many times the scan has been made and averaged

        signal_amplitude (1,1)              % initial signal amplitude
        signal_no_cycles (1,1)              % how many cycles were fired
        signal_central_freq (1,1)           % central frequency of the packet (w or f)

        % Variables related to the transducer
        trans_sampl_freq (1,1)              % transducer sampling frequency
        transducer_operating_freq (1,1)     % transducer maximum frequency
        
        trans_width (1,1) double            % transducer width
        trans_length (1,1) double           % transducer length

        trans_longitudinal (1,1) logical    % if transducer is longitudinal

        % Variables related to the materials that were scanned
        mat_c (1,1)          % material speed
        mat_t (1,1)          % material temperature
        

        % Variables related to the environment in which the scan happened 
        env_X (1,1)                     % transducer centrepoint X
        env_Y (1,1)                     % transducer centrepoint Y
        env_Z (1,1)                     % transducer centrepoint Z

        env_t (1,1)                     % environment temperature 

        env_min_expected (1,1)          % minimum expected depth
        env_max_expected (1,1)          % maximum expected depth
    end
    
    methods

        %%% Constructors
        % Potential inputs
        %   varargin{1} - signal
        %   varargin{2} - trans_sampl_freq
        %   varargin{3} - 
        %   varargin{4} - 
        %   varargin{5} - 
        %   varargin{6} - 
        %   varargin{7} -  
        %   varargin{8} - 
        function obj = A_scan(varargin)

            if nargin > 0
                signal = varargin{1};
    
                obj.signal = signal;
                obj.len = length(signal);    
            end

            if nargin > 1
                % Set sampl freq
                sampl_freq = varargin{2};
                obj.trans_sampl_freq = sampl_freq;

                time_vec = double(1:obj.len);
                time_between_samples = double(1/obj.trans_sampl_freq);
                time_vec = double(time_vec * time_between_samples);
                obj.time_vector = time_vec;
            end

            obj = obj.set_coordinates(-100, -100, -100);
            obj = obj.set_timestamp();
            obj.trans_longitudinal = true;
            obj = obj.set_transducer();
            
%             obj = obj.calculate_tof();
%             obj = obj.calculate_SNR(1);

        end

        % Function to set signal 
        function obj = set_signal(obj, varargin)
            obj.signal = varargin{1};
            obj = obj.calculate_tof();

            SNR = 1;
            if nargin > 2
                SNR = varargin{2};
            end
            obj = obj.calculate_SNR(SNR);
        end
        
        % Function to set timestamp
        function obj = set_timestamp(obj, varargin)
            if nargin == 1
                obj.timestamp = datetime("now");
            elseif nargin == 2
                date = varargin{1};
                obj.timestamp = date;
            end
        end

        % Function to set X, Y, Z coordinates
        function obj = set_coordinates(obj, varargin)
            if nargin < 3 || nargin > 4
                error('Only enter X, Y, (and optionally) Z as parameters');
            end

            obj.env_X = varargin{1};
            obj.env_Y = varargin{2};

            if nargin == 4
                obj.env_Z = varargin{3};
            end
        end

        % Function to set transducer properties
        % Input
        %   @param1 - transducer_width
        %   @param2 - transducer_length
        %   @param3 - transducer_longitudinal T/F (if it's longitudinal or shear)   
        %   @param4 - transducer_operating_freq
        function obj = set_transducer(obj, varargin)
            if nargin > 1
                obj.trans_width = varargin{1};
            end

            if nargin > 2
                obj.trans_length = varargin{2};
            end  

            if nargin > 3
                obj.trans_longitudinal = varargin{3};
            end

            if nargin > 4
                obj.transducer_operating_freq = varargin{4};
            end           
        end

        % Function to set material speed
        function obj = set_c(obj, input_mat_c)
            obj.mat_c = input_mat_c;
        end

        % Function to find out material speed
        function obj = calibrate_c(obj, depth)
            if (obj.tof <= 0)
                obj.tof = calculate_tof(obj);
            end
        
            %   depth = time/2 * c_mat 
            % material speed = depth / (1/sampling freq x time of flight) 
            obj.mat_c = depth / (1/obj.trans_sampl_freq * obj.tof);
        end

        % Function to calculate tof from signal
        function obj = calculate_tof(obj, varargin)
            try
                if isempty(obj.signal)
                    return
                end

                %%% Input questions
                peak_deadzone = 125;
                obj.initial_deadzone = 250;
                max_signal = obj.initial_deadzone + 600 -1;
                if max_signal > length(obj.signal)
                    max_signal = length(obj.signal);
                end
                noiseless = false;
    
                if nargin > 1
                    
                    if varargin{2} > length(obj.signal)
                        warning('peak deadzone larger than the signal. please input a smaller initial_deadzone')
                    else 
                        peak_deadzone = varargin{1};
                    end
                end
                if nargin > 2
                    obj.initial_deadzone = varargin{2};
                    if obj.initial_deadzone <= 0
                        obj.initial_deadzone = 1;
                    end

                    if obj.initial_deadzone > length(obj.signal)
                        error('initial deadzone larger than the signal. please input a smaller initial_deadzone')
                    end
                end
                if nargin > 3
                    max_signal = varargin{3};
                    if max_signal > length(obj.signal)
                        max_signal = length(obj.signal);
                    end
                end


                %%% Step 1: Get short_signal containing first and second peak
                
                % try
                %     temp_short = obj.signal(obj.initial_deadzone:max_signal);  
                %     figure(1);
                %     subplot(2, 1, 1)
                %     plot(temp_short);
                %     temp_short = abs(hilbert(temp_short));
                %     temp_short(temp_short<max(temp_short)*0.2) = 0;
                %     subplot(2, 1, 2)
                %     plot(temp_short);
                % 
                % 
                %     % Find largest peak, take peak_deadzone units before
                    % [obj.short_signal_peaks(1, 1), obj.short_signal_peaks(2, 1)] = max(temp_short);
                %     % Remove amplitude around it
                %     remove_start = obj.short_signal_peaks(2, 1) - peak_deadzone;
                %     remove_end = obj.short_signal_peaks(2, 1) - peak_deadzone;
                %     if remove_start < 1
                %         remove_start = 1;
                %     end
                %     if remove_end > length(temp_short)
                %         remove_end = length(temp_short);
                %     end
                % 
                %     index_start_1 = remove_start;
                %     temp_short(remove_start:remove_end) = 0;
                % 
                %     % Find second largest peak not in the first peak's
                %     % deadzone, get peak_deadzone units after it.
                %     [obj.short_signal_peaks(1, 2), obj.short_signal_peaks(2, 2)] = max(temp_short(1:end));
                %     % if obj.short_signal_peaks(2, 1) > obj.short_signal_peaks(2, 2) % Second signal was found behind the first
                %     %     temp_short(-1);
                %     % end
                %     if obj.short_signal_peaks(2, 2) + peak_deadzone > length(temp_short)
                %         index_end = length(temp_short);
                %     else
                %         index_end = obj.short_signal_peaks(2, 2)+peak_deadzone;
                %     end
                %     % Select that as f to autocorrelate
                %     obj.short_signal = obj.signal(index_start+obj.initial_deadzone:index_end+obj.initial_deadzone);
                % catch ex
                %     fprintf('Second signal was found behind the first (signal)');
                %     rethrow(ex);
                % end
                obj.short_signal = obj.signal(obj.initial_deadzone:max_signal);
                [obj.short_signal_peaks(1, 1), obj.short_signal_peaks(2, 1)] = max(obj.short_signal);
 
                %%% Step 2: Apply cross correlation on short_signal
                [obj.short_xcorr, ~] = xcorr(obj.short_signal);
                obj.short_xcorr = obj.short_xcorr((length(obj.short_xcorr)+1)/2:end);
                hb = abs(hilbert(obj.short_xcorr));

                %%% Step 3: Get max value of cross corr  except 1 and subtract them
                obj.short_xcorr_peaks(1, 1) = hb(1);
                obj.short_xcorr_peaks(2, 1)= 1;
                [obj.short_xcorr_peaks(1, 2), obj.short_xcorr_peaks(2, 2)] = max(hb(peak_deadzone:end-100));

                
                obj.short_signal_peaks(2, 2) = obj.short_signal_peaks(2, 1) + obj.short_xcorr_peaks(2, 2);
                if obj.short_signal_peaks(2, 2) > length(obj.short_signal)
                    obj.short_signal_peaks(2, 2) = obj.short_signal_peaks(2, 1) - obj.short_xcorr_peaks(2, 2);
                end
                obj.short_signal_peaks(1, 2) = obj.short_signal(obj.short_signal_peaks(2, 2));


                %%% Optional: Step 3.5 interpolate peaks.
                try
                    % obj.short_xcorr_peaks(2, 1) = A_scan.interp_gradient(obj.short_xcorr, obj.short_xcorr_peaks(2, 1));
                    tried = A_scan.interp_gradient(hb, obj.short_xcorr_peaks(2, 2));
                    if abs(tried - obj.short_xcorr_peaks(2, 2)) < 1
                        obj.short_xcorr_peaks(2, 2) = tried;
                    end
                catch ex
                    fprintf('Gradient Interpolation failed');
                end

                try
                    % % If not within predicted range of +- 10% 
                    % if abs(obj.peaks(3) - obj.peaks(1)) > 1.1 * 2 * abs(obj.peaks(2) - obj.peaks(1))
                    obj.tof = abs(obj.short_xcorr_peaks(2, 2) - obj.short_xcorr_peaks(2, 1));

                catch ex
                    fprintf('Not even 1 peak was identified here')
                    rethrow(ex);
                end
            catch ex
                obj.tof = -1;
                warning('\nTOF caused an error on this A_scan\n');

                figure(2);
                plot(obj.signal(obj.initial_deadzone:max_signal));
                hold on;
                plot(obj.short_signal_peaks(2, :), obj.short_signal_peaks(1, :));
                hold off;
                title(sprintf('Exception thrown on coordinates: %d %d', obj.env_X, obj.env_Y ));

                fprintf(ex.message);
                rethrow(ex);
            end
        end

        % Function to calculate depth from tof
        function obj = calculate_depth(obj, varargin)
            % Use object c_mat
            if nargin == 1
                if obj.mat_c == 0
                    error('speed of sound in material of the object is 0. Set mat_c to a positive non-zero value or call the function with the desired mat_c as parameter');
                end
                % Calculation 
            % Given c_mat
            elseif nargin >= 2
                obj.mat_c = varargin{1};
            end

            if nargin >= 3
                time_between_samples = varargin{2};
            else
                time_between_samples = 1 / obj.trans_sampl_freq;
            end
        
            if obj.tof == 0
                obj = obj.calculate_tof() ; % Calculated value
            end

            obj.depth = (double(obj.tof/2) * time_between_samples) * obj.mat_c; 
        end

        % TODO
        function obj = calculate_SNR(obj, method)
            if method == 1 % Coherent SNR

            elseif method == 2 % Random SNR

            else 
                error('No such method. Please input 1 for coherent SNR and 2 for random SNR');
            end 
        end

    end

    methods(Static)            
        % Inputs:
        %    signal - signals
        %    i - the intdex at which interpolation happens
        % Outputs:
        %    int_p - interpolated x
        function int_x = interp_gradient(signal, i)

            if i <= 1 || i >= length(signal)-1
                error('Index out of array')
            end
        
            grad1 = gradient(signal);
            % Get the 2 lines
            a1 = signal(i-1) - (i-1) * grad1(i-1);
            a2 = signal(i+1) - (i+1) * grad1(i+1);
            
            % intersection of 2 lines
            int_x = (a2 - a1) / (grad1(i-1) - grad1(i+1));
        end

        % Inputs:
        %    signal - signal
        %    percentage - the percentage under which signal is cut off
        % Outputs:
        %    int_p - noise free signal
        function signal_noise_free = noise_free_array(signal, percentage) 
            if percentage < 0 || percentage > 100
                error('Percentage not percentaging (not a number between 0 and 100)')
            end
        
            % Find max value 
            max_amplitude = max(signal);
            percentage = percentage / 100;

            % Eliminate signal points that are (in absolute value)
            % lower than percentage * max_amplitude
            noise_mask = (abs(signal) > max_amplitude * percentage);
            signal_noise_free = signal .* noise_mask;
        end       
    end
end

