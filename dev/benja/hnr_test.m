%%  Test script for WORLD analysis/synthesis with new waveform generator
% 2018/04/04: First version

%% Add path
addpath('../WORLD')
addpath('../python-PESQ')
addpath(genpath('../covarep'))
%% Script Variables
filenames = ["female_runn", "female_sust", "male_runn", "male_sust"];
% The form is {{objective_measure, [param_min,param_leap,param_max], param_show}}
parameters_values_dictionary =  {{"rd",[0.35,0.05,4], 14}, ...
                                {"f0_multiplier",[0.5,0.1,4], 7}, ...
                                {"rpp_k",[0,0.2,5], 5}, ...
                                {"jitter_amplitude",[0,2,50], 5}, ...
                                {"jitter_frequency",[0,100,22000], 44}, ...
                                {"shimmer_amplitude",[0,2,50], 5}, ...
                                {"shimmer_frequency",[0,100,22000], 44}}; ...
                                
hnr_command = "SMILExtract -C ../SMILExtract/config/emobase/emobase.conf -I ../outputs/tmp/y.wav -O ../outputs/hnr_test/";
%jit_shim_command = "SMILExtract -C ../SMILExtract/config/is09-13/IS10_paraling.conf -I ../outputs/tmp/y.wav -O ../outputs/analysis/";


%% Objective Measures generation
%
tic
for file = filenames
    disp('Reading: ../inputs/' + file + '.wav')
    
    % Audio Inputs
    [x, fs] = audioread('../inputs/' + file + '.wav');
    
    % WORLD Options
    dio_options = struct('f0_floor', 80,...
                         'f0_ceil', 400,...
                         'target_fs', fs,...
                         'frame_period', 2,... % Time in [ms] between 2 filters
                         'channels_in_octave', 24);

    cheaptrick_options = struct('q1', -0.15,...
                                'fft_size', 2048);

    D4C_options = struct('threshold', 0.85,...
                         'fft_size', 4096);

    seeds_options = struct('fft_size', 1024 * 2 ^ ceil(log2(fs / 48000)),...
                         'number_of_samples', 2 ^ ceil(log2(fs / 2)));
    
    % WORLD
    % DIO (F0 Estimation)
    f0_parameter = Dio(x, fs, dio_options);
    f0_parameter.f0 = StoneMask(x, fs, f0_parameter.temporal_positions, f0_parameter.f0);

    % CheapTrick (Spectral Envelope Estimation)
    spectrum_parameter = CheapTrick(x, fs, f0_parameter, cheaptrick_options);

    % D4C (Aperiodicity Estimation)
    source_parameter = D4CRequiem(x, fs, f0_parameter, D4C_options);

    % Seeds
    seeds_signals = GetSeedsSignals(fs, seeds_options);
    

    %y = SynthesisRequiem(source_parameter, spectrum_parameter, seeds_signals, vq_options, 'off');
    
    for parameter_values = parameters_values_dictionary
        param_show_counter = 0;
        param_show = parameter_values{1}{3};
        fprintf('Objective Measure: ' + parameter_values{1}{1} + ' -> ')
            
        % Definition of parameters that need to be defined so that
        % parameter_values{1}{1} works
        vq_options = struct();

        if (parameter_values{1}{1} == 'rpp_k')
            if (file == 'female_runn' || file == 'female_sust')
                vq_options.rd_param = 0.5;
                vq_options.spectrum_filtering_exponential = 12;
                vq_options.spectrum_filtering_samples = 75;
            else
                vq_options.rd_param = 0.35;
                vq_options.spectrum_filtering_exponential = 8.5;
                vq_options.spectrum_filtering_samples = 45;
            end
        elseif (parameter_values{1}{1} == 'rd')
            if (file == 'female_runn' || file == 'female_sust')
                vq_options.spectrum_filtering_exponential = 12;
                vq_options.spectrum_filtering_samples = 75;
            else
                vq_options.spectrum_filtering_exponential = 8.5;
                vq_options.spectrum_filtering_samples = 45;
            end
        elseif (parameter_values{1}{1} == 'jitter_amplitude')
            if (file == 'female_runn' || file == 'female_sust')
                vq_options.jitter_frequency = 50*22;
            else
                vq_options.jitter_frequency = 50*10;
            end
        elseif (parameter_values{1}{1} == 'jitter_frequency')
            if (file == 'female_runn' || file == 'female_sust')
                vq_options.jitter_amplitude = 38;
            else
                vq_options.jitter_amplitude = 38;
            end
        elseif (parameter_values{1}{1} == 'shimmer_amplitude')
            if (file == 'female_runn' || file == 'female_sust')
                vq_options.shimmer_frequency = 22;
            else
                vq_options.shimmer_frequency = 10;
            end
        elseif (parameter_values{1}{1} == 'shimmer_frequency')
            if (file == 'female_runn' || file == 'female_sust')
                vq_options.shimmer_amplitude = 38;
            else
                vq_options.shimmer_amplitude = 38;
            end
        end
        
        % Parameters under test
        param_min = parameter_values{1}{2}(1);
        param_leap = parameter_values{1}{2}(2);
        param_max = parameter_values{1}{2}(3);
        
        %Debug version
        %param_min = parameter_values{1}{2}(1);
        %param_leap = parameter_values{1}{2}(1);
        %param_max = parameter_values{1}{2}(1);

        % Temporary Variables
        %data_out_se = [];
        %data_out_cpp = [];
        %data_out_ps = [];
        %data_out_pesq = [];
        
        for param = param_min:param_leap:param_max
            
            param_show_counter = param_show_counter + 1;
            if mod(param_show_counter, param_show) == 0
               fprintf('%.2f, ',param)
            end
        
            if parameter_values{1}{1} == "rd"
                vq_options.rd_param = param;
            elseif parameter_values{1}{1} == "f0_multiplier"
                vq_options.f0_multiplier = param;
            elseif parameter_values{1}{1} == "rpp_k"
                vq_options.rpp_k = param;
            elseif parameter_values{1}{1} == "jitter_amplitude"
                vq_options.jitter_amplitude = param;
            elseif parameter_values{1}{1} == "jitter_frequency"
                vq_options.jitter_frequency = param;
            elseif parameter_values{1}{1} == "shimmer_amplitude"
                vq_options.shimmer_amplitude = param;
            elseif parameter_values{1}{1} == "shimmer_frequency"
                vq_options.shimmer_frequency = param;
            end

            y = SynthesisRequiem(source_parameter, spectrum_parameter, seeds_signals, vq_options, 'off');

            % Audiowrites
            audiowrite('../outputs/tmp/x.wav',x./max(abs(x)),fs)
            audiowrite('../outputs/tmp/y.wav',y./max(abs(y)),fs)

            % Objective Measurements
            %cpp_y = cpp_fix(y,fs,0,'line',1);
            %ps_y = peakslope(y,fs);
            %pesq_value = pesq('../outputs/tmp/x.wav','../outputs/tmp/y.wav');
            
            %if (file ~= "female_runn" && file ~= "male_runn" && parameter_values{1}{1} ~= "jitter_amplitude" && parameter_values{1}{1} ~= "jitter_frequency" && parameter_values{1}{1} ~= "shimmer_amplitude" && parameter_values{1}{1} ~= "shimmer_frequency" )
            %    se = SpectralEnvelopeH1H2_covarep(y(5000:5000+fs/4),fs);
            %end
            % SmileExtract
            [~,~] = system(hnr_command + parameter_values{1}{1} + "_hnr_" + file + ".arff");
            %[~,~] = system(jit_shim_command + parameter_values{1}{1} + "_jitter_and_shimmer_" + file + ".arff");


            % This line concatenates the generated data in the data dumps lists
            %data_out_cpp = [data_out_cpp mean(cpp_y(:,1))];
            %data_out_ps = [data_out_ps mean(ps_y(:,2))];
            %data_out_pesq = [data_out_pesq pesq_value];
            %if (file ~= "female_runn" && file ~= "male_runn" && parameter_values{1}{1} ~= "jitter_amplitude" && parameter_values{1}{1} ~= "jitter_frequency" && parameter_values{1}{1} ~= "shimmer_amplitude" && parameter_values{1}{1} ~= "shimmer_frequency" )
%                data_out_se = [data_out_se se];
 %           end
            
        end
        
        fprintf('\n')
        
        %csvwrite("../outputs/analysis/" +  parameter_values{1}{1} + "_cpp_" + file + ".csv", data_out_cpp)
        %csvwrite("../outputs/analysis/" +  parameter_values{1}{1} + "_ps_" + file + ".csv", data_out_ps)
        %csvwrite("../outputs/analysis/" +  parameter_values{1}{1} + "_pesq_" + file + ".csv", data_out_pesq)
        %if (file ~= "female_runn" && file ~= "male_runn" && parameter_values{1}{1} ~= "jitter_amplitude" && parameter_values{1}{1} ~= "jitter_frequency" && parameter_values{1}{1} ~= "shimmer_amplitude" && parameter_values{1}{1} ~= "shimmer_frequency" )
%            csvwrite("../outputs/analysis/" +  parameter_values{1}{1} + "_se_" + file + ".csv", data_out_se)
%        end
    end
end
toc
%}