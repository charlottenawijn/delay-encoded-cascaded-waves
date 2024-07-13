% Optimize the gain in SNR based on the Fourier-domain conding function to decode a
% delay-encoded cascaded ultrasound transmission within a given range of
% delays for N pulses within a pulse train.

% Charlotte Nawijn and Guillaume Lajoinie, University of Twente, 2023-2024


%% general
clear
close all
clc

%% settings (stored in struct S)
S.plot_optimization_result = 1;                 % make figure of all results of optimization for each number of pulses and sacrificial_depth
Nshow = 5;                                      % how many of the conding function options to show after the optimization

S.attempts = 20;                                % number of random attempts for which to perform optimization

P.sacrificial_depth = [10, 15, 20]*1e-3;        % sacrificial depth (m); dictates maximum length pulse train


% uncomment which transducer to simulate
S.probe = 'P4-1';
% S.probe = 'C1-6';
% S.probe = 'L12-3v';
% S.probe = 'L22-8v';


%% load 1-way waveform of the Verasonics (stored in struct T)
if strcmp(S.probe, 'C1-6')
    % C1-6: (PulseCode 5 -21 32 42 1; 21 -42 32 21 1) (last index of row is how many times repeated)
    path_transmit = 'single_pulses\C1-6-D\';
    transmit_data = load([path_transmit 'C1_6_D_infoTW.mat'], 'TW');
    transmit_wvfm_temp = transmit_data.TW(1).Wvfm1Wy;
    transmit_wvfm = transmit_wvfm_temp(transmit_data.TW(1).States(1,2)+1:end)'/max(transmit_wvfm_temp);
    T.f_sampling = 250e6;
    T.transmit_time = 0: 1/T.f_sampling: (length(transmit_wvfm)-1)/T.f_sampling;
    
elseif strcmp(S.probe, 'P4-1')
    % P4-1: (PulseCode 15 60 29 -60 1) (last index is how many times repeated)
    path_transmit = 'single_pulses\P4-1\';
    transmit_data = load([path_transmit 'P4_1_infoTW.mat'], 'TW');
    transmit_wvfm_temp = transmit_data.TW.Wvfm1Wy;
    transmit_wvfm = transmit_wvfm_temp(transmit_data.TW.States(1,2)+1:end)'/max(transmit_wvfm_temp);
    T.f_sampling = 250e6;
    T.transmit_time = 0: 1/T.f_sampling: (length(transmit_wvfm)-1)/T.f_sampling;                                        % sampling frequency is 250 MHz of the pulse Verasonics makes

elseif strcmp(S.probe, 'L12-3v')
    % L12-3: (PulseCode 1 -5 7 9 1; 5 -9 7 5 1) (last index of row is how many times repeated)
    path_transmit = 'single_pulses\L12-3v\';
    transmit_data = load([path_transmit 'L12_3v_infoTW.mat'], 'TW');
    transmit_wvfm_temp = transmit_data.TW.Wvfm1Wy;
    transmit_wvfm = transmit_wvfm_temp(transmit_data.TW(1).States(1,2)+1:end)'/max(transmit_wvfm_temp);
    T.f_sampling = 250e6;
    T.transmit_time = 0: 1/T.f_sampling: (length(transmit_wvfm)-1)/T.f_sampling;

elseif strcmp(S.probe, 'L22-8v')
    % L12-3: (PulseCode 1 -5 7 9 1; 5 -9 7 5 1) (last index of row is how many times repeated)
    path_transmit = 'single_pulses\L22-8v\';
    transmit_data = load([path_transmit 'L22_8v_infoTW.mat'], 'TW');
    transmit_wvfm_temp = transmit_data.TW.Wvfm1Wy;
    transmit_wvfm = transmit_wvfm_temp(transmit_data.TW(1).States(1,2)+1:end)'/max(transmit_wvfm_temp);
    T.f_sampling = 250e6;
    T.transmit_time = 0: 1/T.f_sampling: (length(transmit_wvfm)-1)/T.f_sampling;
    
end


%% ultrasound parameters (stored in struct P)
% define start and end frequencies (in Hz) over which to evaluate the conding function
% and the driving frequency (Hz)
if strcmp(S.probe, 'C1-6')            % Trans.Bandwidth = 1.785 to 5.015 MHz
    P.f_ROI(1) = 1e6;
    P.f_ROI(2) = 6e6;

	P.f_driving = 2.2e6;

elseif strcmp(S.probe, 'P4-1')        % Trans.Bandwidth = 1.5 to 3.5 MHz
    P.f_ROI(1) = 1e6;
    P.f_ROI(2) = 4e6;

    P.f_driving = 1.7e6;

elseif strcmp(S.probe, 'L12-3v')      % Trans.Bandwidth = 4.03 to 11.05 MHz
    P.f_ROI(1) = 3e6;
    P.f_ROI(2) = 12e6;

    P.f_driving = 8.9290e6;             

elseif strcmp(S.probe, 'L22-8v')      % Trans.Bandwidth = 4.03 to 11.05 MHz
    P.f_ROI(1) = 8e6;
    P.f_ROI(2) = 22e6;

    P.f_driving = 15.6250e6;            

end

P.f_sampling = 250e6/4;                 % sampling frequency (Hz)
P.time_step = 1/P.f_sampling; 	        % in s

P.N_pulses = 3;%2:8;                       % number of pulses per train
P.dt0 = T.transmit_time(end);           % duration of 1 pulse (s)
P.ds0 = ceil(P.dt0*P.f_sampling);       % duration of 1 pulse (samples)
P.N_cy = P.dt0*P.f_driving;             % number of cycles per pulse in the train

P.c = 1490;


P.end_depth = 12e-2;                    % imaging depth in m

% vectors for the transmit pulse (to plot for reference)
P.time = 0:1/P.f_sampling: 2*P.end_depth/P.c;

% frequency vectors
P.f_full = linspace(0, P.f_sampling, length(P.time));
P.df = mean(diff(P.f_full));
P.f = P.f_ROI(1): P.df: P.f_ROI(2);      % frequency space (Hz), based on the transmit-receive transfer function of the P4-1; at -20 dB level
P.w = 2.*pi.*P.f;                        % angular frequency space (rad/s)
P.w_full = 2.*pi.*P.f_full;


%% Optimize gain in SNR for different sacrificial_depths and number of pulses in one train
SNR_gain = NaN(length(P.sacrificial_depth), length(P.N_pulses), Nshow);
ds_all = cell(length(P.sacrificial_depth), length(P.N_pulses));

for sac_depth_idx = 1: length(P.sacrificial_depth)

    S.sacrificial_depth = P.sacrificial_depth(sac_depth_idx);
    P.sacrificial_max_time = 2*S.sacrificial_depth/P.c;
    P.sacrificial_depth_max_samples = ceil(P.sacrificial_max_time/P.time_step);

    for pulse_no_idx = 1: length(P.N_pulses)

        % optimization parameters
        S.N_pulses = P.N_pulses(pulse_no_idx);           % number of pulses per train

        S.N_var = 3*P.ds0;                               % variation to either side from starting point?
        S.N_iter(pulse_no_idx) = (5*S.N_pulses)^2;       % number of iterations per attempt

        % initialize
        SNR_gain_err_best = cell(1,S.attempts);     % stores the best SNR gain per attempt
        ds_best = cell(1,S.attempts);               % stores the corresponding delays in samples
        iter_best = cell(1,S.attempts);             % stores the index of the best iteration (to check if we have enough iterations)

        %% Optimize for all random attempts/starting points

        for start_idx = 1: S.attempts

            %% Check all pulses would fit in the sacrificial_depth (without delays)
            if S.N_pulses*P.ds0 > P.sacrificial_depth_max_samples

                SNR_gain_err_best{start_idx} = NaN;
                ds_best{start_idx} = NaN(1, S.N_pulses);
                iter_best{start_idx} = NaN;

                continue
            end


            %% Initialize the optimization
            % generate S.N_pulses times a random delay in samples (ds) to initialize the
            % optimization, between P.ds0 and P.sacrificial_depth_max_samples-P.ds0

            ds_init = randi([P.ds0, P.sacrificial_depth_max_samples-P.ds0], 1, S.N_pulses-1);
            ds_init = [1, sort(ds_init)];

            % inter-pulse samples (ips)
            ips_init = diff(ds_init);

            % keep generating new ds until the separation ips is at least the
            % duration of each pulse, and the total train fits within
            % the sacrificial depth
            while min(ips_init) < P.ds0 || max(ds_init) > (P.sacrificial_depth_max_samples - P.ds0)  
                ds_init = randi([0, P.sacrificial_depth_max_samples-P.ds0], ...
                    1, S.N_pulses-1);
                ds_init = [1, sort(ds_init)];

                ips_init = diff(ds_init);
            end


            coding_function = sum(exp(-1i .* (ds_init' .* P.time_step) .* P.w), 1);

            SNR_gain_err = - (P.f_ROI(2) - P.f_ROI(1)) ./...
                ( sum( 1./ (abs(coding_function)).^2 , 'omitnan') .* mean(diff(P.f)) );


            %% Optimization
            % initialize storing the best options
            iter_best_idx = 1;
            SNR_gain_err_start{start_idx} = SNR_gain_err;
            SNR_gain_err_lowest = SNR_gain_err;

            ds_temp = ds_init;

            for iter_idx = 1: S.N_iter(pulse_no_idx)
                % number of total attempts (limited to 1e6)
                N_attempt = min(1e6, 1e3*iter_idx);

                % number of samples with respect to the starting point for all
                % attempts
                ips_variation = randi([-S.N_var,S.N_var], S.N_pulses-1, N_attempt);

                % construct full matrix with delays by adding the delay for the first pulse (=1)
                % stored as: (pulse idx, attempt index)
                N_ini_mat = repmat(ds_temp',[1,N_attempt]);
                N_delays = N_ini_mat;
                N_delays(2:end,:) = N_delays(2:end,:)+ips_variation;
                N_delays = sort(N_delays,1);

                %% Check the inter-pulse samples are feasible
                ips_init = diff(N_delays);

                % 1. remove options for which the inter-pulse samples is less than
                % the duration of a single pulse in samples
                ips_init_min = min(ips_init,[],1);
                N_delays(:, ips_init_min < P.ds0) = [];
                ips_variation(:, ips_init_min < P.ds0) = [];

                % 2. remove options for which the inter-pulse samples is more than
                % the maximum allowed delay
                ips_init_max = max(N_delays);
                N_delays(:, ips_init_max > (P.sacrificial_depth_max_samples-P.ds0)) = [];
                ips_variation(:, ips_init_max > (P.sacrificial_depth_max_samples-P.ds0)) = [];

                % 3. remove options for which the inter-pulse samples are negative
                N_init_min = min(N_delays);
                N_delays(:, N_init_min < 0) = [];
                ips_variation(:, N_init_min < 0) = [];

                % number of attempts after the feasibility-check
                N_attempt_loc = size(ips_variation, 2);


                %% For the remaining attempts, determine the - SNR gain based on the conding function
                decoding_factor_temp = sum(exp(-1i .* P.w' .* permute(N_delays * P.time_step, [3 2 1])), 3);

                SNR_gain_err = - (P.f_ROI(2) - P.f_ROI(1)) ./...
                    ( sum( 1./ (abs(decoding_factor_temp)).^2 , 'omitnan') .* mean(diff(P.f)) );

                % find attempt with lowest error
                [SNR_gain_err_min, SNR_gain_err_min_idx] = min(SNR_gain_err);

                % save the current SNR gain if it is lower than the lowest one saved so far
                if SNR_gain_err_min < SNR_gain_err_lowest
                    SNR_gain_err_lowest = SNR_gain_err_min;

                    % make lowest error new starting point
                    ds_temp = N_delays(:, SNR_gain_err_min_idx)';

                    % save the index of the best iteration
                    iter_best_idx = iter_idx;
                end
            end

            SNR_gain_err_best{start_idx} = SNR_gain_err_lowest;
            ds_best{start_idx} = ds_temp;
            iter_best{start_idx} = iter_best_idx;

            clc
            fprintf('Number of pulses: %d\n', S.N_pulses)
            fprintf('Attempt: %d out of %d \n', start_idx, S.attempts)
            fprintf('SNR gain best: %3.2f dB\n', 10*log10(-SNR_gain_err_lowest))


            % store options resulting in SNR loss as NaN
            if SNR_gain_err_best{start_idx} > 0    
                SNR_gain_err_best{start_idx} = NaN;
                ds_best{start_idx} = NaN(1, S.N_pulses);
                iter_best{start_idx} = NaN;
            end

        end



        %% Show details of the optimization for S.N_pulses
        error_best = cell2mat(SNR_gain_err_best);
        SNR_gain_temp = 10*log10(-error_best);

        % Save and show the conding functions of the Nshow best attempts
        Nshow = min(Nshow, S.attempts);
        [SNR_gain_sel, sel_idx] = maxk(SNR_gain_temp, Nshow);

        ds_best_sel_temp = cell2mat(ds_best');
        ds_best_sel = ds_best_sel_temp(sel_idx, :);

        SNR_gain(sac_depth_idx, pulse_no_idx, :) = SNR_gain_sel;
        ds_all{sac_depth_idx, pulse_no_idx} = ds_best_sel;


        if S.plot_optimization_result

            fig_handle = figure(S.N_pulses*10 + sac_depth_idx);
            clf

            sgtitle(sprintf('Optimization details N = %d', S.N_pulses))

            subplot(2,2,1)
            hold on
            title(sprintf('%d best conding functions', Nshow))
            for n = 1: Nshow
                ds = ds_best_sel(n, :);
                coding_function = sum(exp(-1i .* (ds' .* P.time_step) .* P.w), 1);

                plot(P.f.*1e-6, abs(coding_function));
            end
            grid on
            xlabel('frequency (MHz)')
            ylabel('|conding function|')
            xlim([P.f_ROI]*1e-6)
            set(gca, 'Box', 'on')
            ylim([0 S.N_pulses])

            % SNR gain for all attempt
            subplot(2,2,2)
            hold on
            title('SNR gain')
            plot(squeeze(SNR_gain_temp), '.', 'markersize', 10)
            hold all
            xlabel('attempt')
            ylabel('best SNR gain (dB)');
            grid on
            set(gca, 'Box', 'on')
            ylim([-8 ceil(10*log10(max(P.N_pulses)))])
            yline(10*log10(S.N_pulses), 'r--')
            xlim([1 S.attempts])


            % which iteration returned best SNR gain (if equal to S.N_iter, it
            % might need to be increased)
            subplot(2,2,3)
            hold on
            title('best iteration')
            plot([iter_best{:}], '.', 'markersize', 10)
            xlabel('attempt')
            ylabel('index')
            grid on
            yline(S.N_iter(pulse_no_idx) , 'r--')
            ylim([0 S.N_iter(pulse_no_idx)+1])
            set(gca, 'Box', 'on')
            xlim([1 S.attempts])

            subplot(2,2,4)
            hold on
            title('lowest error')
            plot(error_best, '.', 'markersize', 10)
            hold all
            xlabel('attempt')
            ylabel('error minimized')
            grid on
            set(gca, 'Box', 'on')
            xlim([1 S.attempts])

        end
    end


    %% Plot the results of all S.N_pulses
    % boxplot showing the maximum, 75th percentile, median, 25th percentile,
    % and minimum for each number of pulses per train

    colors = colororder('gem');

    figure(1000+sac_depth_idx-1)
    clf
    hold on

    Npulse_theo = 1: max(P.N_pulses)+1;
    SNR_gain_separate = 10.*log10(Npulse_theo);
    plt_theo_max = plot(Npulse_theo, SNR_gain_separate, '.-', 'color','k');

    plt = cell(length(P.sacrificial_depth), 1);
    % plt_max = NaN(1, length(P.sacrificial_depth));
    for pulse_no_idx = 1:length(P.N_pulses)
        plt{sac_depth_idx} = boxplot(squeeze(SNR_gain(sac_depth_idx, pulse_no_idx, :)),...
            P.N_pulses(pulse_no_idx), 'Positions', P.N_pulses(pulse_no_idx), ...
            'Colors', colors(sac_depth_idx,:), 'BoxStyle', 'filled', 'Symbol', '.');

        plt_max = plot( P.N_pulses(pulse_no_idx), ...
            max(SNR_gain(sac_depth_idx, pulse_no_idx, :)), ...
            'x', 'markersize', 10, 'color', colors(sac_depth_idx,:), 'linewidth', 1.5);
    end

    xlabel('number of pulses')
    ylabel('SNR gain (dB)')
    grid on

    plt_leg = legend([plt_theo_max(1), plt_max], 'max theoretical gain',...
        sprintf('max for sacrificial_depth %2.2f cm', S.sacrificial_depth*1e2));
    plt_leg.Location = 'southeast';

    xticks(Npulse_theo)
    xticklabels(Npulse_theo)
    xlim([1 max(P.N_pulses)+1])
    ylim([0 ceil(max(SNR_gain_separate))])


end







