classdef EI_pool < handle
    % EI_pool A pool of exictatory and inhibitory neurons (potentially on the GPU)

    properties (Access = public)
        % Parameters
        on_GPU logical

        N uint32 {mustBeInteger, mustBePositive}

        T double {mustBeFloat, mustBePositive} % sec
        dt double {mustBeFloat, mustBePositive} = 0.0001 % sec

        N_input uint32 {mustBeInteger, mustBePositive} = 50
        check_input_period double {mustBeFloat, mustBePositive} = 0.010 % sec
        input_stopped logical = false
        input_Hz_in_dt double {mustBeFloat} % Hz in dts

        conn_prob_input2E double {mustBeFloat, mustBePositive} % %
        conn_strength_input2E double {mustBeFloat, mustBeNonnegative}
        conn_prob_E2E double {mustBeFloat, mustBePositive} % %
        conn_strength_E2E double {mustBeFloat, mustBeNonnegative}

        conn_max double {mustBeFloat, mustBePositive} = 1.0 % a.u.

        synapses_glutamate_tau_rise double {mustBeFloat, mustBePositive} = 0.0026 % sec (Human - Hunt et al., Cerebral Cortex, 2023)
        synapses_glutamate_tau_fall double {mustBeFloat, mustBePositive} = 0.0313 % sec (Human - Hunt et al., Cerebral Cortex, 2023)

        Vm_tau double {mustBeFloat, mustBePositive} = 0.025 % sec (Rall, Biophysical Journal, 1969)
        Vm_base double {mustBeFloat} = 0.0 % a.u.
        Vm_thresh double {mustBeFloat} = 1.0 % a.u.
        Vm_min double {mustBeFloat} = -1.0 % a.u.

        axonal_delay_min double {mustBeFloat, mustBePositive} = 0.001 % sec (Human - Lemarechal et al., Brain, 2022)
        axonal_delay_max double {mustBeFloat, mustBePositive} = 0.005 % sec (Human - Lemarechal et al., Brain, 2022)

        refractoriness_absolute double {mustBeFloat, mustBePositive} = 0.003 % sec (Hodgkin-Huzley model - Follman et al., Physical Rev., 2015)
        refractoriness_absolute_in_dt uint32 {mustBeInteger, mustBePositive} % dts
        refractoriness_relative_tau double {mustBeFloat, mustBePositive} = 0.005 % sec (Hodgkin-Huzley model - Follman et al., Physical Rev., 2015)

        % Variables
        %   (Have to have these separate so they can be sparse)
        conn_matrix_input2E (:, :) % can't specify otherwise will loose sparseness
        conn_weights_input2E (:, :) % ''
        conn_delays_input2E (:, :) % ''

        conn_matrix_E2E (:, :) % ''
        conn_weights_E2E (:, :) % ''
        conn_delays_E2E (:, :) % ''

        %   (Have to have these separate so they can be sparse)
        synapses_rise_input2E (:, :) % can't specify otherwise will loose sparseness
        synapses_fall_input2E (:, :) % ''

        synapses_rise_E2E (:, :) % ''
        synapses_fall_E2E (:, :) % ''

        Vm (1, :) double {mustBeFloat} = []

        input_spikes (1, :) logical = []
        spikes (1, :) logical = []

        %   (Have to have these separate because others are separate due to sparseness)
        spikes_delayed_input2E (:, :)
        spikes_delayed_E2E (:, :)
        spikes_times (1, :) uint32 {mustBeInteger, mustBeNonnegative} = []

        refractoriness (1, :) double {mustBeFloat, mustBeNonnegative} = []

        % Plotting
        fontsize double = 14;

        save_to_plot_synapses logical
        save_to_plot_Vm logical
        save_to_plot_spikes logical
        save_to_plot_refractoriness logical
        save_to_plot_weights logical
        save_to_plot_Hz logical
        save_to_plot_stimulation logical
        b_plot_digraph logical

        save_weights_t_period double {mustBeFloat, mustBePositive} = 0.100 % sec
        conn_weights_input2E_t (:, :, :) double {mustBeFloat} = [] % pre, post, t
        conn_weights_E2E_t (:, :, :) double {mustBeFloat} = [] % ''

        synapses_input2E_t (:, :) double {mustBeFloat} = []
        synapses_E2E_t (:, :) double {mustBeFloat} = []

        Vm_t (:, :) double {mustBeFloat} = [] % t, post

        input_spikes_t (:, :) logical = [] % t, post
        spikes_t (:, :) logical = [] % t, post

        refractoriness_t (:, :) double {mustBeFloat} = [] % t, post

        save_Hz_t_period double {mustBeFloat, mustBePositive} = 0.100 % sec
        input_spike_count_individual uint32 {mustBeInteger, mustBeNonnegative} = []
        input_spike_count uint32 {mustBeInteger, mustBeNonnegative} = 0
        input_Hz_t (:, :) double {mustBeFloat, mustBeNonnegative} = [] % t, post
        spike_count_individual uint32 {mustBeInteger, mustBeNonnegative} = []
        spike_count uint32 {mustBeInteger, mustBeNonnegative} = 0
        Hz_t (:, :) double {mustBeFloat, mustBeNonnegative} = [] % t, post

        save_stimulation_t_period double {mustBeFloat, mustBePositive} = 0.100 % sec
        stimulation_t (:, 1) logical = []
    end
    
    methods (Static, Access = private)
        function [weights, matrix, delays] = set_up_weights( ...
                dim1, ...
                dim2, ...
                conn_prob,  ...
                conn_max, ...
                recurrent, ...
                axonal_delay_max, ...
                axonal_delay_min, ...
                dt, ...
                on_GPU ...
                )
            arguments
                dim1 uint32 {mustBeInteger, mustBePositive}
                dim2 uint32 {mustBeInteger, mustBePositive}
                conn_prob double {mustBeFloat, mustBePositive}
                conn_max double {mustBeFloat}
                recurrent logical
                axonal_delay_max double {mustBeFloat, mustBePositive}
                axonal_delay_min double {mustBeFloat, mustBePositive}
                dt double {mustBeFloat, mustBePositive}
                on_GPU logical = false
            end
            % set_up_weights Sets up the weights, matrix, and delays
            if on_GPU
                weights = gpuArray(sparse(dim1, dim2));
                matrix = gpuArray(sparse(dim1, dim2));
            else
                weights = sparse(dim1, dim2);
                matrix = sparse(dim1, dim2);
            end
            % weights
            weights = ...
                rand(dim1, dim2) .* sparse(rand(dim1, dim2) <= conn_prob); % Uniform
            weights = min(weights, conn_max);
            weights = max(weights, 0.0);
            if recurrent
                weights = weights - spdiags( ...
                    diag(weights), ...
                    0, ...
                    dim1, ...
                    dim2 ...
                    );
            end
            % matrix
            matrix = weights > 0;
            % Stabilize the distribution
            for i = 1:1000
                weights = weights + ...
                    ( ...
                    ((weights * 0.2) + 0.01) .* ...
                    (sprandn(matrix) * 0.5) ... % don't go too fast
                    );
                weights = min(weights, conn_max);
                weights = max(weights, 0.0);
            end
            matrix = weights > 0;
            % delays
            if on_GPU
                delays = gpuArray(sparse(dim1, dim2));
            else
                delays = sparse(dim1, dim2);
            end
            delays = ( ...
                randi((axonal_delay_max - axonal_delay_min) / ...
                dt, dim1, dim2 ...
                ) + (axonal_delay_min / dt)) .* matrix;
        end
    end

    methods (Access = protected)
        function plot_synapses_refractoriness_Vm(pool, export)
            arguments
                pool EI_pool
                export logical
            end
            % plot_synapses_refractoriness_Vm Plots the network's synapses, neurons'
            % refractoriness, and neurons' Vm

            tmp_synapses_input2E_t = gather(pool.synapses_input2E_t);
            tmp_refractoriness_t = gather(pool.refractoriness_t);
            tmp_Vm_t = gather(pool.Vm_t);
            if ~isempty(tmp_synapses_input2E_t) || ~isempty(tmp_refractoriness_t) || ~isempty(tmp_Vm_t)
                tmp_synapses_E2E_t = gather(pool.synapses_E2E_t);

                fig = figure();

                if ~isempty(tmp_synapses_E2E_t)
                    nexttile
                    plot(tmp_synapses_E2E_t');
                    title("Synapses")
                    xlabel("Time (s)")
                    xticks(0:0.1 / pool.dt:pool.T / pool.dt)
                    xticklabels(0:0.1:pool.T)
                    ylabel("(a.u.)")
                end
                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/synapses.pdf','ContentType','vector')
                end

                fig = figure();
                tiledlayout(2, 1)

                if ~isempty(tmp_refractoriness_t)
                    nexttile
                    plot(tmp_refractoriness_t(:, 10));
                    title("Refractoriness")
                    xticklabels([])
                    ylabel("(a.u.)")
                end

                if ~isempty(tmp_Vm_t)
                    nexttile
                    plot(tmp_Vm_t(:, 10));
                    title("Vm")
                    xlabel("Time (s)")
                    xticks(0:0.1 / pool.dt:pool.T / pool.dt)
                    xticklabels(0:0.1:pool.T)
                    ylabel("(a.u.)")
                end
                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/Vm.pdf','ContentType','vector')
                end
            end
        end

        function plot_spikes(pool, export)
            arguments
                pool EI_pool
                export logical
            end
            % plot_spikes Plots the network's spikes

            tmp_input_spikes_t = gather(pool.input_spikes_t);
            tmp_spikes_t = gather(pool.spikes_t);
            if ~isempty(tmp_spikes_t)
                fig = figure();

                nexttile
                [t, n] = find(tmp_spikes_t);
                scatter( ...
                    t, ...
                    n, ...
                    ".", ...
                    'red' ...
                    );
                title("Spikes")
                xlabel("Time (s)")
                xticks(0:0.1 / pool.dt:pool.T / pool.dt)
                xticklabels(0:0.1:pool.T)
                ylabel("Neuron")
                legend( ...
                    {"Excitatory"}, "Location", ...
                    "southoutside", ...
                    "NumColumns", ...
                    3 ...
                    )

                nexttile("south")
                [t, n] = find(tmp_input_spikes_t);
                scatter(t, n, ".", 'magenta');
                title("Input spikes")
                xlabel("Time (s)")
                xticks(0:0.1 / pool.dt:pool.T / pool.dt)
                xticklabels(0:0.1:pool.T)
                ylabel("Neuron")
                legend({"Input", "Other"})

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/spikes.pdf','ContentType','vector')
                end
            end
        end

        function plot_weights(pool, export)
            arguments
                pool EI_pool
                export logical
            end
            % plot_weights Plots the network's weights

            tmp_conn_weights_input2E_t = pool.conn_weights_input2E_t;
            tmp_conn_weights_E2E_t = pool.conn_weights_E2E_t;
            if ~isempty(tmp_conn_weights_input2E_t)
                bar_N = 100;
                max_W = 10.0;

                tmp_weights = tmp_conn_weights_input2E_t(:, :, 1);
                tmp_weights = tmp_weights(pool.conn_matrix_input2E);
                tmp_input2E_mean_start = mean(tmp_weights(:));
                tmp_counts_input2E_start = histcounts( ...
                    tmp_weights, ...
                    0:max_W/bar_N:max_W ...
                    );
                tmp_weights = tmp_conn_weights_input2E_t(:, :, end);
                tmp_weights = tmp_weights(pool.conn_matrix_input2E);
                tmp_input2E_mean_end = mean(tmp_weights(:));
                tmp_counts_input2E_end = histcounts( ...
                    tmp_weights, ...
                    0:max_W/bar_N:max_W ...
                    );

                fig = figure();
                tiledlayout(2,1)

                nexttile
                tmp_h3 = bar( ...
                    0.0:max_W / (bar_N - 1):max_W, ...
                    [tmp_counts_input2E_start(:)], ...
                    "stacked" ...
                    );
                set(gca,'YScale','log')
                tmp_h3(1).FaceColor = "magenta";
                title("Start of simulation input➔E weights")
                xlim([-0.01 max_W + 0.01])
                xticks([])
                ylabel("Frequency")

                xline(tmp_input2E_mean_start, "magenta", "LineStyle", "--")

                legend({"input➔E"}, "Location", "southoutside")

                nexttile
                tmp_h4 = bar( ...
                    0.0:max_W / (bar_N - 1):max_W, ...
                    [tmp_counts_input2E_end(:)], ...
                    "stacked" ...
                    );
                set(gca,'YScale','log')
                tmp_h4(1).FaceColor = "magenta";
                title("End of simulation input➔E weights")
                xlim([-0.01 max_W + 0.01])
                ylabel("Frequency")

                xline(tmp_input2E_mean_end, "magenta", "LineStyle", "--")

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/histogram_weights_input2E.pdf','ContentType','vector')
                end

                fig = figure();
                tiledlayout(2,1)

                tmp_weights = sparse(tmp_conn_weights_E2E_t(:, :, 1));
                tmp_weights = tmp_weights(pool.conn_matrix_E2E);
                tmp_E2E_mean_start = mean(tmp_weights(:));

                tmp_weights = sparse(tmp_conn_weights_E2E_t(:, :, end));
                tmp_weights = tmp_weights(pool.conn_matrix_E2E);
                tmp_E2E_mean_end = mean(tmp_weights(:));

                nexttile

                tmp_weights = tmp_conn_weights_E2E_t(:, :, 1);
                tmp_weights = tmp_weights(pool.conn_matrix_E2E);
                tmp_counts_E2E_start = histcounts( ...
                    tmp_weights, ...
                    0:max_W/bar_N:max_W ...
                    );

                tmp_h5 = bar( ...
                    0.0:max_W / (bar_N - 1):max_W, ...
                    [tmp_counts_E2E_start(:)], ...
                    "stacked" ...
                    );
                set(gca,'YScale','log')
                tmp_h5(1).FaceColor = "red";
                title("Start of simulation E➔E weights")
                xlim([-0.01 max_W + 0.01])
                xticks([])
                ylabel("Frequency")

                xline(tmp_E2E_mean_start, "red", "LineStyle", "--")

                legend({ ...
                    "E➔E", ...
                    "E➔E mean", ...
                    }, ...
                    "Location", ...
                    "southoutside", ...
                    "NumColumns", ...
                    3 ...
                    )

                nexttile

                tmp_weights = tmp_conn_weights_E2E_t(:, :, end);
                tmp_weights = tmp_weights(pool.conn_matrix_E2E);
                tmp_counts_E2E_end = histcounts( ...
                    tmp_weights, ...
                    0:max_W/bar_N:max_W ...
                    );


                tmp_h6 = bar( ...
                    0.0:max_W / (bar_N - 1):max_W, ...
                    [tmp_counts_E2E_end(:)], ...
                    "stacked" ...
                    );
                set(gca,'YScale','log')
                tmp_h6(1).FaceColor = "red";
                title("End of simulation E➔E weights")
                xlim([-0.01 max_W + 0.01])
                xlabel("Weight (a.u.)")
                ylabel("Frequency")

                xline(tmp_E2E_mean_end, "red", "LineStyle", "--")
                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/histogram_weights_E2E.pdf','ContentType','vector')
                end

                fig = figure();
                tmp_weights = reshape( ...
                    pool.conn_weights_E2E_t, ...
                    pool.N * pool.N, ...
                    size(pool.conn_weights_E2E_t, 3) ...
                    );
                tmp_weights( ~any(tmp_weights,2), : ) = [];
                plot(tmp_weights(1:100, :)')
                box off
                xlabel("Time (s)")
                xlim([-0.5 pool.T / pool.save_weights_t_period])
                xticks(0:10.0 / pool.save_resource_pools_t_period:pool.T / pool.save_resource_pools_t_period)
                xticklabels(0:10:pool.T)
                % yticklabels([0])
                % yticks([0])
                % ylim([0 1])
                ylabel("$w$", "Interpreter", "latex")

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/weights_E2E_traces.pdf','ContentType','vector')
                end

                fig = figure();
                tmp_weights = tmp_conn_weights_E2E_t(:, :, end);
                tmp_weights = tmp_weights(pool.conn_matrix_E2E);
                histogram(tmp_weights, 30, "FaceColor", "black")
                box off
                xlabel("$w$", "Interpreter", "latex")
                % xticklabels([0])
                % xticks([0])
                % xlim([0 1])
                ylabel("Frequency")
                set(gca, 'YScale', 'log')

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/weights_E2E_histogram.pdf','ContentType','vector')
                end

                fig = figure();
                tmp_mean_diffs = zeros(pool.N, size(tmp_conn_weights_E2E_t, 3));
                for t = 1:size(tmp_conn_weights_E2E_t, 3)
                    for post = 1:pool.N
                        tmp_significant_weights = tmp_conn_weights_E2E_t(:, post, t) > 0;
                        tmp_mean_diffs(post, t) = mean( ...
                            diff(find(tmp_significant_weights(pool.conn_matrix_E2E(:, post)))) ...
                            );
                    end
                end

                plot(mean(tmp_mean_diffs, "omitnan"), "black");
                xlabel("Time (s)")
                xlim([-0.5 pool.T / pool.save_weights_t_period])
                xticks(0:10.0 / pool.save_weights_t_period:pool.T / pool.save_weights_t_period)
                xticklabels(0:10:pool.T)
                ylabel("Distance")
                box off

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/weights_E2E_distance.pdf','ContentType','vector')
                end
            end
        end

        function plot_Hz(pool, export)
            arguments
                pool EI_pool
                export logical
            end
            % plot_Hz Plots the network's firing rate

            tmp_input_Hz_t = gather(pool.input_Hz_t);
            tmp_Hz_t = gather(pool.Hz_t);
            if ~isempty(tmp_input_Hz_t)
                fig = figure();
                % tiledlayout(3, 1)

                % nexttile([2, 1])
                imagesc(tmp_Hz_t')
                axis xy
                xlabel("Time (s)")
                xlim([-0.5 pool.T / pool.save_Hz_t_period])
                xticks(0:100.0 / pool.save_Hz_t_period:pool.T / pool.save_Hz_t_period)
                xticklabels(0:100:pool.T)
                ylim([0 pool.N])
                ylabel("Neuron")
                yticklabels([])
                colormap hot
                cb = colorbar;
                ylabel(cb, "Firing rate (Hz)")
                box off

                tmp_stimulation_t = gather(pool.stimulation_t);
                if ~isempty(tmp_stimulation_t)
                    tmp_stimulation_range = 1:length(tmp_stimulation_t);
                    tmp_stimulation_x = tmp_stimulation_range(tmp_stimulation_t);

                    % hold on
                    % scatter( ...
                    %     tmp_stimulation_x, ...
                    %     zeros(length(tmp_stimulation_x)), ...
                    %     49, ...
                    %     'm', ...
                    %     "square", ...
                    %     "filled", ...
                    %     "HandleVisibility", ...
                    %     "off" ...
                    %     )
                    % hold off
                end

                % nexttile
                % imagesc(tmp_input_Hz_t')
                % title("Input Hz")
                % axis xy
                % xlabel("Time (s)")
                % xlim([-0.5 pool.T / pool.save_Hz_t_period])
                % xticks(0:1.0 / pool.save_Hz_t_period:pool.T / pool.save_Hz_t_period)
                % xticklabels(0:pool.T)
                % ylim([0 pool.N_input])
                % ylabel("Neuron")
                % colormap hot
                % cb = colorbar;
                % ylabel(cb, "Firing rate (Hz)")

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/Hz.pdf','ContentType','vector')
                end

                tmp_Hz_t_E = tmp_Hz_t(:, 1:pool.N);
                tmp_selected = tmp_Hz_t_E > 50;
                [~, tmp_i] = sortrows(tmp_selected', 1:size(tmp_Hz_t_E, 1), "descend");

                fig = figure();
                imagesc(tmp_Hz_t_E(:, tmp_i)')
                title("Sorted excitatory Hz")
                axis xy
                xlabel("Time (s)")
                xlim([-0.5 pool.T / pool.save_Hz_t_period])
                xticks(0:10.0 / pool.save_resource_pools_t_period:pool.T / pool.save_resource_pools_t_period)
                xticklabels(0:10:pool.T)
                ylabel("Excitatory neuron (sorted > 25 Hz)")
                yticks([]);
                colormap hot
                cb = colorbar();
                ylabel(cb, "Firing rate (Hz)")

                if ~isempty(tmp_stimulation_t)
                    hold on
                    scatter( ...
                        tmp_stimulation_x, ...
                        ones(length(tmp_stimulation_x)), ...
                        49, ...
                        'm', ...
                        "square", ...
                        "filled", ...
                        "HandleVisibility", ...
                        "off" ...
                        )
                    hold off
                end
                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/Hz_sorted.pdf','ContentType','vector')
                end

                fig = figure();
                hold on
                shadedErrorBar( ...
                    1:numel(tmp_input_Hz_t(:, 1)), ...
                    tmp_input_Hz_t', ...
                    {@mean, @std}, ...
                    "lineprops", ...
                    '-m' ...
                    )
                shadedErrorBar( ...
                    1:numel(tmp_Hz_t(:, 1)), ...
                    tmp_Hz_t', ...
                    {@mean, @std}, ...
                    "lineprops", ...
                    '-r' ...
                    )
                hold off
                title("Hz")
                xlabel("Time (s)")
                xticks(0:1.0 / pool.save_Hz_t_period:pool.T / pool.save_Hz_t_period)
                xticklabels(0:pool.T)
                ylim([-7 250])
                ylabel("Mean firing rate (Hz)")
                legend({"Input", "Excitatory"}, "location", "northwest")

                if ~isempty(tmp_stimulation_t)
                    hold on
                    scatter( ...
                        tmp_stimulation_x, ...
                        ones(length(tmp_stimulation_x)) * -3.5, ...
                        49, ...
                        'm', ...
                        "square", ...
                        "filled", ...
                        "HandleVisibility", ...
                        "off" ...
                        )
                    hold off
                end
                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/Hz_means.pdf','ContentType','vector')
                end
            end
        end

        function plot_digraph(pool, export)
            arguments
                pool EI_pool
                export logical
            end
            % plot_digraph Plots the network's digraph

            if pool.b_plot_digraph
                tmp_conn_weights_E2E_t = gather(pool.conn_weights_E2E_t);
                tmp_threshold_weights = tmp_conn_weights_E2E_t(:, :, end);
                tmp_threshold_weights = mean(nonzeros(tmp_threshold_weights(:)));

                fig = figure();
                % Graph
                tmp_G_E2E = digraph();
                for pre = 1:pool.N
                    for post = 1:pool.N
                        if tmp_conn_weights_E2E_t(pre, post, end) > tmp_threshold_weights
                            tmp_G_E2E = addedge(tmp_G_E2E, pre, post, tmp_conn_weights_E2E_t(pre, post, end));
                        end
                    end
                end
                % Plot
                isolated_nodes = find(indegree(tmp_G_E2E) + outdegree(tmp_G_E2E) == 0);
                tmp_G_E2E = rmnode(tmp_G_E2E, isolated_nodes);
                h = plot(tmp_G_E2E);%, "Layout", "circle");
                box off
                axis off
                % title("> mean of E➔E weights at end of simulation")
                tmp_G_E2E.Edges.LWidths = 0.00001 + tmp_G_E2E.Edges.Weight/max(tmp_G_E2E.Edges.Weight) * 5;
                h.LineWidth = tmp_G_E2E.Edges.LWidths;
                h.EdgeAlpha = 0.2;

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/digraph.pdf','ContentType','vector')
                end

                tmp_threshold_weights = 0;

                tmp_in_degree = zeros(length(tmp_conn_weights_E2E_t(1, 1, :)), 1);
                tmp_out_degree = zeros(length(tmp_conn_weights_E2E_t(1, 1, :)), 1);
                tmp_in_closeness = zeros(length(tmp_conn_weights_E2E_t(1, 1, :)), 1);
                tmp_out_closeness = zeros(length(tmp_conn_weights_E2E_t(1, 1, :)), 1);
                tmp_betweenness = zeros(length(tmp_conn_weights_E2E_t(1, 1, :)), 1);
                tmp_pagerank = zeros(length(tmp_conn_weights_E2E_t(1, 1, :)), 1);
                tmp_hubs = zeros(length(tmp_conn_weights_E2E_t(1, 1, :)), 1);
                tmp_authorities = zeros(length(tmp_conn_weights_E2E_t(1, 1, :)), 1);
                for t = 1:length(tmp_conn_weights_E2E_t(1, 1, :))
                    % Graph
                    tmp_G_E2E = digraph();
                    for pre = 1:pool.N
                        for post = 1:pool.N
                            if tmp_conn_weights_E2E_t(pre, post, t) > tmp_threshold_weights
                                tmp_G_E2E = addedge(tmp_G_E2E, pre, post, tmp_conn_weights_E2E_t(pre, post, t));
                            end
                        end
                    end
                    isolated_nodes = find(indegree(tmp_G_E2E) + outdegree(tmp_G_E2E) == 0);
                    tmp_G_E2E = rmnode(tmp_G_E2E, isolated_nodes);
                    tmp_functional_weights = tmp_conn_weights_E2E_t(:, :, t);
                    tmp_functional_weights = tmp_functional_weights( ...
                        tmp_functional_weights > tmp_threshold_weights ...
                        );
                    tmp_functional_delays = pool.conn_delays_E2E(:, :);
                    tmp_functional_delays = tmp_functional_delays( ...
                        tmp_conn_weights_E2E_t(:, :, t) > tmp_threshold_weights ...
                        );
                    % Degree
                    tmp_in_degree(t) = mean( ...
                        centrality( ...
                        tmp_G_E2E, ...
                        "indegree", ...
                        "Importance", ...
                        tmp_functional_weights ...
                        ));
                    tmp_out_degree(t) = mean( ...
                        centrality( ...
                        tmp_G_E2E, ...
                        "outdegree", ...
                        "Importance", ...
                        tmp_functional_weights ...
                        ));
                    % Closeness
                    tmp_in_closeness(t) = mean( ...
                        centrality( ...
                        tmp_G_E2E, ...
                        "incloseness", ...
                        "Cost", ...
                        tmp_functional_delays ...
                        ));
                    tmp_out_closeness(t) = mean( ...
                        centrality( ...
                        tmp_G_E2E, ...
                        "outcloseness", ...
                        "Cost", ...
                        tmp_functional_delays ...
                        ));
                    % Betweenness
                    tmp_betweenness(t) = mean( ...
                        centrality( ...
                        tmp_G_E2E, ...
                        "betweenness", ...
                        "Cost", ...
                        tmp_functional_delays ...
                        ));
                    % Pagerank
                    tmp_pagerank(t) = mean( ...
                        centrality( ...
                        tmp_G_E2E, ...
                        "pagerank", ...
                        "Importance", ...
                        tmp_functional_weights ...
                        ));
                    % Hubs
                    tmp_hubs(t) = mean( ...
                        centrality( ...
                        tmp_G_E2E, ...
                        "hubs", ...
                        "Importance", ...
                        tmp_functional_weights, ...
                        "MaxIterations", ...
                        5000 ...
                        ));
                    % Authorities
                    tmp_authorities(t) = mean( ...
                        centrality( ...
                        tmp_G_E2E, ...
                        "authorities", ...
                        "Importance", ...
                        tmp_functional_weights, ...
                        "MaxIterations", ...
                        5000 ...
                        ));
                end

                fig = figure();
                tiledlayout(2, 2)

                % title(tiledlayout, "Graph measures (> mean of E➔E weights at end of simulation)")

                nexttile
                plot(tmp_in_degree)
                box off
                xlabel("Time")
                xticks([])
                ylabel("In degree")

                nexttile
                plot(tmp_out_degree)
                box off
                xlabel("Time")
                xticks([])
                ylabel("Out degree")

                nexttile
                plot(tmp_in_closeness)
                box off
                xlabel("Time")
                xticks([])
                ylabel("In closeness")

                nexttile
                plot(tmp_out_closeness)
                box off
                xlabel("Time")
                xticks([])
                ylabel("Out closeness")

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/digraph_measures_1.pdf','ContentType','vector')
                end

                fig = figure();
                tiledlayout(2, 2)

                title(tiledlayout, "Graph measures (> mean of E➔E weights at end of simulation)")

                nexttile
                plot(tmp_betweenness)
                box off
                xlabel("Time")
                xticks([])
                ylabel("Betweeness")

                nexttile
                plot(tmp_pagerank)
                box off
                xlabel("Time")
                xticks([])
                ylabel("Pagerank")

                nexttile
                plot(tmp_hubs)
                box off
                xlabel("Time")
                xticks([])
                ylabel("Hubs")

                nexttile
                plot(tmp_authorities)
                box off
                xlabel("Time")
                xticks([])
                ylabel("Authorities")

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/digraph_measures_2.pdf','ContentType','vector')
                end
            end
        end    
    end

    methods (Access = protected, Abstract = true)
        check_input(pool, t)
    end

    methods
        function pool = EI_pool( ...
                T, ...
                N, ...
                on_GPU, ...
                conn_prob_input2E, ...
                conn_strength_input2E, ...
                conn_prob_E2E, ...
                conn_strength_E2E, ...
                save_to_plot_synapses, ...
                save_to_plot_Vm, ...
                save_to_plot_spikes, ...
                save_to_plot_refractoriness, ...
                save_to_plot_weights, ...
                save_to_plot_Hz, ...
                save_to_plot_stimulation, ...
                b_plot_digraph ...
                )
            arguments
                T double {mustBeFloat, mustBePositive}
                N uint32 {mustBeInteger, mustBePositive}
                on_GPU logical = false
                conn_prob_input2E double {mustBeFloat, mustBeNonnegative} = 0.0 % %
                conn_strength_input2E double {mustBeFloat, mustBeNonnegative} = 0.0 % a.u.
                conn_prob_E2E double {mustBeFloat, mustBeNonnegative} = 0.0 % %
                conn_strength_E2E double {mustBeFloat, mustBeNonnegative} = 0.0 % a.u.
                save_to_plot_synapses logical = false
                save_to_plot_Vm logical = false
                save_to_plot_spikes logical = false
                save_to_plot_refractoriness logical = false
                save_to_plot_weights logical = false
                save_to_plot_Hz logical = false
                save_to_plot_stimulation logical = false
                b_plot_digraph logical = false
            end
            % EI_pool Construct an instance of an EI pool (on the GPU)
            pool = pool@handle();

            pool.on_GPU = on_GPU;

            pool.conn_prob_input2E = conn_prob_input2E;
            pool.conn_strength_input2E = conn_strength_input2E;
            pool.conn_prob_E2E = conn_prob_E2E;
            pool.conn_strength_E2E = conn_strength_E2E;

            % Plotting
            pool.save_to_plot_synapses = save_to_plot_synapses;
            pool.save_to_plot_Vm = save_to_plot_Vm;
            pool.save_to_plot_spikes = save_to_plot_spikes;
            pool.save_to_plot_refractoriness = save_to_plot_refractoriness;
            pool.save_to_plot_weights = save_to_plot_weights;
            pool.save_to_plot_Hz = save_to_plot_Hz;
            pool.save_to_plot_stimulation = save_to_plot_stimulation;
            pool.b_plot_digraph = b_plot_digraph;

            % General
            pool.T = T;

            pool.N = N;

            % Input
            % none

            % Weights
            %   input➔E
            [weights, matrix, delays] = pool.set_up_weights( ...
                pool.N_input, ...
                pool.N, ...
                pool.conn_prob_input2E, ...
                pool.conn_max, ...
                false, ...
                pool.axonal_delay_max, ...
                pool.axonal_delay_min, ...
                pool.dt, ...
                pool.on_GPU ...
                );
            pool.conn_weights_input2E = weights;
            pool.conn_matrix_input2E = matrix;
            pool.conn_delays_input2E = delays;
            %   E➔E
            [weights, matrix, delays] = pool.set_up_weights( ...
                pool.N, ...
                pool.N, ...
                pool.conn_prob_E2E, ...
                pool.conn_max, ...
                true, ...
                pool.axonal_delay_max, ...
                pool.axonal_delay_min, ...
                pool.dt, ...
                pool.on_GPU ...
                );
            pool.conn_weights_E2E = weights;
            pool.conn_matrix_E2E = matrix;
            pool.conn_delays_E2E = delays;

            % Neurons
            pool.refractoriness_absolute_in_dt = pool.refractoriness_absolute / pool.dt;

            % States            
            if pool.on_GPU
                %   Input
                pool.synapses_rise_input2E = gpuArray(sparse(pool.N_input, pool.N)); % double
                pool.synapses_fall_input2E = gpuArray(sparse(pool.N_input, pool.N)); % double
                %   E➔E
                pool.synapses_rise_E2E = gpuArray(sparse(pool.N, pool.N)); % double
                pool.synapses_fall_E2E = gpuArray(sparse(pool.N, pool.N)); % double
            else
                %   Input
                pool.synapses_rise_input2E = sparse(pool.N_input, pool.N); % double
                pool.synapses_fall_input2E = sparse(pool.N_input, pool.N); % double
                %   E➔E
                pool.synapses_rise_E2E = sparse(pool.N, pool.N); % double
                pool.synapses_fall_E2E = sparse(pool.N, pool.N); % double
            end

            if pool.on_GPU
                pool.refractoriness = ones(1, pool.N, "gpuArray"); % a.u.
                pool.Vm = zeros(1, pool.N, "gpuArray"); % a.u.
            else
                pool.refractoriness = ones(1, pool.N); % a.u.
                pool.Vm = zeros(1, pool.N); % a.u.
            end

            if pool.on_GPU
                pool.input_spikes = false(1, pool.N_input, "gpuArray");
                pool.spikes = false(1, pool.N, "gpuArray");
                pool.spikes_delayed_input2E = sparse( ...
                    pool.N_input * pool.N, ... % pre & post
                    pool.axonal_delay_max / pool.dt, ... % delay
                    "gpuArray" ...
                    );
                pool.spikes_delayed_E2E = sparse( ...
                    pool.N * pool.N, ... % pre & post
                    pool.axonal_delay_max / pool.dt, ... % delay
                    "gpuArray" ...
                    );
                pool.spikes_times = zeros(1, pool.N, "gpuArray");
            else
                pool.input_spikes = false(1, pool.N_input);
                pool.spikes = false(1, pool.N);
                pool.spikes_delayed_input2E = sparse( ...
                    pool.N_input * pool.N, ... % pre and post
                    pool.axonal_delay_max / pool.dt ...
                    );
                pool.spikes_delayed_E2E = sparse( ...
                    pool.N * pool.N, ... % pre and post
                    pool.axonal_delay_max / pool.dt ...
                    );
                pool.spikes_times = zeros(1, pool.N);
            end

            % Plotting
            if pool.save_to_plot_synapses
                if pool.on_GPU
                    pool.synapses_input2E_t = zeros(pool.N_input, T / pool.dt, "gpuArray");
                    pool.synapses_E2E_t = zeros(pool.N, T / pool.dt, "gpuArray");
                else
                    pool.synapses_input2E_t = zeros(pool.N_input, T / pool.dt);
                    pool.synapses_E2E_t = zeros(pool.N, T / pool.dt);
                end
            end
            if pool.save_to_plot_Vm
                if pool.on_GPU
                    pool.Vm_t = zeros(T / pool.dt, pool.N, "gpuArray");
                else
                    pool.Vm_t = zeros(T / pool.dt, pool.N);
                end
            end
            if pool.save_to_plot_spikes
                if pool.on_GPU
                    pool.input_spikes_t = false(T / pool.dt, pool.N_input, "gpuArray");
                    pool.spikes_t = false(T / pool.dt, pool.N, "gpuArray");
                else
                    pool.input_spikes_t = false(T / pool.dt, pool.N_input);
                    pool.spikes_t = false(T / pool.dt, pool.N);
                end
            end
            if pool.save_to_plot_refractoriness
                if pool.on_GPU
                    pool.refractoriness_t = ones(T / pool.dt, pool.N, "gpuArray");
                else
                    pool.refractoriness_t = ones(T / pool.dt, pool.N);
                end
            end
            if pool.save_to_plot_weights
                if pool.on_GPU
                    pool.conn_weights_input2E_t = zeros( ...
                        pool.N_input, ...
                        pool.N, ...
                        (T / pool.dt) / (pool.save_weights_t_period / pool.dt), ...
                        "gpuArray" ...
                        );
                    pool.conn_weights_E2E_t = zeros( ...
                        pool.N, ...
                        pool.N, ...
                        (T / pool.dt) / (pool.save_weights_t_period / pool.dt), ...
                        "gpuArray" ...
                        );
                else
                    pool.conn_weights_input2E_t = zeros( ...
                        pool.N_input, ...
                        pool.N, ...
                        (T / pool.dt) / (pool.save_weights_t_period / pool.dt) ...
                        );
                    pool.conn_weights_E2E_t = zeros( ...
                        pool.N, ...
                        pool.N, ...
                        (T / pool.dt) / (pool.save_weights_t_period / pool.dt) ...
                        );
                end
            end
            if pool.save_to_plot_Hz
                if pool.on_GPU
                    pool.input_spike_count_individual = zeros(1, pool.N_input, "uint32", "gpuArray");
                    pool.input_Hz_t = zeros( ...
                        (T / pool.dt) / (pool.save_Hz_t_period / pool.dt), ...
                        pool.N_input, ...
                        "gpuArray" ...
                        );
                    pool.spike_count_individual = zeros(1, pool.N, "uint32", "gpuArray");
                    pool.Hz_t = zeros( ...
                        (T / pool.dt) / (pool.save_Hz_t_period / pool.dt), ...
                        pool.N, ...
                        "gpuArray" ...
                        );
                else
                    pool.input_spike_count_individual = zeros(1, pool.N_input, "uint32");
                    pool.input_Hz_t = zeros( ...
                        (T / pool.dt) / (pool.save_Hz_t_period / pool.dt), ...
                        pool.N_input ...
                        );
                    pool.spike_count_individual = zeros(1, pool.N, "uint32");
                    pool.Hz_t = zeros( ...
                        (T / pool.dt) / (pool.save_Hz_t_period / pool.dt), ...
                        pool.N ...
                        );
                end
            end
            if pool.save_to_plot_stimulation
                if pool.on_GPU
                    pool.stimulation_t = false( ...
                        1, ...
                        (T / pool.dt) / (pool.save_stimulation_t_period / pool.dt), ...
                        "gpuArray" ...
                        );
                else
                    pool.stimulation_t = false( ...
                        1, ...
                        (T / pool.dt) / (pool.save_stimulation_t_period / pool.dt) ...
                        );
                end
            end
        end

        function set_input( ...
                pool, ...
                input_Hz ...
                )
            arguments
                pool EI_pool
                input_Hz double {mustBeFloat, mustBeNonnegative}
            end
            % set_input Sets the input firing rates

            % Input
            pool.input_Hz_in_dt = input_Hz * pool.dt;
        end

        function step( ...
                pool, ...
                t ...
                )
            arguments
                pool EI_pool
                t uint32
            end
            % step Runs a single step of the simulation

            % Input
            %   stimulation
            pool.check_input(t)
            if ~pool.input_stopped
                if pool.on_GPU
                    pool.input_spikes = rand(1, pool.N_input, "gpuArray") <= pool.input_Hz_in_dt;
                else
                    pool.input_spikes = rand(1, pool.N_input) <= pool.input_Hz_in_dt;
                end
            else
                if pool.on_GPU
                    pool.input_spikes = zeros(1, pool.N_input, "gpuArray");
                else
                    pool.input_spikes = zeros(1, pool.N_input);
                end
            end

            % Synapses
            %   Input
            pool.synapses_rise_input2E = pool.synapses_rise_input2E + ( ...
                ((-pool.synapses_rise_input2E + reshape(pool.spikes_delayed_input2E(:, 1), pool.N_input, pool.N)) ...
                / pool.synapses_glutamate_tau_rise) * pool.dt);
            pool.synapses_fall_input2E = pool.synapses_fall_input2E + ( ...
                ((-pool.synapses_fall_input2E + pool.synapses_rise_input2E) / pool.synapses_glutamate_tau_fall) * pool.dt);
            %   E➔E
            pool.synapses_rise_E2E = pool.synapses_rise_E2E + ( ...
                ((-pool.synapses_rise_E2E + reshape(pool.spikes_delayed_E2E(:, 1), pool.N, pool.N)) ...
                / pool.synapses_glutamate_tau_rise) * pool.dt);
            pool.synapses_fall_E2E = pool.synapses_fall_E2E + ( ...
                ((-pool.synapses_fall_E2E + pool.synapses_rise_E2E) / pool.synapses_glutamate_tau_fall) * pool.dt);

            % Integrate Vm
            if  pool.on_GPU
                warning("Needs reworking, see below CPU version.")
            else
                %   Excitatory
                pool.Vm(1:pool.N) = pool.Vm(1:pool.N) + ((-pool.Vm(1:pool.N) / pool.Vm_tau) * pool.dt) + ...
                    ( ...
                    ... %   input➔E
                    ((sum(pool.conn_weights_input2E .* pool.synapses_fall_input2E, 1) * ...
                    pool.conn_strength_input2E) .* pool.refractoriness) + ...
                    ... %   E➔E
                    ((sum(pool.conn_weights_E2E .* pool.synapses_fall_E2E, 1) * ...
                    pool.conn_strength_E2E) .* pool.refractoriness) ...
                    );
            end

            % Threshold
            pool.spikes = pool.Vm >= pool.Vm_thresh;
            pool.Vm(pool.spikes) = pool.Vm_base;
            %   (Gather the spikes for loops below)
            %   input➔E
            pool.spikes_delayed_input2E(:, 1) = false; % we've done this dt's spikes
            pool.spikes_delayed_input2E = circshift(pool.spikes_delayed_input2E, -1, 2); % step forward in time
            for pre = 1:pool.N_input % can't vectorize as is not all row or all column
                if pool.input_spikes(pre)
                    for post = find(pool.conn_matrix_input2E(pre, :))
                        pool.spikes_delayed_input2E( ...
                            sub2ind([pool.N_input pool.N], pre, post), ...
                            round(pool.conn_delays_input2E(pre, post)) ...
                            ) = true;
                    end
                end
            end
            %   E➔E
            pool.spikes_delayed_E2E(:, 1) = false; % ''
            pool.spikes_delayed_E2E = circshift(pool.spikes_delayed_E2E, -1, 2); % ''
            for pre = 1:pool.N % ''
                if pool.spikes(pre)
                    for post = find(pool.conn_matrix_E2E(pre, :))
                        pool.spikes_delayed_E2E( ...
                            sub2ind([pool.N pool.N], pre, post), ...
                            round(pool.conn_delays_E2E(pre, post)) ...
                            ) = true;
                    end
                end
            end

            % Refractory period
            %   absolute
            pool.spikes_times(pool.spikes) = t; % previous dt's
            pool.Vm(t <= pool.spikes_times + pool.refractoriness_absolute_in_dt) = pool.Vm_base;
            %   relative - only those that are out of the absolute period
            pool.refractoriness(pool.spikes) = 0.0; % previous dt's
            tmp_mask_refractoriness = t > pool.spikes_times + pool.refractoriness_absolute_in_dt;
            pool.refractoriness(tmp_mask_refractoriness) = pool.refractoriness(tmp_mask_refractoriness) + ...
                ((1.0 - pool.refractoriness(tmp_mask_refractoriness)) / pool.refractoriness_relative_tau) * pool.dt;

            % Update history of states (for plotting purposes)
            if pool.save_to_plot_synapses
                pool.synapses_input2E_t(:, t) = pool.synapses_fall_input2E(:, 1);
                pool.synapses_E2E_t(:, t) = pool.synapses_fall_E2E(:, 1);
            end
            if pool.save_to_plot_refractoriness
                pool.refractoriness_t(t, :) = pool.refractoriness;
            end
            if pool.save_to_plot_Vm
                pool.Vm_t(t, :) = pool.Vm;
            end
            if pool.save_to_plot_spikes
                pool.input_spikes_t(t, :) = pool.input_spikes;
                pool.spikes_t(t, :) = pool.spikes;
            end
            if pool.save_to_plot_weights
                if ~mod(t, pool.save_weights_t_period / pool.dt)
                    tmp_t = t / (pool.save_weights_t_period / pool.dt);
                    pool.conn_weights_input2E_t(:, :, tmp_t) = ...
                        pool.conn_weights_input2E;
                    pool.conn_weights_E2E_t(:, :, tmp_t) = ...
                        pool.conn_weights_E2E;
                end
            end
            if pool.save_to_plot_Hz
                % Update spike count (used to calculate Hz)
                pool.input_spike_count_individual = pool.input_spike_count_individual ...
                    + cast(pool.input_spikes, "uint32");
                pool.spike_count_individual = pool.spike_count_individual + cast(pool.spikes, "uint32");
                if ~mod(t, pool.save_Hz_t_period / pool.dt)
                    % Update Hz
                    tmp_t = t / (pool.save_Hz_t_period / pool.dt);
                    if pool.on_GPU
                        pool.input_Hz_t(tmp_t, :) = ...
                            gpuArray(cast(pool.input_spike_count_individual)) ...
                            / pool.save_Hz_t_period;
                        pool.Hz_t(tmp_t, :) = ...
                            gpuArray(cast(pool.spike_count_individual)) ...
                            / pool.save_Hz_t_period;
                    else
                        pool.input_Hz_t(tmp_t, :) = ...
                            cast(pool.input_spike_count_individual, "double") ...
                            / pool.save_Hz_t_period;
                        pool.Hz_t(tmp_t, :) = ...
                            cast(pool.spike_count_individual, "double") ...
                            / pool.save_Hz_t_period;
                    end
                    pool.spike_count_individual(:) = 0;
                    pool.input_spike_count_individual(:) = 0;
                end
            end
            if pool.save_to_plot_stimulation
                % Update spike count (used to calculate Hz)
                if ~mod(t, pool.save_stimulation_t_period / pool.dt)
                    % Update stimulation
                    pool.check_input(t);
                    pool.stimulation_t(t / (pool.save_stimulation_t_period / pool.dt), 1) = ...
                        ~pool.input_stopped;
                end
            end
        end

        function run(...
                pool, ...
                plot_during ...
                )
            arguments
                pool EI_pool
                plot_during logical
            end
            % run Runs the simulation

            % Run
            for t = 1:pool.T / pool.dt
                if ~mod(t, 1 / pool.dt)
                    % Plot the progress
                    if plot_during
                        pool.plot()
                    end
                end
                pool.step(t)
            end
        end

        function plot(pool, export)
            arguments
                pool EI_pool
                export logical
            end
            % plot Plot various aspects and statistics of the pool and
            % simulation

            pool.plot_synapses_refractoriness_Vm(export)
            pool.plot_spikes(export)
            pool.plot_weights(export)
            pool.plot_Hz(export)
            pool.plot_digraph(export)
        end

        function finish(pool)
            arguments
                pool EI_pool
            end
            % finish Finish running the pool and gather resources

            % Variables - gather either because distributed or on the GPU
            pool.conn_matrix_input2E = gather(pool.conn_matrix_input2E);
            pool.conn_weights_input2E = gather(pool.conn_weights_input2E);
            pool.conn_delays_input2E = gather(pool.conn_delays_input2E);

            pool.conn_matrix_E2E = gather(pool.conn_matrix_E2E);
            pool.conn_weights_E2E = gather(pool.conn_weights_E2E);
            pool.conn_delays_E2E = gather(pool.conn_delays_E2E);

            pool.synapses_rise_input2E = gather(pool.synapses_rise_input2E);
            pool.synapses_fall_input2E = gather(pool.synapses_fall_input2E);

            pool.synapses_rise_E2E = gather(pool.synapses_rise_E2E);
            pool.synapses_fall_E2E = gather(pool.synapses_fall_E2E);

            pool.Vm = gather(pool.Vm);

            pool.input_spikes = gather(pool.input_spikes);
            pool.spikes = gather(pool.spikes);

            pool.spikes_delayed_input2E = gather(pool.spikes_delayed_input2E);
            pool.spikes_delayed_E2E = gather(pool.spikes_delayed_E2E);
            pool.spikes_times = gather(pool.spikes_times);

            pool.refractoriness = gather(pool.refractoriness);

            % Plotting
            pool.conn_weights_input2E_t = gather(pool.conn_weights_input2E_t);
            pool.conn_weights_E2E_t = gather(pool.conn_weights_E2E_t);

            pool.synapses_input2E_t = gather(pool.synapses_input2E_t);
            pool.synapses_E2E_t = gather(pool.synapses_E2E_t);

            pool.Vm_t = gather(pool.Vm_t);

            pool.refractoriness_t = gather(pool.refractoriness_t);

            pool.input_spikes_t = gather(pool.input_spikes_t);
            pool.spikes_t = gather(pool.spikes_t);

            pool.input_spike_count_individual = gather(pool.input_spike_count_individual);
            pool.input_Hz_t = gather(pool.input_Hz_t);
            pool.spike_count_individual = gather(pool.spike_count_individual);
            pool.Hz_t = gather(pool.Hz_t);

            pool.stimulation_t = gather(pool.stimulation_t);
        end
    end
end
