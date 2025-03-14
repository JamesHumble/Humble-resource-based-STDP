classdef RB_STDP_EI_pool < EI_pool
    % RB_STDP_EI_pool A pool of exictatory and inhibitory neurons with all plastic
    % connections (potentially on the GPU)

    properties (Access = public)
        % Parameters
        input_stop_T double {mustBeFloat, mustBePositive} = 20.0 % sec
        
        STDP_tau double {mustBeFloat, mustBePositive} = 0.020 % sec

        AD logical = false

        % Variables
        plasticity_type string

        STDP_traces (1, :) double {mustBeFloat, mustBeNonnegative} = []

        resource_neighbor_sharing_profile (1, :) double {mustBeFloat, mustBeNonnegative} = []
        resource_neighbor_sharing_max double {mustBeFloat, mustBePositive} = 3 % IF CHANGED, CHANGE NEIGHBORS ARRAY BELOW
        resource_pools (1, :) {mustBeFloat, mustBeNonnegative} = [] % a.u.
        resource_pools_tau double {mustBeFloat, mustBeNonnegative} = 10.0 % sec

        % Plotting
        save_to_plot_STDP_traces logical
        save_to_plot_resource_pools logical
        save_to_plot_non_potentiable logical
        save_to_plot_dependence logical
        save_to_plot_AD logical

        STDP_traces_t (:, :) double {mustBeFloat, mustBeNonnegative} = [] % t, post

        save_resource_pools_t_period double {mustBeFloat, mustBePositive} = 0.100 % sec
        resource_pools_t (:, :) double {mustBeFloat, mustBeNonnegative} = [] % t, post (excitatory)

        non_potentiable_weight_t (:, 1) double {mustBeFloat, mustBeNonnegative} = []
        non_potentiable_required_t (:, 1) double {mustBeFloat, mustBeNonnegative} = []
        non_potentiable_got_t (:, 1) double {mustBeFloat, mustBeNonnegative} = []
        non_potentiable_i uint32 {mustBeNonnegative} = 1
        non_potentiable_N uint32 {mustBeNonnegative} = 20000000
        non_potentiable_n uint32 {mustBeNonnegative} = 1
        save_potentiation_N_period double {mustBeFloat, mustBePositive} = 1000 % potentiation events

        dependence_weight_t (:, 1) double {mustBeFloat, mustBeNonnegative} = []
        dependence_t (:, 1) double {mustBeFloat, mustBeNonpositive} = []
        dependence_i uint32 {mustBeNonnegative} = 1
        dependence_N uint32 {mustBeNonnegative} = 20000000
        dependence_n uint32 {mustBeNonnegative} = 1
        save_dependence_N_period double {mustBeFloat, mustBePositive} = 1000 % potentiation events

        AD_delete_period double {mustBeFloat, mustBePositive} = 1.000 % sec
        AD_delete_amount double {mustBeFloat, mustBePositive} = 1 % number of synapses
        conn_matrix_E2E_t (:, :, :) double {mustBeFloat} = [] % pre, post, t
    end

    methods (Access = protected)
        function plot_STDP_traces(pool, export)
            arguments
                pool EI_pool
                export logical
            end
            % plot_STDP_traces Plots the network's excitatory STDP traces

            tmp_STDP_traces_t = gather(pool.STDP_traces_t);
            if ~isempty(tmp_STDP_traces_t)

                fig = figure();
                plot(tmp_STDP_traces_t);
                title("STDP traces")
                xlabel("Time (s)")
                xticks(0:0.1 / pool.dt:pool.T / pool.dt)
                xticklabels(0:0.1:pool.T)
                ylabel("(a.u.)")

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/STDP_traces.pdf','ContentType','vector')
                end
            end
        end

        function plot_resource_pools(pool, export)
            arguments
                pool EI_pool
                export logical
            end
            % plot_resource_pools Plots the network's resource pool

            tmp_resource_pools_t = gather(pool.resource_pools_t);
            if ~isempty(tmp_resource_pools_t)

                fig = figure();

                plot(tmp_resource_pools_t, "LineWidth", 2);
                box off
                xlabel("Time (s)")
                xticks(0:10.0 / pool.save_resource_pools_t_period:pool.T / pool.save_resource_pools_t_period)
                xticklabels(0:10:pool.T)
                yticklabels([0])
                yticks([0])
                ylabel("$p$", "Interpreter", "latex")
                % set(gca,'YScale','log')

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/resource_pools.pdf','ContentType','vector')
                end
            end
        end

        function plot_non_potentiable(pool, export)
            arguments
                pool EI_pool
                export logical
            end
            % plot_non_potentiable Plots the STDP's non-potentiability

            tmp_non_potentiable_weight_t = gather(pool.non_potentiable_weight_t);
            if ~isempty(tmp_non_potentiable_weight_t)
                tmp_non_potentiable_required_t = gather(pool.non_potentiable_required_t);
                tmp_non_potentiable_got_t = gather(pool.non_potentiable_got_t);

                fig = figure();
                scatter( ...
                    tmp_non_potentiable_weight_t(1:pool.non_potentiable_i - 1), ...
                    (tmp_non_potentiable_got_t(1:pool.non_potentiable_i - 1) ...
                    ./ tmp_non_potentiable_required_t(1:pool.non_potentiable_i - 1)) * 100, ...
                    '.' ...
                    );
                box off
                xlabel("$w$", "Interpreter", "latex")
                ylabel("$f(\tau)$ / Resources available (\%)", "Interpreter", "latex")
                ylim([0 100])

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/non_potentiable.pdf','ContentType','vector')
                end

                fig = figure();
                histogram( ...
                    (tmp_non_potentiable_got_t(1:pool.non_potentiable_i - 1) ...
                    ./ tmp_non_potentiable_required_t(1:pool.non_potentiable_i - 1)) * 100, ...
                    30, ...
                    "FaceColor", ...
                    "black" ...
                    );
                box off
                set(gca,'YScale','log')
                xlabel("$f(\tau)$ / Resources available (\%)", "Interpreter", "latex")
                ylabel("Frequency")

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/non_potentiable_histogram.pdf','ContentType','vector')
                end

                fig = figure();
                scatter( ...
                    tmp_non_potentiable_required_t(1:pool.non_potentiable_i - 1), ...
                    tmp_non_potentiable_got_t(1: pool.non_potentiable_i - 1), ...
                    '.' ...
                    );
                box off
                xlabel("Required")
                ylabel("Got")

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/non_potentiable_required_vs_got.pdf','ContentType','vector')
                end

                fig = figure();
                scatter( ...
                    tmp_non_potentiable_weight_t(1:pool.non_potentiable_i - 1), ...
                    (tmp_non_potentiable_got_t(1:pool.non_potentiable_i - 1) ...
                    ./ tmp_non_potentiable_weight_t(1:pool.non_potentiable_i - 1)) * 100, ...
                    '.' ...
                    );
                box off
                xlabel("$w$", "Interpreter", "latex")
                set(gca,'XScale','log')
                xlim([0.1 max(tmp_non_potentiable_weight_t(1:pool.non_potentiable_i - 1))])
                ylabel("$\Delta w$ (\%)", "Interpreter", "latex")

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/non_potentiable_weight_vs_got.pdf','ContentType','vector')
                end
            end
        end

        function plot_dependence(pool, export)
            arguments
                pool EI_pool
                export logical
            end
            % plot_dependence Plots the STDP's depression weight dependence

            tmp_dependence_weight_t = gather(pool.dependence_weight_t);
            if ~isempty(tmp_dependence_weight_t)
                tmp_dependence_t = gather(pool.dependence_t);

                fig = figure();
                scatter( ...
                    tmp_dependence_weight_t(1:10:pool.dependence_i - 1), ...
                    (tmp_dependence_t(1:10:pool.dependence_i - 1) ...
                    ./ tmp_dependence_weight_t(1:10:pool.dependence_i - 1)) * 100, ...
                    "filled", ...
                    "blue" ...
                    );

                tmp_non_potentiable_weight_t = gather(pool.non_potentiable_weight_t);
                if ~isempty(tmp_non_potentiable_weight_t)
                    tmp_non_potentiable_got_t = gather(pool.non_potentiable_got_t);
                    hold on
                    scatter( ...
                        tmp_non_potentiable_weight_t(1:10:pool.non_potentiable_i - 1), ...
                        (tmp_non_potentiable_got_t(1:10:pool.non_potentiable_i - 1) ...
                        ./ tmp_non_potentiable_weight_t(1:10:pool.non_potentiable_i - 1)) * 100, ...
                        "filled", ...
                        "blue" ...
                        );
                    hold off
                end

                box off
                xlabel("$w$", "Interpreter", "latex")
                set(gca,'XScale','log')
                xlim([0.5 max(tmp_non_potentiable_weight_t(1:pool.non_potentiable_i - 1))])
                ylim([-0.5 2.5])
                yline(0, "--")
                ylabel("$\Delta w$ (\%)", "Interpreter", "latex")
                legend off

                axes("Position", [0.7 0.75 0.2 0.15]);
                if ~isempty(tmp_non_potentiable_weight_t)
                    tmp_non_potentiable_got_t = gather(pool.non_potentiable_got_t);
                    hold on
                    scatter( ...
                        tmp_non_potentiable_weight_t(1:100:pool.non_potentiable_i - 1), ...
                        (tmp_non_potentiable_got_t(1:100:pool.non_potentiable_i - 1) ...
                        ./ tmp_non_potentiable_weight_t(1:100:pool.non_potentiable_i - 1)) * 100, ...
                        ".", ...
                        "blue" ...
                        );
                    hold off
                end
                box off
                set(gca,'XScale','log')
                xticks([0.001, 0.01, 1])
                yticks([1 10 100])
                yline(0, "--")
                set(gca,'YScale','log')

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/dependence.pdf','ContentType','vector')
                end
            end
        end

        function plot_AD(pool, export)
            arguments
                pool EI_pool
                export logical
            end
            % plot_AD Plots the network Alzhemeier's disease components

            tmp_conn_matrix_E2E_t = gather(pool.conn_matrix_E2E_t);
            if ~isempty(tmp_conn_matrix_E2E_t)

                tmp_connection_sum = sum(reshape(tmp_conn_matrix_E2E_t, pool.N * pool.N, []), 1);

                fig = figure();
                plot( ...
                    (tmp_connection_sum / max(tmp_connection_sum)) * 100, ...
                    "black", ...
                    "LineWidth", ...
                    2 ...
                    );

                box off
                xlabel("Time (s)")
                xticks(0:10.0 / pool.save_resource_pools_t_period:pool.T / pool.save_resource_pools_t_period)
                xticklabels(0:10:pool.T)
                ylabel("% of pre learning connections")

                xline(528, "--")
                yline(17, "--", "17 %")

                xline(784, "--")
                yline(72, "--", "72 %")

                fontsize(fig, pool.fontsize, "points")
                drawnow
                if export
                    exportgraphics(gcf,'figures/connections.pdf','ContentType','vector')
                end
            end
        end

        function check_input(pool, t)
            arguments
                pool EI_pool
                t uint32
            end
            % check_input Checks whether stimulation is finished or not

            if t * pool.dt == pool.input_stop_T
                pool.input_stopped = true;
            end
        end

        function check_AD(pool, t)
            arguments
                pool EI_pool
                t uint32
            end
            % check_AD Checks whether stimulation is Alzheimer's disease or not

            if t * pool.dt >= 40 && ~mod(t, pool.AD_delete_period / pool.dt)
                % Delete connections
                for post = 1:pool.N
                    tmp_current_conns = find(pool.conn_matrix_E2E(:, post));
                    if ~isempty(tmp_current_conns)
                        tmp_subset_to_delete = ...
                            randsample(tmp_current_conns, pool.AD_delete_amount);
                        % pool.resource_pools(post) = pool.resource_pools(post) + ...
                        %     sum(pool.conn_weights_E2E(tmp_subset_to_delete, post));
                        pool.conn_weights_E2E(tmp_subset_to_delete, post) = 0.0;
                        pool.conn_matrix_E2E(tmp_subset_to_delete, post) = false;
                    end
                end
            end
        end
    end

    methods
        function pool = RB_STDP_EI_pool( ...
                T, ...
                N, ...
                on_GPU, ...
                conn_prob_input2E, ...
                conn_strength_input2E, ...
                conn_prob_E2E, ...
                conn_strength_E2E, ...
                plasticity_type, ...
                save_to_plot_synapses, ...
                save_to_plot_Vm, ...
                save_to_plot_spikes, ...
                save_to_plot_refractoriness, ...
                save_to_plot_weights, ...
                save_to_plot_Hz, ...
                save_to_plot_stimulation, ...
                save_to_plot_STDP_traces, ...
                save_to_plot_resource_pools, ...
                b_plot_digraph, ...
                save_to_plot_non_potentiable, ...
                save_to_plot_dependence, ...
                AD, ...
                save_to_plot_AD ...
                )
            arguments
                T double {mustBeFloat, mustBePositive}
                N uint32 {mustBeInteger, mustBePositive}
                on_GPU logical = true
                conn_prob_input2E double {mustBeFloat, mustBeNonnegative} = 0.0 % %
                conn_strength_input2E double {mustBeFloat, mustBeNonnegative} = 0.0 % a.u.
                conn_prob_E2E double {mustBeFloat, mustBeNonnegative} = 0.0 % %
                conn_strength_E2E double {mustBeFloat, mustBeNonnegative} = 0.0 % a.u.
                plasticity_type string = ""
                save_to_plot_synapses logical = false
                save_to_plot_Vm logical = false
                save_to_plot_spikes logical = false
                save_to_plot_refractoriness logical = false
                save_to_plot_weights logical = false
                save_to_plot_Hz logical = false
                save_to_plot_stimulation logical = false
                save_to_plot_STDP_traces logical = false
                save_to_plot_resource_pools logical = false
                b_plot_digraph logical = false
                save_to_plot_non_potentiable logical = false
                save_to_plot_dependence logical = false
                AD logical = false
                save_to_plot_AD logical = false
            end
            % RB_STDP_EI_pool Construct an instance of an EI pool (on the GPU)
            pool = pool@EI_pool( ...
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
                );

            pool.conn_strength_input2E = conn_strength_input2E;
            pool.AD = AD;
            pool.plasticity_type = plasticity_type;

            % Plotting
            pool.save_to_plot_STDP_traces = save_to_plot_STDP_traces;
            pool.save_to_plot_resource_pools = save_to_plot_resource_pools;
            pool.save_to_plot_non_potentiable = save_to_plot_non_potentiable;
            pool.save_to_plot_dependence = save_to_plot_dependence;
            pool.save_to_plot_AD = save_to_plot_AD;

            % General
            % none

            % Input
            % none

            % Weights
            %   reset weights to Uniform
            %   input➔E
            pool.conn_weights_input2E = ...
                rand(pool.N_input, pool.N) .* sparse(rand(pool.N_input, pool.N) <= pool.conn_prob_input2E);
            pool.conn_weights_input2E = min(pool.conn_weights_input2E, pool.conn_max);
            pool.conn_weights_input2E = max(pool.conn_weights_input2E, 0.0);
            pool.conn_matrix_input2E = pool.conn_weights_input2E > 0;
            pool.conn_delays_input2E = ( ...
                randi((pool.axonal_delay_max - pool.axonal_delay_min) / ...
                pool.dt, pool.N_input, pool.N ...
                ) + (pool.axonal_delay_min / pool.dt)) .* pool.conn_matrix_input2E;
            %   E➔E
            pool.conn_weights_E2E = 0.0;

            if pool.plasticity_type == "rbSTDP"
                % Resources
                pool.resource_neighbor_sharing_profile = zeros(1, length(pool.resource_neighbor_sharing_max));
                for n = 1:pool.resource_neighbor_sharing_max
                    pool.resource_neighbor_sharing_profile(n) = exp(-n);
                end
                pool.resource_neighbor_sharing_profile = ...
                    pool.resource_neighbor_sharing_profile / sum(pool.resource_neighbor_sharing_profile) / 2;
            end

            % Neurons
            % none

            % States
            if pool.on_GPU
                pool.STDP_traces = sparse(1, pool.N, "gpuArray"); % so keeps STDP update sparse
                if pool.plasticity_type == "rbSTDP"
                    % Comment/uncomment for different initial distributions
                    % pool.resource_pools = rand(1, pool.N) * 4;
                    pool.resource_pools = exp(randn(1, pool.N)) * 2;
                end
            else
                pool.STDP_traces = sparse(1, pool.N); % ''
                if pool.plasticity_type == "rbSTDP"
                    % Comment/uncomment for different initial distributions
                    % pool.resource_pools = rand(1, pool.N) * 4;
                    pool.resource_pools = exp(randn(1, pool.N)) * 2;
                end
            end

            % Plotting
            if pool.save_to_plot_STDP_traces
                if pool.on_GPU
                    pool.STDP_traces_t = zeros(T / pool.dt, pool.N, "gpuArray");
                else
                    pool.STDP_traces_t = zeros(T / pool.dt, pool.N);
                end
            end
            if pool.plasticity_type == "rbSTDP" && pool.save_to_plot_resource_pools
                if pool.on_GPU
                    pool.resource_pools_t = zeros( ...
                        (T / pool.dt) / (pool.save_resource_pools_t_period / pool.dt), ...
                        pool.N, ...
                        "gpuArray" ...
                        );
                else
                    pool.resource_pools_t = zeros( ...
                        (T / pool.dt) / (pool.save_resource_pools_t_period / pool.dt), ...
                        pool.N ...
                        );
                end
            end
            if pool.save_to_plot_non_potentiable
                if pool.on_GPU
                    pool.non_potentiable_weight_t = zeros(1, pool.non_potentiable_N, "gpuArray");
                    pool.non_potentiable_required_t = zeros(1, pool.non_potentiable_N, "gpuArray");
                    pool.non_potentiable_got_t = zeros(1, pool.non_potentiable_N, "gpuArray");
                else
                    pool.non_potentiable_weight_t = zeros(1, pool.non_potentiable_N);
                    pool.non_potentiable_required_t = zeros(1, pool.non_potentiable_N);
                    pool.non_potentiable_got_t = zeros(1, pool.non_potentiable_N);
                end
            end
            if pool.save_to_plot_non_potentiable
                if pool.on_GPU
                    pool.dependence_weight_t = zeros(1, pool.dependence_N, "gpuArray");
                    pool.dependence_t = zeros(1, pool.dependence_N, "gpuArray");
                else
                    pool.dependence_weight_t = zeros(1, pool.dependence_N);
                    pool.dependence_t = zeros(1, pool.dependence_N);
                end
            end
            if pool.save_to_plot_AD
                if pool.on_GPU
                    pool.conn_matrix_E2E_t = zeros( ...
                        pool.N, ...
                        pool.N, ...
                        (T / pool.dt) / (pool.save_weights_t_period / pool.dt), ...
                        "gpuArray" ...
                        );
                else
                    pool.conn_matrix_E2E_t = zeros( ...
                        pool.N, ...
                        pool.N, ...
                        (T / pool.dt) / (pool.save_weights_t_period / pool.dt) ...
                        );
                end
            end

            % Plot
            % none
        end

        function step( ...
                pool, ...
                t, ...
                STDP_learning_rate_E2E ...
                )
            arguments
                pool RB_STDP_EI_pool
                t uint32
                STDP_learning_rate_E2E double {mustBeFloat, mustBeNonnegative} % a.u.
            end
            % step Runs a single step of the simulation

            % Check if we're running with Alzheimer's disease
            if pool.AD
                pool.check_AD(t);
            end

            % Plasticity
            %   E➔E - eSTDP
            if STDP_learning_rate_E2E > 0
                switch pool.plasticity_type
                    case "rbSTDP"
                        % LTP - if post spike, ...
                        tmp_post_ind = find(pool.spikes(1:pool.N));
                        for post = 1:length(tmp_post_ind)
                            % ... take resources from:
                            %   1) nearby pre connections (assuming they are equidistant) and
                            %   2) the pool (if necessary)
                            tmp_pre_ind = find(pool.conn_matrix_E2E(1:pool.N, tmp_post_ind(post)));
                            tmp_pre_ind_N = length(tmp_pre_ind);
                            for pre = 1:length(tmp_pre_ind)
                                tmp_amount_req = pool.STDP_traces(tmp_pre_ind(pre)) * STDP_learning_rate_E2E;
                                tmp_amount_got = 0.0;

                                % 1) If there are enough resources in neighbors, take from them ...
                                for neighbor_offset = [-3, -2, -1, 1, 2, 3]
                                    neighbor = pre + neighbor_offset;
                                    % Check each neighbor - only if a valid neighbor
                                    if neighbor >= 1 && neighbor <= tmp_pre_ind_N
                                        tmp_adj_amount_req = tmp_amount_req ...
                                            * pool.resource_neighbor_sharing_profile(abs(neighbor_offset));
                                        if pool.conn_weights_E2E(tmp_pre_ind(neighbor), tmp_post_ind(post)) ...
                                                >= tmp_adj_amount_req
                                            % If there is enough, take it ...
                                            pool.conn_weights_E2E(tmp_pre_ind(pre), tmp_post_ind(post)) = ...
                                                pool.conn_weights_E2E(tmp_pre_ind(pre), tmp_post_ind(post)) + ...
                                                tmp_adj_amount_req;
                                            pool.conn_weights_E2E(tmp_pre_ind(neighbor), tmp_post_ind(post)) = ...
                                                pool.conn_weights_E2E(tmp_pre_ind(neighbor), tmp_post_ind(post)) - ...
                                                tmp_adj_amount_req;
                                            tmp_amount_got = tmp_amount_got + tmp_adj_amount_req;
                                        else
                                            % ... if not, take what there is
                                            pool.conn_weights_E2E(tmp_pre_ind(pre), tmp_post_ind(post)) = ...
                                                pool.conn_weights_E2E(tmp_pre_ind(pre), tmp_post_ind(post)) + ...
                                                pool.conn_weights_E2E(tmp_pre_ind(neighbor), tmp_post_ind(post));
                                            tmp_amount_got = tmp_amount_got + ...
                                                pool.conn_weights_E2E(tmp_pre_ind(neighbor), tmp_post_ind(post));
                                            pool.conn_weights_E2E(tmp_pre_ind(neighbor), tmp_post_ind(post)) = 0.0;
                                        end
                                    end
                                end

                                % 2) If there are enough resources in the pool take them (if needed) ...
                                if pool.resource_pools(tmp_post_ind(post)) > 0 && tmp_amount_got < tmp_amount_req
                                    tmp_amount_left = tmp_amount_req - tmp_amount_got;
                                    if pool.resource_pools(tmp_post_ind(post)) >= tmp_amount_left
                                        pool.conn_weights_E2E(tmp_pre_ind(pre), tmp_post_ind(post)) = ...
                                            pool.conn_weights_E2E(tmp_pre_ind(pre), tmp_post_ind(post)) ...
                                            + tmp_amount_left;
                                        tmp_amount_got = tmp_amount_got + tmp_amount_left;
                                        pool.resource_pools(tmp_post_ind(post)) = ...
                                            pool.resource_pools(tmp_post_ind(post)) - tmp_amount_left;
                                    else
                                        % ... else, take what there is
                                        pool.conn_weights_E2E(tmp_pre_ind(pre), tmp_post_ind(post)) = ...
                                            pool.conn_weights_E2E(tmp_pre_ind(pre), tmp_post_ind(post)) ...
                                            + pool.resource_pools(tmp_post_ind(post));
                                        tmp_amount_got = tmp_amount_got + ...
                                            pool.resource_pools(tmp_post_ind(post));
                                        pool.resource_pools(tmp_post_ind(post)) = 0.0;
                                    end
                                end

                                % Save non-potentiation
                                if pool.save_to_plot_non_potentiable 
                                    if ~mod(pool.non_potentiable_n, pool.save_potentiation_N_period) && ...
                                            tmp_amount_req > 0 % && tmp_amount_req ~= tmp_amount_got
                                        pool.non_potentiable_weight_t(pool.non_potentiable_i) = ...
                                            pool.conn_weights_E2E(tmp_pre_ind(pre), tmp_post_ind(post));
                                        pool.non_potentiable_required_t(pool.non_potentiable_i) = tmp_amount_req;
                                        pool.non_potentiable_got_t(pool.non_potentiable_i) = tmp_amount_got;
                                        pool.non_potentiable_i = pool.non_potentiable_i + 1;
                                        if pool.non_potentiable_i > pool.non_potentiable_N
                                            warning("Overflow non_potentiable arrays.")
                                        end
                                        pool.non_potentiable_n = 1;
                                    end
                                    pool.non_potentiable_n = pool.non_potentiable_n + 1;
                                end
                            end
                        end

                        pool.conn_weights_E2E = pool.conn_weights_E2E .* pool.conn_matrix_E2E;

                        % LTD - normal additive LTD plus resources go in to the pool
                        tmp_weights_before = pool.conn_weights_E2E;
                        pool.conn_weights_E2E = pool.conn_weights_E2E + ...
                            ( ...
                            ( ...
                            - ( ...
                            ... % post-pre LTD - if pre spike, update with the post STDP traces
                            reshape(pool.spikes_delayed_E2E(:, 1), pool.N, pool.N) ...
                            .* repmat(pool.STDP_traces(1:pool.N), pool.N, 1)' ...
                            .* pool.conn_weights_E2E ...
                            ) ...
                            ) * STDP_learning_rate_E2E * 0.18 ...
                            );

                        pool.conn_weights_E2E = pool.conn_weights_E2E .* pool.conn_matrix_E2E;

                        % Limit weights
                        pool.conn_weights_E2E = max(pool.conn_weights_E2E, 0.0);

                        % Save weight dependence
                        if pool.save_to_plot_dependence
                            if ~mod(pool.dependence_n, pool.save_dependence_N_period)
                                for pre = 1:pool.N
                                    for post = 1:pool.N
                                        if pool.conn_matrix_E2E(pre, post) && ...
                                                pool.conn_weights_E2E(pre, post) ~= tmp_weights_before(pre, post)
                                            pool.dependence_weight_t(pool.dependence_i) = ...
                                                tmp_weights_before(pre, post);
                                            pool.dependence_t(pool.dependence_i) = ...
                                                pool.conn_weights_E2E(pre, post) - tmp_weights_before(pre, post);
                                            pool.dependence_i = pool.dependence_i + 1;
                                            if pool.dependence_i > pool.dependence_N
                                                warning("Overflow dependence array.")
                                            end
                                            pool.dependence_n = 1;
                                        end
                                    end
                                end
                            end
                            pool.dependence_n = pool.dependence_n + 1;
                        end

                        % Update resource pools with any LTD-based additional resources
                        pool.resource_pools = pool.resource_pools + sum(tmp_weights_before - pool.conn_weights_E2E, 1);
                    case "addSTDP"
                        % STDP
                        pool.conn_weights_E2E = pool.conn_weights_E2E + ...
                            ( ...
                            ( ...
                            + ( ...
                            ... % pre-post LTP - if post spike, update with the pre STDP traces
                            repmat(pool.spikes .* pool.STDP_traces, pool.N, 1)' ...
                            ) ...
                            - ( ...
                            ... % post-pre LTD - if pre spike, update with the post STDP traces
                            reshape(pool.spikes_delayed_E2E(:, 1), pool.N, pool.N) ...
                            .* repmat(pool.STDP_traces, pool.N, 1)' ...
                            ) ...
                            ) * STDP_learning_rate_E2E ...
                            );
                        pool.conn_weights_E2E = pool.conn_weights_E2E .* pool.conn_matrix_E2E;

                        % Limit weights
                        pool.conn_weights_E2E = min(pool.conn_weights_E2E, pool.conn_max);
                        pool.conn_weights_E2E = max(pool.conn_weights_E2E, 0.0);
                    case "mulSTDP"
                        % STDP
                        pool.conn_weights_E2E = pool.conn_weights_E2E + ...
                            ( ...
                            ( ...
                            + ( ...
                            ... % pre-post LTP - if post spike, update with the pre STDP traces
                            repmat(pool.spikes .* pool.STDP_traces, pool.N, 1)' ...
                            ) ...
                            - ( ...
                            ... % post-pre LTD - if pre spike, update with the post STDP traces
                            reshape(pool.spikes_delayed_E2E(:, 1), pool.N, pool.N) ...
                            .* repmat(pool.STDP_traces, pool.N, 1)' ...
                            .* (pool.conn_weights_E2E) ...
                            ) ...
                            ) * STDP_learning_rate_E2E ...
                            );
                        pool.conn_weights_E2E = pool.conn_weights_E2E .* pool.conn_matrix_E2E;

                        % Limit weights
                        pool.conn_weights_E2E = min(pool.conn_weights_E2E, pool.conn_max);
                        pool.conn_weights_E2E = max(pool.conn_weights_E2E, 0.0);
                    otherwise
                        disp("Unknown plasticity type")
                end
            end

            % Background spikes
            % Input spikes
            % Synapses
            % Integrame Vm
            % Limit inhibition
            % Threshold
            % Refractory period
            % Update history of states (for plotting purposes)
            step@EI_pool(pool, t);

            % STDP traces
            pool.STDP_traces = pool.STDP_traces + ...
                (((-pool.STDP_traces + pool.spikes) / pool.STDP_tau) * pool.dt);

            % Resource pool decay
            if pool.plasticity_type == "rbSTDP"
                pool.resource_pools = pool.resource_pools + ((-pool.resource_pools / pool.resource_pools_tau) * pool.dt);
            end

            % Update history of states (for plotting purposes)
            if pool.save_to_plot_STDP_traces
                pool.STDP_traces_t(t, :) = pool.STDP_traces(1, :);
            end
            if pool.plasticity_type == "rbSTDP" && pool.save_to_plot_resource_pools
                if ~mod(t, pool.save_resource_pools_t_period / pool.dt)
                    pool.resource_pools_t(t / (pool.save_resource_pools_t_period / pool.dt), :) =  pool.resource_pools;
                end
            end
            if pool.save_to_plot_AD
                if ~mod(t, pool.save_weights_t_period / pool.dt)
                    pool.conn_matrix_E2E_t(:, :, t / (pool.save_weights_t_period / pool.dt)) = ...
                        pool.conn_matrix_E2E;
                end
            end
        end

        function run(...
                pool, ...
                plot_during, ...
                STDP_learning_rate_E2E ...
                )
            arguments
                pool RB_STDP_EI_pool
                plot_during logical
                STDP_learning_rate_E2E double {mustBeFloat, mustBeNonnegative} % a.u.
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
                % step@EI_pool inside
                pool.step( ...
                    t, ...
                    STDP_learning_rate_E2E ...
                    )
            end
        end

        function plot(pool, export)
            arguments
                pool RB_STDP_EI_pool
                export logical
            end
            % plot Plot various aspects and statistics of the pool and
            % simulation

            plot@EI_pool(pool, export);

            pool.plot_STDP_traces(export)
            pool.plot_resource_pools(export)
            pool.plot_non_potentiable(export)
            pool.plot_dependence(export)
            pool.plot_AD(export)
        end

        function finish(pool)
            arguments
                pool RB_STDP_EI_pool
            end
            % finish Finish running the pool and gather resources
            finish@EI_pool(pool);

            % Variables - gather either because distributed or on the GPU
            pool.STDP_traces = gather(pool.STDP_traces);

            pool.resource_pools = gather(pool.resource_pools);

            pool.non_potentiable_required_t = gather(pool.non_potentiable_required_t);
            pool.non_potentiable_got_t = gather(pool.non_potentiable_got_t);

            pool.dependence_weight_t = gather(pool.dependence_weight_t);
            pool.dependence_t = gather(pool.dependence_t);

            % Plotting
            pool.STDP_traces_t = gather(pool.STDP_traces_t);
            pool.resource_pools_t = gather(pool.resource_pools_t);
            pool.conn_matrix_E2E_t = gather(pool.conn_matrix_E2E_t);
        end
    end
end
