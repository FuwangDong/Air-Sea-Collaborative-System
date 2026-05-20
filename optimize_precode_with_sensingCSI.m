function [w_comm, w_sens, power_vec_c, power_vec_s] = optimize_precode_with_sensingCSI(N, K, s, jk, b, com_lb, sen_lb)
    % 参数设置
  %  com_lb = 13;
  %  sen_lb = 15;
    M = 4;  % 天线数
    p_max = 20;  % 最大功率
    N0 = 1e-14;  % 通信噪声
    N1 = 1e-14;  % 感知噪声
    acc = 2;  % 常数，控制干扰项权重
    shijianxishu = 80; %时间系数
    com_lb_g = 2^(com_lb) - 1;  % 通信SINR门限gamma

    %% 信道构造
    h_c = cell(1, N);
    H_c = cell(1, N);
    for n = 1:N
        h_c{n} = generate_fixed_position_channel_c(M, s, b(:,n));
        H_c{n} = h_c{n} * h_c{n}';
    end

    h_s_k = cell(K, N);
    H_s_k = cell(K, N);
    for k = 1:K
        for n = 1:N
            h_s_k{k, n} = generate_fixed_position_channel_s(M, s, jk(:, k));
            H_s_k{k, n} = h_s_k{k, n} * h_s_k{k, n}';
        end
    end

    %% SCA迭代初始化
    MaxSCA = 80;  % 最大SCA迭代次数
    tol = 1e-3;  % 收敛容忍度

    % 初始化参考点
    mu_ref_sca = 1e-3 * ones(K, N);
    phi_ref_sca = 1e-3 * ones(K, N);

    % 松弛变量（用于约束的松弛）
    lambda_dc  = 100;  % DC松弛项权重，较小值
    lambda_sen = 100;  % 感知松弛项权重
    lambda_com = 100;  % 通信松弛项权重
    lambda_pow = 100;  % 功率松弛项权重

    prev_obj_value = Inf;  % 记录前一次的目标函数值

    for iter = 1:MaxSCA
        %fprintf("SCA iteration %d\n", iter);

        % 清洗参考点，避免非法值
        mu_ref_sca = max(real(double(mu_ref_sca)), 1e-6);
        phi_ref_sca = max(real(double(phi_ref_sca)), 1e-6);

        % 进入CVX之前更新参考点
        cvx_begin sdp quiet
        cvx_solver mosek

        % CVX变量
        variable Wc(M, M, N) hermitian semidefinite  % 通信波束
        variable Ws(M, M, K, N) hermitian semidefinite  % 感知波束
        variable mu_sca(K, N)  % 感知信号的信噪比
        variable phi_sca(K, N)  % 干扰项的限制
        variable s_dc(K, N)  % DC松弛变量
        variable s_sen(K)  % 感知松弛变量
        variable s_com(N)  % 通信松弛变量
        variable s_pow(N)  % 功率松弛变量

        % 目标函数：最小化总功率 + 所有松弛变量的惩罚项
        expression total_power
        total_power = 0;
        for n = 1:N
            total_power = total_power + trace(Wc(:, :, n));
            for k = 1:K
                total_power = total_power + trace(Ws(:, :, k, n));
            end
        end

        minimize(total_power + lambda_dc * sum(s_dc(:)) + lambda_sen * sum(s_sen) + lambda_com * sum(s_com) + lambda_pow * sum(s_pow))

        subject to
            % 变量域限制
           % mu_sca >= 0;
           % phi_sca >= 0;
           % s_com >= 0;
           % s_pow >= 0;
           % for k = 1:K
           %     s_dc(k, :) >= 0;
           %     s_sen(k) >= 0;
           % end
            for n = 1:N
                s_com(n) >= 0;
                s_pow(n) >= 0;
                for k = 1:K
                    mu_sca(k,n)  >= 0;
                    phi_sca(k,n) >= 0;
                    s_dc(k,n)    >= 0;
                end
            end
            for k = 1:K
                s_sen(k) >= 0;
            end

            %% 通信SINR约束 + 松弛
            for n = 1:N
                expression interf
                interf = N0;
                for k = 1:K
                    interf = interf + acc * trace(H_c{n} * Ws(:, :, k, n));
                end
                acc * trace(H_c{n} * Wc(:, :, n)) + s_com(n) >= com_lb_g * interf;
            end

            %% 累计感知约束 + 松弛
            for k = 1:K
                sum(mu_sca(k, :)) + s_sen(k) >= sen_lb;
            end

            %% 干扰项 + 松弛
            for k = 1:K
                for n = 1:N
                    expression interf_s
                    interf_s = N1;
                    for i = 1:K
                        if i ~= k
                            interf_s = interf_s + acc * trace(H_s_k{k, n} * Ws(:, :, i, n));
                        end
                    end
                    interf_s = interf_s + acc * trace(H_s_k{k, n} * Wc(:, :, n));
                    interf_s <= phi_sca(k, n);
                end
            end

            %% DC约束 + 松弛
            for k = 1:K
                for n = 1:N
                    expression nu
                    nu = 0.5 * square_pos(mu_sca(k, n) + phi_sca(k, n)) ...
                        - mu_ref_sca(k, n) * mu_sca(k, n) ...
                        - phi_ref_sca(k, n) * phi_sca(k, n) ...
                        + 0.5 * (mu_ref_sca(k, n)^2 + phi_ref_sca(k, n)^2);

                    trace(H_s_k{k, n} * Ws(:, :, k, n)) + s_dc(k, n) >= nu;
                end
            end

            %% 功率约束 + 松弛
            for n = 1:N
                expression p
                p = trace(Wc(:, :, n));
                for k = 1:K
                    p = p + trace(Ws(:, :, k, n));
                end
                p <= p_max + s_pow(n);
                p >= 0.1;
            end

        cvx_end

        % 输出每次迭代的优化状态
        %disp('Iteration Status:');
        %disp(['Iteration ' num2str(iter)]);
        %disp(['Objective Value: ' num2str(cvx_optval)]);
        %disp('----------------------------');

        % 收敛检测
        mu_new = max(real(double(mu_sca)), 1e-6);
        phi_new = max(real(double(phi_sca)), 1e-6);

        % 稳定的参考点更新
        mu_ref_sca = mu_ref_sca * 0.5 + mu_new * 0.5;
        phi_ref_sca = phi_ref_sca * 0.5 + phi_new * 0.5;

           % 收敛检测（基于目标函数的变化）
        if abs(prev_obj_value - cvx_optval) < tol
        %    disp('Convergence achieved!');
            break;
        end

        % 更新目标函数值
        prev_obj_value = cvx_optval;
    end

    %% 输出结果
    w_comm = zeros(M, N);
    w_sens = zeros(M, K, N);
    power_vec_c = zeros(1, N);
    power_vec_s = zeros(1, N);

    %% 分解通信波束形成向量并判断秩  不要改动!!
    num_random_trials = 1000; % 随机采样次数
    for n = 1:N
        W = Wc(:, :, n);
        [V, D] = eig(W);
        [~, idx] = max(diag(D));
        w_comm(:, n) = sqrt(D(idx, idx)) * V(:, idx);
        
        r = rank(W, 1e-6);
        if r ~= 1
            is_communication_rank_one = false;
            [M, ~] = size(W);
            target_power = trace(W); % 目标功率
            
            % 多次随机采样选择最佳
            best_error = inf;
            for trial = 1:num_random_trials
                w_trial = (randn(M, 1) + 1i*randn(M, 1))/sqrt(2);
                w_trial = w_trial * sqrt(target_power)/norm(w_trial); % 功率归一化
                
                % 计算与W的匹配误差
                error = norm(w_trial*w_trial' - W, 'fro');
                if error < best_error
                    best_error = error;
                    w_comm(:, n) = w_trial;
                end
            end
        end
        power_vec_c(n) = shijianxishu*norm(w_comm(:, n))^2;
    end

    % 分解感知波束形成向量并判断秩
    for k = 1:K
        for n = 1:N
            W = Ws(:, :, k, n);
            [V, D] = eig(W);
            [~, idx] = max(diag(D));
            w_sens(:, k, n) = sqrt(D(idx, idx)) * V(:, idx);
            
            r = rank(W, 1e-6);
            if r ~= 1
                is_sensing_rank_one = false;
                [M, ~] = size(W);
                target_power = trace(W);
                
                best_error = inf;
                for trial = 1:num_random_trials
                    w_trial = (randn(M, 1) + 1i*randn(M, 1))/sqrt(2);
                    w_trial = w_trial * sqrt(target_power)/norm(w_trial);
                    
                    error = norm(w_trial*w_trial' - W, 'fro');
                    if error < best_error
                        best_error = error;
                        w_sens(:, k, n) = w_trial;
                    end
                end
            end
            power_vec_s1(k, n) = norm(w_sens(:, k, n))^2;
        end
    end
    for n = 1:N
       for k = 1:K
           power_vec_s(n)=shijianxishu*power_vec_s1(k, n)+power_vec_s(n);
       end
    end

end


%{
function [w_comm, w_sens, power_vec_c, power_vec_s] = ...
    optimize_precode_with_sensingCSI(N, K, s, jk, b, com_lb, sen_lb)
    % 参数设置
    com_lb = 13;
    sen_lb = 15;
    M = 4;  % 天线数
    p_max = 20;  % 最大功率
    N0 = 1e-14;  % 通信噪声
    N1 = 1e-14;  % 感知噪声
    acc = 0.5;  % 常数，控制干扰项权重

    com_lb_g = 2^(com_lb) - 1;  % 通信SINR门限gamma

    %% 信道构造
    h_c = cell(1, N);
    H_c = cell(1, N);
    for n = 1:N
        h_c{n} = generate_fixed_position_channel_c(M, s, b(:,n));
        H_c{n} = h_c{n} * h_c{n}';
    end

    h_s_k = cell(K, N);
    H_s_k = cell(K, N);
    for k = 1:K
        for n = 1:N
            h_s_k{k, n} = generate_fixed_position_channel_s(M, s, jk(:, k));
            H_s_k{k, n} = h_s_k{k, n} * h_s_k{k, n}';
        end
    end

    %% SCA迭代初始化
    MaxSCA = 80;  % 最大SCA迭代次数
    tol = 1e-3;  % 收敛容忍度

    mu_ref_sca = 1e-3 * ones(K, N);
    phi_ref_sca = 1e-3 * ones(K, N);

    % 松弛变量（用于约束的松弛）
    lambda_dc  = 10;  % DC松弛项权重，较小值
    lambda_sen = 10;  % 感知松弛项权重
    lambda_com = 10;  % 通信松弛项权重
    lambda_pow = 10;  % 功率松弛项权重

    for iter = 1:MaxSCA
        fprintf("SCA iteration %d\n", iter);

        mu_ref_sca = max(real(double(mu_ref_sca)), 1e-6);
        phi_ref_sca = max(real(double(phi_ref_sca)), 1e-6);

        % 进入CVX之前更新参考点
        cvx_begin sdp quiet
        cvx_solver mosek

        % CVX变量
        variable Wc(M, M, N) hermitian semidefinite  % 通信波束
        variable Ws(M, M, K, N) hermitian semidefinite  % 感知波束
        variable mu_sca(K, N)  % 感知信号的信噪比
        variable phi_sca(K, N)  % 干扰项的限制
        variable s_dc(K, N)  % DC松弛变量
        variable s_sen(K)  % 感知松弛变量
        variable s_com(N)  % 通信松弛变量
        variable s_pow(N)  % 功率松弛变量

        % 目标函数：最小化总功率 + 所有松弛变量的惩罚项
        expression total_power
        total_power = 0;
        for n = 1:N
            total_power = total_power + trace(Wc(:, :, n));
            for k = 1:K
                total_power = total_power + trace(Ws(:, :, k, n));
            end
        end

        minimize(total_power + lambda_dc * sum(s_dc(:)) + lambda_sen * sum(s_sen) + lambda_com * sum(s_com) + lambda_pow * sum(s_pow))

        subject to
            % 变量域限制
            mu_sca >= 0;
            phi_sca >= 0;
            s_com >= 0;
            s_pow >= 0;
            for k = 1:K
                s_dc(k, :) >= 0;
                s_sen(k) >= 0;
            end

            %% 通信SINR约束 + 松弛
            for n = 1:N
                expression interf
                interf = N0;
                for k = 1:K
                    interf = interf + acc * trace(H_c{n} * Ws(:, :, k, n));
                end
                acc * trace(H_c{n} * Wc(:, :, n)) + s_com(n) >= com_lb_g * interf;
            end

            %% 累计感知约束 + 松弛
            for k = 1:K
                sum(mu_sca(k, :)) + s_sen(k) >= sen_lb;
            end

            %% 干扰项 + 松弛
            for k = 1:K
                for n = 1:N
                    expression interf_s
                    interf_s = N1;
                    for i = 1:K
                        if i ~= k
                            interf_s = interf_s + acc * trace(H_s_k{k, n} * Ws(:, :, i, n));
                        end
                    end
                    interf_s = interf_s + acc * trace(H_s_k{k, n} * Wc(:, :, n));
                    interf_s <= phi_sca(k, n);
                end
            end

            %% DC约束 + 松弛
            for k = 1:K
                for n = 1:N
                    expression nu
                    nu = 0.5 * square_pos(mu_sca(k, n) + phi_sca(k, n)) ...
                        - mu_ref_sca(k, n) * mu_sca(k, n) ...
                        - phi_ref_sca(k, n) * phi_sca(k, n) ...
                        + 0.5 * (mu_ref_sca(k, n)^2 + phi_ref_sca(k, n)^2);

                    trace(H_s_k{k, n} * Ws(:, :, k, n)) + s_dc(k, n) >= nu;
                end
            end

            %% 功率约束 + 松弛
            for n = 1:N
                expression p
                p = trace(Wc(:, :, n));
                for k = 1:K
                    p = p + trace(Ws(:, :, k, n));
                end
                p <= p_max + s_pow(n);
            end

        cvx_end

        % 输出每次迭代的优化状态
        disp('Iteration Status:');
        disp(['Iteration ' num2str(iter)]);
        disp(['Objective Value: ' num2str(cvx_optval)]);
        disp('----------------------------');
        
        % 更新参考点
        mu_new = max(real(double(mu_sca)), 1e-6);
        phi_new = max(real(double(phi_sca)), 1e-6);

        % 收敛检测
        if max(abs(mu_new(:) - mu_ref_sca(:))) < tol && max(abs(phi_new(:) - phi_ref_sca(:))) < tol
            break;
        end

        mu_ref_sca = mu_new;
        phi_ref_sca = phi_new;
    end

    %% 输出结果
    w_comm = zeros(M, N);
    w_sens = zeros(M, K, N);
    power_vec_c = zeros(1, N);
    power_vec_s = zeros(1, N);

    num_random_trials = 500;

    % 通信波束
    for n = 1:N
        W = Wc(:, :, n);
        [V, D] = eig(W);
        [~, idx] = max(diag(D));
        w_comm(:, n) = sqrt(D(idx, idx)) * V(:, idx);
        power_vec_c(n) = norm(w_comm(:, n))^2;
    end
    
    % 感知波束
    for k = 1:K
        for n = 1:N
            W = Ws(:, :, k, n);
            [V, D] = eig(W);
            [~, idx] = max(diag(D));
            w_sens(:, k, n) = sqrt(D(idx, idx)) * V(:, idx);
            power_vec_s(n) = power_vec_s(n) + norm(w_sens(:, k, n))^2;
        end
    end
end
%}