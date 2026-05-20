function [w_comm, w_sens, power_vec_c, power_vec_s] = ...
    optimize_precode_with_sensingCSI1(N, K, s, jk, b, com_lb, sen_lb)
    
    %% 参数设置
    M = 4;
    p_max = 20;
    N0 = 1e-14;
    N1 = 1e-14;
    acc = 0.5;
    
    com_lb_g = 2^(com_lb) - 1;
    
    fprintf('\n======= 极端信道条件优化 =======\n');
    fprintf('N=%d, K=%d, M=%d, p_max=%.1f\n', N, K, M, p_max);
    fprintf('com_lb=%.2f (SINR=%.2f), sen_lb=%.2f\n', com_lb, com_lb_g, sen_lb);
    
    %% 生成信道矩阵
    H_c = cell(1, N);
    H_s_k = cell(K, N);
    
    for n = 1:N
        h_c = generate_fixed_position_channel_c(M, s, b(:, n));
        H_c{n} = h_c * h_c';
    end
    
    for k = 1:K
        for n = 1:N
            h_s = generate_fixed_position_channel_s(M, s, jk(:, k));
            H_s_k{k, n} = h_s * h_s';
        end
    end
    
    %% 诊断信道条件
    Hc_norms = cellfun(@(H) norm(H, 'fro'), H_c);
    Hs_norms = cellfun(@(H) norm(H, 'fro'), H_s_k);
    
    fprintf('信道范数诊断:\n');
    fprintf('  H_c: 平均值=%.2e, 最小值=%.2e, 最大值=%.2e\n', ...
        mean(Hc_norms), min(Hc_norms), max(Hc_norms));
    fprintf('  H_s: 平均值=%.2e, 最小值=%.2e, 最大值=%.2e\n', ...
        mean(Hs_norms(:)), min(Hs_norms(:)), max(Hs_norms(:)));
    
    %% 分析可行性
    fprintf('\n可行性分析:\n');
    
    % 估计最佳可能的SINR
    best_sinr_est = p_max * max(Hc_norms) / N0;
    fprintf('  理论最大SINR: %.2f (需求: %.2f)\n', best_sinr_est, com_lb_g);
    
    if best_sinr_est < com_lb_g
        fprintf('  ⚠ 通信需求可能过高，降低目标值...\n');
        % 自动调整通信需求
        com_lb_g = min(com_lb_g, best_sinr_est * 0.9);
        com_lb = log2(1 + com_lb_g);
        fprintf('  调整后: com_lb=%.2f (SINR=%.2f)\n', com_lb, com_lb_g);
    end
    
    % 估计感知可行性
    best_snr_est = p_max * max(Hs_norms(:)) / N1;
    fprintf('  理论最大感知SNR: %.2f (需求: %.2f)\n', best_snr_est, sen_lb);
    
    if best_snr_est * N < sen_lb
        fprintf('  ⚠ 感知需求可能过高，降低目标值...\n');
        sen_lb = min(sen_lb, best_snr_est * N * 0.9);
        fprintf('  调整后: sen_lb=%.2f\n', sen_lb);
    end
    
    %% 方案1: 尝试简化的可行性问题
    fprintf('\n尝试简化问题...\n');
    
    % 先只求解通信问题
    [Wc_feas, comm_feas] = solve_comm_only(N, M, p_max, H_c, com_lb_g, N0, acc);
    
    if ~comm_feas
        fprintf('  ✗ 通信问题不可行\n');
        % 返回零解
        [w_comm, w_sens, power_vec_c, power_vec_s] = get_zero_solution(M, N, K);
        return;
    end
    
    fprintf('  ✓ 通信问题可行\n');
    
    %% 方案2: 逐步构建可行解
    fprintf('\n逐步构建可行解...\n');
    
    % 初始化
    mu_ref = ones(K, N) * 0.1;
    phi_ref = ones(K, N) * 1.0;
    
    MaxSCA = 30;
    tol = 1e-4;
    
    for iter = 1:MaxSCA
        fprintf('  SCA迭代 %d/%d: ', iter, MaxSCA);
        
        % 求解子问题
        [success, Wc, Ws, mu, phi, obj_val] = solve_extreme_subproblem(...
            N, K, M, p_max, H_c, H_s_k, com_lb_g, sen_lb, ...
            N0, N1, acc, mu_ref, phi_ref);
        
        if ~success
            fprintf('失败，尝试重启\n');
            % 重启机制
            mu_ref = mu_ref * 0.5;
            phi_ref = phi_ref * 2.0;
            continue;
        end
        
        fprintf('目标值=%.6f\n', obj_val);
        
        % 检查解的质量
        total_power = 0;
        for n = 1:N
            total_power = total_power + trace(Wc(:,:,n));
            for k = 1:K
                total_power = total_power + trace(Ws(:,:,k,n));
            end
        end
        
        fprintf('    总功率=%.6f, 平均mu=%.6f\n', total_power, mean(mu(:)));
        
        % 更新参考点
        mu_ref = 0.3 * mu_ref + 0.7 * mu;
        phi_ref = 0.3 * phi_ref + 0.7 * phi;
        
        % 检查收敛
        if iter > 1
            mu_change = norm(mu(:) - mu_prev(:)) / norm(mu_prev(:) + eps);
            if mu_change < tol
                fprintf('    ✓ 收敛\n');
                break;
            end
        end
        
        mu_prev = mu;
        
        if iter == MaxSCA
            fprintf('    ⚠ 达到最大迭代\n');
        end
    end
    
    %% 提取波束形成向量
    fprintf('\n提取波束形成向量...\n');
    
    [w_comm, w_sens, power_vec_c, power_vec_s] = extract_solution(Wc, Ws, M, N, K);
    
    %% 验证解
    fprintf('\n验证解:\n');
    
    [comm_satisfied, sen_satisfied] = verify_solution_extreme(...
        w_comm, w_sens, H_c, H_s_k, com_lb_g, sen_lb, N0, N1, acc, p_max);
    
    if comm_satisfied && sen_satisfied
        fprintf('  ✓ 解满足所有约束\n');
    else
        if ~comm_satisfied
            fprintf('  ⚠ 通信约束可能不满足\n');
        end
        if ~sen_satisfied
            fprintf('  ⚠ 感知约束可能不满足\n');
        end
    end
    
    fprintf('总通信功率: %.6f W\n', sum(power_vec_c));
    fprintf('总感知功率: %.6f W\n', sum(power_vec_s));
end

%% 辅助函数1: 只求解通信问题
function [Wc, feasible] = solve_comm_only(N, M, p_max, H_c, com_lb_g, N0, acc)
    feasible = false;
    Wc = [];
    
    try
        cvx_begin sdp quiet
        cvx_solver mosek
        
        variable Wc_var(M,M,N) hermitian semidefinite
        
        minimize(0)  % 可行性问题
        
        subject to
        for n = 1:N
            trace(H_c{n} * Wc_var(:,:,n)) >= com_lb_g * N0;
            trace(Wc_var(:,:,n)) <= p_max;
        end
        
        cvx_end
        
        if strcmp(cvx_status, 'Solved')
            feasible = true;
            Wc = Wc_var;
        end
        
    catch
        feasible = false;
    end
end

%% 辅助函数2: 极端条件下的子问题求解
function [success, Wc, Ws, mu, phi, obj_val] = solve_extreme_subproblem(...
    N, K, M, p_max, H_c, H_s_k, com_lb_g, sen_lb, N0, N1, acc, mu_ref, phi_ref)
    
    success = false;
    Wc = [];
    Ws = [];
    mu = [];
    phi = [];
    obj_val = Inf;
    
    try
        cvx_begin sdp quiet
        cvx_solver mosek
        cvx_precision high
        
        % 变量
        variable Wc_var(M,M,N) hermitian semidefinite
        variable Ws_var(M,M,K,N) hermitian semidefinite
        variable mu_var(K,N)
        variable phi_var(K,N)
        
        % 松弛变量
        variable slack_comm(N) nonnegative
        variable slack_sen(K) nonnegative
        variable slack_dc(K,N) nonnegative
        
        % 目标函数
        expression total_power
        total_power = 0;
        for n = 1:N
            total_power = total_power + trace(Wc_var(:,:,n));
            for k = 1:K
                total_power = total_power + trace(Ws_var(:,:,k,n));
            end
        end
        
        % 添加正则化防止退化解
        reg_term = 1e-6 * (sum_square_abs(vec(Wc_var)) + sum_square_abs(vec(Ws_var)));
        
        minimize(total_power + reg_term + ...
            100 * sum(slack_comm) + 100 * sum(slack_sen) + 10 * sum(slack_dc(:)))
        
        subject to
        % 变量范围
        mu_var >= 1e-6;
        phi_var >= 1e-6;
        
        %% 通信约束（宽松版本）
        for n = 1:N
            expression interf
            interf = N0;
            for k = 1:K
                interf = interf + 0.1 * trace(H_c{n} * Ws_var(:,:,k,n));
            end
            
            trace(H_c{n} * Wc_var(:,:,n)) + slack_comm(n) >= com_lb_g * interf;
        end
        
        %% 感知约束（宽松版本）
        for k = 1:K
            sum(mu_var(k,:)) + slack_sen(k) >= sen_lb;
        end
        
        %% 干扰约束
        for k = 1:K
            for n = 1:N
                expression interf_s
                interf_s = N1;
                
                for i = 1:K
                    if i ~= k
                        interf_s = interf_s + 0.1 * trace(H_s_k{k,n} * Ws_var(:,:,i,n));
                    end
                end
                
                interf_s = interf_s + 0.1 * trace(H_s_k{k,n} * Wc_var(:,:,n));
                interf_s <= phi_var(k,n) + 1e-6;
            end
        end
        
        %% 简化的DC约束
        for k = 1:K
            for n = 1:N
                % 线性近似
                nu_approx = mu_ref(k,n) * phi_var(k,n) + phi_ref(k,n) * mu_var(k,n) ...
                    - mu_ref(k,n) * phi_ref(k,n);
                
                trace(H_s_k{k,n} * Ws_var(:,:,k,n)) + slack_dc(k,n) >= nu_approx;
            end
        end
        
        %% 功率约束
        for n = 1:N
            expression p
            p = trace(Wc_var(:,:,n));
            for k = 1:K
                p = p + trace(Ws_var(:,:,k,n));
            end
            p <= p_max;
            
            % 确保最小功率
            p >= 1e-3;
        end
        
        %% 额外稳定性约束
        for n = 1:N
            trace(Wc_var(:,:,n)) >= 1e-3;
            for k = 1:K
                trace(Ws_var(:,:,k,n)) >= 1e-3;
            end
        end
        
        cvx_end
        
        if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
            success = true;
            Wc = Wc_var;
            Ws = Ws_var;
            mu = mu_var;
            phi = phi_var;
            obj_val = cvx_optval;
        end
        
    catch ME
        fprintf('求解异常: %s\n', ME.message);
    end
end

%% 辅助函数3: 提取解
function [w_comm, w_sens, power_vec_c, power_vec_s] = extract_solution(Wc, Ws, M, N, K)
    w_comm = zeros(M, N);
    w_sens = zeros(M, K, N);
    power_vec_c = zeros(1, N);
    power_vec_s = zeros(1, N);
    
    % 通信波束
    for n = 1:N
        W = double(Wc(:,:,n));
        [V, D] = eig(W);
        [~, idx] = max(diag(D));
        w_comm(:, n) = sqrt(abs(D(idx, idx))) * V(:, idx);
        power_vec_c(n) = norm(w_comm(:, n))^2;
    end
    
    % 感知波束
    for k = 1:K
        for n = 1:N
            W = double(Ws(:,:,k,n));
            [V, D] = eig(W);
            [~, idx] = max(diag(D));
            w_sens(:, k, n) = sqrt(abs(D(idx, idx))) * V(:, idx);
            power_vec_s(n) = power_vec_s(n) + norm(w_sens(:, k, n))^2;
        end
    end
end

%% 辅助函数4: 验证解
function [comm_ok, sen_ok] = verify_solution_extreme(...
    w_comm, w_sens, H_c, H_s_k, com_lb_g, sen_lb, N0, N1, acc, p_max)
    
    [K, N] = size(H_s_k);
    
    comm_ok = true;
    sen_ok = true;
    
    tol = 1e-2;
    
    % 验证通信约束
    for n = 1:N
        interf = N0;
        for k = 1:K
            w_s = w_sens(:, k, n);
            interf = interf + acc * real(w_s' * H_c{n} * w_s);
        end
        
        w_c = w_comm(:, n);
        signal = acc * real(w_c' * H_c{n} * w_c);
        
        sinr = signal / (interf + eps);
        
        if sinr < com_lb_g - tol
            fprintf('    时隙%d: SINR=%.2f < %.2f\n', n, sinr, com_lb_g);
            comm_ok = false;
        end
    end
    
    % 验证感知约束
    for k = 1:K
        total_snr = 0;
        
        for n = 1:N
            w_s = w_sens(:, k, n);
            signal = real(w_s' * H_s_k{k,n} * w_s);
            
            interf = N1;
            
            for i = 1:K
                if i ~= k
                    w_i = w_sens(:, i, n);
                    interf = interf + acc * real(w_i' * H_s_k{k,n} * w_i);
                end
            end
            
            w_c = w_comm(:, n);
            interf = interf + acc * real(w_c' * H_s_k{k,n} * w_c);
            
            if interf > 0
                snr = signal / interf;
                total_snr = total_snr + snr;
            end
        end
        
        if total_snr < sen_lb - tol
            fprintf('    目标%d: 总SNR=%.2f < %.2f\n', k, total_snr, sen_lb);
            sen_ok = false;
        end
    end
end

%% 辅助函数5: 返回零解
function [w_comm, w_sens, power_vec_c, power_vec_s] = get_zero_solution(M, N, K)
    w_comm = zeros(M, N);
    w_sens = zeros(M, K, N);
    power_vec_c = zeros(1, N);
    power_vec_s = zeros(1, N);
end