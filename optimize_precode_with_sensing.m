function [w_comm, w_sens, power_vec_c,power_vec_s] = optimize_precode_with_sensing(N, K, s, jk, b,com_lb,sen_lb)
    %% 参数设置
    M = 4;
    p_max = 25;
    ZZ = 1e-14;
    N0 = ZZ; % 通信的噪声
    N1 = ZZ;    % 感知的噪声
    acc = 0.5;
    ruo_bc = 1;
    com_lb_g = 2^(com_lb/1) - 1;
    % 通信信道 h_c 和 H_c
    h_c = cell(1, N);
    H_c = cell(1, N);
    for n = 1:N
        h_c{n} = generate_fixed_position_channel_c(M,s,b(:, n));
        H_c{n} = h_c{n} * h_c{n}';
    end
    % 感知信道 h_s_k 和 H_s_k
    h_s_k = cell(K, N);
    H_s_k = cell(K, N);
    for k = 1:K
        for n = 1:N
            h_s_k{k, n} = generate_fixed_position_channel_s(M,s,jk(:, k));
            H_s_k{k, n} = h_s_k{k, n} * h_s_k{k, n}';
        end
    end
    %% SDR求解
    %power_s = zeros(N);
    cvx_precision high   
    cvx_begin sdp quiet
        variable  Wc_h(M, M, N) hermitian semidefinite
        variable  Ws_k_h(M, M, K, N) hermitian semidefinite
        variables Sk(K)                % 功率松弛变量
        variables Cs(N)                % 功率松弛变量

        variables epsilonc(N)
        variables power_s(N)
        variables epsilons

        variables epsilon
        variables epsilons1

        variable u_curv(N-2) nonnegative
        

        lambda_curv = 0;  % 可调权重（越大越平滑）

        % 目标函数
        objective = 0;
        
        for n = 1:N
            epsilonc(n) >=0;
            objective = objective + trace(Wc_h(:, :, n));
            for k = 1:K
                objective = objective + trace(Ws_k_h(:, :, k, n));
            end
        end

        objective = objective +  1e3 * sum(epsilonc);  % lambda 是惩罚权重
        %objective = objective + 0*sum(epsilonc)+ 0*(epsilons)+lambda_curv*sum(u_curv);  % lambda下的惩罚项 1e3  * sum(epsilonc)+ 1e5  * (epsilons)
        %objective = objective;
        minimize(objective)

        subject to
        for n = 1:N
            sinr_constraint_s = 0;
            sinr_constraint = trace(H_c{n} * Wc_h(:, :, n));
            for k = 1:K
                sinr_constraint_s = sinr_constraint_s + trace(H_c{n} * Ws_k_h(:, :, k, n));
            end
            %acc*sinr_constraint- acc*com_lb_g*sinr_constraint_s>= com_lb_g * N0;
            %Cs(n) == acc*sinr_constraint- acc*com_lb_g*sinr_constraint_s;
            acc*sinr_constraint>= com_lb_g * N0;
            Cs(n) == acc*sinr_constraint;
        end
        % 感知约束
        for k = 1:K
            interference_sum = 0;
            for n = 1:N
                interference_sum = interference_sum + trace(H_s_k{k, n} * Ws_k_h(:, :, k, n));
            end
            acc*interference_sum / N1>= sen_lb ;
            Sk(k) == acc*interference_sum;
        end
        % 功率限制
        for n = 1:N
            power_sum = trace(Wc_h(:, :, n));
            power_sums = 0;
            for k = 1:K
                power_sum = power_sum + trace(Ws_k_h(:, :, k, n));
                power_sums = power_sums + trace(Ws_k_h(:, :, k, n));
            end
            power_sum <=  p_max+epsilonc(n); 
            %power_sum  <=  p_max;
            epsilonc(n) >= 0;
            %0.5<=power_sums;
            power_s(n) == power_sums; 
            epsilons>=0;
        end

        %for n = 2:N-1
            % 二阶差分 p(n+1)-2p(n)+p(n-1)
        %    power_s(n+1) - 2*power_s(n) + power_s(n-1) <= u_curv(n-1);
        %    -(power_s(n+1) - 2*power_s(n) + power_s(n-1)) <= u_curv(n-1);
       % end
        delta = 0.1;  % 允许的相邻时隙最大功率变化（单位同 p_s）
        for n = 1:N-1
            -delta <= power_s(n+1) - power_s(n) <= delta;
        end

    cvx_end
    %% 输出功率向量
    power_vec_c = zeros(1, N);
    power_vec_s = zeros(1, N);
    %for n = 1:N
    %    power_vec_s(n) = 0;
    %    power_vec_c(n) = 0;
    %    power_vec_c(n) = trace(Wc_h(:, :, n));
    %    for k = 1:K
     %       power_vec_s(n) = power_vec_s(n) + trace(Ws_k_h(:, :, k, n))/ruo_bc;
     %   end
    %end
    
    % 分解通信波束形成向量
    %for n = 1:N
        % 获取协方差矩阵
    %    W = Wc_h(:, :, n);
        % 特征值分解
    %    [V, D] = eig(W);
        % 取最大特征值对应的特征向量
    %    [~, idx] = max(diag(D));
     %   w_comm(:, n) = sqrt(D(idx, idx)) * V(:, idx);
    %end
    % 分解感知波束形成向量
    %for k = 1:K
    %    for n = 1:N
            % 获取协方差矩阵
    %        W = Ws_k_h(:, :, k, n);
            % 特征值分解
    %        [V, D] = eig(W);
            % 取最大特征值对应的特征向量
    %        [~, idx] = max(diag(D));
     %       w_sens(:, k, n) = sqrt(D(idx, idx)) * V(:, idx);
     %   end
    %end

    % 分解通信波束形成向量并判断秩
    %for n = 1:N
        % 获取协方差矩阵
    %    W = Wc_h(:, :, n);
        % 特征值分解
    %    [V, D] = eig(W);
        % 取最大特征值对应的特征向量
    %    [~, idx] = max(diag(D));
    %    w_comm(:, n) = sqrt(D(idx, idx)) * V(:, idx);
        % 判断秩是否为1
    %    r = rank(W, 1e-6); % 可设置容差
    %    fprintf('Wc_h(:,:, %d) 的秩为 %d\n', n, r);
    %    if r == 1
    %        disp('通信矩阵秩为1 ✅');
    %    else
    %        disp('通信矩阵秩大于1 ⚠️');
    %    end
    %end
    
    % 分解感知波束形成向量并判断秩
    %for k = 1:K
    %    for n = 1:N
            % 获取协方差矩阵
    %        W = Ws_k_h(:, :, k, n);
            % 特征值分解
    %        [V, D] = eig(W);
            % 取最大特征值对应的特征向量
    %        [~, idx] = max(diag(D));
    %        w_sens(:, k, n) = sqrt(D(idx, idx)) * V(:, idx);
            % 判断秩是否为1
     %       r = rank(W, 1e-6);
     %       fprintf('Ws_k_h(:,:, %d, %d) 的秩为 %d\n', k, n, r);
    %        if r == 1
     %           disp('感知矩阵秩为1 ✅');
    %        else
    %            disp('感知矩阵秩大于1 ⚠️');
    %        end
     %   end
    %end


    % 10 28 
    %{
    % 初始化标志变量
    is_communication_rank_one = true;
    is_sensing_rank_one = true;
    
    % 分解通信波束形成向量并判断秩
    for n = 1:N
        % 获取协方差矩阵
        W = Wc_h(:, :, n);
        % 特征值分解
        [V, D] = eig(W);
        % 取最大特征值对应的特征向量
        [~, idx] = max(diag(D));
        w_comm(:, n) = sqrt(D(idx, idx)) * V(:, idx);
        % 判断秩是否为1
        r = rank(W, 1e-6); % 可设置容差
        if r ~= 1
            is_communication_rank_one = false;
            % 使用SVD作为替代提取方法
            [U, S, V_svd] = svd(W);
            w_comm(:, n) = U(:, 1) * S(1, 1); % 取最大的奇异值对应的左奇异向量
        end
    end
    
    % 分解感知波束形成向量并判断秩
    for k = 1:K
        for n = 1:N
            % 获取协方差矩阵
            W = Ws_k_h(:, :, k, n);
            % 特征值分解
            [V, D] = eig(W);
            % 取最大特征值对应的特征向量
            [~, idx] = max(diag(D));
            w_sens(:, k, n) = sqrt(D(idx, idx)) * V(:, idx);
            % 判断秩是否为1
            r = rank(W, 1e-6);
            if r ~= 1
                is_sensing_rank_one = false;
                % 使用SVD作为替代提取方法
                [U, S, V_svd] = svd(W);
                w_sens(:, k, n) = U(:, 1) * S(1, 1); % 取最大的奇异值对应的左奇异向量
            end
        end
    end
    
    % 最后输出结果
    if is_communication_rank_one
        %disp('所有通信矩阵的秩为1 ✅');
    else
        disp('存在通信矩阵的秩大于1 ⚠️,切换到SVD作为替代提取方法');
    end
    
    if is_sensing_rank_one
        %disp('所有感知矩阵的秩为1 ✅');
    else
        disp('存在感知矩阵的秩大于1 ⚠️，切换到SVD作为替代提取方法');
    end
    %}

    % 分解通信波束形成向量并判断秩
num_random_trials = 500; % 随机采样次数
for n = 1:N
    W = Wc_h(:, :, n);
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
    power_vec_c(n) = norm(w_comm(:, n))^2;
end

% 分解感知波束形成向量并判断秩
for k = 1:K
    for n = 1:N
        W = Ws_k_h(:, :, k, n);
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
       power_vec_s(n)=power_vec_s1(k, n)+power_vec_s(n);
   end
end

end
