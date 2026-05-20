function [s_uav, s_usv, total_energy, uav_energy, usv_energy] = optimize_trajectory_with_fly(D_max,Nh,H,N_f,Q,v_uav_max,v_usv_max,obs_pos,r1,delta,P_t,gamma_th,com_energy, s_uav, s_usv, start_uav, end_uav, start_usv, end_usv)
%function [s_uav, s_usv, total_energy, uav_energy, usv_energy] = optimize_trajectory_with_fly()
%D_max =D_max+100;
% ===== 参数设置 =====
E_uav_max = 3e4;
E_usv_max = 5e4;

num_obs = length(obs_pos);    
% 无人船参数 USV
rho_usv = 1000; % 水流密度
C_usv = 0.001; % 阻力因子
A_usv = 25; % 船的面积
alpha = 0.5*rho_usv*C_usv*A_usv; % USV速度惩罚系数
N0 = 1e-14 ;     %  下行链路噪声功率（图中标记为N₀^{DL}），通常设为1e-15
% ==== SCA主循环 ====
tol = 0.01;        
prev_obj = inf;
total_energy = 0; % 初始化总能耗
water_speed_max = 0; % 最大水流速度
% 水流速度函数
flow_direction = 1;
U0   = flow_direction * 0.8 * water_speed_max;   % 主流强度
A    = 40.6 * water_speed_max;                   % 波动幅值
kx   = 0.06;                                     % x方向波数（越大波越密）
ky   = 0.03;                                     % y方向波数
phi  = 0;                                        % 相位（可改成随时间变化的phi = omega*t 实现行波）
% 定义水流速度函数，可用于动态仿真
water_vel_func = @(x,y) [
  (  U0 - A * ky * sin(kx * x + phi) .* sin(ky * y) );   % X方向速度
   ( -A * kx * cos(kx * x + phi) .* cos(ky * y) )        % Y方向速度
];
% 预计算水流速度
water_vel_usv = zeros(2, N_f);
for n = 1:N_f
    water_vel = water_vel_func(s_usv(1,n), s_usv(2,n));
    water_vel_usv(1,n) = water_vel(1);
    water_vel_usv(2,n) = water_vel(2);
end
% 自适应松弛参数
slack_weight = 100;  % 初始松弛惩罚权重
max_slack_weight = 1e6; % 最大松弛惩罚权重
min_slack_weight = 1;  % 最小松弛惩罚权重
converged = false;
%% 迭代开始
for iter = 1:10
    fprintf('迭代 %d... ', iter);
    % 计算当前梯度和距离 处理障碍物
    grad_usv_list = zeros(3, N_f, num_obs);
    dist_usv_list = zeros(N_f, num_obs);
    for obs = 1:num_obs
        for n = 1:N_f   
            % USV梯度计算 (仅XY平面)
            usv_to_obs = s_usv(1:2,n) - obs_pos(1:2,obs);
            dist_usv = norm(usv_to_obs);
            dist_usv_list(n, obs) = dist_usv;
            if dist_usv > 1e-3
                grad_usv_list(:, n, obs) = [usv_to_obs/dist_usv; 0];
            else
                grad_usv_list(:, n, obs) = [1; 0; 0];
            end
        end
    end
   %MRT 参数
   Kc = 0.00001;
   gc = 0.00001;
   B = 1;
   rho = 0.03; % 参考距离信道增益
   rho0 = rho*( ( (Kc/(Kc+1))^0.5+gc*(1/(Kc+1))^0.5 )^0.5  )^2;

    for i = 1:N_f
        rho0 = rho*(  (Kc/(Kc+1))^0.5  +  gc*(1/(Kc+1))^0.5  )^2;
        numerator = P_t(i) * Q * rho0;
        denominator =(2^(gamma_th/1)-1)*N0;

        %numerator = P_t(i) * Q * rho0;
        %denominator =0.1*(2^(gamma_th/B)-1)*N0;
        % 2. 应用推导公式: d^4 <= P_t Q rho0^2 / (gamma_th N0)
        d_3d_max(i) = real((numerator / denominator)^(1/2));
    end
    %disp(d_3d_max);

    %% CVX
    d_0 = 0.6;
    P0 = 80;
    P1 = 88.6;
    v0 =4.03;
    Acanshu = 0.503;
    s_uavcanshu = 0.05;
    rou = 1.225;
    oumu = 300;
    r_uav = 0.4;
    % 处理无人机能耗 
    lamda_B = zeros(1,N_f);
    sudu_B = zeros(1,N_f);
    for n = 2:N_f
        v_uav = norm(s_uav(:,n) - s_uav(:,n-1)) / delta;
        sudu_B(n) = v_uav;
        lamda_B(n) = (sqrt(1 + (v_uav^4)/(4*v0^4)) - 0.5*(v_uav^2)/(v0^2))^(1/2);
    end
    cvx_precision high
    cvx_begin quiet
    cvx_solver mosek
        variables s_uav_new(3, N_f) s_usv_new(3, N_f)
        variable lambda_B1(N_f) % 松弛变量 λ_B[n]
        variable distance_songchi(N_f-2) % 松弛变量 λ_B[n]
        variable phi_diff(N_f-1)  % 角度差（作为辅助变量）
        variables vrel_mag(N_f-1)   % 每段相对速度的模长（>=0）
        expression obj_uav
        expression obj_usv
        obj_uav = 0;
        obj_usv = 0;
        obj = 0;
        obj1  = 0;
        % UAV能耗模型
        for n = 2:N_f
            v_uav = norm(s_uav_new(:,n) - s_uav_new(:,n-1)) / delta;
            thm1 = P0*(1+(3/(oumu*r_uav)^2)*pow_pos(v_uav,2 ));
            thm2 = 0.5*d_0*rou*s_uavcanshu*Acanshu*pow_pos(v_uav,3);
            obj = obj + thm1 + thm2+P1*lambda_B1(n);
            obj_uav = obj_uav + thm1 + thm2 + P1*lambda_B1(n);
        end
        % USV运动惩罚 (矢量计算)
        %for n = 2:N_f
            % 仅处理XY平面
        %    usv_displacement = s_usv_new(1:2,n) - s_usv_new(1:2,n-1);
        %    v_usv_ground = usv_displacement / delta;
        %    v_water = water_vel_usv(1:2, n-1);
         %   v_relative = v_usv_ground - v_water;
        %    v_rel_mag = norm(v_relative);
        %    obj = obj + alpha * square_pos(v_rel_mag);
        %end

        for n = 2:N_f
            usv_displacement = s_usv_new(1:2,n) - s_usv_new(1:2,n-1);
            v_usv_ground = usv_displacement / delta;
            v_water = water_vel_usv(1:2, n-1);
            v_relative = v_usv_ground - v_water;
        
            vrel_mag(n-1) >= norm(v_relative);   % SOC 约束
            vrel_mag(n-1) >= 0;
        
            obj = obj + alpha * pow_pos(vrel_mag(n-1), 3);   % 三次方
            obj_usv = obj_usv + alpha * pow_pos(vrel_mag(n-1), 3);   % 三次方
        end

        obj1 = 10000*obj_uav/E_uav_max+  10000*obj_usv/E_usv_max+ slack_weight * sum(distance_songchi);

        obj = obj + com_energy;
        minimize(obj1 )
   %% 约束处理
        subject to
        % 边界条件
        s_uav_new(:,1) == start_uav;
        s_uav_new(:,end) == end_uav;
        s_usv_new(:,1) == start_usv;
        s_usv_new(:,end) == end_usv;
        for n = 1:N_f
           %USV高度限制
           s_uav_new(3, n) ==100;
           %s_uav_new(1, n) <=300;
           %0<=s_uav_new(1, n); 
           %s_uav_new(2, n) <=300;
           %0<=s_uav_new(2, n); 


           s_usv_new(3, n) ==0;
           %s_usv_new(1, n) <=300;
           %0<=s_usv_new(1, n); 
           %s_usv_new(2, n) <=300;
           %0<=s_usv_new(2, n); 
        end
        % 避障约束
        for obs = 1:num_obs
            for n = 1:N_f                
                if dist_usv_list(n, obs) > 1e-3
                    dist_usv_list(n, obs) + grad_usv_list(1:2, n, obs)' * (s_usv_new(1:2,n) - s_usv(1:2,n)) >= r1;
                end
            end
        end
        % 速度约束
        for n = 1:N_f-1
            norm(s_uav_new(:,n+1) - s_uav_new(:,n)) <= v_uav_max * delta;
            norm(s_usv_new(:,n+1) - s_usv_new(:,n)) <= v_usv_max * delta;
        end

        % 通信约束
        for n = 2:N_f-1
             % 位置坐标
             uav_pos = s_uav_new(:,n);
             usv_pos = s_usv_new(:,n);
             % 距离约束 (二阶锥形式)
             norm(uav_pos - usv_pos) <= d_3d_max(n)+distance_songchi(n-1);
             distance_songchi(n-1)>=0;
             distance_songchi(n-1) <= 1000; % 限制松弛变量范围
        end

        % UAV能耗模型
        for n = 2:N_f
            v_uav = (s_uav_new(:,n) - s_uav_new(:,n-1)) / delta;
            v_uav = real(v_uav);
            pow_pos(inv_pos(  (lambda_B1(n))),2 )<= ...
            2*lamda_B(n)*(  (lambda_B1(n))- lamda_B(n) ) + ...
            (lamda_B(n))^2 + ...
            (sudu_B(n)/v0)^2 + ...
            2/v0^2 * (  (s_uav(:,n) - s_uav(:,n-1)) / delta   )'*(   (s_uav_new(:,n) - s_uav_new(:,n-1)) / delta - (s_uav(:,n) - s_uav(:,n-1)) / delta    );
        end
    cvx_end
    % 更新轨迹
    if strcmp(cvx_status, 'Solved')
        s_uav = s_uav_new;
        s_usv = s_usv_new;
        fprintf('目标值 = %.2f\n', cvx_optval);
        val_slack_sen = sum(double(distance_songchi));
        fprintf('惩罚值 = %.2f\n', val_slack_sen);
        % 计算总能耗
        %total_energy = cvx_optval;
 

        % 更新水流速度计算
        for n = 1:N_f
            water_vel = water_vel_func(s_usv(1,n), s_usv(2,n));
            water_vel_usv(1,n) = water_vel(1);
            water_vel_usv(2,n) = water_vel(2);
        end

        % 检查收敛
        if  abs(cvx_optval - prev_obj)/prev_obj < tol
            fprintf('===== 收敛于%d次迭代 =====\n', iter);
            converged = true;
            break;
        end
        prev_obj = cvx_optval;
            % 成功时逐步收紧松弛惩罚
        slack_weight = min(slack_weight * 2, max_slack_weight);
    else
        % 失败时放松约束并继续迭代
        slack_weight = max(slack_weight / 5, min_slack_weight);
        %fprintf('松弛权重调整为%.1f\n',  slack_weight);
    end
end
         % 如果迭代结束但未收敛，计算当前轨迹的能耗
    %if ~converged
        % 计算能耗 
        %total_energy = 0;
        uav_energy = 0;
        usv_energy = 0;
        % 计算UAV能耗
        for n = 2:N_f
            v_uav = norm(s_uav(:,n) - s_uav(:,n-1)) / delta;
            lamba_B(n) = (sqrt(1 + (v_uav^4)/(v0^4)) - 0.5*(v_uav^2)/(v0^2))^(1/2);
            thm1 = P0*(1 + (3/(oumu*r_uav)^2)*v_uav^2);
            thm2 = 0.5*d_0*rou*s_uavcanshu*Acanshu*v_uav^3;
            uav_energy = uav_energy + thm1 + thm2 + P1*lamba_B(n);
        end
        
        % 计算USV能耗
        for n = 2:N_f
            usv_displacement = s_usv(1:2,n) - s_usv(1:2,n-1);
            v_usv_ground = usv_displacement / delta;
        % 注意：这里必须用最终轨迹位置重新算水流
            v_water = water_vel_func(s_usv(1,n-1), s_usv(2,n-1));
            v_relative = v_usv_ground - v_water;
            v_rel_mag = norm(v_relative);
            usv_energy = usv_energy + alpha * (v_rel_mag)^3;
        end
        
        total_energy = usv_energy+uav_energy + com_energy;
        %total_energy = 1000;
        %fprintf('求解Over: %s\n', cvx_status);
        %fprintf('Orignal\n');
    %end
end