%% 2025.10.14 代码重构 完整版
%% 2026.02.28 重构-
clear; clc; close all;
cvx_solver mosek
%% 0.文件存储名称
tic;  % 开始计时
clust_xieru = 'clusterceshi.txt'; 
hover_xietu =  'uav_sKceshi.txt';
uav_xieru = 'uav_pathKceshi.txt';
usv_xieru = 'usv_pathKceshi.txt';
pc_xieru =       'Pc_Kceshi.txt';
ps_xieru =       'Ps_Kceshi.txt'; 
%% 1. 共用参数
uav_start =   [0, 0,100];       % UAV起始点
uav_end   =   [300, 0,100];   % UAV终止点
usv_start =   [0, 0,0];         % USV起始点
usv_end =     [300, 0,0];         % USV终止点
Nh = 1;               % 预备参数，忽略
delta = 1;                    % 时隙间隔
r1 = 0;                       % 和障碍物的安全距离
obs_pos = [
            65, 35, 0;    % 地面障碍物2
           120, 190, 0;    % 高空障碍物
           210, 100, 0;   % 地面障碍物1
           ]';   
E_uav_max = 3e4;
E_usv_max = 3e4; 
%% 2.确定目标分布 
num_points = 40;  % 目标数量
area_size = 300;
num_cols = 8;
num_rows = ceil(num_points / num_cols);
rng(42);
GT_positions = zeros(num_points, 2);
count = 0;
for i = 1:num_rows
    for j = 1:num_cols
        if count >= num_points, break; end
        count = count + 1;
        cell_width = area_size / num_cols;
        cell_height = area_size / num_rows;
        x0 = (j-1) * cell_width;
        y0 = (i-1) * cell_height;
        x = x0 + 0.5 * cell_width + rand() * 0.1 * cell_width;
        y = y0 + 0.2* cell_height + rand() * 0.6 * cell_height;
        %x = x0 + 0.5 * cell_width ;
        %y = y0 + 0.5* cell_height ;
        %x = x0 + 2.5 * cell_width ;
        %y = y0 + 0* cell_height ;
        GT_positions(count, :) = [x, y];
    end
end
std_deviation_target =40; % 目标标准差   10 - 50
mean_point = mean(GT_positions, 1);
% 计算当前点集的标准差
distances = sqrt(sum((GT_positions - mean_point).^2, 2));
current_std_deviation = std(distances);
disp(['Current Standard Deviation: ', num2str(current_std_deviation)]);
% 如果当前标准差与目标标准差不符，进行调整
adjust_factor = std_deviation_target / current_std_deviation;
% 调整点的分布，使其标准差接近目标标准差
%GT_positions_adjusted = (GT_positions - mean_point) * adjust_factor + mean_point;
% 计算调整后的标准差
%distances_adjusted = sqrt(sum((GT_positions_adjusted - mean_point).^2, 2));
%adjusted_std_deviation = std(distances_adjusted);
%disp(['Adjusted Standard Deviation: ', num2str(adjusted_std_deviation)]);
%GT_positions = GT_positions_adjusted;
%disp(GT_positions);

% 调整标准差
GT_positions_adjusted = (GT_positions - mean_point) * adjust_factor + mean_point;

% ---- 保证所有 y 坐标大于 0 ----
min_y = min(GT_positions_adjusted(:,2));
if min_y < 0
    % 平移所有点，使最小 y = 0.1（或其他安全边界）
    GT_positions_adjusted(:,2) = GT_positions_adjusted(:,2) - min_y + 0.1;
end

% 重新计算标准差
distances_adjusted = sqrt(sum((GT_positions_adjusted - mean_point).^2, 2));
adjusted_std_deviation = std(distances_adjusted);
disp(['Adjusted Standard Deviation: ', num2str(adjusted_std_deviation)]);
GT_positions = GT_positions_adjusted;
disp(GT_positions);

%% 3.悬停点优化参数
SNR_k   =  10;        % 感知累计功率 
SNR_one = 3;        % 感知瞬时功率 
gamma_value = 13;  % 飞行模式下通信速率要求
M = 4;               % 天线数量
Q = M;               % 复用参数，无用
P =  2;              % 预设感知功率
PC = 2;              % 预设感知功率
H = 100;             % 无人机高度
sensing_nums = 8;    % 设置每个悬停点的感知目标上限
N0_UL = 1e-14;       % 感知干扰功率 -110dBm
N0 = 1e-14;          % 通信干扰功率 -110dBm
acc = 0.5;           % 脉冲时间  tp/(tp+to)
neta = 0.1;          % 对目标的  RCS m^2 
beta = 0.03;         % 1m时候的感知信道增益 14.8dBm
Kc = 0.00001;        % 通信信道参数 选值来自参考文献[23]
gc = 0.00001;        % 通信信道参数 选值来自参考文献[23]
beta_0 = 0.03;       % 通信信道增益
T_K= 10;             % AO迭代上限 如果收敛 提前跳出 
threshold = 1e-3;    % 设定一个阈值来判断迭代结果的变化程度
% 通信信道距离转换 
rho0 = beta_0*(  (Kc/(Kc+1))^0.5  +  gc*(1/(Kc+1))^0.5  )^2;
numerator = PC * M * rho0;
denominator =(2^(gamma_value/1)-1)*N0;
% 2. 应用推导公式: d^2 <= P_t Q rho0^2 / (gamma_th N0)
fly_com_d = real((numerator / denominator)^(1/2));
if (fly_com_d^2-H^2)<0
    disp('通信距离少于下限');
end
fly_com_d = sqrt(fly_com_d^2-H^2);
disp('在悬停点优化阶段的通信距离:');
disp(fly_com_d);

% 运动参数
Vmax = 20;                       % UAV 最大速度 (m/s)
Vmax_usv = 10;                   % USV 最大速度 (m/s) 5-10
v_uav_max = Vmax;
v_usv_max = Vmax_usv;
GT_positions11 = GT_positions ;  % 后面画图用
% 定义起点和终点
s_start = uav_start(1:2);        % 起始点
f_end = uav_end(1:2);            % 终止点
com_lb = gamma_value;            % 悬停模式下通信速率要求
sen_lb = SNR_k;                  % 悬停模式下感知信噪比要求
%% 4.虚拟基站放置 
num_points = size(GT_positions, 1);
% 找最小满足距离限制的聚类数量
for k = 1:num_points
    [idx, C] = kmeans(GT_positions, k, 'Replicates', 100);
    valid = true;
    for i = 1:k
        cluster_points = GT_positions(idx == i, :);
        if isempty(cluster_points)
            continue;
        end
        % 如果超过 就调整 
        % 如果某个簇的目标点数量超过sensing_nums，进行调整
        if size(cluster_points, 1) > sensing_nums
            % 调整聚类，移除多余的点并将其分配到其他聚类
            idx = redistribute_points(GT_positions, idx, C, i, sensing_nums);
        end
        % 还超过就直接+1 
        if size(cluster_points, 1) > sensing_nums
            valid = false;
            break;
        end
        % (P/xgz) 均分感知功率 
        xgz = size(cluster_points, 1);
        Upsilon1 = (P/xgz)*M*acc*(neta*beta^2)*((16*pi)^(-1))*(N0_UL^(-1));
        coverage_radius = (( Upsilon1/(SNR_one))^(0.5) -H^2  )^(0.5); % UAV 覆盖半径 D
        max_dist_threshold = coverage_radius;
        distances = sqrt(sum((cluster_points - C(i, :)).^2, 2));
        if any(distances > coverage_radius)
            valid = false;
            break;
        end
    end
    if valid
        break;
    end
end
% 构造 clusters 和 VBS_positions
clusters = cell(k, 1);
VBS_positions = zeros(k, 2);
for i = 1:k
    clusters{i} = find(idx == i);           % 索引集合
    VBS_positions(i, :) = C(i, :);          % 簇中心
end
writematrix(VBS_positions, 'VBS_positions.txt', 'Delimiter', 'tab');
disp(VBS_positions);
G = size(VBS_positions, 1); % 虚拟基站数量
%% 5.TSP对虚拟基站排序
num_VBS = G;
% 计算距离矩阵（包括起点和终点）
all_points = [s_start; VBS_positions; f_end]; % 起点 + 所有VBS + 终点
num_all = size(all_points, 1);
dist_matrix = zeros(num_all);
for i = 1:num_all
    for j = 1:num_all
        dist_matrix(i,j) = norm(all_points(i,:) - all_points(j,:));
    end
end
% 贪心法路径排序（考虑起点和终点）
unvisited = 2:num_all-1;  % 跳过起点(1)和终点(end)
path = 1;                 % 从起点开始
current_node = 1;
while ~isempty(unvisited)
    dists = dist_matrix(current_node, unvisited);
    [~, idx] = min(dists);
    next_node = unvisited(idx);
    path = [path, next_node];
    current_node = next_node;
    unvisited(unvisited == next_node) = [];
end
path = [path, num_all];   % 添加终点
% 提取排序后的所有点顺序
sorted_all_points = all_points(path, :);
disp(length(sorted_all_points));
% 区分不同类型的位置
sorted_s_start = sorted_all_points(1, :);                 % 起点
sorted_f_end = sorted_all_points(end, :);                 % 终点
% 用 Gurobi 求解
sorted_all_points = TSP(sorted_all_points);
disp( sorted_all_points );
sorted_VBS_positions = sorted_all_points(2:end-1, :);     % 中间VBS位置
% 调整clusters顺序
[~, idx] = ismember(sorted_VBS_positions, VBS_positions, 'rows');
clusters = clusters(idx);
disp("每个虚拟基站的聚簇个数:");
disp(clusters);
% 运行算法
% 创建已覆盖GT的全局记录
covered_GTs = [];
% 直接修改聚簇，删除重复的GT（保持原始顺序）
for i = 1:length(clusters)
    % 获取当前聚簇
    current_cluster = clusters{i};
    % 筛选出还未被覆盖的GT
    to_keep = [];
    for j = 1:length(current_cluster)
        gt = current_cluster(j);
        if ~ismember(gt, covered_GTs)
            to_keep = [to_keep, gt];
            covered_GTs = [covered_GTs, gt];
        end
    end
    % 直接修改当前聚簇（保持原始顺序）
    clusters{i} = to_keep;
end
G = size(sorted_VBS_positions, 1); % 更新VBS数量
fid = fopen(clust_xieru, 'w');
for i = 1:length(clusters)
    fprintf(fid, '%s\n', strjoin(string(clusters{i}), ','));
end
fclose(fid);
%% 6.CVX 优化 悬停点位置、USV位置、悬停时间、起始时间

% **** 重要
cvx_clear  
%cvx_solver mosek
s_current = ones(G, 2);

v_a_f_fz = 10*ones(G+1);

v_s_f_fz = 10*ones(G+1);
v_s_h_fz = 10*ones(G);

t_move_fz = 10*ones(G+1);
t_in_fz = 10*ones(G+1);

grad_t = zeros(G+1, 1);
grad_v = zeros(G+1, 1);
f_current = zeros(G+1, 1);

% 无人机参数 UAV
d0_uav_c = 0.6;
P0_uav_c = 80; 
P1_uav_c = 88.6;
v0_uav_c =4.03;
A_uav_c = 0.503;
fai_uav_c = 0.05;
rou_uav_c = 1.225;
Utip_uav_c = 120;
% 无人船参数 USV
rho_usv = 1000; % 水流密度
C_usv = 0.001; % 阻力因子
A_usv = 25; % 船的面积
a_usv_c = 0.5*rho_usv*C_usv*A_usv;

% 初始化信任域
trust_t = 10*ones(G+1, 1);
trust_v = 10*ones(G+1, 1);

cvx_clear  
%cvx_solver mosek

prev_total_energy = Inf; % 初始化前一次的total_energy为无穷大
for kk = 1:100
    fprintf('悬停点求解，第%d次SCA  ', kk); 
    % 计算当前目标函数值和梯度
    for n = 1:G+1
        lamda_B(n)= (   (1 + (v_a_f_fz(n))^4 *( (4*v0_uav_c^4)^(-1) )  )^(1/2) ...
                        -1/2*( (v_a_f_fz(n))^2 )*( (2*v0_uav_c^2)^(-1) )  )^(1/2);
        % 添加数值稳定性检查
        if t_move_fz(n) < 1e-6
            t_move_fz(n) = 1e-6;
        end
        if v_a_f_fz(n) < 1e-6
            v_a_f_fz(n) = 1e-6;
        end
        
        % 计算梯度和当前函数值
        grad_t(n) = P0_uav_c + (3*P0_uav_c/(Utip_uav_c^2))*v_a_f_fz(n)^2 + 0.5*d0_uav_c*rou_uav_c*fai_uav_c*A_uav_c*v_a_f_fz(n)^3;
        grad_v(n) = (3*P0_uav_c/(Utip_uav_c^2))*2*v_a_f_fz(n)*t_move_fz(n) + 0.5*d0_uav_c*rou_uav_c*fai_uav_c*A_uav_c*3*t_move_fz(n)*v_a_f_fz(n)^2;
        f_current(n) = t_move_fz(n)*P0_uav_c*(1 + 3*v_a_f_fz(n)^2/(Utip_uav_c^2)) +0.5*d0_uav_c*rou_uav_c*fai_uav_c*A_uav_c*...
        t_move_fz(n)*v_a_f_fz(n)^3;
        
        % 添加正则化项
        reg_factor = 0.01;
        grad_t(n) = grad_t(n) + reg_factor * (t_move_fz(n) - 10);
        grad_v(n) = grad_v(n) + reg_factor * (v_a_f_fz(n) - 10);
    end
    % 使用MOSEK求解器（如果可用）  
    % 或调整CVX精度
    
    % 在进入CVX之前计算位移参考和梯度
    d_current = cell(G+1,1);
    d_current{1} = s_current(1,:) - sorted_s_start;
    for n=2:G
        d_current{n} = s_current(n,:) - s_current(n-1,:);
    end
    d_current{G+1} = sorted_f_end - s_current(G,:);
    
    f0 = zeros(G+1,1);

    cvx_precision high
    cvx_begin quiet
        % 优化变量 每个变量的意思 
        % s(G,2):UAV的G个悬停点 f(G,2)和s(G,2)是一样的 后面我加了约束 s(G,2) = f(G,2)
        % t_in(G,1) G个悬停点的悬停时间
        % t_move(G+1,1) G+1个飞行模式的时间
        % usv_s(G,2) USV在第G个阶段下飞行模式的终点 
        % usv_f(G,2) USV在第G个阶段下悬停模式的终点 
        % t_s(G,1) =  t_in(G,1)     t_m(G+1,1) =  t_move(G+1,1)
        % 当时写的时候为了松弛有很多垃圾变量没删（重复定义的 没用到的） 最后我统一删除 
        variables s(G,2) f(G,2) 
        variables t_in(G,1)       
        variables t_move(G+1,1)    
        variables Rs_time(G,1)
        variables usv_s(G,2) usv_f(G,2) t_s(G,1) t_m(G+1,1) 
        % 无人机的
        variables v_a_f(G+1)     %无人机在第G+1个飞行模式下的平均速度
        variables lambda_B1(G+1) %处理无人机能耗第三项引入的辅助变量
        % 无人船的
        % 飞行
        variables v_s_f(G+1,2)     %无人船在第G+1个飞行模式下的平均速度
        % 悬停
        variables v_s_h(G,2)       %无人船在第G个悬停模式下的平均速度
        expression expr % 定义表达式

        %expression total_energy
        
        expression total_energy
        expression UAV_energy
        expression USV_energy

        subject to
            % 所有约束条件...
            % 1.起点到第一个悬停点的移动约束
            norm(s(1,:) - sorted_s_start) <= Vmax * t_move(1);
            norm(usv_s(1,:) - sorted_s_start) <= Vmax_usv * t_m(1);
            
            % 2.悬停点之间的移动约束
            for g = 1:G-1
                norm(s(g+1,:) - f(g,:)) <= Vmax * t_move(g+1);
                norm(usv_s(g+1,:) - usv_f(g,:)) <= Vmax_usv * t_m(g+1);
            end
        
            % 3.对于UAV的悬停点形成
            for g = 1:G
                norm(s(g,:)- f(g,:)) <=0;
            end
        
            % 5.UAV悬停时，USV的轨迹动向
            for g = 1:G
                norm(usv_f(g,:) - usv_s(g,:)) <= Vmax_usv * t_s(g);
            end
            
            % 6.最后一个VBS到终点的移动约束
            norm(sorted_f_end - f(G,:)) <= Vmax * t_move(G+1);
            norm(sorted_f_end - usv_f(G,:)) <= Vmax_usv * t_m(G+1);
        
            % 7. 覆盖半径约束  
           for g = 1:G
                covered_GTs = clusters{g};
                for k_idx = 1:length(covered_GTs)
                    k = covered_GTs(k_idx);
                    GT_pos = GT_positions(k,:);
                
                    % 约束：悬停点必须在覆盖半径内
                    %norm(s(g,:) - GT_pos) <=  coverage_radius;
                    %Rs_time(g) >= 1*(length(covered_GTs))*Upsilon*pow_pos((pow_pos(H, 2)+pow_pos(norm(s(g,:) - GT_pos), 2)),2);
                
                    Hove = (length(covered_GTs));
                    UpsilonH = (P/Hove)*M*acc*(neta*beta^2)*((16*pi)^(-1))*(N0_UL^(-1));
                    coverage_radius = (( UpsilonH/(SNR_one))^(0.5) -H^2  )^(0.5); % UAV 覆盖半径 D
                    Upsilon = SNR_k/(UpsilonH); 

                    norm(s(g,:) - GT_pos) <=  coverage_radius;
                    %norm(s(g,:) - GT_pos) <=  0;  % 这里是二维度坐标 
                    Rs_time(g) >= 1*Upsilon*pow_pos((pow_pos(H, 2)+pow_pos(norm(s(g,:) - GT_pos), 2)),2);
                end
           end

            % 8.UAV悬停时间
            for g = 1:G
                t_in(g) >= Rs_time(g);
                norm(f(g,:) - s(g,:)) <= Vmax * t_in(g);
            end
            
            % 9.通信距离约束
            for g = 1:G
                norm(s(g,:) - usv_s(g,:)) <= fly_com_d;
                norm(s(g,:) - usv_f(g,:)) <= fly_com_d;
            end
            
            % 10.时间一致性约束
            for g = 1:G
                t_s(g,:) == t_in(g,:);
            end

            for g = 1:G+1
                t_m(g,:) == t_move(g,:);
            end
        
            %11.UAV的速度、位置、时间 三者关系处理

            %abs(s(1,:)) <= expr;
            %norm(s(1,:)) >= expr;

            norm(s(1,:)) <= t_move(1) * v_a_f_fz(1) + t_move_fz(1) * v_a_f(1) - t_move_fz(1) * v_a_f_fz(1);

            for n = 2:G
                norm(s(n,:)-s(n-1,:)) <= t_move(n) * v_a_f_fz(n) + t_move_fz(n) * v_a_f(n) - t_move_fz(n) * v_a_f_fz(n);
            end
            norm(sorted_f_end - s(G,:)) <=  t_move(G+1) * v_a_f_fz(G+1) + t_move_fz(G+1) * v_a_f(G+1) - t_move_fz(G+1) * v_a_f_fz(G+1);

            % UAV 能耗第三项
            for n=1:G+1
                pow_pos(inv_pos(lambda_B1(n)),2) <= ...
                2*lamda_B(n)*( lambda_B1(n)-lamda_B(n) ) + ...
                ( lamda_B(n) )^2 + ...
                ( v_a_f_fz(n)/v0_uav_c )^2 + ...
                2/v0_uav_c^2 *(   v_a_f_fz(n) * ( v_a_f(n) - v_a_f_fz(n) ) );
            end

            % 飞行模式下
            norm( usv_s(1,:) ) <= t_move(1) * v_s_f_fz(1) + t_move_fz(1) * v_s_f(1) - t_move_fz(1) * v_s_f_fz(1);
            for n = 2:G
                norm(usv_s(n,:)-usv_f(n-1,:)) <= t_move(n) * v_s_f_fz(n) + t_move_fz(n) * v_s_f(n) - t_move_fz(n) * v_s_f_fz(n);
            end
            norm(sorted_f_end - usv_f(G,:)) <= t_move(G+1) * v_s_f_fz(G+1) + t_move_fz(G+1) * v_s_f(G+1) - t_move_fz(G+1) * v_s_f_fz(G+1);

            % 悬停模式下
            for n = 1:G
                norm(usv_f(n,:)-usv_s(n,:)) <= t_in(n) * v_s_h_fz(n) + t_in_fz(n) * v_s_f(n) - t_in_fz(n) * v_s_h_fz(n);
            end

            %14. 变量范围约束
            for n = 1:G
                t_move(n) <= 200;
                t_move(n) >= 0.1;
                v_a_f(n) >= 0.1;
                v_a_f(n) <= 50;
                % 信任域约束
                abs(t_move(n) - t_move_fz(n)) <= trust_t(n);
                abs(v_a_f(n) - v_a_f_fz(n)) <= trust_v(n);
                t_in(n) <= 100;
                t_in(n) >= 0.1;
                v_s_f(n) >= 0.1;
                v_s_f(n) <= 50;
                v_s_h(n) >= 0.1;
                v_s_h(n) <= 50;
                % 信任域约束
                abs(t_move(n) - t_move_fz(n)) <= trust_t(n);
                abs(v_a_f(n) - v_a_f_fz(n)) <= trust_v(n);
                abs(t_in(n) - t_in_fz(n)) <= trust_t(n);
                abs(v_s_f(n) - v_s_f_fz(n)) <= trust_v(n);
                abs(v_s_h(n) - v_s_h_fz(n)) <= trust_v(n);
            end
           n = G+1;
                t_move(n) <= 100;
                t_move(n) >= 0.1;
                v_a_f(n) >= 0.1;
                v_a_f(n) <= 50;
                % 信任域约束
                abs(t_move(n) - t_move_fz(n)) <= trust_t(n);
                abs(v_a_f(n) - v_a_f_fz(n)) <= trust_v(n);
                v_s_f(n) >= 0.1;
                v_s_f(n) <= 50;
                abs(v_a_f(n) - v_a_f_fz(n)) <= trust_v(n);
            %15. 目标函数
            total_energy = 0;
            UAV_energy = 0;
            USV_energy = 0;
            for n = 1:G
                term1 = f_current(n);
                term2 = grad_t(n) * (t_move(n) - t_move_fz(n));
                term3 = grad_v(n) * (v_a_f(n) - v_a_f_fz(n));
                %term4 = t_move(n)*v_s_f_fz(n)^2+t_move_fz(n)*v_s_f(n)^2+v_s_f_fz(n)^2+2*t_move_fz(n)*v_s_f(n);
                %term5 = t_in(n)*v_s_h_fz(n)^2+t_in_fz(n)*v_s_h(n)^2+v_s_h_fz(n)^2+2*t_in_fz(n)*v_s_h(n);
                %term4 = t_move_fz(n)*v_s_f_fz(n)^2+v_s_f_fz(n)^2*(t_move(n)-t_move_fz(n))+2*t_move_fz(n)*v_s_f_fz(n)*(v_s_f(n)-v_s_f_fz(n));
                term4 = t_move_fz(n)*v_s_f_fz(n)^3 + v_s_f_fz(n)^3*(t_move(n)-t_move_fz(n)) + 3*t_move_fz(n)*v_s_f_fz(n)^2*(v_s_f(n)-v_s_f_fz(n));
                
                %term5 = t_in_fz(n)*v_s_h_fz(n)^2+v_s_h_fz(n)^2*(t_in(n)-t_in_fz(n))+2*t_in_fz(n)*v_s_h_fz(n)*(v_s_h(n)-v_s_h_fz(n));
                term5 = t_in_fz(n)*v_s_h_fz(n)^3 + v_s_h_fz(n)^3*(t_in(n)-t_in_fz(n)) + 3*t_in_fz(n)*v_s_h_fz(n)^2*(v_s_h(n)-v_s_h_fz(n));
                
                term6 = t_in(n)*(P1_uav_c+P0_uav_c);
                term7 =t_move_fz(n)*P1_uav_c*lamda_B(n)+P1_uav_c*lamda_B(n)*(t_move(n)-t_move_fz(n))+P1_uav_c*t_move_fz(n)*(lambda_B1(n)-lamda_B(n)); % 能耗第三项
                % 添加小常数避免数值问题
                epsilon = 1e-6;
                total_energy = total_energy + term1 + term2 + term3 + a_usv_c*term4 + a_usv_c*term5 +term6+term7+ epsilon;
                UAV_energy = UAV_energy + (term1 + term2 + term3 + term6 + term7 + epsilon);
                USV_energy = USV_energy + (a_usv_c*term4 + a_usv_c*term5);
            end
            n = G+1;
            term1 = f_current(n);
            term2 = grad_t(n) * (t_move(n) - t_move_fz(n));
            term3 = grad_v(n) * (v_a_f(n) - v_a_f_fz(n));
            %term4 = t_move(n)*v_s_f_fz(n)^2+t_move_fz(n)*v_s_f(n)^2+v_s_f_fz(n)^2+2*t_move_fz(n)*v_s_f(n);
            %term4 = t_move_fz(n)*v_s_f_fz(n)^2+v_s_f_fz(n)^2*(t_move(n)-t_move_fz(n))+2*t_move_fz(n)*v_s_f_fz(n)*(v_s_f(n)-v_s_f_fz(n));
            term4 = t_move_fz(n)*v_s_f_fz(n)^3 + v_s_f_fz(n)^3*(t_move(n)-t_move_fz(n)) + 3*t_move_fz(n)*v_s_f_fz(n)^2*(v_s_f(n)-v_s_f_fz(n));
            
            term7 =t_move_fz(n)*P1_uav_c*lamda_B(n)+P1_uav_c*lamda_B(n)*(t_move(n)-t_move_fz(n))+P1_uav_c*t_move_fz(n)*(lambda_B1(n)-lamda_B(n)); % 能耗第三项
            % 添加小常数避免数值问题
            epsilon = 1e-6;
            total_energy = total_energy + term1 + term2 + term3 + a_usv_c*term4+term7+ epsilon;
            UAV_energy = UAV_energy + (term1 + term2 + term3 + term7 + epsilon);
            USV_energy = USV_energy + (a_usv_c*term4);
            minimize(10000*UAV_energy/E_uav_max+10000*USV_energy/E_usv_max)
            %minimize(total_energy)
    cvx_end

    % 检查求解状态
    if strcmp(cvx_status, 'Solved') && isfinite(cvx_optval)
        % 更新参考点
        s_current = s;
        v_a_f_fz = v_a_f;
        v_s_f_fz = v_s_f;
        v_s_h_fz = v_s_h;
        t_move_fz = t_move;
        t_in_fz = t_in;
        % 自适应调整信任域
        for n = 1:G+1
            % 计算实际函数值与线性近似的差异
            actual_f = t_move(n)*P0_uav_c*(1 + 3*v_a_f(n)^2/(Utip_uav_c^2)) + t_move(n)*v_a_f(n)^3;
            approx_f = f_current(n) + grad_t(n)*(t_move(n) - t_move_fz(n)) + grad_v(n)*(v_a_f(n) - v_a_f_fz(n));
            approx_error = abs(actual_f - approx_f);
            
            % 根据近似误差调整信任域
            if approx_error > 0.1
                trust_t(n) = trust_t(n) * 0.8;
                trust_v(n) = trust_v(n) * 0.8;
            else
                trust_t(n) = min(trust_t(n) * 1.2, 10);
                trust_v(n) = min(trust_v(n) * 1.2, 5);
            end
        end
        
        disp(['迭代 ', num2str(kk), ': 目标函数值 = ', num2str(cvx_optval)]);

    if abs( (total_energy - prev_total_energy)/prev_total_energy ) < threshold
        break;
    end
    % 更新前一次的total_energy
    prev_total_energy = total_energy;

    else
        disp(['迭代 ', num2str(kk), ': 求解失败，状态 = ', cvx_status]);
        break;
    end
end
%disp(v_a_f_fz);
%disp(v_a_f_shiliang_fz);
disp('第一阶段能耗');
disp(total_energy);
fprintf('UAV能耗=%.6f, USV能耗=%.6f\n', double(UAV_energy), double(USV_energy));

uav_s = s(:, 1:2);  % 取前两列
uav_s(:, 3) = H;
uav_s = uav_s';
usv_s = usv_s(:, 1:2);  % 取前两列
usv_s(:, 3) = 0;
usv_s = usv_s';
usv_f = usv_f(:, 1:2);  % 取前两列
usv_f(:, 3) = 0;
usv_f = usv_f';
Flight_time = ceil(t_move);
Hover_time = ceil(t_in);
disp("每个飞行模式的持续时间:");
disp(Flight_time');
disp("每个悬停模式的悬停时间:");
disp(Hover_time');
disp("系统运行总时间:");
disp(sum(t_in) + sum(t_move));
leng = length(Flight_time);


%% 7.时间离散化求解波束和速度--分界线
energy_sum = 0;      % 总功耗
energy_sum_uav_p = 0; % UAV的推进能耗
energy_sum_uav_s = 0; % UAV的通感能耗
energy_sum_uav_c = 0; % UAV的通感能耗
energy_sum_usv = 0 ; %USV的推进能耗
%% 8.生成初始路径 UAV和USV的初始路径

% UAV的初始轨迹点 只生成第一阶段的
start_point = uav_start;
end_point = uav_s(:,1);
num_points = Flight_time(1);
disp(num_points);
path_uav = generatePath(start_point, end_point, num_points);
path_uav = path_uav';
% USV的初始轨迹点
start_point = usv_start;
end_point = usv_s(:,1);
num_points = Flight_time(1);
path_usv = generatePath(start_point, end_point, num_points);
path_usv = path_usv';

% 提取初始位置点
N_f1 = num_points+1;
s_uav = path_uav;
start_uav = path_uav(:, 1);
end_uav = path_uav(:, num_points+1);
s_usv = path_usv;
start_usv = path_usv(:, 1);
end_usv = path_usv(:, num_points+1);


threshold = 1e-2; % 设定一个阈值来判断迭代结果的变化程度
prev_total_energy = Inf; % 初始化前一次的total_energy为无穷大

% 飞行模式下的交替优化算法
for kk = 1:1:T_K
    fprintf('飞行模式第一阶段，第%d次AO', kk); 
    disp(' ');
    P_t = zeros(1,N_f1);     % 存储功率结果
    C_rate = zeros(1,N_f1);  % 存储功率结果
    D_max = zeros(1,N_f1);   % 存储最远距离 
    for i = 1:N_f1
        D_max(i) = norm(s_uav( :,i) - s_usv( :,i));
    end
    for i = 1:N_f1
            uav_pos11 = s_uav( :,i);
            uav_pos11 = uav_pos11';
            usv_pos11 = s_usv( :,i);
            usv_pos11 = usv_pos11';
            [P_t(i),C_rate(i)] = optimize_precode_with_communication(M,uav_pos11, usv_pos11, gamma_value);
    end
    %disp('通信速率：');
    %disp(C_rate);
    %disp('分配功率：');
    %disp(P_t);
    com_energy = sum(P_t);
    
    P_t = real(P_t);
    Pc_set = P_t;      % 存储通信功率
    Ps_set = zeros(1,N_f1);  % 储存感知功率
    [ s_uav, s_usv, total_energy, UAV_energy, USV_energy ] = optimize_trajectory_with_fly(D_max,Nh,H,N_f1,Q,...
                    v_uav_max,v_usv_max,obs_pos,r1,delta,...
                    P_t,gamma_value,com_energy,s_uav,s_usv,start_uav, end_uav, start_usv, end_usv);
    %disp(total_energy);
    fprintf('通信能耗=%.6f, UAV能耗=%.6f, USV能耗=%.6f\n', double(com_energy), double(UAV_energy), double(USV_energy));

    if abs( (total_energy - prev_total_energy)/prev_total_energy ) < threshold
        break;
    end
    % 更新前一次的total_energy
    prev_total_energy = total_energy;
end
end_path_uav = s_uav;
end_path_usv = s_usv;

energy_sum_uav_c = energy_sum_uav_c+com_energy;
energy_sum_uav_p = energy_sum_uav_p+UAV_energy;
energy_sum_usv = energy_sum_usv+USV_energy;
energy_sum = energy_sum+total_energy;


fprintf('悬停模式第一阶段，第%d次AO', kk); 
num_points = Hover_time(1);
N_f1 = num_points+1;
targets = GT_positions(:, 1:2);  % 取前两列
targets(:, 3) = 0;
targets = targets';
covered_GTs = clusters{1};
GT_positions = zeros(3, length(covered_GTs));
for k_idx = 1:length(covered_GTs)
    k = covered_GTs(k_idx);
    GT_pos = targets(:,k);
    GT_positions(:, k_idx) = GT_pos; % 填入结果矩阵
end



% 重复生成悬停点的轨迹
v = uav_s(:,1);
path_uav = repmat(v, 1, num_points+1);

start_point = usv_s(:,1);
end_point = usv_f(:,1);
path_usv = generatePath(start_point, end_point, num_points);
path_usv = path_usv';

s_uav = path_uav;
start_uav = path_uav(:, 1);
end_uav = path_uav(:, num_points+1);

s_usv = path_usv;
start_usv = path_usv(:, 1);
end_usv = path_usv(:, num_points+1);

prev_total_energy = Inf; % 初始化前一次的total_energy为无穷大
for kk=1:1:T_K
    fprintf('悬停模式第一阶段，第%d次AO', kk);
    D_max = zeros(1,N_f1);  % 存储距离
    for i = 1:N_f1
        D_max(i) = norm(s_uav( :,i) - s_usv( :,i));
    end
    % 无干扰
    %[Wc_h, Ws_k_h, power_vec,power_vec_s] = optimize_precode_with_sensing(N_f1, length(covered_GTs), start_uav, GT_positions, s_usv,com_lb,sen_lb);
    % 有干扰
    [Wc_h, Ws_k_h, power_vec,power_vec_s] = optimize_precode_with_sensingCSI(N_f1, length(covered_GTs), start_uav, GT_positions, s_usv,com_lb,sen_lb);
    
    P_t1 = power_vec;
    P_t2 = power_vec_s;
    com_energy = sum(P_t1);
    com_energy = com_energy+sum(power_vec_s);
    %disp(P_t1);
    %disp(P_t2);
    %disp(com_energy);

    [ s_uav, s_usv, total_energy, UAV_energy, USV_energy ] = optimize_trajectory_with_hover(D_max,Nh,H,N_f1,Q,...
                    v_uav_max,v_usv_max,obs_pos,r1,delta,...
                    P_t1,com_lb,com_energy,s_uav,s_usv,start_uav, end_uav, start_usv, end_usv,length(covered_GTs),Wc_h,Ws_k_h);
    if abs( (total_energy - prev_total_energy)/prev_total_energy ) < threshold
        break;
    end

    fprintf('通感能耗=%.6f, UAV能耗=%.6f, USV能耗=%.6f\n', double(com_energy), double(UAV_energy), double(USV_energy));

    % 更新前一次的total_energy
    prev_total_energy = total_energy;
end
Pc_set =[Pc_set,P_t1(2:end)];
Ps_set =[Ps_set,power_vec_s(2:end)];

% 轨迹合并
energy_sum_uav_c = energy_sum_uav_c+com_energy;
energy_sum_uav_p = energy_sum_uav_p+UAV_energy;
energy_sum_usv = energy_sum_usv+USV_energy;
energy_sum = energy_sum+total_energy;

end_path_uav1 = s_uav;
end_path_usv1 = s_usv;
end_path_uav = [end_path_uav,end_path_uav1(:, 2:end)];
end_path_usv = [end_path_usv,end_path_usv1(:, 2:end)];

%% 9. 中间阶段 
for i = 1:1:(leng-2)
    
    %飞行模式
    num_points = Flight_time(i+1);

    start_point = uav_s(:,i);
    end_point = uav_s(:,i+1);
    path_uav = generatePath(start_point, end_point, num_points);
    path_uav = path_uav';
    
    start_point = usv_f(:,i);
    end_point = usv_s(:,i+1);
    path_usv = generatePath(start_point, end_point, num_points);
    path_usv = path_usv';
    
    
    N_f1 = num_points+1;
    P_t1 = zeros(1, N_f1) + 20;
    
    s_uav = path_uav;
    start_uav = path_uav(:, 1);
    end_uav = path_uav(:, num_points+1);
    
    s_usv = path_usv;
    start_usv = path_usv(:, 1);
    end_usv = path_usv(:, num_points+1);

    P_t = zeros(1,N_f1);  % 存储功率结果  
prev_total_energy = Inf; % 初始化前一次的total_energy为无穷大
for kk = 1:1:T_K
    fprintf('当前阶段：%d,飞行模式的第%d次AO', i,kk);
    %fprintf('飞行模式：%d', kk);
    D_max = zeros(1,N_f1);  % 存储功率结果

    for yy = 1:N_f1
        D_max(yy) = norm(s_uav( :,yy) - s_usv( :,yy));
    end
    for i1 = 1:N_f1
            uav_pos11 = s_uav( :,i1);
            uav_pos11 = uav_pos11';
            usv_pos11 = s_usv( :,i1);
            usv_pos11 = usv_pos11';
            %disp(uav_pos11);
            % %disp(usv_pos11);
            P_t(i1) = optimize_precode_with_communication(Q,uav_pos11, usv_pos11, gamma_value);
            P_t(i1) = real( P_t(i1)  );
    end
    %disp(P_t);
    com_energy = sum(P_t);
    P_t = real(P_t);
    
    [ s_uav, s_usv, total_energy, UAV_energy, USV_energy ] = optimize_trajectory_with_fly(D_max,Nh,H,N_f1,Q,...
                    v_uav_max,v_usv_max,obs_pos,r1,delta,...
                    P_t,gamma_value,com_energy,s_uav,s_usv,start_uav, end_uav, start_usv, end_usv);
    if abs( (total_energy - prev_total_energy)/prev_total_energy ) < threshold
        break;
    end
    % 更新前一次的total_energy
    prev_total_energy = total_energy;
end
Pc_set =[Pc_set,P_t(2:end)];

zeros_array = zeros(1, N_f1);
Ps_set =[Ps_set,zeros_array(2:end)];

end_path_uav3= s_uav;
end_path_usv3= s_usv;
energy_sum_uav_c = energy_sum_uav_c+com_energy;
energy_sum_uav_p = energy_sum_uav_p+UAV_energy;
energy_sum_usv = energy_sum_usv+USV_energy;
energy_sum = energy_sum+total_energy;

    % 悬停模式
    num_points = Hover_time(i+1);
    N_f1 = num_points+1;

    
    v = uav_s(:,i+1);
    path_uav = repmat(v, 1, num_points+1);
    
    start_point = usv_s(:,i+1);
    end_point = usv_f(:,i+1);
    path_usv = generatePath(start_point, end_point, num_points);
    path_usv = path_usv';

    s_uav = path_uav;
    start_uav = path_uav(:, 1);
    end_uav = path_uav(:, num_points+1);
    
    s_usv = path_usv;
    start_usv = path_usv(:, 1);
    end_usv = path_usv(:, num_points+1);
    prev_total_energy = Inf; % 初始化前一次的total_energy为无穷大
    for kk = 1:1:T_K
        fprintf('当前阶段：%d,悬停模式的第%d次AO', i,kk);
        D_max = zeros(1,N_f1);  % 存储功率结果
        for yy = 1:N_f1
            D_max(yy) = norm(s_uav( :,yy) - s_usv( :,yy));
        end
        [Wc_h, Ws_k_h, power_vec,power_vec_s] = optimize_precode_with_sensingCSI(N_f1, length(covered_GTs), start_uav, GT_positions, s_usv,com_lb,sen_lb);
        P_t1 = real(power_vec);
        %disp(power_vec);
        %disp(power_vec_s);
        com_energy = sum(P_t1);
        com_energy = com_energy+sum(power_vec_s);
        [ s_uav, s_usv, total_energy, UAV_energy, USV_energy ] = optimize_trajectory_with_hover(D_max,Nh,H,N_f1,Q,...
                    v_uav_max,v_usv_max,obs_pos,r1,delta,...
                    P_t1,com_lb,com_energy,s_uav,s_usv,start_uav, end_uav, start_usv, end_usv,length(covered_GTs),Wc_h,Ws_k_h);
         s_uav = s_uav;
         s_usv = s_usv;
        if abs( (total_energy - prev_total_energy)/prev_total_energy ) < threshold
         break;
        end
        % 更新前一次的total_energy
        prev_total_energy = total_energy;
    end
    Pc_set =[Pc_set,P_t1(2:end)];
    Ps_set =[Ps_set,power_vec_s(2:end)];
    end_path_uav4 = s_uav;
    end_path_usv4 = s_usv;
    energy_sum_uav_c = energy_sum_uav_c+com_energy;
energy_sum_uav_p = energy_sum_uav_p+UAV_energy;
energy_sum_usv = energy_sum_usv+USV_energy;
energy_sum = energy_sum+total_energy;

    end_path_uav = [end_path_uav,end_path_uav3(:, 2:end),end_path_uav4(:, 2:end)];
    end_path_usv = [end_path_usv,end_path_usv3(:, 2:end),end_path_usv4(:, 2:end)];
end

%% 10.最后一段 只有飞行模式
start_point = uav_s(:,leng-1);
end_point = uav_end;
num_points = Flight_time(leng);
path_uav = generatePath(start_point, end_point, num_points);
path_uav = path_uav';

start_point = usv_f(:,leng-1);
end_point = usv_end;
path_usv = generatePath(start_point, end_point, num_points);
path_usv = path_usv';

N_f1 = num_points+1;
s_uav = path_uav;
start_uav = path_uav(:, 1);
end_uav = path_uav(:, num_points+1);

s_usv = path_usv;
start_usv = path_usv(:, 1);
end_usv = path_usv(:, num_points+1);

    
P_t = zeros(1,N_f1);  % 存储功率结果  
prev_total_energy = Inf; % 初始化前一次的total_energy为无穷大
for kk = 1:1:1
    fprintf('最后一阶段的飞行模式：%d', kk);
    for i1 = 1:N_f1
            uav_pos11 = s_uav( :,i1);
            uav_pos11 = uav_pos11';
            usv_pos11 = s_usv( :,i1);
            usv_pos11 = usv_pos11';
            P_t(i1) = optimize_precode_with_communication(Q,uav_pos11, usv_pos11, gamma_value);
            P_t(i1) = real( P_t(i1)  );
    end
    for yy = 1:N_f1
        D_max(yy) = norm(s_uav( :,yy) - s_usv( :,yy));
    end
    [s_uav, s_usv, total_energy, UAV_energy, USV_energy ] = optimize_trajectory_with_fly(D_max,Nh,H,N_f1,Q,...
                    v_uav_max,v_usv_max,obs_pos,r1,delta,...
                    P_t,gamma_value,com_energy,s_uav,s_usv,start_uav, end_uav, start_usv, end_usv);

    if abs(total_energy - prev_total_energy) < 1e-1
            break;
    end
        % 更新前一次的total_energy
    prev_total_energy = total_energy;
end
%disp(P_t);
P_t = real(P_t);
Pc_set =[Pc_set,P_t(2:end)];

zeros_array = zeros(1, N_f1);
Ps_set =[Ps_set,zeros_array(2:end)];


end_path_uav3 = s_uav;
end_path_usv3 = s_usv;
energy_sum_uav_c = energy_sum_uav_c+com_energy;
energy_sum_uav_p = energy_sum_uav_p+UAV_energy;
energy_sum_usv = energy_sum_usv+USV_energy;
energy_sum = energy_sum+total_energy;

end_path_uav = [end_path_uav,end_path_uav3(:, 2:end)];
end_path_usv = [end_path_usv,end_path_usv3(:, 2:end)];
% 合并
path_uav = end_path_uav;

%path_usv = [end_path_usv,path_usv_lianxu(:, 2:end)];
path_usv = end_path_usv;

%disp(path_uav);
window_size = 2.5;  % 设置平滑窗口大小
% 对 UAV 路径的每一维进行平滑
path_uav_smooth = zeros(size(path_uav));  % 创建与原始路径相同的大小
for i = 1:3
    path_uav_smooth(i, :) = smooth(path_uav(i, :), window_size);  % 平滑每一维（x, y, z）
end

% 对 USV 路径的每一维进行平滑
path_usv_smooth = zeros(size(path_usv));  % 创建与原始路径相同的大小
for i = 1:3
    path_usv_smooth(i, :) = smooth(path_usv(i, :), window_size);  % 平滑每一维（x, y, z）
end
path_uav = path_uav_smooth';
path_usv = path_usv_smooth';


%% ====== 计算 UAV 速度变化惩罚能耗 E_tr ======
%{
alpha_uav = 1;     % 你论文里的效率系数（无量纲）
m_uav     = 3.0;     % UAV 质量(kg)，按你设定改
eta_uav   = 0.5 * alpha_uav * m_uav;   % = 1/2 * alpha * m

% 速度向量 v[n] = (p[n+1]-p[n]) / delta, 大小: (N-1)×3
v_uav_vec = diff(path_uav,1,1) / delta;

% Δv[n] = v[n+1]-v[n], 大小: (N-2)×3
dv_uav_vec = diff(v_uav_vec,1,1);

% 论文里的惩罚：eta * sum(||Δv||^2) * delta
E_tr_uav = eta_uav * sum( sum(dv_uav_vec.^2, 2) ) * delta;

disp(['UAV 惩罚能耗 = ', num2str(E_tr_uav), ' (J if consistent)']);
%}
alpha_uav = 1;      
m_uav     = 3.0;    
eta_uav   = 0.5 * alpha_uav * m_uav;   % = 1/2 * alpha * m

% 速度向量 v[n] = (p[n+1]-p[n]) / delta
v_uav_vec = diff(path_uav,1,1) / delta;    % (N-1)×3

% 每个时刻速度平方 ||v[n]||^2
v_sq = sum(v_uav_vec.^2,2);                 % (N-1)×1

% 相邻时刻速度平方之差
dv_sq = diff(v_sq);                         % (N-2)×1

% 取绝对值再累加
E_tr_uav = eta_uav * sum(abs(dv_sq)) * delta;

disp(['UAV 惩罚能耗 = ', num2str(E_tr_uav), ' J']);

% 将 UAV 路径保存为 uav_path.txt


writematrix(path_uav, uav_xieru, 'Delimiter', 'tab');
writematrix(Pc_set, pc_xieru, 'Delimiter', 'tab');
writematrix(Ps_set, ps_xieru, 'Delimiter', 'tab');

% 将 USV 路径保存为 usv_path.txt
writematrix(path_usv, usv_xieru, 'Delimiter', 'tab');

disp("总能耗");
disp(energy_sum);
fprintf('通感能耗=%.6f, UAV能耗=%.6f, USV能耗=%.6f\n', double(energy_sum_uav_c), double(energy_sum_uav_p), double(energy_sum_usv));
disp("总能耗");
disp(double(energy_sum_uav_c)+double(energy_sum_uav_p)+double(energy_sum_usv)+E_tr_uav);
%% 绘图
% 显示路径点
%disp('路径点集合:');
%disp(path_uav);
%disp(path_usv);
% 绘制路径
figure;
hold on;
grid on;
view(3); % 3D视图
scatter3(obs_pos(1,:), obs_pos(2,:), obs_pos(3,:), 150, 'k', 'x', 'LineWidth',2.5);
% 悬停点
scatter3(uav_s(1,:), uav_s(2,:), uav_s(3,:), 100, 'p', '^', 'LineWidth', 3, ...
         'MarkerFaceColor', [1, 1, 1]); % 浅红色填充
% 将 UAV 路径保存为 uav_path.txt
uav_s1 = uav_s';

writematrix(uav_s1, hover_xietu, 'Delimiter', 'tab');

%scatter3(usv_s(1,:), usv_s(2,:), usv_s(3,:), 150, 'p', '^', 'LineWidth', 3, ...
%         'MarkerFaceColor', [1, 1, 1]); % 浅红色填充
%scatter3(usv_f(1,:), usv_f(2,:), usv_f(3,:), 150, 'p', '^', 'LineWidth', 3, ...
%         'MarkerFaceColor', [1, 1, 1]); % 浅红色填充

targets = GT_positions11(:, 1:2);  % 取前两列
targets(:, 3) = 0;
targets = targets';

scatter3(VBS_positions(:,1), VBS_positions(:,2),0, 100, 'r', '^', 'LineWidth', 3, ...
         'MarkerFaceColor', [1, 1, 1]); % VBS 位置
% 画虚拟基站服务圈
for g = 1:size(VBS_positions,1)
   viscircles(VBS_positions(g,:), coverage_radius, 'Color', 'r', 'LineStyle', '--');
end

% 感知目标
%scatter3(targets(1,:), targets(2,:), targets(3,:), 150, 'k', '^', 'LineWidth', 3, ...
%         'MarkerFaceColor', [1, 1, 1]); % 浅灰色填充
%scatter3(targets(1,:), targets(2,:), targets(3,:), 150, 'k', '^', 'LineWidth', 0.1, ...
%         'MarkerFaceColor', [0, 0.5, 0]); % 使用绿色填充
scatter3(targets(1,:), targets(2,:), targets(3,:), ...
    150, ...  % 标记大小保持不变
    '^', ...  % 三角形标记
    'MarkerFaceColor', [0.0, 0.5, 0.0], ...  % 柔和的浅绿色
    'MarkerEdgeColor', 'none', ...  % 移除边框线
    'LineWidth', 0.1);  % 进一步确保无边框线


plot3(0, 0, 100, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % 红色实心圆点
plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % 红色实心圆点

plot3(300, 0, 100, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % 红色实心圆点
plot3(300, 0, 0, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % 红色实心圆点

% 创建安全区域 (灰色虚线圆圈)
theta = linspace(0, 2*pi, 50);
for i = 1:size(obs_pos,2)
    % XY平面安全圆
    x_circle = obs_pos(1,i) + r1 * cos(theta);
    y_circle = obs_pos(2,i) + r1 * sin(theta);
    z_circle = ones(size(theta)) * obs_pos(3,i);
    plot3(x_circle, y_circle, z_circle, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 2);
end

water_speed_max = 30;
     % 水流速度场设置
% 创建10x10网格
Width = 300;
grid_size = Width/20;
x_grid = 0:grid_size:Width; 
y_grid = 0:grid_size:Width; 
[X_grid, Y_grid] = meshgrid(x_grid, y_grid);
water_level = 0.1; % 略高于0的水面高度
% 水流速度模型
water_vel_x = zeros(size(X_grid));
water_vel_y = zeros(size(Y_grid));
for i = 1:size(X_grid,1)
    for j = 1:size(X_grid,2)
        x = X_grid(i,j);
        y = Y_grid(i,j);
        
        % 主流方向（x方向）
        main_flow = water_speed_max * (0.8 + 0.2*cos(0.05*x));
        
        % 涡旋分量（y方向）
        vortex = 0.4 * water_speed_max * sin(0.05*x) .* cos(0.05*y);
        
        water_vel_x(i,j) = main_flow;
        water_vel_y(i,j) = vortex;
    end
end

% 水流速度函数
water_vel_func = @(x,y) [water_speed_max * (1 + 0.2*cos(0.05*x)); 
                         water_speed_max * 0.4 * sin(0.05*x) .* cos(0.05*y)];


% 绘制水流场 (更小的箭头)
quiver_scale = 0.5; % 缩小箭头尺寸
arrow_size = 0.2; % 箭头头大小
% 在三维空间绘制水流场
%quiver3(X_grid, Y_grid, water_level * ones(size(X_grid)),...
%        water_vel_x, water_vel_y, zeros(size(water_vel_x)), quiver_scale, ...
%       'Color', [0.2, 0.6, 0.8], 'LineWidth', 1.2, 'MaxHeadSize', arrow_size);
hold on;

% 添加流速标签
text_offset_x = -3;
text_offset_y = -3;
for i = 1:size(X_grid,1)
    for j = 1:size(X_grid,2)
        % 计算速度大小
        v = norm([water_vel_x(i,j), water_vel_y(i,j)]);
        
        % 在箭头旁边添加流速文本
       % text(X_grid(i,j) + text_offset_x, Y_grid(i,j) + text_offset_y, ...
        %     sprintf('%.1f', v), 'FontSize', 9, 'Color', 'b');
    end
end
%% 速度差异
% 计算UAV的速度（欧氏距离）
vel_uav = [0; vecnorm(diff(path_uav), 2, 2)];  % 在路径上的每个点的速度
% 计算USV的速度（欧氏距离）
vel_usv = [0; vecnorm(diff(path_usv), 2, 2)];  % 在路径上的每个点的速度

% 设置颜色映射
colormap jet;  % 使用Jet色图
colorbar;      % 显示颜色条
caxis([min([vel_uav; vel_usv]), max([vel_uav; vel_usv])]);  % 设置颜色条范围，统一UAV和USV的速度范围

% 获取颜色映射的RGB值
cmap = colormap;  % 获取颜色图
ncol = size(cmap, 1);  % 获取颜色图的长度

% 绘制UAV路径的直线，速度决定颜色
for i = 1:length(path_uav)-1
    % 获取当前路径段的速度
    color_val = vel_uav(i);  % 获取每一段路径的速度
    % 将速度值映射到颜色图的索引
    idx = round((color_val - min([vel_uav; vel_usv])) / (max([vel_uav; vel_usv]) - min([vel_uav; vel_usv])) * (ncol - 1)) + 1;
    % 获取该速度值对应的颜色
    color = cmap(idx, :);  % 获取对应的颜色
    plot3([path_uav(i,1), path_uav(i+1,1)], ...
          [path_uav(i,2), path_uav(i+1,2)], ...
          [path_uav(i,3), path_uav(i+1,3)], ...
          'LineWidth', 3, 'Color', color);
end

% 绘制USV路径的直线，速度决定颜色
for i = 1:length(path_usv)-1
    % 获取当前路径段的速度
    color_val = vel_usv(i);  % 获取每一段路径的速度
    % 将速度值映射到颜色图的索引
    idx = round((color_val - min([vel_uav; vel_usv])) / (max([vel_uav; vel_usv]) - min([vel_uav; vel_usv])) * (ncol - 1)) + 1;
    % 获取该速度值对应的颜色
    color = cmap(idx, :);  % 获取对应的颜色
    plot3([path_usv(i,1), path_usv(i+1,1)], ...
          [path_usv(i,2), path_usv(i+1,2)], ...
          [path_usv(i,3), path_usv(i+1,3)], ...
          'LineWidth', 3, 'Color', color);
end


% 绘制路径线
%plot3(path_uav(:,1), path_uav(:,2), path_uav(:,3), 'r-', 'LineWidth', 3, 'MarkerSize', 8);
%plot3(path_uav(1,1), path_uav(1,2), path_uav(1,3), 'r', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % 起点
%plot3(path_uav(end,1), path_uav(end,2), path_uav(end,3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); % 终点

%hold on;


% 绘制路径线
%plot3(path_usv(:,1), path_usv(:,2), path_usv(:,3), 'b-', 'LineWidth', 3, 'MarkerSize', 8);
%plot3(path_usv(1,1), path_usv(1,2), path_usv(1,3), 'r', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % 起点
%plot3(path_usv(end,1), path_usv(end,2), path_usv(end,3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); % 终点

%hold on;
%set(gca, 'Color', [0.5 0.5 0.5]);
%set(gca, 'GridColor', [0.5 0.5 0.5]); % 浅灰网格线
% 添加标签和标题
xlabel('X');
ylabel('Y');
zlabel('Z');
%title('从(1,2,5)到(3,6,5)的直角路径');
%legend('路径', '起点 (1,2,5)', '终点 (3,6,5)', '转折点 (1,6,5)', 'Location', 'best');
hold off;
% 设置三维坐标轴范围
xlim([0 300]);   % x轴范围从0到300
ylim([0 300]);   % y轴范围从0到300
zlim([0 100]);   % z轴范围从0到300


% 绘制速度曲线
% 每两点间时间间隔为 1 秒

% 计算速度（欧氏距离）
vel_uav = vecnorm(diff(path_uav), 2, 2);  % 长度 N-1
vel_usv = vecnorm(diff(path_usv), 2, 2);

% 时间轴（每秒一个速度）
time_uav = 1:length(vel_uav);
time_usv = 1:length(vel_usv);

% 绘图
figure;
plot(time_uav, vel_uav, 'r-', 'LineWidth', 1.5); 
hold on;
plot(time_usv, vel_usv, 'b-', 'LineWidth', 1.5);


xlabel('Time (s)');
ylabel('Speed (m/s)');
title('UAV and USV Speed vs. Time');
legend('UAV Speed', 'USV Speed');
grid on;

elapsed_time = toc;  % 结束计时并获取耗时（单位：秒）

fprintf('代码运行时间为 %.4f 分钟\n', elapsed_time/60);
% 折线
function path = generatePath123(start_point, end_point, num_points)
    num_points = num_points-1;
    % 输入验证
    validateattributes(start_point, {'numeric'}, {'vector', 'numel', 3});
    validateattributes(end_point, {'numeric'}, {'vector', 'numel', 3});
    
    % 计算三维方向增量
    dx = end_point(1) - start_point(1);
    dy = end_point(2) - start_point(2);
    dz = end_point(3) - start_point(3);
    
    % 计算绝对位移（避免除零错误）
    abs_dx = abs(dx);
    abs_dy = abs(dy);
    abs_dz = abs(dz);
    total_disp = abs_dx + abs_dy + abs_dz;
    
    % 处理零位移情况
    if total_disp == 0
        path = repmat(start_point, num_points+2, 1);
        return;
    end
    
    % 基于位移比例分配点数
    total_points = num_points + 2;
    y_points = max(2, ceil(total_points * abs_dy / total_disp));
    x_points = max(2, ceil(total_points * abs_dx / total_disp));
    z_points = max(2, ceil(total_points * abs_dz / total_disp));
    
    % 路径方向决策逻辑
    if dy ~= 0  % y方向有变化
        % 先沿y方向移动（无论增加还是减少）
        y_values = linspace(start_point(2), end_point(2), y_points);
        z_values_y = linspace(start_point(3), end_point(3), y_points);
        y_path = [start_point(1)*ones(1,y_points); y_values; z_values_y]';
        
        % 再沿x方向移动
        x_values = linspace(start_point(1), end_point(1), x_points);
        z_values_x = linspace(y_path(end,3), end_point(3), x_points);
        x_path = [x_values; end_point(2)*ones(1,x_points); z_values_x]';
        
        % 合并路径并删除重复点
        path = [y_path; x_path(2:end,:)];
    else % y方向无变化
        % 直接沿x方向移动
        x_values = linspace(start_point(1), end_point(1), total_points);
        z_values = linspace(start_point(3), end_point(3), total_points);
        path = [x_values' repmat(start_point(2), total_points, 1) z_values'];
    end
    
    % 确保输出点数正确
    if size(path, 1) > total_points
        path = path(1:total_points, :);
    elseif size(path, 1) < total_points
        % 线性插值补充点数
        final_path = interp1(1:size(path,1), path, linspace(1, size(path,1), total_points));
        path = final_path;
    end
end
% 直线
function path = generatePath(start_point, end_point, num_points)
    num_points =num_points +1;
    % 验证输入
    validateattributes(start_point, {'numeric'}, {'vector', 'numel', 3});
    validateattributes(end_point, {'numeric'}, {'vector', 'numel', 3});
    validateattributes(num_points, {'numeric'}, {'scalar', 'integer', '>=', 2});

    % 强制转为行向量
    start_point = reshape(start_point, 1, 3);
    end_point = reshape(end_point, 1, 3);

    % 插值比例 t（num_points × 1）
    t = linspace(0, 1, num_points)';

    % 将起点和终点都扩展为 num_points × 3 大小
    start_mat = repmat(start_point, num_points, 1);    % 每一行是起点
    delta_mat = repmat(end_point - start_point, num_points, 1);  % 差值向量

    % 路径计算
    path = start_mat + delta_mat .* t;  % 广播兼容：num_points × 3
end

% 调整聚类函数：redistribute_points
function idx = redistribute_points(GT_positions, idx, C, cluster_idx, max_size)
    % 重新分配聚类中的多余目标点
    % 获取当前聚类的所有目标点
    cluster_points = GT_positions(idx == cluster_idx, :);
    
    % 找到多余的目标点
    extra_points = cluster_points(max_size+1:end, :);
    
    % 遍历多余目标点，重新分配到其他聚类
    for i = 1:size(extra_points, 1)
        % 计算每个目标点到各个簇中心的距离
        distances = sqrt(sum((C - extra_points(i, :)).^2, 2));  % 计算到各个聚类中心的距离
        [~, new_cluster] = min(distances);  % 找到最近的聚类
        
        % 将目标点分配给最近的聚类
        idx(find(ismember(GT_positions, extra_points(i, :), 'rows'))) = new_cluster;
    end
end