% fixed_path_tsp_gurobi.m
% 固定起点和终点的路径TSP (MTZ模型, Gurobi MATLAB API)

function [path] = TSP(GT_positions)


%% 参数
n = length(GT_positions);             % 节点数量
disp(n);
V = 1:n;
s = 1;              % 起点
t = n;              % 终点

GT_positions = GT_positions;
Xcoord = GT_positions(:,1); % 提取所有行的第1列
Ycoord = GT_positions(:,2); % 提取所有行的第2列

%% 水流参数设置
water_speed_max = 0; % 最大水流速度 (m/s)  2
flow_direction = 1;    % 水流方向（1表示正方向，-1表示反方向）
U0 = flow_direction * 0.8 * water_speed_max;  % 主流强度
A = 40.6 * water_speed_max;                   % 波动幅值
kx = 0.06;                                    % x方向波数
ky = 0.03;                                    % y方向波数
phi = 0;                                      % 相位
v_uav = 5;                                    % UAV速度 (m/s)
rho_usv = 1000; % 水流密度
C_usv = 0.001; % 阻力因子
A_usv = 25; % 船的面积
alpha = 0.5*rho_usv*C_usv*A_usv; % USV速度惩罚系数

% 定义水流速度函数
water_vel_func = @(x,y) [
   ( U0 - A * ky * sin(kx * x + phi) .* sin(ky * y) );   % X方向速度
   ( -A * kx * cos(kx * x + phi) .* cos(ky * y) )       % Y方向速度
];

%% 计算代价矩阵（考虑水流影响，累加每个小段）
c = zeros(n,n);
for i=1:n
    for j=1:n
        if i~=j
            % 计算位移向量
            dx = Xcoord(j)-Xcoord(i);
            dy = Ycoord(j)-Ycoord(i);
            distance = sqrt(dx^2 + dy^2);
            
            % 计算航向向量（单位向量）
            heading_vector = [dx, dy]/distance;
            
            % 根据分辨率计算分段数
            resolution = 3;
            num_segments = max(ceil(distance / resolution), 1);

            % 初始化总能耗
            total_energy = 0;
            
            % 将路径分割成多个小段
            for seg = 1:num_segments
                % 计算当前小段的起点和终点
                start_frac = (seg-1)/num_segments;
                end_frac = seg/num_segments;
                
                start_x = Xcoord(i) + start_frac * dx;
                start_y = Ycoord(i) + start_frac * dy;
                end_x = Xcoord(i) + end_frac * dx;
                end_y = Ycoord(i) + end_frac * dy;
                
                % 计算小段中点（用于水流速度计算）
                mid_x = (start_x + end_x)/2;
                mid_y = (start_y + end_y)/2;
                
                % 获取中点处的水流速度
                water_vel = water_vel_func(mid_x, mid_y);
                
                % 计算小段位移向量
                seg_dx = end_x - start_x;
                seg_dy = end_y - start_y;
                seg_distance = sqrt(seg_dx^2 + seg_dy^2);
                
                % 计算UAV对地速度向量
                v_uav_ground = v_uav * heading_vector;
                
                % 计算相对速度向量（UAV相对于水的速度）
                v_relative = v_uav_ground - water_vel';
                
                % 计算相对速度大小
                v_rel_mag = norm(v_relative);
                
                % 计算小段航行时间（秒）
                time_seconds = seg_distance / norm(v_uav_ground);

                d_0 = 0.6;
                P0 = 80;
                P1 = 88.6;
                v0 =4.03;
                Acanshu = 0.503;
                s_uavcanshu = 0.05;
                rou = 1.225;
                oumu = 300;
                r_uav = 0.4;
                v_uav1 = 3;
                lamda_Bw= (  (1 +(v_uav1)^4/(1*(v0)^4))^(1/2)-1/2*(  ( v_uav1^2)/(1*(v0^2))  )  )^(1/2);
                thm1 = P0*(1+( 3/(oumu*r_uav)^2)*(v_uav1^2 ) );
                thm2 = 0.5*d_0*rou*s_uavcanshu*Acanshu*(v_uav1^3);
                total_energyw = thm1 + thm2+P1*lamda_Bw;
                % 累加小段能耗（与时间成正比，与相对速度立方成正比）
                total_energy = total_energy + time_seconds*alpha * (v_rel_mag^3) + time_seconds*total_energyw;
                % 累加小段能耗（与时间成正比，与相对速度立方成正比）
                % = total_energy + 1 * (v_rel_mag^2);
            end
            
            c(i,j) = total_energy;
        end
    end
end





%% 索引映射
Map = zeros(n);
cnt = 0;
for i=1:n
    for j=1:n
        if i~=j
            cnt = cnt + 1;
            Map(i,j) = cnt;
        end
    end
end
num_x = n*(n-1);
num_u = n;
N = num_x + num_u;
disp(Map);
%% 初始化 Gurobi 模型
model.modelname = 'PathTSP';
model.modelsense = 'min';

%% 变量定义
% x(i,j)二进制变量
for i=1:n
    for j=1:n
        if i~=j
            idx = Map(i,j);
            model.varnames{idx} = sprintf('x_%d_%d', i, j);
            model.vtype(idx) = 'B'; % 二进制
            model.obj(idx) = c(i,j);  % 目标函数系数为距离
            model.lb(idx) = 0;
            model.ub(idx) = 1;
        end
    end
end

% u(i)连续变量
for i=1:n
    idx = num_x + i;
    model.varnames{idx} = sprintf('u_%d', i);
    model.vtype(idx) = 'C'; % 连续
    model.obj(idx) = 0;
    if i==s
        model.lb(idx) = 0;
        model.ub(idx) = 0;  % u_s=0
    else
        model.lb(idx) = 0;
        model.ub(idx) = n-1;
    end
end

%% 约束矩阵
A = sparse(0,N);
rhs = [];
sense = [];
rownames = {};

% (1) 起点约束
row = sparse(1,N);
for j=1:n
    if j~=s
        row(1, Map(s,j)) = 1;
    end
end
A(end+1,:) = row;
rhs(end+1) = 1;
sense(end+1) = '=';
rownames{end+1} = 'StartOut';

row = sparse(1,N);
for i=1:n
    if i~=s
        row(1, Map(i,s)) = 1;
    end
end
A(end+1,:) = row;
rhs(end+1) = 0;
sense(end+1) = '=';
rownames{end+1} = 'StartIn';

% (2) 终点约束
row = sparse(1,N);
for i=1:n
    if i~=t
        row(1, Map(i,t)) = 1;
    end
end
A(end+1,:) = row;
rhs(end+1) = 1;
sense(end+1) = '=';
rownames{end+1} = 'EndIn';

row = sparse(1,N);
for j=1:n
    if j~=t
        row(1, Map(t,j)) = 1;
    end
end
A(end+1,:) = row;
rhs(end+1) = 0;
sense(end+1) = '=';
rownames{end+1} = 'EndOut';

% (3) 中间节点约束
for k=1:n
    if k~=s && k~=t
        % 入度=1
        row_in = sparse(1,N);
        for i=1:n
            if i~=k
                row_in(1, Map(i,k)) = 1;
            end
        end
        A(end+1,:) = row_in;
        rhs(end+1) = 1;
        sense(end+1) = '=';
        rownames{end+1} = sprintf('In_%d',k);

        % 出度=1
        row_out = sparse(1,N);
        for j=1:n
            if j~=k
                row_out(1, Map(k,j)) = 1;
            end
        end
        A(end+1,:) = row_out;
        rhs(end+1) = 1;
        sense(end+1) = '=';
        rownames{end+1} = sprintf('Out_%d',k);
    end
end

% (4) MTZ约束
for i=1:n
    for j=1:n
        if i~=j && i~=t && j~=s
            row = sparse(1,N);
            row(1, num_x+i) = 1;
            row(1, num_x+j) = -1;
            row(1, Map(i,j)) = n-1;
            A(end+1,:) = row;
            rhs(end+1) = n-2;
            sense(end+1) = '<';
            rownames{end+1} = sprintf('MTZ_%d_%d',i,j);
        end
    end
end

%% 加入模型
model.A = A;
model.rhs = rhs;
model.sense = char(sense);   % ⭐修正: 转成char array
model.rownames = rownames;

%% 求解
params.outputflag = 1;
result = gurobi(model, params);

%% 输出结果
%fprintf('\n最小总代价: %.2f\n', result.objval);
%fprintf('路径:\n');

% 重建X矩阵
X = zeros(n,n);
for i=1:n
    for j=1:n
        if i~=j
            idx = Map(i,j);
            if result.x(idx)>0.5
                X(i,j)=1;
            end
        end
    end
end

% 输出路径
current = s;
path = current;
while current ~= t
    next = find(X(current,:)==1);
    path = [path, next];
    current = next;
end
%disp(path);



    % 重建路径
    X = zeros(n,n);
    for i=1:n
        for j=1:n
            if i~=j
                idx = Map(i,j);
                if result.x(idx) > 0.5
                    X(i,j) = 1;
                end
            end
        end
    end

    % 构建路径序列
    current = s;
    path = current;
    while current ~= t
        next = find(X(current,:) == 1);
        path = [path, next];
        current = next;
    end

    % 返回按路径排序后的坐标点集合
    path = GT_positions(path, :);

end
