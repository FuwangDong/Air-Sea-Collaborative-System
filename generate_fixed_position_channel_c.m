function h_channel = generate_fixed_position_channel_c(Q,uav_pos,usv_pos)
    f_c = 28e9;          % 载频 28GHz (来自图片描述)
    c = 3e8;             % 光速
    lambda = c / f_c;    % 波长
    d = lambda/2;        % 天线间距 (半波长)
    K_c = 0.00001;            % Rician K因子 (图片中设定为常数)
    g_c = 0.00001; % CN(0,1)

    dx = uav_pos(1) - usv_pos(1);
    dy = uav_pos(2) - usv_pos(2);
    dz = uav_pos(3) - usv_pos(3);
    d_c = sqrt(dx^2 + dy^2 + dz^2);
    phi = acos(dz / d_c);  % 从图片中推导: φ = arccos((H-0)/d_c)
    array_vector = zeros(Q, 1);
    for q = 0:Q-1
        array_vector(q+1) = exp(1j * 2 * pi * q * d * cos(phi) / lambda);
    end

    K_term = sqrt(K_c/(K_c+1)) + sqrt(1/(K_c+1)) * g_c;
    rho0 = 0.03;
    Upsilon = K_term * sqrt(rho0) / d_c;
    h_channel = Upsilon* array_vector;
end