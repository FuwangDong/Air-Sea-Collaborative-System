function h_channel = generate_fixed_position_channel_s(Q,uav_pos,usv_pos)
    % 参数设置
    neta = 0.1;
    beta = 0.03;         % 1m时候的感知信道增益 14.8dBm
    f_c = 28e9;          % 载频 28GHz (来自图片描述)
    c = 3e8;             % 光速
    lambda = c / f_c;    % 波长
    d = lambda/2;        % 天线间距 (半波长)

    dx = uav_pos(1) - usv_pos(1);
    dy = uav_pos(2) - usv_pos(2);
    dz = uav_pos(3) - usv_pos(3);
    d_c = sqrt(dx^2 + dy^2 + dz^2);
    phi = acos(dz / d_c);  
    array_vector = zeros(Q, 1);
    for q = 0:Q-1
        array_vector(q+1) = exp(1j * 2 * pi * q * d * cos(phi) / lambda);
    end
    Upsilon = 0.25*beta*((neta)^(0.5) )*(pi)^(-0.5);
    Upsilon =  Upsilon*(d_c^(-2) );
    % ========== 最终信道向量 ========== (公式2)
    h_channel = Q^(0.5)*Upsilon* array_vector;
end