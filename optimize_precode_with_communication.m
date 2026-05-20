function [power,rate] = optimize_precode_with_communication(Q,uav_pos,usv_pos,gamma_vector)
%% Simulation parameters
N              = Q;       % Number of transmit antennas
K              = 1;       % Number of users
%gamma_vector   = 15 % range of SINR requirement (in dB)
sigma_2        = 1e-14;       % noise power
ch_realization = 1;     % channel realizations
%disp(gamma_vector);
%% Generate Channles - Normalized Rayleigh fading channel
xzx =generate_fixed_position_channel_c(Q,uav_pos,usv_pos);
H = repmat(xzx, [1, K, ch_realization]); % 正确扩展为4×K×ch_realization
%% Variables for save the results 
P_sdp = zeros(ch_realization, length(gamma_vector));
rate_sdp = zeros(ch_realization, length(gamma_vector));
    %% Optimization
    for i=1:ch_realization
        B =1;
        h = H(:,:,i);
        for j=1:length(gamma_vector)
            gamma = (2^(gamma_vector/B)-1);
            [minimun_power_sdp,g] = sdp1(h,gamma,sigma_2,N,K);
            P_sdp(i,j) = minimun_power_sdp;
            rate_sdp(i, j) = B*log2(1+g/ sigma_2);
        end
    end
power = (mean(P_sdp,1));
rate = mean(rate_sdp, 1);          % Return average rate for each gamma value
%disp(rate);
end

function [minimun_power_sdp,g1] = sdp1(h,gamma,sigma_2,N,K)
% 利用SDP求解功率最小化问题
    cvx_precision high
    cvx_begin quiet
    variable B(N,N,K) hermitian
    variable x(K,1)
    variables xi(K)                % 功率松弛变量
    variables g1(K)                % 功率松弛变量
    obj = 0;
    for k=1:K
        obj = obj + trace(B(:,:,k))+10000*sum(x)+10000*sum(xi);
    end
    minimize(obj)
    subject to
        for i=1:K
            c = 0;
            for j=1:K
                if j~=i
                    c = c + real(trace(h(:,i)*h(:,i)'*B(:,:,j)));
                end
            end
            real(trace(h(:,i)*h(:,i)'*B(:,:,i))) -x(i,1) >= gamma*sigma_2;
            g1(i)  == real(trace(h(:,i)*h(:,i)'*B(:,:,i)));
            x(i,1) >= 0;
            xi(i) >=0;
            trace(B(:,:,i)) <= 20+xi(i);  % 每个波束成形向量的功率不超过30
            trace(B(:,:,i)) >= 0;  % 每个波束成形向量的功率不超过30
            B(:,:,i) == hermitian_semidefinite(N);
        end
    cvx_end
    %B_full = full(B);  % 添加这行转换
    minimun_power_sdp = 0;
    %disp( trace(B) );
    for k=1:K
        %disp(B(:,:,k));
        minimun_power_sdp = minimun_power_sdp + trace( B(:,:,k) );
    end
    %disp(minimun_power_sdp);
end
