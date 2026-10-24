function prob_out = compute_KDE(t, data, timewindow, tt)
% COMPUTE_KDE 计算 0-1 时间序列中每一时刻的高斯核密度估计概率
% 输入：
%   t          - 原始时间戳向量
%   data       - 离散数据向量
%   timewindow - 高斯窗口有效跨度 (与时间戳单位一致)
%   tt         - (可选) 目标时间戳向量，未传入时会在原时间戳 t 上进行计算
% 输出：
%   prob_out   - 计算出的每一时刻的概率 (对应传入的 tt 或 t)

    if nargin < 4
        tt = [];
    end

    % 根据时间戳计算平均采样间隔
    dt = mean(diff(t));
    if isnan(dt) || dt == 0
        warning('计算出的采样间隔异常，返回全 NaN。');
        if isempty(tt)
            prob_out = NaN(size(data));
        else
            prob_out = NaN(size(tt));
        end
        return;
    end
    
    % 标准差计算 (窗口跨度 / 6可以覆盖 99.7% 的高斯分布)
    sigma_width = timewindow / 6; 
    
    % 预处理数据：将 NaN 设为 0，将大于 0 的值设为 1
    data(isnan(data)) = 0;
    data(data > 0) = 1;

    if isempty(tt)
        % 未传入 tt，执行原本的 conv 计算
        % 构建高斯核的时间轴范围
        t_kernel = -timewindow/2 : dt : timewindow/2;
        
        % 生成高斯分布概率密度函数并归一化
        gaussian_kernel = exp(-(t_kernel.^2) / (2 * sigma_width^2));
        gaussian_kernel = gaussian_kernel / sum(gaussian_kernel);
        
        % 计算每一时刻的概率 (卷积计算高斯平滑)
        prob_out = conv(data, gaussian_kernel, 'same');
    else
        % 传入了 tt，在目标时间点上分别精确计算 KDE 
        % （与之前的 conv 采用完全一样的归一化常数与窗口截断）
        t_kernel = -timewindow/2 : dt : timewindow/2;
        Z = sum(exp(-(t_kernel.^2) / (2 * sigma_width^2)));
        
        prob_out = zeros(length(tt), 1);
        for i = 1:length(tt)
            % 限定截断窗口为 [-timewindow/2, timewindow/2]
            dist = t - tt(i);
            in_window = abs(dist) <= (timewindow/2 + 1e-9);
            
            % 提取有重叠部分的数据与对应的权重
            w = exp(-(dist(in_window).^2) / (2 * sigma_width^2));
            d_in = data(in_window);
            
            % 为确保不论行列排布都能正确点乘求和
            prob_out(i) = sum(d_in(:) .* w(:)) / Z;
        end
        
        % 保持与输入 tt 维度一致(如果是行向量则转置回去)
        if isrow(tt)
            prob_out = prob_out';
        end
    end
end
