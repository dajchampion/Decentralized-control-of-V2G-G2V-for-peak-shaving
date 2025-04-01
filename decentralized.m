clc; clear; close all;
rng(0);

%% 参数设定
numEV = 20; T = 24; Pmax = 7; Pmin = -7;
baseLoad = [150, 140, 130, 125, 120, 130, 180, 250, 350, 400, ...
            420, 450, 470, 460, 440, 430, 450, 500, 550, 520, ...
            480, 400, 300, 200];
lambda = 0.1 + 0.002 * baseLoad;
% 扩展 lambda 为 20x24 的矩阵
lambda = repmat(lambda, numEV, 1);  % 将 1x24 的 lambda 向量重复 20 次

SOC_max = 50 * ones(numEV, 1); 
SOC_max_expanded=repmat(SOC_max, 1, T); % 20x24

SOC_min = 10 * ones(numEV, 1);
SOC_min_expanded=repmat(SOC_min, 1, T); % 20x24

SOC_init = randi([15, 40], numEV, 1);
SOC_init_expanded = repmat(SOC_init, 1, T);  % 将 SOC_init 重复 24 次，形成 20x24 矩阵

SOC_target = randi([30, 45], numEV, 1);

EV_arrival = randi([1 12], numEV, 1);
EV_departure = EV_arrival + randi([6 12], numEV, 1);

EV_status = zeros(numEV, T);
for i = 1:numEV
    EV_status(i, EV_arrival(i):min(EV_departure(i), T)) = 1;
end

%% 对 EV_arrival 进行升序排序，并重新排列所有 EV 相关信息
[EV_arrival, sortIdx] = sort(EV_arrival);

SOC_max = SOC_max(sortIdx);
SOC_max_expanded = SOC_max_expanded(sortIdx, :);

SOC_min = SOC_min(sortIdx);
SOC_min_expanded = SOC_min_expanded(sortIdx, :);

SOC_init = SOC_init(sortIdx);
SOC_init_expanded = SOC_init_expanded(sortIdx, :);

SOC_target = SOC_target(sortIdx);

EV_departure = EV_departure(sortIdx);
EV_status = EV_status(sortIdx, :);

%%  获取唯一的到达时间及其分组索引
[unique_arrival_times, ~, groupIdx] = unique(EV_arrival, 'stable');

% 按到达时间分组
numGroups = length(unique_arrival_times);
EV_groups = cell(numGroups, 1); % 用 cell 存储不同组的索引

for g = 1:numGroups
    EV_groups{g} = find(groupIdx == g); % 获取属于该组的 EV 索引
end

% 显示分组结果
disp('EV 分组结果:');
for g = 1:numGroups
    fprintf('到达时间 %d: EV 索引 -> ', unique_arrival_times(g));
    disp(EV_groups{g}');
end



%% 生成约束条件中的上三角矩阵
n=T;
A=ones(n);
A=triu(A);

%% decentralized control (updated baseload, first come first optimize)

numGroups = length(EV_groups); % 组数
P_ev_c=[];
SOC_c=[];
for g = 1:numGroups
    EV_idx = EV_groups{g}; % 当前组的EV索引
    fprintf('优化第 %d 组，EV 数量: %d\n', g, length(EV_idx));

    % 提取当前组的EV信息
    SOC_init_expanded_group = SOC_init_expanded(EV_idx,:);
    SOC_min_expanded_group = SOC_min_expanded(EV_idx,:);
    SOC_max_expanded_group = SOC_max_expanded(EV_idx,:);
    EV_status_group = EV_status(EV_idx, :);
    numEV_g=length(EV_idx);

    % 执行凸优化（示例：使用quadprog）

    cvx_begin
    variables P_ev(numEV_g, T)
    variable SOC(numEV_g, T)

    % 目标函数：最小化充电成本和负荷波动
    minimize(sum(sum(lambda(1:numEV_g,:) .* P_ev)) + sum(sum((sum(P_ev, 1) + baseLoad - mean(baseLoad)).^2)));

    % 约束条件
    subject to
        Pmin * EV_status_group <= P_ev <= Pmax * EV_status_group;
        % SOC(:,1) == SOC_init + P_ev(:,1);
        % for t = 2:T
        %     SOC(:,t) == SOC(:,t-1) + P_ev(:,t);
        % end
        SOC == SOC_init_expanded_group + P_ev * A;
        SOC_min_expanded_group <= SOC <= SOC_max_expanded_group;
        SOC(:, EV_departure(i)) >= SOC_target(i);
    cvx_end
    % 处理优化结果
    fval = cvx_optval;
    fprintf('第 %d 组优化完成，目标值: %f\n', g, fval);

    baseLoad = baseLoad+sum(P_ev, 1);
    P_ev_c=[P_ev_c;P_ev];
    SOC_c=[SOC_c;SOC];
     
end

save P_ev_c.txt P_ev_c -ascii;
%% 绘图
load P_ev.txt;
load P_ev_a.txt;
load P_ev_b.txt;
load P_ev_c.txt;

colors = lines(7); 
figure;
subplot(3,1,1);
plot(1:T, sum(P_ev_c,1) + baseLoad, 'Color', colors(1,:), 'LineWidth', 2); hold on;  % results of decentralized method, groups
plot(1:T, sum(P_ev_b,1) + baseLoad, 'Color', colors(2,:), 'LineWidth', 2); hold on;  % results of decentralized method, first come first optimize
plot(1:T, sum(P_ev_a,1) + baseLoad, 'Color', colors(3,:), 'LineWidth', 2); hold on; % results of decentralized method ignoring the order of precedence
plot(1:T, sum(P_ev,1) + baseLoad, 'Color', colors(4,:), 'LineWidth', 2);hold on;  % results of centralized method
plot(1:T, baseLoad, 'k', 'LineWidth', 2); 

xlabel('时间 (小时)'); ylabel('负荷 (kW)');
title('优化后总负荷'); legend('计及分组及时序的分散优化后','计及时序的分散优化后','分散优化后', '集中优化后', '基准负荷');

subplot(3,1,2);
imagesc(P_ev_c); colormap(jet); colorbar;
xlabel('时间 (小时)'); ylabel('EV编号');
title('EV 充放电策略');

subplot(3,1,3);
plot(1:T, SOC_c', 'k', 'LineWidth', 2);
xlabel('时间（小时）'); ylabel('电池能量');
title('EV能量变化');



