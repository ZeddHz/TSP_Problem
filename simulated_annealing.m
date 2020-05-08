function f = simulated_annealing(D, mute, MAX_ITER)
% 模拟退火解决中国省会城市TSP问题的简单实现
% 包含台北，香港，澳门

    %% 默认参数设置
    clear;
    if ~exist('mute', 'var')
        mute = 0; % 是否显示各种提示信息
    end
    if ~exist('MAX_ITER', 'var')
        MAX_ITER = 50; % 每个温度内的最大状态交换尝试次数
                       % 这个参数很玄乎没什么固定规律的样子
    end
    if ~exist('D', 'var')
        data=xlsread('省会经纬坐标.xlsx');
        C=data;
        D=zeros(35,35);% D为城市间的距离矩阵
        % 形成两两之间对应的矩阵（对称阵，可以只看上三角或下三角）
        [LA1,LA2]=meshgrid(C(:,2));
        [LO1,LO2]=meshgrid(C(:,1));

        % 计算两两之间的距离，单位为公里
        R = distance(LA1,LO1,LA2,LO2,almanac('earth','wgs84'));
        D = num2str(R,'%10.2f');
        %disp(D);
    end
    
     %% 模拟退火算法参数设置
    rng(0);
    n = size(D, 1);
    
    T_range_factor = exp(0:-0.1:-5); % 温度的范围系数

    solution = [1, randperm(n-1) + 1]; % 生成一个解，假定从1开始
    if ~mute
        fprintf('邻域内总共有%d个解。\n', numel(P));
    end
    
    f = TSP_distance(D, solution); % 求出当前解的总距离
    
    if ~mute
        disp('初始的路径为：');   
        disp(solution);
        disp('路径长度为：');
        disp(f);
    end
    
    P = generate_neighbors(solution); % 生成初始解的邻域P
    
    TMAX = f;
    % 温度范围动态地随初始解的好坏而变化（初始解越差，要求探索的需求就越大）
    T_range = T_range_factor * TMAX; % 尽量保证开始的温度足够高 结束的足够低
                                     
    first = zeros(1, MAX_ITER);
    final = zeros(1, MAX_ITER);
    for t = T_range  % 外循环
        for i = 1:MAX_ITER  % 内循环
            index = randi(numel(P), 1); % 产生一个1~|P|之间的随机整数
            neighbor = P{index}; % 在P中随机取一个解
            f_neighbor = TSP_distance(D, neighbor); % 计算这个解的距离
            Pt = exp(-(f_neighbor - f)/t); % 转移概率
            if Pt > 1, Pt = 1; end % 说明f_neighbor<f
            if Pt > rand
                f = f_neighbor;
                solution = neighbor;
                P = generate_neighbors(solution);
            end
            % 记录外循环第一轮的转移概率
            if ~mute && t == T_range(1)
                first(i) = Pt;
            end
            % 记录外循环最后一轮的转移概率
            if ~mute && t == T_range(end)
                final(i) = Pt;
            end
        end
        if ~mute && t == T_range(1)
            fprintf('初始温度下的转移概率中位数为%f\n', median(first));
        end
        if ~mute && t == T_range(end)
            fprintf('最后温度下的转移概率中位数为%f\n', median(final));
        end
    end
    if ~mute
    fprintf('最后搜索得到的最优路径为：\n');
    disp(solution);
    disp('路径长度为：');
    disp(f);
    end

