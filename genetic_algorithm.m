function f = genetic_algorithm(D, mute, Pm)
% 遗传算法解决中国省会城市TSP问题的简单实现
% 包含台北，香港，澳门

    %% 默认参数设置
    clear;
    if ~exist('mute', 'var')
        mute = 0; % 是否显示各种提示信息
    end
    if ~exist('Pm', 'var')
        Pm = 0.2; % 变异概率，越大收敛越慢但是解一般越好
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
    
    %% 遗传算法参数设置
    rng(1); % 控制随机数生成（种子）
    n = size(D, 1);
    N = 100; % 群体规模（必须设为偶数）
    TOL = 20; % 最大容忍次数(连续TOL次rate不上升，或找不到更优解，则停止迭代)
    
    %% 初始可行解
    solutions = zeros(N, n); % 由于最后一个城市默认是返回解的第一个城市，故可行解里没有必要包含最后访问的城市
    fs = zeros(N, 1);   % 适应度
    for i = 1:N
        solutions(i, :) = [1, randperm(n-1) + 1]; % 生成N个解，假定从1开始
        fs(i) = TSP_distance(D, solutions(i, :));
    end
    Pu = max(fs) - fs + 1; % 加1为了防止max（fs）对应的那一位为0
    P = Pu/sum(Pu); % 由此，适应度越大的那一项可能选取的概率越小
    cumP = cumsum(P); % 计算累积概率
    best = min(fs); % 计算最短路径
    avg = mean(fs);
    rate = best/avg;
    if ~mute
    disp('初始解的群体中最短的路径长度为：');
    disp(best);
    disp('初始解的群体中平均路径长度为：');
    disp(avg);
    end
    
    %% 进化
    tol = 0;
    count = 0;
    while 1
        count = count + 1;
        if ~mute
            fprintf('当前第%d次迭代\n', count);
        end
        
         % 使用轮盘赌的方式选出父代的染色体
        parents = zeros(size(solutions));
        for i = 1:N
            index = sum(cumP <= rand) + 1; % 随机数rand落在的区间编号
            parents(i, :) = solutions(index, :);
        end
        
        new_solutions = zeros(size(solutions));
        assert(mod(N, 2) == 0); % N必须为偶数，否则出错
        
        % 交配操作；这里默认N为偶数，每两个父代一起产生两个子代
        for i = 1:N/2 
            % 产生的子代1取父代1的前一半染色体，后一半则由父代2提供；同理于子代2
            p1 = parents(2*i-1, :);
            p2 = parents(2*i, :);
            middle = ceil(n/2);
            s1 = p1(1:middle);
            res1 = setdiff(p2, s1, 'stable'); % 找到p2中含有的而s1不中含有的基因，相当于p2\s1
            s1 = [s1, res1];
            s2 = p2(1:middle);
            res2 = setdiff(p1, s2, 'stable');
            s2 = [s2, res2];
            new_solutions(2*i-1, :) = s1;
            new_solutions(2*i, :) = s2;
        end
        
        % 变异操作；变异的方式为随机取个城市与第一个城市交换
        for i = 1:N 
            if rand < Pm
                temp = randperm(n-1) + 1;
                k = temp(1);
                new_solutions(i, [1, k]) = new_solutions(i, [k, 1]);
            end
        end
        
        % 至此，新的种群已经生成完毕，代替旧种群后开始新一轮的计算
        solutions = new_solutions;
        for i = 1:N
            fs(i) = TSP_distance(D, solutions(i, :));
        end
        Pu = max(fs) - fs + 1;
        P = Pu/sum(Pu);
        cumP = cumsum(P);
        best_new = min(fs);
        avg = mean(fs);
        rate_new = best_new/avg;
        if ~mute
            disp('最短的路径长度为：');
            disp(best_new);
            disp('平均的路径长度为：');
            disp(avg);
        end
        tol = tol + 1;
        % 如果本轮更新，最优解有变好，则tol归0
        if best_new < best || rate_new > rate  
            best = best_new;
            rate = rate_new;
            tol = 0;
        end
        if tol >= TOL
            break
        end
        if count > 5000
            break
        end
    end 
    
    %% 输出最优解
    [f, index] = min(fs);
    solution = solutions(index, :);
    if ~mute
    fprintf('最后搜索得到的最优路径为：\n');
    disp(solution); % 最后还要回到起点1
    disp('路径长度为：');
    disp(f);
    end

