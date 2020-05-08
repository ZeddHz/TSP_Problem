function P = generate_neighbors(solution)
    % 本函数根据 交换任意两个城市 的原则，产生模拟退火算法中solution的邻域
    n = numel(solution);
    P = cell(1, n*(n-1)/2); % 使用元组来存放solution的领域
    t = 1;
    for i = 1:n
        for j = i+1:n
            solution_new = solution;
            solution_new([i, j]) = solution_new([j, i]); % 交换两个元素
            P{t} = solution_new;
            t = t + 1;
        end
    end