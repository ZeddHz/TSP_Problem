function f = TSP_distance(D, solution)
% 本函数计算给定solution的距离，其中城市之间的距离由D给出。
    n = numel(solution);
    sum = 0;
    for i = 1:n-1
        sum = sum + D(solution(i), solution(i+1));
    end
    sum = sum + D(solution(n), solution(1));
    f = sum;