function f = TSP_distance(D, solution)
% �������������solution�ľ��룬���г���֮��ľ�����D������
    n = numel(solution);
    sum = 0;
    for i = 1:n-1
        sum = sum + D(solution(i), solution(i+1));
    end
    sum = sum + D(solution(n), solution(1));
    f = sum;