function P = generate_neighbors(solution)
    % ���������� ���������������� ��ԭ�򣬲���ģ���˻��㷨��solution������
    n = numel(solution);
    P = cell(1, n*(n-1)/2); % ʹ��Ԫ�������solution������
    t = 1;
    for i = 1:n
        for j = i+1:n
            solution_new = solution;
            solution_new([i, j]) = solution_new([j, i]); % ��������Ԫ��
            P{t} = solution_new;
            t = t + 1;
        end
    end