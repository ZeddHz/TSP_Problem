function f = simulated_annealing(D, mute, MAX_ITER)
% ģ���˻����й�ʡ�����TSP����ļ�ʵ��
% ����̨������ۣ�����

    %% Ĭ�ϲ�������
    clear;
    if ~exist('mute', 'var')
        mute = 0; % �Ƿ���ʾ������ʾ��Ϣ
    end
    if ~exist('MAX_ITER', 'var')
        MAX_ITER = 50; % ÿ���¶��ڵ����״̬�������Դ���
                       % �������������ûʲô�̶����ɵ�����
    end
    if ~exist('D', 'var')
        data=xlsread('ʡ�ᾭγ����.xlsx');
        C=data;
        D=zeros(35,35);% DΪ���м�ľ������
        % �γ�����֮���Ӧ�ľ��󣨶Գ��󣬿���ֻ�������ǻ������ǣ�
        [LA1,LA2]=meshgrid(C(:,2));
        [LO1,LO2]=meshgrid(C(:,1));

        % ��������֮��ľ��룬��λΪ����
        R = distance(LA1,LO1,LA2,LO2,almanac('earth','wgs84'));
        D = num2str(R,'%10.2f');
        %disp(D);
    end
    
     %% ģ���˻��㷨��������
    rng(0);
    n = size(D, 1);
    
    T_range_factor = exp(0:-0.1:-5); % �¶ȵķ�Χϵ��

    solution = [1, randperm(n-1) + 1]; % ����һ���⣬�ٶ���1��ʼ
    if ~mute
        fprintf('�������ܹ���%d���⡣\n', numel(P));
    end
    
    f = TSP_distance(D, solution); % �����ǰ����ܾ���
    
    if ~mute
        disp('��ʼ��·��Ϊ��');   
        disp(solution);
        disp('·������Ϊ��');
        disp(f);
    end
    
    P = generate_neighbors(solution); % ���ɳ�ʼ�������P
    
    TMAX = f;
    % �¶ȷ�Χ��̬�����ʼ��ĺû����仯����ʼ��Խ�Ҫ��̽���������Խ��
    T_range = T_range_factor * TMAX; % ������֤��ʼ���¶��㹻�� �������㹻��
                                     
    first = zeros(1, MAX_ITER);
    final = zeros(1, MAX_ITER);
    for t = T_range  % ��ѭ��
        for i = 1:MAX_ITER  % ��ѭ��
            index = randi(numel(P), 1); % ����һ��1~|P|֮����������
            neighbor = P{index}; % ��P�����ȡһ����
            f_neighbor = TSP_distance(D, neighbor); % ���������ľ���
            Pt = exp(-(f_neighbor - f)/t); % ת�Ƹ���
            if Pt > 1, Pt = 1; end % ˵��f_neighbor<f
            if Pt > rand
                f = f_neighbor;
                solution = neighbor;
                P = generate_neighbors(solution);
            end
            % ��¼��ѭ����һ�ֵ�ת�Ƹ���
            if ~mute && t == T_range(1)
                first(i) = Pt;
            end
            % ��¼��ѭ�����һ�ֵ�ת�Ƹ���
            if ~mute && t == T_range(end)
                final(i) = Pt;
            end
        end
        if ~mute && t == T_range(1)
            fprintf('��ʼ�¶��µ�ת�Ƹ�����λ��Ϊ%f\n', median(first));
        end
        if ~mute && t == T_range(end)
            fprintf('����¶��µ�ת�Ƹ�����λ��Ϊ%f\n', median(final));
        end
    end
    if ~mute
    fprintf('��������õ�������·��Ϊ��\n');
    disp(solution);
    disp('·������Ϊ��');
    disp(f);
    end

