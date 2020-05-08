function f = genetic_algorithm(D, mute, Pm)
% �Ŵ��㷨����й�ʡ�����TSP����ļ�ʵ��
% ����̨������ۣ�����

    %% Ĭ�ϲ�������
    clear;
    if ~exist('mute', 'var')
        mute = 0; % �Ƿ���ʾ������ʾ��Ϣ
    end
    if ~exist('Pm', 'var')
        Pm = 0.2; % ������ʣ�Խ������Խ�����ǽ�һ��Խ��
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
    
    %% �Ŵ��㷨��������
    rng(1); % ������������ɣ����ӣ�
    n = size(D, 1);
    N = 100; % Ⱥ���ģ��������Ϊż����
    TOL = 20; % ������̴���(����TOL��rate�����������Ҳ������Ž⣬��ֹͣ����)
    
    %% ��ʼ���н�
    solutions = zeros(N, n); % �������һ������Ĭ���Ƿ��ؽ�ĵ�һ�����У��ʿ��н���û�б�Ҫ���������ʵĳ���
    fs = zeros(N, 1);   % ��Ӧ��
    for i = 1:N
        solutions(i, :) = [1, randperm(n-1) + 1]; % ����N���⣬�ٶ���1��ʼ
        fs(i) = TSP_distance(D, solutions(i, :));
    end
    Pu = max(fs) - fs + 1; % ��1Ϊ�˷�ֹmax��fs����Ӧ����һλΪ0
    P = Pu/sum(Pu); % �ɴˣ���Ӧ��Խ�����һ�����ѡȡ�ĸ���ԽС
    cumP = cumsum(P); % �����ۻ�����
    best = min(fs); % �������·��
    avg = mean(fs);
    rate = best/avg;
    if ~mute
    disp('��ʼ���Ⱥ������̵�·������Ϊ��');
    disp(best);
    disp('��ʼ���Ⱥ����ƽ��·������Ϊ��');
    disp(avg);
    end
    
    %% ����
    tol = 0;
    count = 0;
    while 1
        count = count + 1;
        if ~mute
            fprintf('��ǰ��%d�ε���\n', count);
        end
        
         % ʹ�����̶ĵķ�ʽѡ��������Ⱦɫ��
        parents = zeros(size(solutions));
        for i = 1:N
            index = sum(cumP <= rand) + 1; % �����rand���ڵ�������
            parents(i, :) = solutions(index, :);
        end
        
        new_solutions = zeros(size(solutions));
        assert(mod(N, 2) == 0); % N����Ϊż�����������
        
        % �������������Ĭ��NΪż����ÿ��������һ����������Ӵ�
        for i = 1:N/2 
            % �������Ӵ�1ȡ����1��ǰһ��Ⱦɫ�壬��һ�����ɸ���2�ṩ��ͬ�����Ӵ�2
            p1 = parents(2*i-1, :);
            p2 = parents(2*i, :);
            middle = ceil(n/2);
            s1 = p1(1:middle);
            res1 = setdiff(p2, s1, 'stable'); % �ҵ�p2�к��еĶ�s1���к��еĻ����൱��p2\s1
            s1 = [s1, res1];
            s2 = p2(1:middle);
            res2 = setdiff(p1, s2, 'stable');
            s2 = [s2, res2];
            new_solutions(2*i-1, :) = s1;
            new_solutions(2*i, :) = s2;
        end
        
        % �������������ķ�ʽΪ���ȡ���������һ�����н���
        for i = 1:N 
            if rand < Pm
                temp = randperm(n-1) + 1;
                k = temp(1);
                new_solutions(i, [1, k]) = new_solutions(i, [k, 1]);
            end
        end
        
        % ���ˣ��µ���Ⱥ�Ѿ�������ϣ��������Ⱥ��ʼ��һ�ֵļ���
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
            disp('��̵�·������Ϊ��');
            disp(best_new);
            disp('ƽ����·������Ϊ��');
            disp(avg);
        end
        tol = tol + 1;
        % ������ָ��£����Ž��б�ã���tol��0
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
    
    %% ������Ž�
    [f, index] = min(fs);
    solution = solutions(index, :);
    if ~mute
    fprintf('��������õ�������·��Ϊ��\n');
    disp(solution); % ���Ҫ�ص����1
    disp('·������Ϊ��');
    disp(f);
    end

