clear, clc

%delta_T = [0.14, 0.12, 0.10, 0.08, 0.06, 0.04, 0.02];
delta_T = [0.10];

for num=1:5
    %delta_t = delta_T(num);
    delta_t = delta_T(1);
    '初始化种群---------------';
    %'适应度函数定义';
    %f= @(x)x .* sin(x) .* cos(2 * x) - 2 * x .* sin(3 * x); % 函数表达式（适应度函数）
    N = 100;                         % 种群个数
    %d = 1;                          % 空间维数
    d = 2;                          % 空间维数
    ger = 100;                      % 最大迭代次数
    %limit = [0, 100000];                % 设置位置参数限制  d x 2 (d维，上下限)
    limit = [0, 100000; 0, 100000];                % 设置位置参数限制  d x 2 (d维，上下限)
    %vlimit = [-1, 1];               % 设置速度限制  d x 2 (d维，上下限)
    vlimit = [-1, 1; -1, 1];               % 设置速度限制  d x 2 (d维，上下限)
    w = 0.8;                        % 惯性权重
    c1 = 0.5;                       % 自我学习因子
    c2 = 0.5;
    % 群体学习因子
    for i = 1:d
        %x(:,i) = limit(i, 1) + (limit(i, 2) - limit(i, 1)) * rand(N, 1);%初始种群的位置
        x(:,i) = linspace(0, 100000, N)';%初始种群的位置
    end
    v = rand(N, d);                  % 初始种群的速度
    xm = x;                          % 每个个体的历史最佳位置
    ym = zeros(1, d);                % 种群的历史最佳位置
    fxm = zeros(N, 1);               % 每个个体的历史最佳适应度
    fym = -inf;                      % 种群历史最佳适应度（寻max则初始化为-inf，寻min则初始化为inf）
    '-------------------------';

    '群体更新';
    record = zeros(ger, 1);          % 记录器(记录每次迭代的群体最佳适应度)
    for iter = 1:ger
        fx = zeros(N, 1);
        for i = 1:N
            %fx(i) = P(x(i), 0, delta_t); % 个体当前适应度
            fx(i) = P4(x(i, 1), x(i, 2), delta_t); % 个体当前适应度
        end
        '更新个体历史最佳适应度和个体历史最佳位置';
        for i = 1:N      
            if fxm(i) < fx(i) % < / > -----------------
                fxm(i) = fx(i);     % 更新个体历史最佳适应度
                xm(i,:) = x(i,:);   % 更新个体历史最佳位置
            end 
        end

        '更新群体历史最佳适应度和群体历史最佳位置';
        if fym < max(fxm) % < / > ----------------------------
            [fym, nmax] = max(fxm);   % 更新群体历史最佳适应度---
            ym = xm(nmax, :);      % 更新群体历史最佳位置
        end

        '速度更新';
        v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x);% 速度更新
        % 边界速度处理-----------------------------------------
        new_v = zeros(N, d);
        for i = 1:d
            vi = v(:,i);
            vi(vi > vlimit(i,2)) = vlimit(i,2);
            vi(vi < vlimit(i,1)) = vlimit(i,1);
            new_v(:,i) = vi;
        end
        v = new_v;

        '位置更新';
        x = x + v; % 位置更新
        % 边界位置处理-------------------------------------------
        new_x = zeros(N, d);
        for i = 1:d
            xi = x(:,i);
            xi(xi > limit(i,2)) = limit(i,2);
            xi(xi < limit(i,1)) = limit(i,1);
            new_x(:,i) = xi;
        end
        x = new_x;
        record(iter) = fym; % 最佳适应度记录
    end

    figure(3);plot(record);title('收敛过程')
    % f = fopen('2_1.csv', 'a');
    f = fopen('4_true.csv', 'a');
    fprintf(f,'delta_t:%.2f, 最大值:%.4f, 变量取值:%.4f, %.2f\n', delta_t, fym, ym(1), ym(2));
    fclose(f);
    disp(['最大值：',num2str(fym)]);
    disp(['变量取值：',num2str(ym(1)), ' ', num2str(ym(2))]);
end