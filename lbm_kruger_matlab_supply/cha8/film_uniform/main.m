% 创建对象并运行模拟
tic
% filepath: /d:/BaiduSyncdisk/VScode/Books/lbm_kruger/lbm_principles_practice/code-master/chapter8/run_film_uniform.m
% 创建 FilmUniform 类的实例
film = FilmUniform();

% 运行模拟
film.run_simulation();

toc 
%
%%
   % 设置文件路径和文件名格式
    filePattern = 'film_uniformfilm_uniform%07d.dat';
    numFrames = 37; % 根据实际生成的文件数量设置
    NX = 160; % 根据实际情况设置
    NY = 1600; % 根据实际情况设置
    
    % 创建一个figure窗口
    figure;
    
    % 循环读取文件并绘制动画
    for k = 0:numFrames
        % 生成文件名
        filename = sprintf(filePattern, k * 100);
        
        % 检查文件是否存在
        if ~isfile(filename)
            continue;
        end
        
        % 读取文件中的数据
        data = load(filename);
        
        % 将数据重塑为矩阵
        rho = reshape(data, [NX, NY])';
        
        % 绘制密度场
        imagesc(rho);
        colorbar;
        title(sprintf('Density Field at Step %d', k * 100));
        xlabel('X');
        ylabel('Y');
        axis equal tight;
        
        % 暂停以创建动画效果
        pause(0.1);
    end
    %
