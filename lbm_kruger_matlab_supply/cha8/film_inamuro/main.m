%% done
tic
% 创建对象并运行模拟
simulation = FilmInamuro();
simulation.run_simulation();
toc
%
%%
% 设置文件路径和文件名格式
filePattern = 'film_inamurofilm_inamuro%07d.dat';
numFrames = 168; % 根据实际生成的文件数量设置
NX = 40; % 根据实际情况设置
NY = 400; % 根据实际情况设置

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
    pause(0.01);
end