% % 创建对象并运行模拟
tic
simulation = Gaussian1DBGK();
simulation.run_simulation();
toc
%% 绘图

filePath = 'gaussian_1d_bgk';
filePattern = 'gaussian_1d_bgk%06d.dat';

% 设置动画参数
numFrames = 21; % 总帧数
pauseTime = 0.5; % 每帧之间的暂停时间

% 创建图形窗口
figure;
hold on;
xlabel('Position');
ylabel('Phase');
title('Gaussian 1D BGK Animation');

% 读取并绘制每一帧
for frame = 0:numFrames-1
    % 构建文件名
    fileName = sprintf([filePath, filePattern], frame);
    
    % 检查文件是否存在
    if exist(fileName, 'file')
        % 读取数据
        data = load(fileName);
        
        % 清除当前图形
        cla;
        
        % 绘制数据
        plot(data, 'LineWidth', 2);
        
        % 更新图形
        drawnow;
        
        % 暂停
        pause(pauseTime);
    else
        disp(['File not found: ', fileName]);
    end
end

hold off;