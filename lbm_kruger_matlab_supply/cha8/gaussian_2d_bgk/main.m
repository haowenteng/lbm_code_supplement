%% 创建对象并运行模拟
tic
% 创建 Gaussian2DBGK 类的实例
gaussian = Gaussian2DBGK();

% 运行模拟 
gaussian.run_simulation();
toc
%
%%
 % 设置文件路径和文件名格式
    filePath = '';
    filePattern = 'gaussian_2d_bgkgaussian_2d_bgk%06d.dat';
    
    % 设置动画参数
    numFrames = 21; % 总帧数
    pauseTime = 0.2; % 每帧之间的暂停时间
    
    % 创建图形窗口
    figure;
    colormap('jet');
    colorbar;
    xlabel('X');
    ylabel('Y');
    title('Gaussian 2D BGK Animation');
    
    % 读取并绘制每一帧
    for frame = 0:numFrames-1
        % 构建文件名
        fileName = sprintf([filePath, filePattern], frame*10);
        
        % 检查文件是否存在
        if exist(fileName, 'file')
            % 读取数据
            data = load(fileName);
            
            % 清除当前图形
            cla;
            
            % 绘制数据
            imagesc(data);
            axis equal tight;
            colorbar;
            
            % 更新图形
            drawnow;
            
            % 暂停
            pause(pauseTime);
        else
            disp(['File not found: ', fileName]);
        end
    end
    %