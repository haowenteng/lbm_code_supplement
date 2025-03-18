%% 创建对象并运行模拟
tic
gaussian = Gaussian2DTRT();

% 运行模拟
gaussian.run_simulation();
toc
%
%% 设置文件路径和文件名格式
    filePath = '';
    filePattern = 'gaussian_2d_trtgaussian_2d_trt%06d.dat';
    
    % 设置动画参数
    numFrames = 210; % 总帧数
    pauseTime = 0.1; % 每帧之间的暂停时间
     
    % 创建图形窗口
    figure;
    colormap('jet');
    colorbar;
    xlabel('X');
    ylabel('Y');
    title('Gaussian 2D TRT Animation');
    
    % 读取并绘制每一帧
    for frame = 0:10:numFrames-1
        % 构建文件名
        fileName = sprintf([filePath, filePattern], frame);
        
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