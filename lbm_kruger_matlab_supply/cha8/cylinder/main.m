%% done
tic
flow = cylinder2();  
flow.run_simulation();      % 运行模拟，耗时700cpus，cpp耗时40cpus
toc
%
%% 
% 设置文件路径和文件名格式
file_path = 'cylinder/';
file_format = 'cylinder_cylinder%07d.dat';
start_frame = 200;
end_frame = 6000;
frame_step = 200;

% 读取第一个文件以获取网格大小
first_file = sprintf([file_path, file_format], start_frame);
data = load(first_file);
[NY, NX] = size(data);

% 创建一个 figure 窗口
figure;

% 创建视频对象
video = VideoWriter('cylinder_simulation.avi');
open(video);

% 循环读取文件并绘制动画
for frame = start_frame:frame_step:end_frame
    % 构建文件名
    file_name = sprintf([file_path, file_format], frame);
     
    % 读取数据
    data = load(file_name);
    
    % 绘制数据
    imagesc(data);
    colorbar;
    title(['Frame: ', num2str(frame)]);
    xlabel('X');
    ylabel('Y');
    
    % 获取当前帧并写入视频
    frame_data = getframe(gcf);
    writeVideo(video, frame_data);
    
    % 暂停以控制动画速度
    pause(0.01);
end

% 关闭视频对象
close(video);

disp('动画生成完毕');