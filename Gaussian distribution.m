%clear;
clear all; close all; clc;

% 第一组数据
mu1 = [0 0 ];  %均值
S1 = [.1 0 ; 0 .1];  %协方差
data1 = mvnrnd(mu1, S1, 100);   %产生高斯分布数据

%第二组数据
mu2 = [1.25 1.25 ];
S2 = [.1 0 ; 0 .1];
data2 = mvnrnd(mu2, S2, 100);

% 第三组数据
mu3 = [-1.25 1.25 ];
S3 = [.1 0 ; 0 .1];
data3 = mvnrnd(mu3, S3, 100);

% 显示数据
plot(data1(:,1), data1(:,2), 'b+');
hold on;
plot(data2(:,1), data2(:,2), 'r+');
plot(data3(:,1), data3(:,2), 'g+');
grid on;

% 三类数据合成一个不带标号的数据类
data = [data1; data2; data3]; 

% 使用肘部法则确定聚类数量
max_clusters = 10;
N_range = 1:max_clusters;
WCSS = zeros(1, max_clusters);

for k = N_range
    [pattern, ~] = kmeans(data, k);
    WCSS(k) = sum(sum((data - pattern).^2));
end

% 绘制WCSS随聚类数量变化的图表
figure;
plot(N_range, WCSS, 'o-');
xlabel('Number of Clusters');
ylabel('WCSS');
title('Elbow Method to Determine Optimal Number of Clusters');
grid on;

% 使用 kmeans 函数进行聚类
[centers, cluster_idx] = kmeans(data, optimal_k);

% 检查 centers 的维度
if size(centers, 2) ~= size(data, 2)
    error('centers does not have the same number of columns as data');
end

% 检查 data 是否至少有两列
if size(data, 2) < 2
    error('data must have at least two columns for plotting');
end

% 如果没有错误，继续绘制
figure;
colors = {'r', 'g', 'b', 'y', 'm', 'c', 'w', 'k'}; % 添加更多颜色，如果需要
hold on;
for i = 1:optimal_k
    cluster_points = data(cluster_idx == i, :);
    plot(cluster_points(:, 1), cluster_points(:, 2), strrep(colors{i}, 'o', '*'));
    plot(centers(i, 1), centers(i, 2), 'ko'); % 绘制聚类中心
end
grid on;
hold off;
