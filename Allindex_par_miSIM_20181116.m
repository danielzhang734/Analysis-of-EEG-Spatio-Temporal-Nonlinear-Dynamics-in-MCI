% 输入：实验数据中的每个mat文件
% 输出：每个mat文件都有一个矩阵，30*30*【（1000-（4-1）*7-1=1000-22=4969】个miSIM,代表每个导联对，都有4969个miSIM
% 目的：为了后面给这4969个miSIM聚类做准备。
% 注意：要保证原始数据的格式，没有其他多余的内容，否则因为变量不一致的原因，会导致计算结果出现错误，例如之前我在原始数据中加入了一个测试数据，导致一直在算个测试数据。因为计算的对象变量名称是测试的，但是忘记改了。
clear;
RootPath= 'Z:\Sansa_documenetation\D_Resting_Each_Rythm_EP_Extract30ChannelData\ALPHA';
SavePath= 'Z:\Sansa_documenetation\R_mci_resting_nonlinear\ALPHA';
AllmatPath = fullfile(RootPath,'*.mat');
dirRoot_Path = dir(AllmatPath);
DataNumber = length(dirRoot_Path);
m = 14;  tau = 20;  keiler=50;    knn=6;  StateSpaceLength = 5000-(m-1)*tau-1;

for i_data = 1:DataNumber
    DataName = dirRoot_Path(i_data).name;
    prefix=DataName(1:end-4);
    Index_S_name = [prefix,'S.mat'];     Index_H_name = [prefix,'H.mat'];   Index_N_name = [prefix,'N.mat'];    Index_Corr_name = [prefix,'Corr.mat'];
    DataPath = fullfile(RootPath,DataName);
    Data = cell2mat(struct2cell(load(DataPath)));
    Data_permute = permute(Data,[3,2,1]);
    [trial,point,channel] = size(Data_permute);
    if tiral>20
        Data_permute = Data_permute(1:20,:,:)
        [trial,point,channel] = size(Data_permute);
    end
    data_Sxy_unit = zeros(trial,StateSpaceLength);% 98*4969 trial 要平均起来，相当于每对channel，每个时刻只有一个值，一共4969个时刻
    data_Syx_unit = zeros(trial,StateSpaceLength);
    data_Hxy_unit = zeros(trial,StateSpaceLength);
    data_Hyx_unit = zeros(trial,StateSpaceLength);
    data_Nxy_unit = zeros(trial,StateSpaceLength);
    data_Nyx_unit = zeros(trial,StateSpaceLength);
    data_CorrR_unit = zeros(trial,StateSpaceLength);
    Index_S = zeros(30,30,StateSpaceLength); Index_H = zeros(30,30,StateSpaceLength); Index_N = zeros(30,30,StateSpaceLength);Index_Corr = zeros(30,30);
    %             channel = 3;trial = 2;
            display(['the ',num2str(i_data),'file is processing']);
    for k = 1:30 % X导联循环
        
        for j = k+1:channel % Y导联循环
            if k~=j, % 矩a阵对角线上的值为
                
                data_k_j = Data_permute(:,:,[k;j]); %  data_k_j 98*1000*2
                parfor i = 1:trial
                    test_biv = data_k_j(i,:,:); % test_biv 1000*2
                    x = squeeze(test_biv);
                    synchro_results = synchro_miSIM(x,m,tau,keiler,knn);%(x,11,3,16,6);
                    data_Sxy_unit(i,:) = synchro_results(1,:); % syncho_miSIM 函数的输出是：data_Sxy_unit 98*4969
                    data_Syx_unit(i,:) = synchro_results(2,:);
                    data_Hxy_unit(i,:) = synchro_results(3,:);
                    data_Hyx_unit(i,:)= synchro_results(4,:);
                    data_Nxy_unit(i,:) = synchro_results(5,:);
                    data_Nyx_unit(i,:) = synchro_results(6,:);
                    data_CorrR_unit(i,:) = synchro_results(7,:);
                end
                Index_S(k,j,:) = sum(data_Sxy_unit)/trial;
                Index_S(j,k,:) = sum(data_Syx_unit)/trial;
                Index_H(k,j,:) = sum(data_Hxy_unit)/trial;
                Index_H(j,k,:) = sum(data_Hyx_unit)/trial;
                Index_N(k,j,:) = sum(data_Nxy_unit)/trial;
                Index_N(j,k,:) = sum(data_Nyx_unit)/trial;
                Index_Corr(k,j) = mean(data_CorrR_unit(:));
            end
        end
        display(['the ',num2str(k),'channel has finished processing']);
    end
    % saving data
    sname = fullfile(SavePath,Index_S_name);
    save(sname,'Index_S');   
    sname = fullfile(SavePath,Index_H_name);
    save(sname,'Index_H');
    sname = fullfile(SavePath,Index_N_name);
    save(sname,'Index_N');
    sname = fullfile(SavePath,Index_Corr_name);
    save(sname,'Index_Corr');
    
    display([sname,'        is finished processing']);
end


