function out=synchro_miSIM(x,m,tau,theiler,nn)
%	Given a bivariate input data (x) calculates the cross-correlation and the
%   non-linear interpendences S,H and N.
%	Results are in the vector out in the following order:
%	S(X|Y) S(Y|X) H(X|Y) H(Y|X) N(X|Y) N(Y|X) cross-correlation.
%	Parameters to be set are: m (embedding dimension), tau (time delay),
%	theiler (theiler correction in # data points) and NN (# of nearest neighbors).
%
%m=10;                   %embedding dimension
%tau=2;                  %time lag
%theiler=50;             %theiler correction
%nn = 10;                %number of nearest neighbors
%输出：ndatos个点  ---> 5000-(11-1)*3_1=5000-31=ndatos  5000-(4-1)*8-1=5000-25=ndatos

% 2018年11月22日15:03:17  xss 改成miSIM，所有index的

ndatos = 5000-(m-1)*tau-1;  % length（x）是取行数
% xn=zscore(x);    % 将这个导联的数据做标准化
cross = xcorr(x(:,1),x(:,2),0,'biased'); % c = xcorr(x,y,maxlags,'option') 同时指定maxlags和option的互相关计算.
out=zeros(4,ndatos); % 初始化
sxy = zeros(ndatos,1); syx = zeros(ndatos,1); hxy = zeros(ndatos,1); hyx = zeros(ndatos,1); nxy = zeros(ndatos,1); nyx = zeros(ndatos,1);% 初始化
auxx = zeros(1,nn);                         % 后面的代码就是和这个auxx（k）进行比较，找到这nn个最小的临近点。
indexx = zeros(1,nn);
auxy = zeros(1,nn);
indexy = zeros(1,nn);
for i = 1:ndatos; % 遍历每一个元素 现在元素的个数是1200,1200-（10-1）*2-1=1200-18-1=1200-19=ndatos
    for k=1:nn                                      %INICIALIZE AUX
        auxx(k) = 100000000;                         % 后面的代码就是和这个auxx（k）进行比较，找到这nn个最小的临近点。
        indexx(k) = 100000000;
        auxy(k) = 100000000;
        indexy(k) = 100000000;
    end
    auxx(nn+1) = 0;
    auxy(nn+1) = 0;
    indexx(nn+1) = 100000000;
    indexy(nn+1) = 100000000;
    rrx = 0; rry = 0;
    distx = zeros(1,(ndatos-(m-1)*tau-1));
    disty = zeros(1,(ndatos-(m-1)*tau-1));
    for j = 1:ndatos-(m-1)*tau-1      %第j个时间点              %*****************计算距离
        distx(j) = 0;
        disty(j) = 0;
        %DISTANCES
        %           ii = 1;
        %           jj = 1;
        %           a_inspect_i = 0;
        %           a_inspect_j = 0;
        for k=0:m-1     % k = 0 1 2 3 4 5 6 7 8 9 一共10个维度。公式中m是从1开始，这里是从0开始，是一样的道理。
            %             a_inspect_i(ii) = i+k*tau;ii = ii + 1;
            %             a_inspect_j(jj) = j+k*tau;jj = jj + 1;
            distx(j) = distx(j)+(x(i+k*tau,1)-x(j+k*tau,1)).^2;   % 迭代相加.
            disty(j) = disty(j)+(x(i+k*tau,2)-x(j+k*tau,2)).^2;   % 第2列它的最近邻点
        end
        if ((abs(i-j)) > theiler)
            % 对于所有j>50的距离，加入到最近10个点的计算中。j<50则不加入最近距离的计算了。
            % 遍历这10个最近的距离，找到一个更小的，就放到它对应的位置。因为auxx是按照从
            % 大到小的顺序来排列的。因为最后一个数字是0，所以距离不可能比0还小，所以就可
            % 以挡住了，就能找到最小的那个值了。
            if (distx(j) < auxx(1))
                flagx=0;
                for k=1:nn+1
                    if (distx(j) < auxx(k))
                        auxx(k) = auxx(k+1);
                        indexx(k) = indexx(k+1);
                    else
                        auxx(k-1) = distx(j);
                        indexx(k-1) = j;
                        flagx=1;
                    end
                    if flagx==1;break;end
                end
            end
            if (disty(j) < auxy(1))
                flagy=0;
                for k=1:nn+1
                    if (disty(j) < auxy(k))
                        auxy(k) = auxy(k+1);
                        indexy(k) = indexy(k+1);
                    else
                        auxy(k-1) = disty(j);
                        indexy(k-1) = j;
                        flagy=1;
                    end
                    if flagy==1,break,end
                end
            end
        end
        rrx = rrx + distx(j);        %SIZE OF THE ATTRACTORS  吸引子的大小
        % rrx：分子。当n=1时，将它和其它所有位移的距离相加起来了。
        rry = rry + disty(j);         % 不管是不是j<50,最终的rrx还是把所有的距离都加起来了。
    end
    %*****************
    rxx = 0; ryy = 0; rxy = 0; ryx = 0;
    for k=1:nn  % %  nn: number of nearest neighbors
        rxx = auxx(k) + rxx; % auxx是 距离x1 点最近的k个点
        ryy = auxy(k) + ryy; % auxy是 距离x1 的与k个最近的映射点
        rxy = distx(indexy(k)) + rxy;
        ryx = disty(indexx(k)) + ryx;
        % rxy: y的最近k个点的索引。也就是将y最近k个临近点的时间点j，x按照这个时间点的那个窗和x的距离有多远。
        % 它是比较随时间变化的模式是否一致，所以是和自己的值进行比较。就好像大学里的图书馆和自习室的模式其实类似的。
        % 和rxx不一样。rxx知识将最近的k个临近点加起来了。
    end
    rxx = rxx/nn;
    ryy = ryy/nn;
    rxy = rxy/nn;
    ryx = ryx/nn;
    % 把它们都变成矩阵 （ndatos*1）ryx(i))表示第i个sample的值。
    sxy(i) = rxx/rxy;
    syx(i) = ryy/ryx;
    hxy(i) = log(rrx/(ndatos-(m-1)*tau-2) / rxy); % H指标分子换了。
    hyx(i) = log(rry/(ndatos-(m-1)*tau-2) / ryx);
    nxy(i) = 1-rxy/rrx;
    nyx(i) = 1-ryx/rry;     
end

out(1,:)=sxy';
out(2,:)=syx';
out(3,:)=hxy';
out(4,:)=hyx';
out(5,:)=nxy';
out(6,:)=nyx';
out(7,:)=cross;

