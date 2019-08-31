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
%�����ndatos����  ---> 5000-(11-1)*3_1=5000-31=ndatos  5000-(4-1)*8-1=5000-25=ndatos

% 2018��11��22��15:03:17  xss �ĳ�miSIM������index��

ndatos = 5000-(m-1)*tau-1;  % length��x����ȡ����
% xn=zscore(x);    % �������������������׼��
cross = xcorr(x(:,1),x(:,2),0,'biased'); % c = xcorr(x,y,maxlags,'option') ͬʱָ��maxlags��option�Ļ���ؼ���.
out=zeros(4,ndatos); % ��ʼ��
sxy = zeros(ndatos,1); syx = zeros(ndatos,1); hxy = zeros(ndatos,1); hyx = zeros(ndatos,1); nxy = zeros(ndatos,1); nyx = zeros(ndatos,1);% ��ʼ��
auxx = zeros(1,nn);                         % ����Ĵ�����Ǻ����auxx��k�����бȽϣ��ҵ���nn����С���ٽ��㡣
indexx = zeros(1,nn);
auxy = zeros(1,nn);
indexy = zeros(1,nn);
for i = 1:ndatos; % ����ÿһ��Ԫ�� ����Ԫ�صĸ�����1200,1200-��10-1��*2-1=1200-18-1=1200-19=ndatos
    for k=1:nn                                      %INICIALIZE AUX
        auxx(k) = 100000000;                         % ����Ĵ�����Ǻ����auxx��k�����бȽϣ��ҵ���nn����С���ٽ��㡣
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
    for j = 1:ndatos-(m-1)*tau-1      %��j��ʱ���              %*****************�������
        distx(j) = 0;
        disty(j) = 0;
        %DISTANCES
        %           ii = 1;
        %           jj = 1;
        %           a_inspect_i = 0;
        %           a_inspect_j = 0;
        for k=0:m-1     % k = 0 1 2 3 4 5 6 7 8 9 һ��10��ά�ȡ���ʽ��m�Ǵ�1��ʼ�������Ǵ�0��ʼ����һ���ĵ���
            %             a_inspect_i(ii) = i+k*tau;ii = ii + 1;
            %             a_inspect_j(jj) = j+k*tau;jj = jj + 1;
            distx(j) = distx(j)+(x(i+k*tau,1)-x(j+k*tau,1)).^2;   % �������.
            disty(j) = disty(j)+(x(i+k*tau,2)-x(j+k*tau,2)).^2;   % ��2����������ڵ�
        end
        if ((abs(i-j)) > theiler)
            % ��������j>50�ľ��룬���뵽���10����ļ����С�j<50�򲻼����������ļ����ˡ�
            % ������10������ľ��룬�ҵ�һ����С�ģ��ͷŵ�����Ӧ��λ�á���Ϊauxx�ǰ��մ�
            % ��С��˳�������еġ���Ϊ���һ��������0�����Ծ��벻���ܱ�0��С�����ԾͿ�
            % �Ե�ס�ˣ������ҵ���С���Ǹ�ֵ�ˡ�
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
        rrx = rrx + distx(j);        %SIZE OF THE ATTRACTORS  �����ӵĴ�С
        % rrx�����ӡ���n=1ʱ����������������λ�Ƶľ�����������ˡ�
        rry = rry + disty(j);         % �����ǲ���j<50,���յ�rrx���ǰ����еľ��붼�������ˡ�
    end
    %*****************
    rxx = 0; ryy = 0; rxy = 0; ryx = 0;
    for k=1:nn  % %  nn: number of nearest neighbors
        rxx = auxx(k) + rxx; % auxx�� ����x1 �������k����
        ryy = auxy(k) + ryy; % auxy�� ����x1 ����k�������ӳ���
        rxy = distx(indexy(k)) + rxy;
        ryx = disty(indexx(k)) + ryx;
        % rxy: y�����k�����������Ҳ���ǽ�y���k���ٽ����ʱ���j��x�������ʱ�����Ǹ�����x�ľ����ж�Զ��
        % ���ǱȽ���ʱ��仯��ģʽ�Ƿ�һ�£������Ǻ��Լ���ֵ���бȽϡ��ͺ����ѧ���ͼ��ݺ���ϰ�ҵ�ģʽ��ʵ���Ƶġ�
        % ��rxx��һ����rxx֪ʶ�������k���ٽ���������ˡ�
    end
    rxx = rxx/nn;
    ryy = ryy/nn;
    rxy = rxy/nn;
    ryx = ryx/nn;
    % �����Ƕ���ɾ��� ��ndatos*1��ryx(i))��ʾ��i��sample��ֵ��
    sxy(i) = rxx/rxy;
    syx(i) = ryy/ryx;
    hxy(i) = log(rrx/(ndatos-(m-1)*tau-2) / rxy); % Hָ����ӻ��ˡ�
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

