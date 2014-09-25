close all;
if exist('initialized','var')==0 || initialized==false
    %变量初始化
    
    disp('flag of initialized not found,initializing...');
    inputdata = load('Entry001.csv');
    X = inputdata(:,1)';
    Y = inputdata(:,2)';
    plot(X,Y);
    kmax=100;%最大迭代次数
    len=length(Y);
    kerlen=150;
    alpha=5000/15;
    gamma=50000/15;
    delta=0.5e8/15;
    beta=0.5e9/15;
    omega=5000/15;
    
    channelN=3%通道数
    ker=cell(1,channelN);%每个卷积核长度不相同,所以用了cell
    ker{1}=[Gauss_1D(22),zeros(1,10)];%第一个卷积核,符合能量守恒
    ker{2}=Gauss_1D(18);
    ker{3}=[zeros(1,30),Gauss_1D(5),zeros(1,10),Gauss_1D(10)]/2;
    g=zeros(len*channelN,1);
    Gi=cell(channelN,1);
    u=zeros(len,1);
    for i=1:channelN
        g((i-1)*len+1:i*len)=conv(Y,ker{i},'same')';%+rand(len,1)/20;%退化得到的第i个频谱
        u=u+g((i-1)*len+1:i*len);
        Gi{i}=generateCM(g((i-1)*len+1:i*len),kerlen,len);%对应的卷积矩阵
    end
    u=u/channelN;%u初始值为各个输入的平均
    %     g1=conv(Y,ker1,'same')';%退化得到的第一个频谱
    %     G1=generateCM(g1,kerlen,len);%对应的卷积矩阵
    %     g2=conv(Y,ker2,'same')';%退化得到的第二个频谱
    %     G2=generateCM(g2,kerlen,len);
    k=0;
    G=zeros(len*channelN*(channelN-1)/2,kerlen*channelN*(channelN-1)/2);
    Delta=generateCM([-1,2,-1],len,len);%Laplacian operator 
    %Dx=eye(len,len);%as no noise added, 
    for i=1:channelN-1
        for j=i+1:channelN
            G(k*len+1:(k+1)*len,(i-1)*kerlen+1:(i)*kerlen)=-Delta*Gi{j};
            G(k*len+1:(k+1)*len,(j-1)*kerlen+1:(j)*kerlen)=Delta*Gi{i};
            k=k+1;
        end
    end
    %fig2=plot(X,g1,X,g2);
    %g=[g1',g2']';
    % u=0.5*(g1+g2);
    h=zeros(kerlen*channelN,1);
    for i=1:channelN
        h(round((i-0.5)*kerlen))=1;
    end
    H=zeros(channelN*len,len);
    %     H1=generateCM(h1,len,len);
    %     H2=H1;
    %     H=[H1,H2];
    Dx=generateCM([1,-1],len,len);%微分
    G_Delta=G'*G;
    uhistory=cell(100,1);
    hhistory=cell(100,1);
    initialized=true;
end
tick1=clock();
%%循环迭代开始
kmax=20; %
for k=1:kmax %控制迭代次数
    %%求解u
    for i=1:channelN
        H((i-1)*len+1:i*len,:)=generateCM(h((i-1)*kerlen+1:i*kerlen),len,len);
    end
    %     H1=generateCM(h1,len,len);
    %     H2=generateCM(h2,len,len);
    %     H=[H1',H2']';
    v=zeros(len,1);
    a=zeros(len,1);
    w=zeros(len,1);
    b=zeros(len,1);
    hUsed=h;
    for j=1:15 %控制循环次数
        ulast=u;
        u=linsolve(H'*H+alpha/gamma*(Dx'*Dx)+omega/gamma*eye(len),H'*g+alpha/gamma*Dx'*(v+a))+omega/gamma*(w+b);
        s=abs(Dx*u-a);
        v=(Dx*u-a)./s.*max(s-1/alpha,0);
        a=a-Dx*u+v;
        w=max(u-b-1/omega,0);%保证u>0
        b=b-u+w;
        norm(ulast-u)/norm(u)
    end
    %     useeds=zeros(len,100);
    %     for n=1:100
    %         useeds(:,n)=u;
    %     end
    %     for n=1:20
    %         for m=1:100
    %             useeds(:,m)=useeds(:,m)*(rand(len,1)*0.1+0.95);
    %             useeds(:,m)=0.5*( useeds(:,m)+abs( useeds(:,m)));
    %             uenergy(m)=gamma/2*norm(H*useeds(:,m)-g)^2+sum(abs(Dx*useeds(:,m)));
    %         end
    %         [AAA,uorder]=sort(energy);
    %         ubuff=useeds;
    %         for m=1:100
    %             useeds(:,m)=ubuff(:,mod(m,20));
    %         end
    %     end
    fig_u=plot(X,u);
    %%求解h
    w=zeros(channelN*kerlen,1);
    b=w;
    Uc=generateCM(u,kerlen,len);
    U=zeros(len*channelN,kerlen*channelN);
    for i=1:channelN
        U((i-1)*len+1:i*len,(i-1)*kerlen+1:i*kerlen)=Uc;
    end
    %     U=[[Uc,zeros(len,kerlen)];[zeros(len,kerlen),Uc]];

    history=zeros(5000,1);
    i=1;
    exit=false;
    while  (exit~=true)%控制循环次数
        hlast=h;
        h=linsolve(U'*U+delta/gamma*G_Delta+beta/gamma*eye(channelN*kerlen),U'*g+beta/gamma*(w+b));
        w=max(h-b-1/beta,0);
        b=b-h+w;
        history(i)=norm(hlast-h)/norm(h);
        if history(i)<2e-4
            exit=true;
        end
        i=i+1;
    end
    disp('h solver, time looped:');
    i-1
    fig_h=plot(h);
    uhistory{k}=u;
    hhistory{k}=h;
    k
%     if k>1 %判断解是否劣化
%         if (gamma/2*norm(H*0.5*(abs(u)+u)-g)^2+sum(abs(Dx*0.5*(abs(u)+u)))>...
%                 gamma/2*norm(H*0.5*(abs(uhistory{k-1})+uhistory{k-1})-g)^2+sum(abs(Dx*0.5*(abs(uhistory{k-1})+uhistory{k-1}))))
%             
%             disp('exit criterion satisfied, iteration stops.')
%             break;
%         end
%     end
end
tick2=clock();
disp('core evaluation ends, time used:');
timeused=tick2-tick1;
timeused(4)*3600+timeused(5)*60+timeused(6) %所用时间,
figure;plot(G*h)
figure;plot([H(1:len,:)*u,g(1:len)])
figure;plot([H(len+1:2*len,:)*u,g(1+len:2*len)])
figure;plot([H(2*len+1:3*len,:)*u,g(1+2*len:3*len)])
figure;plot(u)
for i=1:channelN
    H((i-1)*len+1:i*len,:)=generateCM(h((i-1)*kerlen+1:i*kerlen),len,len);
end
gamma/2*norm(H*u-g)^2+sum(abs(Dx*u))+delta/2*h'*G_Delta*h
