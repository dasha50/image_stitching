%-------------------------------------------------------------------------
%                   Parallax_tolerant Image Stitching 
%-------------------------------------------------------------------------

close all;
clear all;
clc;

%-------------------------
% User defined parameters.
%-------------------------
clear global;
global fitfn resfn degenfn psize numpar
fitfn = 'homography_fit';   %fit function
resfn = 'homography_res';   
degenfn = 'homography_degen';
psize   = 4;
numpar  = 9;

M     = 500;  % Number of hypotheses for RANSAC.
thr   = 0.1;  % RANSAC threshold.

scale = 0.25;    % Scale of input images (maybe for large images you would like to use a smaller scale).

%-------
% Paths.
%-------
addpath('modelspecific');
addpath('mexfiles');
addpath('multigs');

%-------------------
% Compile Mex files.
%-------------------
cd multigs;
if exist('computeIntersection','file')~=3
    mex computeIntersection.c; % <-- for multigs
end
cd ..;

cd mexfiles;
if exist('imagewarping','file')~=3
    mex ../imagewarping.cpp; 
end
if exist('wsvd','file')~=3
    mex ../wsvd.cpp; % We make use of eigen3's SVD in this file.
end
cd ..;

%----------------------
% Setup VLFeat toolbox.
%----------------------
cd vlfeat-0.9.14/toolbox;
feval('vl_setup');
cd ../..;

%---------------------------------------------
% Check if we are already running in parallel.
%---------------------------------------------
poolsize = matlabpool('size');
if poolsize == 0 %if not, we attempt to do it:
    matlabpool open;
end

%------------------
% Images to stitch.
%------------------
path1 = 'images/005/005.JPG';
path2 = 'images/005/006.JPG';

%-------------
% Read images.
%-------------
fprintf('Read images and SIFT matching\n');tic;
fprintf('> Reading images...');tic;
img1 = imresize(imread(sprintf('%s',path1)),scale);
img2 = imresize(imread(sprintf('%s',path2)),scale);
fprintf('done (%fs)\n',toc);

%四个角点坐标
C1=[1;1;1];
C2=[size(img2,2);1;1];
C3=[1;size(img2,1);1];
C4=[size(img2,2);size(img2,1);1];
%--------------------------------------
% SIFT keypoint detection and matching.
%--------------------------------------
fprintf('  Keypoint detection and matching...');tic;
[ kp1,ds1 ] = sift(single(rgb2gray(img1)));
[ kp2,ds2 ] = sift(single(rgb2gray(img2)));
%***************************
%****去掉重复特征点！！！****
%***************************
kp1=kp1(:,1:2);%取前两行特征点坐标
[kp1,i,~]=unique(kp1,'rows','stable');
ds1=ds1(i,:);
kp2=kp2(:,1:2);%取前两行特征点坐标
[kp2,i,~]=unique(kp2,'rows','stable');
ds2=ds2(i,:);
matches = match(ds1,ds2);
%kp:K-by-4,(y, x, scale, orientation)
kp1=kp1';
kp2=kp2';
%{
[ kp1,ds1 ] = vl_sift(single(rgb2gray(img1)),'PeakThresh', 0,'edgethresh',3);
[ kp2,ds2 ] = vl_sift(single(rgb2gray(img2)),'PeakThresh', 0,'edgethresh',3);

matches   = vl_ubcmatch(ds1,ds2);
[matches, scores]= vl_ubcmatch(ds1,ds2);
numBestPoints = 100;
[~,indices] = sort(scores);
%// Get the numBestPoints best matched features
bestMatches = matches(:,indices(1:numBestPoints));
%}

fprintf('done (%fs)\n',toc);

% Normalise point distribution.
fprintf('  Normalising point distribution...');tic;

%data_orig = [ kp1(1:2,matches(1,:)) ; ones(1,size(matches,2)) ; kp2(1:2,matches(2,:)) ; ones(1,size(matches,2)) ];
data_orig = [ kp1(2,matches(1,:)) ; 
    kp1(1,matches(1,:));
    ones(1,size(matches,2)) ; 
    kp2(2,matches(2,:)) ; 
    kp2(1,matches(2,:)) ;
    ones(1,size(matches,2)) ];
%***************************
%****去掉重复特征点！！！****
%***************************
data_orig=data_orig';
data_orig=unique(data_orig,'rows','stable');
data_orig=data_orig';
%data_orig表示所有特征点中匹配好的特征点的齐次坐标
[ dat_norm_img1,T1 ] = normalise2dpts(data_orig(1:3,:));
%dat_norm_img1=T1*data_orig(1:3,:)
[ dat_norm_img2,T2 ] = normalise2dpts(data_orig(4:6,:));
data_norm = [ dat_norm_img1 ; dat_norm_img2 ];
fprintf('done (%fs)\n',toc);

%-----------------
% Outlier removal.
%-----------------
fprintf('Outlier removal\n');tic;
% Multi-GS
rng(0);
[ ~,res,~,~ ] = multigsSampling(100,data_norm,M,10);
con = sum(res<=thr);
[ ~, maxinx ] = max(con);
inliers = find(res(:,maxinx)<=thr);
%inliers=inliers(1:2:end);%使特征点分散，不要太密集


if size(img1,1) == size(img2,1)
    % Show results of RANSAC.
    fprintf('  Showing results of RANSAC...');tic;
    figure;
    imshow([img1 img2]);
    hold on;
    %在左图上对IMG1上所有匹配好的点划红圈
    %plot(data_orig(1,:),data_orig(2,:),'ro','LineWidth',2);
    %在右图上对IMG2上所有匹配好的点划红圈
    %plot(data_orig(4,:)+size(img1,2),data_orig(5,:),'ro','LineWidth',2);
    for i=1:length(inliers)
        %在左图上对IMG1上第i个匹配好的inlier点(内点)划绿圈
        plot(data_orig(1,inliers(i)),data_orig(2,inliers(i)),'go','LineWidth',2);
        %在右图上对IMG2上第i个匹配好的inlier点(内点)划绿圈
        plot(data_orig(4,inliers(i))+size(img1,2),data_orig(5,inliers(i)),'go','LineWidth',2);
        %把匹配的inlier点连线
        %plot([data_orig(1,inliers(i)) data_orig(4,inliers(i))+size(img1,2)],[data_orig(2,inliers(i)) data_orig(5,inliers(i))],'g-');
    end
    title('Ransac''s results');
    fprintf('done (%fs)\n',toc);
end


%所有可以用来匹配的特征点
keypoint_img1=data_orig(1:2,inliers);
keypoint_img2=data_orig(4:5,inliers);
keypoint_norm=data_norm(:,inliers);

%-------------------------------------------------------------
% Step2.Finding Best-fitting Global homography (H).
%-------------------------------------------------------------
penalty = zeros(length(inliers),1); %初始化所有特征点的penalty值
flag=zeros(length(inliers),1); %初始化所有特征点的flag值，表示有没有被选为seed点，0为没有被选过
sum_threshold=0.01; %距离之和sum的阈值

edge_cost=1;
Hg_stack=[1,1,1];
edge_cost_stack=1;
while edge_cost>0.05
if mean(flag)==1
    break;
end    
%随机选取一个seed特征点
rng('shuffle');
seednum=randi(length(inliers));     %seednum为特征点索引
%一个特征点被选择为有效的seed点时
%它之前应该没有被选为seed点过
%而且它的penalty值应该低于所有特征点penalty的均值
if flag(seednum)==1
   continue;
else
   flag(seednum)=1;   %设置该seednum点被已被选为seed点
end   
if penalty(seednum)>mean(penalty)
    continue;
end

%搜索seed点的最近邻K个特征点%
seedfpt=data_orig(4:5,inliers(seednum));    %设置seed feature point的坐标
index=knnsearch(keypoint_img2',seedfpt','K',length(inliers));    
%构建特征点的集合的索引idx_featrues_group
featrues_group_idx=index(1:5); %初始化特征点集合，即只含有seed点和最近邻的3点组成的4个点的集合
%penalty(index(2:5))=penalty(index(2:5))+1;

for i=6:35
    %-------------
    %--更新group--
    %-------------
    featrues_group_idx=[featrues_group_idx,index(i)]; %特征点集合中每次增加一个特征点
    %--------
    %--求Hg--
    %--------
    [ h,A,D1,D2 ] = feval(fitfn,keypoint_norm(:,featrues_group_idx));
    Hg = T2\(reshape(h,3,3)*T1);    %Hg（全局单应性矩阵H）     
    %--------
    %--求Hs--
    %--------
    C1bar=Hg*C1;%左上
    C1bar=round([ C1bar(1)/C1bar(3) ; C1bar(2)/C1bar(3) ]);
    C2bar=Hg*C2;%右上
    C2bar=round([ C2bar(1)/C2bar(3) ; C2bar(2)/C2bar(3) ]);
    C3bar=Hg*C3;%左下
    C3bar=round([ C3bar(1)/C3bar(3) ; C3bar(2)/C3bar(3) ]);
    C4bar=Hg*C4;%右上
    C4bar=round([ C4bar(1)/C4bar(3) ; C4bar(2)/C4bar(3) ]);        
    save('corners.mat','C1bar','C2bar','C3bar','C4bar','img2','img1')
    %求与Hg最符合的相似变换Hs
    %用线性最小二乘法求解||AX-b||^2=0
    %一阶必要条件：A'AX-A'b=0
    %X=(A'A)\A'b
    A_hs=[C1(1),-C1(2),1,0;
          C1(2),C1(1),0,1;
          C2(1),-C2(2),1,0;
          C2(2),C2(1),0,1;
          C3(1),-C3(2),1,0;
          C3(2),C3(1),0,1;
          C4(1),-C4(2),1,0;
          C4(2),C4(1),0,1];
    b_hs=[C1bar(1);
          C1bar(2);
          C2bar(1);
          C2bar(2);
          C3bar(1);
          C3bar(2);
          C4bar(1);
          C4bar(2);];
    X_hs=(A_hs'*A_hs)\A_hs'*b_hs;
    dist_hs=A_hs*X_hs-b_hs;      
    Hs=[X_hs(1),-X_hs(2),X_hs(3);X_hs(2),X_hs(1),X_hs(4)];  
    %-----------------------
    %--求Hg和Hs距离之和sum--
    %-----------------------
    %求经过Hs变换的四角坐标
    C1s=Hs*[1;1;1];   %左上
    C2s=Hs*[size(img1,2);1;1];    %右上
    C3s=Hs*[1;size(img1,1);1];    %左下
    C4s=Hs*[size(img1,2);size(img1,1);1];     %右下   
    sum=norm(C1bar-C1s)+norm(C2bar-C2s)+norm(C3bar-C3s)+norm(C4bar-C4s);
    %距离之和归一化
    img1_size=size(img1,1)*size(img1,2);
    norm_sum=sum/img1_size;
    norm_sum
    %判断H是否满足预先设置的阈值
    if norm_sum>sum_threshold 
       featrues_group_idx=featrues_group_idx(:,1:i-1);
       break;
    end
end

%查看右图是否翻转
TL = Hg\[1;1;1];    %左上角
TL = round([ TL(1)/TL(3) ; TL(2)/TL(3) ]);
BL = Hg\[1;size(img2,1);1];     %左下角
BL = round([ BL(1)/BL(3) ; BL(2)/BL(3) ]);
TR = Hg\[size(img2,2);1;1];     %右上角
TR = round([ TR(1)/TR(3) ; TR(2)/TR(3) ]);
BR = Hg\[size(img2,2);size(img2,1);1];      %右下角
BR = round([ BR(1)/BR(3) ; BR(2)/BR(3) ]);
if TL(1)>TR(1) || TL(2)>BL(2) || BR(1)<BL(1) || BR(2)<TR(2)%如果在Hg下img2发生翻转，则放弃
    continue;
end
%特征点组包含的的特征点的penalty值+1
penalty(featrues_group_idx)=penalty(featrues_group_idx)+1;

%----------------------------------------------------
% Obtaining size of canvas (using global Homography).
%----------------------------------------------------
fprintf('Canvas size and offset (using global Homography)\n');
fprintf('> Getting canvas size...');tic;
 
% Canvas size.(包络)
cw = max([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) - min([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) + 1; %画面宽度
ch = max([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) - min([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) + 1; %画面高度
fprintf('done (%fs)\n',toc);

% Offset for left image.（img1左上角点在canvas上坐标）
fprintf('> Getting offset...');tic;
off = [ 1 - min([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) + 1 ; 1 - min([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) + 1 ];
fprintf('done (%fs)\n',toc);

% Warping source image with global homography 
fprintf('Image stitching with global homography (H) and linear blending\n');
fprintf('> Warping images by global homography...');tic;
warped_img1 = uint8(zeros(ch,cw,3));
warped_img1(off(2):(off(2)+size(img1,1)-1),off(1):(off(1)+size(img1,2)-1),:) = img1;
warped_img2 = imagewarping(double(ch),double(cw),double(img2),Hg,double(off));
warped_img2 = reshape(uint8(warped_img2),size(warped_img2,1),size(warped_img2,2)/3,3);
fprintf('done (%fs)\n',toc);


%-------------------------------------------------------------
% Step3.Evaluate the alignment quality of H
%-------------------------------------------------------------
%******************************************
%在warped_img下采样前，计算重叠区特征点分布
%******************************************
%计算重叠区域mask
w1 = imfill(im2bw(uint8(warped_img1), 0),'holes');
w2 = imfill(im2bw(uint8(warped_img2), 0),'holes');
mask = w1 & w2;%掩膜，重叠区域为1，其他区域为0
%寻找重叠区域的矩形包络rect
[gy, gx] = find(mask); %寻找mask不为零的点
miny = min(gy);
minx = min(gx);
maxy = max(gy);
maxx = max(gx);
%-------------------------------
%求重叠区域mask内的所有特征点
%-------------------------------
%初始化所有特征点
kp_img=zeros(ch,cw);
for i=1:length(index)
    %逐个特征点遍历
    point1=keypoint_img1(:,index(i)); %原图像img1中的特征点
    point2=keypoint_img2(:,index(i)); %img2坐标中的特征点
    %计算point2经过Hg变换后的img1坐标系中的点
    point2=pt_after_H(point2,Hg);
    %从img1坐标系到canvas坐标系，X=x+off(1)-1;Y=y+off(2)-1;
    %有特征点的位置设为1
    kp_img(uint8(point1(2))+off(2)-1,uint8(point1(1))+off(1)-1)=1;
    kp_img(point2(2)+off(2)-1,point2(1)+off(1)-1)=1;
end
%去掉重叠区域以外的特征点
kp_img=mask&kp_img;
kp_img=kp_img(miny:maxy,minx:maxx);%从canvas大小缩小到包络大小
[kpy,kpx]=find(kp_img);%kpx,kpy是纵列
%{
%**********
%**下采样**
%**********
kpy=kpy*(360/ch);
kpx=kpx*(480/cw);
%}
%重叠区域特征点集合
kp=[kpy,kpx];%n*2
%{
%*************************************
%warped_img下采样后
%**************************************

%下采样
warped_img1_ds=imresize(warped_img1,[360 480]);
warped_img2_ds=imresize(warped_img2,[360 480]);

%计算重叠区域mask
w1 = imfill(im2bw(uint8(warped_img1_ds), 0),'holes');
w2 = imfill(im2bw(uint8(warped_img2_ds), 0),'holes');
mask = w1 & w2;%掩膜，重叠区域为1，其他区域为0


%寻找重叠区域的矩形包络rect
[gy, gx] = find(mask); %寻找mask不为零的点
miny = min(gy);
minx = min(gx);
maxy = max(gy);
maxx = max(gx);
%}
mask=mask(miny:maxy,minx:maxx);
%包络区域长和宽
wd=maxx-minx+1;
ht=maxy-miny+1;
%包络矩形区域为1，其余区域为0
rect=zeros(ht,wd);

%计算函数fc图
allone=ones(ht,wd);
[i,j]=find(allone);
ps=[i j];
fc_dim1=fc(ps,kp);%1*m
fc_img=reshape(fc_dim1,ht,wd);
fc_img=fc_img.*mask;

%计算差值图Ed
%canny检测边缘，得到边缘图
gray_warped_img1=rgb2gray(warped_img1);
gray_warped_img2=rgb2gray(warped_img2);
warped_img1_edge = edge(gray_warped_img1,'canny'); 
warped_img2_edge = edge(gray_warped_img2,'canny'); 
Ed=bitxor(warped_img1_edge,warped_img2_edge);
Ed=Ed(miny:maxy,minx:maxx);
%figure;
%imshow(Ed);

%计算fc*Ed图
fced_img=fc_img.*Ed;

%非重叠区代价为无穷大
infinity=ones(ht,wd)*1e10;

%计算(x,y)与(x+1,y)的代价图cost_x
cost_x=fced_img+rowmove(fced_img,-1);
cost_x=cost_x.*mask+rowmove(infinity,-1).*(~mask);  %重叠区代价+非重叠区代价
cost_x(:,wd)=0;
%计算(x,y)与(x,y+1)的代价图cost_y
cost_y=fced_img+colmove(fced_img,-1);
cost_y=cost_y.*mask+colmove(infinity,-1).*(~mask);
cost_y(ht,:)=0;

%------------------------------------
%构建graph，主要参数：dc,nb
%------------------------------------
%———----------------------
%设置source、sink与节点相连的代价dc
%———----------------------
%rect左侧边以及包络内非重叠部分的的像素点与source相连
%rect右侧边的像素与sink相连
%与source相连的点
source=~mask & rect;%包络内非重叠部分为1，其余为0   
source(1:ht,1)=1;
source(1:ht,wd)=0;
source=source*1e10;   %代价无穷大
%与sink相连的点
loi=zeros(ht,wd);
loi(1:ht,wd)=1;  %rect右侧边上的点设为1
sink=loi&mask;   %只连接右侧边上位于重叠区域的点
sink=sink*1e10;   %面积ht*wd，代价无穷大
%构建dc
source=source';
sink=sink';
dc=[source(:)';sink(:)'];%2*(wd*ht)的矩阵

%———------------------------------------
%设置各节点之间的代价nb
%———------------------------------------
%roi=roi(miny:maxy,minx:maxx);
%nb = sparse(wd*ht,wd*ht);
Ix=(1:wd*ht-1)';
Jx=(2:wd*ht)';
Sx=cost_x';
Sx=Sx(:);
Sx=Sx(1:wd*ht-1);

Iy=(1:wd*(ht-1))';
Jy=(1+wd:wd*ht)';
Sy=cost_y';
Sy=Sy(:);
Sy=Sy(1:wd*(ht-1));

I=[Ix;Iy];
J=[Jx;Jy];
S=[Sx;Sy];
nb=sparse(I,J,S,wd*ht,wd*ht);

%{
I=zeros(2*wd*ht,1);
J=zeros(2*wd*ht,1);
S=zeros(2*wd*ht,1);
count=0;
for y=1:ht 
    tic;
    %连接edge，连接像素点与邻域内4点
    for x=1:wd
        if roi(y,x)==1  %如果在重叠区域mask以外，e(s,t)设为无穷大
            %把(x,y)与分别与(x+1,y)和(x,y+1)连接起来，edge cost为无穷大
            if (x < wd)
                count=count+1;
                I(count)=x*ht+y;
                J(count)=(x-1)*ht+y;
                S(count)=1e10;
                %nb((x-1)*ht+y,x*ht+y) = 1e10;
            end
            if (y < ht)
                count=count+1;
                I(count)=(x-1)*ht+y+1;
                J(count)=(x-1)*ht+y;
                S(count)=1e10;
                %nb((x-1)*ht+y,(x-1)*ht+y+1) = 1e10; 
            end
        else   %如果在重叠区域mask以内，e(s,t)按照公式来
            ps=[x+minx-1;y+miny-1];%当前点ps坐标(canvas坐标系)            
            if (x < wd)
                pt=[x+1;y];%pt(x+1,y)
                count=count+1;
                I(count)=x*ht+y;
                J(count)=(x-1)*ht+y;
                S(count)=fc(ps,kp)*Ed(y,x)+fc(pt,kp)*Ed(y,x+1);
                %nb((x-1)*ht+y,x*ht+y) = fc(ps,kp)*Ed(x+minx-1,y+miny-1)+fc(pt,kp)*Ed(x+minx,y+miny-1); 
            end
            if (y < ht)
                pt=[x;y+1];%pt(x,y+1)
                count=count+1;
                I(count)=(x-1)*ht+y+1;
                J(count)=(x-1)*ht+y;
                S(count)=fc(ps,kp)*Ed(y,x)+fc(pt,kp)*Ed(y+1,x);
                %nb((x-1)*ht+y,(x-1)*ht+y+1) = fc(ps,kp)*Ed(x+minx-1,y+miny-1)+fc(pt,kp)*Ed(x+minx-1,y+miny); 
            end            
        end
    end
    fprintf('done (%fs)\n',toc);
end
nb=sparse(I,J,S);
%}

%----------------------
%用graph-cut求最小代价
%----------------------
hinc = BK_Create(wd*ht,2*wd*ht);
BK_SetUnary(hinc,dc);
BK_SetNeighbors(hinc,nb);
edge_cost=BK_Minimize(hinc);

Hg_stack=[Hg_stack;Hg];
edge_cost_stack=[edge_cost_stack;edge_cost];

%显示featrues group
if size(img1,1) == size(img2,1)
    fprintf('  Showing results of features group');tic;
    figure;
    imshow([img1 img2]);
    hold on;    
    %{
    %在左图上对IMG1上所有匹配好的点划红圈
    plot(data_orig(1,:),data_orig(2,:),'ro','LineWidth',2);
    %在右图上对IMG2上所有匹配好的点划红圈
    plot(data_orig(4,:)+size(img1,2),data_orig(5,:),'ro','LineWidth',2);
    %}
    
    %在左图上对IMG1上所有匹配好的点划红圈
    plot(keypoint_img1(1,index),keypoint_img1(2,index),'ro','LineWidth',2);
    %在右图上对IMG2上所有匹配好的点划红圈
    plot(keypoint_img2(1,index)+size(img1,2),keypoint_img2(2,index),'ro','LineWidth',2);
   
    %在左图上对IMG1上第i个匹配好的inlier点(内点)划绿圈
    plot(keypoint_img1(1,featrues_group_idx),keypoint_img1(2,featrues_group_idx),'go','LineWidth',2);
    %在右图上对IMG2上第i个匹配好的inlier点(内点)划绿圈
    plot(keypoint_img2(1,featrues_group_idx)+size(img1,2),keypoint_img2(2,featrues_group_idx),'go','LineWidth',2);
   
    title('features group''s results');
    fprintf('done (%fs)\n',toc);
end
end
Hg_stack=Hg_stack(2:end,:);
edge_cost_stack=edge_cost_stack(2:end);
[~,min_index]=min(edge_cost_stack);
Hg=Hg_stack(3*min_index-2:3*min_index,:);



%{
%显示featrues group
if size(img1,1) == size(img2,1)
    fprintf('  Showing results of features group');tic;
    figure;
    imshow([img1 img2]);
    hold on;    
    %{
    %在左图上对IMG1上所有匹配好的点划红圈
    plot(data_orig(1,:),data_orig(2,:),'ro','LineWidth',2);
    %在右图上对IMG2上所有匹配好的点划红圈
    plot(data_orig(4,:)+size(img1,2),data_orig(5,:),'ro','LineWidth',2);
    %}
    
    %在左图上对IMG1上所有匹配好的点划红圈
    plot(keypoint_img1(1,index),keypoint_img1(2,index),'ro','LineWidth',2);
    %在右图上对IMG2上所有匹配好的点划红圈
    plot(keypoint_img2(1,index)+size(img1,2),keypoint_img2(2,index),'ro','LineWidth',2);
   
    %在左图上对IMG1上第i个匹配好的inlier点(内点)划绿圈
    plot(keypoint_img1(1,featrues_group_idx),keypoint_img1(2,featrues_group_idx),'go','LineWidth',2);
    %在右图上对IMG2上第i个匹配好的inlier点(内点)划绿圈
    plot(keypoint_img2(1,featrues_group_idx)+size(img1,2),keypoint_img2(2,featrues_group_idx),'go','LineWidth',2);
   
    title('features group''s results');
    fprintf('done (%fs)\n',toc);
end
%}
TL = Hg\[1;1;1];    %左上角
TL = round([ TL(1)/TL(3) ; TL(2)/TL(3) ]);
BL = Hg\[1;size(img2,1);1];     %左下角
BL = round([ BL(1)/BL(3) ; BL(2)/BL(3) ]);
TR = Hg\[size(img2,2);1;1];     %右上角
TR = round([ TR(1)/TR(3) ; TR(2)/TR(3) ]);
BR = Hg\[size(img2,2);size(img2,1);1];      %右下角
BR = round([ BR(1)/BR(3) ; BR(2)/BR(3) ]);
cw = max([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) - min([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) + 1; %画面宽度
ch = max([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) - min([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) + 1; %画面高度
off = [ 1 - min([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) + 1 ; 1 - min([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) + 1 ];
warped_img1 = uint8(zeros(ch,cw,3));
warped_img1(off(2):(off(2)+size(img1,1)-1),off(1):(off(1)+size(img1,2)-1),:) = img1;
warped_img2 = imagewarping(double(ch),double(cw),double(img2),Hg,double(off));
warped_img2 = reshape(uint8(warped_img2),size(warped_img2,1),size(warped_img2,2)/3,3);
%显示warped_img1、warped_img2
fprintf('  Showing results of warped imgs');tic;
figure;
imshow([warped_img1 warped_img2]);




%-------------------------------------------------------------
% Step4.Employ the optimal homography to pre-align images
%-------------------------------------------------------------

%-------------------------------------------------------------
% Step5.Use content-preserving warping to refine the alignment
%-------------------------------------------------------------

% Blending images by simple average (linear blending)
fprintf('  Homography linear image blending (averaging)...');tic;
linear_hom = imageblending(warped_img1,warped_img2);
fprintf('done (%fs)\n',toc);
figure;
imshow(linear_hom);
title('Image Stitching with global homography');

%{
figure;
label=BK_GetLabeling(hinc);
result=reshape(label,[wd,ht]);
imagesc(result');
drawnow;
%}


