function linear_hom=main(ima1,ima2)
%-------------------------------------------------------------------------
%                   Parallax_tolerant Image Stitching 
%-------------------------------------------------------------------------

% close all;
% clear all;
% clc;

%-------------------------
% User defined parameters.
%-------------------------
%clear global;
%定义全局变量
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
%MEX文件是一种可在matlab环境中调用的C（或fortran）语言衍生程序，
%MEX文件的后缀名按32位/64位分别为 .mexw32/.mexw64。
%MEX文件是由C或Fortran语言编写的源代码，经matlab编译器处理而生成的二进制文件。
%它是可以被matlab解释器自动装载并执行的动态链接程序，类似windows下的dll文件。
%-------------------
cd multigs;
if exist('computeIntersection','file')~=3  %exist用来判断变量或函数是否存在：返回3 表示是 MEX-file on MATLAB's search path
    mex computeIntersection.c; % <-- for multigs 如编译链接C语言的MEX文件源程序，在MATLAB的控制窗口中输入：mex computeIntersection.c生成一个名为computeIntersection.mexw64的MEX文件
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
% Setup VLFeat toolbox.  VLFeat是图像处理库，在matlab中配置VLFeat
%----------------------
cd vlfeat-0.9.14/toolbox;
feval('vl_setup');
cd ../..;

%---------------------------------------------
% Check if we are already running in parallel.  判断并行运算环境是否启动
%---------------------------------------------
poolsize = matlabpool('size');
if poolsize == 0 %if not, we attempt to do it:
    matlabpool open;
end

%------------------
% Images to stitch.
%------------------
% path1 = 'images/005/005.JPG';
% path2 = 'images/005/006.JPG';

%-------------
% Read images.
%-------------
% 当你需要计算一组Matlab操作的运行时间时，可以使用tic和toc函数。
% tic函数启动一个秒表，表示计时开始；
% toc则停止这个秒表，表示计时结束，并计算出所经历的时间（单位为秒）。
fprintf('Read images and SIFT matching\n');tic;
fprintf('> Reading images...');tic;
img1 = imresize(imread(ima1),scale);  %图像缩放到image*scale
img2 = imresize(imread(ima2),scale);
fprintf('done (%fs)\n',toc);


% size（）：获取矩阵的行数和列数
% size(A,n)如果在size函数的输入参数中再添加一项n，并用1或2为n赋值，则 size将返回矩阵的行数或列数。
% 其中r=size(A,1)该语句返回的时矩阵A的行数， c=size(A,2) 该语句返回的时矩阵A的列数。

%四个角点坐标    
C1=[1;1;1]; 
C2=[size(img2,2);1;1];
C3=[1;size(img2,1);1];
C4=[size(img2,2);size(img2,1);1];
%--------------------------------------
% SIFT keypoint detection and matching.
%--------------------------------------


fprintf('  Keypoint detection and matching...');tic;
% [image, descriptors, locs] = sift(imageFile)
%     descriptors: a K-by-128 matrix, where each row gives an invariant
%         descriptor for one of the K keypoints.  The descriptor is a vector
%         of 128 values normalized to unit length.
%     locs: K-by-4 matrix, in which each row行 has the 4 values for a
%         keypoint location (row, column, scale, orientation).  The 
%         orientation is in the range [-PI, PI] radians.
%这里kp1是loc即位置
%ds1是描述子


[ kp1,ds1 ] = sift(single(rgb2gray(img1))); %kp1表示特征点location Kx4矩阵，K表示K个特征点，ds1，Kx128矩阵，特征点描述子
[ kp2,ds2 ] = sift(single(rgb2gray(img2)));


%***************************
%****去掉重复特征点！！！****
%***************************
%A(:,k:k+m)表示取A矩阵第k~k+m列的全部元素
% [C,IA,IC] = unique(A,'rows','stable') returns the rows of C in the same
%  order that they appear in A, C=A[IA],A=C[IC] 即C是由A中按照IA抽取

kp1=kp1(:,1:2);     %取前两列特征点坐标y,x
[kp1,i,~]=unique(kp1,'rows','stable');  %stable按照原来顺序排，找kp1中非重复的行，i保存的是剩余的行数
ds1=ds1(i,:);       %按顺序取ds内的值
kp2=kp2(:,1:2);     %取前两列特征点坐标
[kp2,i,~]=unique(kp2,'rows','stable');
ds2=ds2(i,:);

matches = match(ds1,ds2); %特征点匹配，matches 2行K列，K是匹配特征点数，2行分别是图一和图二中匹配特征点索引
%kp:K-by-4,(y, x, scale, orientation)
kp1=kp1';   %矩阵转置
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
%data_orig前两行存储的是图一匹配特征点坐标，4,5行存储的是图二的特征点坐标，每一列图一图二的特征点匹配索引
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
data_orig=unique(data_orig,'rows','stable');%按行unique，stable排序顺序不变
data_orig=data_orig';
%data_orig表示所有特征点中匹配好的特征点的齐次坐标

%normalise2dpts齐次坐标归一化，便于求单应性矩阵
%齐次坐标就是将一个原本是n维的向量用一个n+1维向量来表示。
% 比如齐次坐标(8,4,2)、(4,2,1)表示的都是二维点(4,2)。
% 给出点的齐次表达式[X Y H]，就可求得其二维笛卡尔坐标，即
% [X Y H]→  = [x y 1]， 这个过程称为正常化处理。

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
% M     = 500;  % Number of hypotheses假定 for RANSAC.
% thr   = 0.1;  % RANSAC threshold.

[ ~,res,~,~ ] = multigsSampling(100,data_norm,M,10);  %M次ransac，结果记录到res
con = sum(res<=thr);   %res和thr比较后得到的是0,1组成的矩阵，然后按列求和得到con
[ ~, maxinx ] = max(con);
inliers = find(res(:,maxinx)<=thr);
%inliers得到的是M*1的矩阵，该矩阵记录了内点索引

%inliers=inliers(1:2:end);%使特征点分散，不要太密集


if size(img1,1) == size(img2,1)
    % Show results of RANSAC.
    fprintf('  Showing results of RANSAC...');tic;
%     figure;
%     imshow([img1 img2]);
%     hold on;
    %在左图上对IMG1上所有匹配好的点划红圈
    %plot(data_orig(1,:),data_orig(2,:),'ro','LineWidth',2);
    %在右图上对IMG2上所有匹配好的点划红圈
    %plot(data_orig(4,:)+size(img1,2),data_orig(5,:),'ro','LineWidth',2);
    for i=1:length(inliers)
        %在左图上对IMG1上第i个匹配好的inlier点(内点)划绿圈，g：green o:circle linewidth:线宽2
     %   plot(data_orig(1,inliers(i)),data_orig(2,inliers(i)),'go','LineWidth',2);
        %在右图上对IMG2上第i个匹配好的inlier点(内点)划绿圈
     %   plot(data_orig(4,inliers(i))+size(img1,2),data_orig(5,inliers(i)),'go','LineWidth',2);
        %把匹配的inlier点连线
        %plot([data_orig(1,inliers(i)) data_orig(4,inliers(i))+size(img1,2)],[data_orig(2,inliers(i)) data_orig(5,inliers(i))],'g-');
    end
   % title('Ransac''s results');
    fprintf('done (%fs)\n',toc);
end


%所有可以用来匹配的特征点,提取图一，二的内点inliers，作为特征点
keypoint_img1=data_orig(1:2,inliers);
keypoint_img2=data_orig(4:5,inliers);
keypoint_norm=data_norm(:,inliers);

%-------------------------------------------------------------
% Step2.Finding Best-fitting Global homography (H).   best-fitting最佳拟合，最相符
%-------------------------------------------------------------
penalty = zeros(length(inliers),1); %初始化所有特征点的penalty值，zeros(m,n),产生m行n列0矩阵，每次选取的seed点，该点的penalty值加一
flag=zeros(length(inliers),1); %初始化所有特征点的flag标志值，表示有没有被选为seed点，0为没有被选过
%Hg_fitness_threshold=2
Hg_fitness_threshold=3;
sum_threshold=0.01;             %距离之和sum的阈值

edge_cost=1;
Hg_stack=[];           %存储找到的单应性矩阵
Hg_fitness_stack=[];   %存储符合条件单应性矩阵
edge_cost_stack=[];    %存储edge_cost
while edge_cost>0.05
    if mean(flag)==1   %所有点均被选取
        break;
    end    
    %随机选取一个seed特征点
    rng('shuffle');                     %基于时间产生seed，每次都不相同
    seednum=randi(length(inliers));     %seednum为特征点索引，任意1-length的数
    
    %一个特征点被选择为有效的seed点时
    %它之前应该没有被选为seed点过
    %而且它的penalty值应该低于所有特征点penalty的均值？
    if flag(seednum)==1
       continue;
    else
       flag(seednum)=1;   %设置该seednum点被已被选为seed点
    end   
    if penalty(seednum)>mean(penalty)
        continue;
    end

    %搜索seed点的最近邻K个特征点%
    seedfpt=data_orig(4:5,inliers(seednum));    %设置seed feature point的坐标，以img2为标准
    index=knnsearch(keypoint_img2',seedfpt','K',length(inliers));   %index保存的是离seednum最近的特征点点的顺序
    %构建特征点的集合的索引idx_featrues_group
    featrues_group_idx=index(1:9); %初始化特征点集合，即只含有seed点和最近邻的3点组成的4个点的集合
    
    %寻找特征点组
    is_Hg_exist=0;
    for i=10:length(index)
        %-------------
        %--更新group--
        %-------------
        featrues_group_idx=[featrues_group_idx,index(i)]; %特征点集合中每次增加一个特征点
        %--------
        %--求Hg--  ？？？不知道怎么求出来的
        %--------
        [ h,A,D1,D2 ] = feval(fitfn,keypoint_norm(:,featrues_group_idx));%FEVAL Execute the specified function.fitfn = 'homography_fit';  
        Hg = T2\(reshape(h,3,3)*T1);    %Hg（全局单应性矩阵H）    %reshape重新调整矩阵的行数、列数、维数
                                                                 %函数对原数组的抽取是按照列抽取的（对原数组按列抽取，抽取的元素填充为新数组的列） 
        Hg_fitness=fitness(A,h);
        Hg_fitness;    %Hg的拟合度
        %判断H是否满足拟合阈值
        if Hg_fitness>Hg_fitness_threshold 
            if i==10  
                is_Hg_exist=0;
                break;
            else
                %把最后满足拟合阈值条件的Hg作为选中的Hg进行step3
                featrues_group_idx=featrues_group_idx(:,1:i-1);
                [ h,A,D1,D2 ] = feval(fitfn,keypoint_norm(:,featrues_group_idx));
                Hg = T2\(reshape(h,3,3)*T1);    %Hg（全局单应性矩阵H）    
                Hg_fitness=fitness(A,h);
                is_Hg_exist=1;
                break;
            end
        end
        if i==length(index)
            is_Hg_exist=1;
        end
    end
    
    if is_Hg_exist==0
        continue;
    end
    
    %特征点组包含的的特征点的penalty值+1，penalty矩阵中把含featrues_group_idx的序列加一
    penalty(featrues_group_idx)=penalty(featrues_group_idx)+1;
     
    %查看右图是否翻转  ？？？？
    TL = Hg\[1;1;1];    %左上角
    %TL(1)指矩阵中第一个元素，矩阵排序，先列后行，若是TL(1，3)指第一行第三列的元素
    %这里对TL进行归一化，TL是3x1，归一化后是2x1
    TL = round([ TL(1)/TL(3) ; TL(2)/TL(3) ]);
    BL = Hg\[1;size(img2,1);1];     %左下角
    BL = round([ BL(1)/BL(3) ; BL(2)/BL(3) ]);
    TR = Hg\[size(img2,2);1;1];     %右上角
    TR = round([ TR(1)/TR(3) ; TR(2)/TR(3) ]);
    BR = Hg\[size(img2,2);size(img2,1);1];      %右下角
    BR = round([ BR(1)/BR(3) ; BR(2)/BR(3) ]);
    %{
    if TL(1)>TR(1) || TL(2)>BL(2) || BR(1)<BL(1) || BR(2)<TR(2)%如果在Hg下img2发生翻转，则放弃
        continue;
    end
    %}
    % Canvas size.(包络)
    cw = max([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) - min([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) + 1; %画面宽度
    ch = max([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) - min([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) + 1; %画面高度

    % Offset for left image.（img1左上角点在canvas上坐标）
    off = [ 1 - min([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) + 1 ; 1 - min([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) + 1 ];
    % Warping source image with global homography 
    warped_img1 = uint8(zeros(ch,cw,3));
    warped_img1(off(2):(off(2)+size(img1,1)-1),off(1):(off(1)+size(img1,2)-1),:) = img1;    %将img1赋值到画布上的位置
    warped_img2 = imagewarping(double(ch),double(cw),double(img2),Hg,double(off));
    warped_img2 = reshape(uint8(warped_img2),size(warped_img2,1),size(warped_img2,2)/3,3);  %？/3
    
%     imwrite(warped_img1,'warped_img1.png','png');
%     imwrite(warped_img2,'warped_img2.png','png');
%     
    
    
      
    
    %-------------------------------------------------------------
    % Step3.Evaluate the alignment quality of H
    %-------------------------------------------------------------
    %--------------------------
    %对H质量评价前进行预筛选，去掉一些与相似变换Hs差别太大的H
    %--------------------------
    %--------
    %--求Hs-- Homography Screening   
    %--------
    %round 四舍五入
    %四个角点坐标   
    % C1=[1;1;1]; 
    % C2=[size(img2,2);1;1];
    % C3=[1;size(img2,1);1];
    % C4=[size(img2,2);size(img2,1);1];
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
    Sum=norm(C1bar-C1s)+norm(C2bar-C2s)+norm(C3bar-C3s)+norm(C4bar-C4s);
    %距离之和归一化
    img1_size=size(img1,1)*size(img1,2);
    norm_sum=Sum/img1_size;
    %判断H是否满足预先设置的阈值
    if norm_sum>sum_threshold 
       continue;
    end
    
    %-------------------------------------
    %对Hg进行质量评价
    %-------------------------------------
    %BW2 = imfill(BW,'holes')
    %填充二值图像中的空洞区域。 如， 黑色的背景上有个白色的圆圈。 则这个圆圈内区域将被填充。
    %matlab中函数im2bw使用阈值（threshold）变换法把灰度图像（grayscale image）转换成二值图像
    %一般意义上是指只有纯黑（0）、纯白（255）两种颜色的图像。 当然， 也可以是其他任意两种颜色的组合。
    
    %计算重叠区域mask
    w1 = imfill(im2bw(uint8(warped_img1), 0),'holes');
    w2 = imfill(im2bw(uint8(warped_img2), 0),'holes');
    
    mask = w1 & w2;            %掩膜，重叠区域为1，其他区域为0
    
    %重叠区边界一分为二，分别与source和sink相连
    % getborder获取掩膜的边界
    border1=getborder(w1);
    border2=getborder(w2);
    border=border1&border2;
    border1=border1-border;   
    source=border2&mask;%右图的边缘在重叠区的部分连接source
    sink=border1&mask;  %左图的边缘在重叠区的部分连接sink
    %寻找重叠区域的矩形包络rect
    %[i,j]=find(A) 返回矩阵A中非零元素所在的行i,列j
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
        point1=keypoint_img1(:,index(i)); %img1中的特征点坐标
        point2=keypoint_img2(:,index(i)); %img2特征点坐标
        %计算point2经过Hg变换到img1坐标系上的点的坐标
        point2=pt_after_H(point2,Hg);
        %从img1坐标系到canvas坐标系，X=x+off(1)-1;Y=y+off(2)-1; ？为何这么转换  canvas帆布
        %将img1和img2上的特征点转到canvas后的位置设为1
        kp_img(uint8(point1(2))+off(2)-1,uint8(point1(1))+off(1)-1)=1;
        kp_img(point2(2)+off(2)-1,point2(1)+off(1)-1)=1;
    end
    %去掉重叠区域以外的特征点
    kp_img=mask&kp_img;
    kp_img=kp_img(miny:maxy,minx:maxx);%从canvas大小缩小到包络大小
    [kpy,kpx]=find(kp_img);%kpx,kpy是纵列
    %重叠区域特征点集合
    kp=[kpy,kpx];%n*2,得到重叠区内所有特征点
    
    mask=mask(miny:maxy,minx:maxx);%mask缩小到包络大小
    
    wd=maxx-minx+1;%包络的宽
    ht=maxy-miny+1;%包络的高

    %计算函数fc（置信度）图
    all_one=ones(ht,wd);
    [i,j]=find(all_one);
    ps=[i j];
    fc_dim1=fc(ps,kp);%1*m
    fc_img=reshape(fc_dim1,ht,wd);
    fc_img=fc_img.*mask;%只计算重叠区域内的点的fc
    
  %  imwrite(fc_img,'fc_img.png','png');
    
    %计算差值图Ed
    %canny检测边缘，得到边缘图
    gray_warped_img1=rgb2gray(warped_img1);           %rgb2gray是一种函数，功能是将真彩色图像转换为灰度图像。
    gray_warped_img2=rgb2gray(warped_img2);
    warped_img1_edge = edge(gray_warped_img1,'canny'); 
    warped_img2_edge = edge(gray_warped_img2,'canny'); %使用canny方法检测图像边缘，在边缘的地方矩阵值为1
    Ed=bitxor(warped_img1_edge,warped_img2_edge);      %按位异或，相同为0，不同为1，为了得到
    Ed=Ed(miny:maxy,minx:maxx);
    
%     figure;
%     imshow(Ed);



    %计算fc*Ed图
    fced_img=fc_img.*Ed;
    %计算(x,y)与(x+1,y)的代价图cost_x
    cost_x=fced_img+rowmove(fced_img,-1)+1e-6;       % x>0,右移x x<0,左移x
    infinity=ones(ht,wd)*1e10;                  %infinity无穷大
    cost_x=cost_x.*mask+rowmove(infinity,-1).*(~mask);  %包络内重叠区代价（fced）+包络内非重叠区代价（无穷大）
    cost_x(:,wd)=0;
 
    %计算(x,y)与(x,y+1)的代价图cost_y
    cost_y=fced_img+colmove(fced_img,-1)+1e-6;
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
    %{ 
    source=~mask;%包络内非重叠部分为1，其余为0   
    source(1:ht,1)=1;%左侧边为1
    source(1:ht,wd)=0;%右侧边为0
    source=source*1e10;   %代价无穷大
    %与sink相连的点
    sink=zeros(ht,wd);
    sink(1:ht,wd)=1;  %右侧边上的点设为1
    %sink=sink&mask;   %只连接右侧边上位于重叠区域的点
    sink=sink*1e10;   %代价无穷大
    %}
    
    %????
    source=source(miny:maxy,minx:maxx);
    sink=sink(miny:maxy,minx:maxx);
    source=source+~mask;%source与包络内非重叠部分相连
    source=source*1e10;   %代价无穷大，代价即两节点的权值，边界代价无穷大
    sink=sink*1e10;   %代价无穷大
    %构建dc，datacost数据项
    source=source';  %进行转置和（：）'是为了得到一行权重，该行的前后顺序是按图像矩阵行编
    sink=sink';
    dc=[source(:)';sink(:)'];%2*(wd*ht)的矩阵  A=[3 4 2;1 5 3;4 7 1];A(:),变成一列，314457231，再转置变成行
     
    %———------------------------------------
    %设置各节点之间的代价nb   ？？？？不是sc吗？smoothcost平滑项
    %———------------------------------------
    %roi=roi(miny:maxy,minx:maxx);
    %nb = sparse(wd*ht,wd*ht);
    
    %???????
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
    nb=sparse(I,J,S,wd*ht,wd*ht);   %稀疏矩阵，I行J列，值为S
    %稀疏矩阵，就是看起来很松散的，也就是说，在这个矩阵中，绝大多数元素是零元素。
    %----------------------
    %用graph-cut求最小代价，图割
    %----------------------
    hinc = BK_Create(wd*ht,2*wd*ht);   %创建目标体
    BK_SetUnary(hinc,dc);              %dc在这里用到,添加数据项
    BK_SetNeighbors(hinc,nb);          %nb在这里用到,dc,添加平滑项
    edge_cost=BK_Minimize(hinc);       %使能量最小化函数
    
    
    L=BK_GetLabeling(hinc);            %得到哪个点属于哪个类
    maskone=zeros(wd*ht,1);
    maskone=im2uint8(maskone)+1;
    B=L-maskone;    
    C=reshape(B,wd,ht);
    
    C=C';                   %source右图 0,sink左图 1,这是重叠区域的mask,原mask经过了缩放到现在的大小
    Cc=C.*255;
%     figure;
%     imshow(Cc);
    D=uint8(~C);
    
    %imwrite(C,'warped_img1.png','png');
    
    
    Hg_stack=[Hg_stack;Hg];
    edge_cost_stack=[edge_cost_stack;edge_cost];
    Hg_fitness_stack=[Hg_fitness_stack;Hg_fitness];
    %显示featrues group
    if size(img1,1) == size(img2,1)
        fprintf('  Showing results of features group');tic;
       % figure;
      %  imshow([img1 img2]);
       % hold on;   
        %在左图上对IMG1上所有匹配好的点划红圈
     %   plot(keypoint_img1(1,index),keypoint_img1(2,index),'ro','LineWidth',2);
        %在右图上对IMG2上所有匹配好的点划红圈
     %   plot(keypoint_img2(1,index)+size(img1,2),keypoint_img2(2,index),'ro','LineWidth',2);

        %在左图上对IMG1上第i个匹配好的inlier点(内点)划绿圈
     %   plot(keypoint_img1(1,featrues_group_idx),keypoint_img1(2,featrues_group_idx),'go','LineWidth',2);
        %在右图上对IMG2上第i个匹配好的inlier点(内点)划绿圈
      %  plot(keypoint_img2(1,featrues_group_idx)+size(img1,2),keypoint_img2(2,featrues_group_idx),'go','LineWidth',2);

      %  title('features group''s results');
        fprintf('done (%fs)\n',toc);
    end
end
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
%mask=mask(miny:maxy,minx:maxx)
warped_img1(miny:maxy,minx:maxx,1)=warped_img1(miny:maxy,minx:maxx,1).*C;
warped_img1(miny:maxy,minx:maxx,2)=warped_img1(miny:maxy,minx:maxx,2).*C;
warped_img1(miny:maxy,minx:maxx,3)=warped_img1(miny:maxy,minx:maxx,3).*C;

 
 
warped_img2 = imagewarping(double(ch),double(cw),double(img2),Hg,double(off));
warped_img2 = reshape(uint8(warped_img2),size(warped_img2,1),size(warped_img2,2)/3,3);
% 
% 
warped_img2(miny:maxy,minx:maxx,1)=warped_img2(miny:maxy,minx:maxx,1).*D;
warped_img2(miny:maxy,minx:maxx,2)=warped_img2(miny:maxy,minx:maxx,2).*D;
warped_img2(miny:maxy,minx:maxx,3)=warped_img2(miny:maxy,minx:maxx,3).*D;


%显示warped_img1、warped_img2
% fprintf('  Showing results of warped imgs');tic;
% figure;
% imshow([warped_img1 warped_img2]);




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
% figure;
% imshow(linear_hom);
% title('Image Stitching with global homography');
imwrite(linear_hom,'GraphCutBlend.png','png');
%{
figure;
label=BK_GetLabeling(hinc);
result=reshape(label,[wd,ht]);
imagesc(result');
drawnow;
%}


