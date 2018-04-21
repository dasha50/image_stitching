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
%����ȫ�ֱ���
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
%MEX�ļ���һ�ֿ���matlab�����е��õ�C����fortran��������������
%MEX�ļ��ĺ�׺����32λ/64λ�ֱ�Ϊ .mexw32/.mexw64��
%MEX�ļ�����C��Fortran���Ա�д��Դ���룬��matlab��������������ɵĶ������ļ���
%���ǿ��Ա�matlab�������Զ�װ�ز�ִ�еĶ�̬���ӳ�������windows�µ�dll�ļ���
%-------------------
cd multigs;
if exist('computeIntersection','file')~=3  %exist�����жϱ��������Ƿ���ڣ�����3 ��ʾ�� MEX-file on MATLAB's search path
    mex computeIntersection.c; % <-- for multigs ���������C���Ե�MEX�ļ�Դ������MATLAB�Ŀ��ƴ��������룺mex computeIntersection.c����һ����ΪcomputeIntersection.mexw64��MEX�ļ�
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
% Setup VLFeat toolbox.  VLFeat��ͼ����⣬��matlab������VLFeat
%----------------------
cd vlfeat-0.9.14/toolbox;
feval('vl_setup');
cd ../..;

%---------------------------------------------
% Check if we are already running in parallel.  �жϲ������㻷���Ƿ�����
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
% ������Ҫ����һ��Matlab����������ʱ��ʱ������ʹ��tic��toc������
% tic��������һ�������ʾ��ʱ��ʼ��
% toc��ֹͣ��������ʾ��ʱ���������������������ʱ�䣨��λΪ�룩��
fprintf('Read images and SIFT matching\n');tic;
fprintf('> Reading images...');tic;
img1 = imresize(imread(ima1),scale);  %ͼ�����ŵ�image*scale
img2 = imresize(imread(ima2),scale);
fprintf('done (%fs)\n',toc);


% size��������ȡ���������������
% size(A,n)�����size��������������������һ��n������1��2Ϊn��ֵ���� size�����ؾ����������������
% ����r=size(A,1)����䷵�ص�ʱ����A�������� c=size(A,2) ����䷵�ص�ʱ����A��������

%�ĸ��ǵ�����    
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
%     locs: K-by-4 matrix, in which each row�� has the 4 values for a
%         keypoint location (row, column, scale, orientation).  The 
%         orientation is in the range [-PI, PI] radians.
%����kp1��loc��λ��
%ds1��������


[ kp1,ds1 ] = sift(single(rgb2gray(img1))); %kp1��ʾ������location Kx4����K��ʾK�������㣬ds1��Kx128����������������
[ kp2,ds2 ] = sift(single(rgb2gray(img2)));


%***************************
%****ȥ���ظ������㣡����****
%***************************
%A(:,k:k+m)��ʾȡA�����k~k+m�е�ȫ��Ԫ��
% [C,IA,IC] = unique(A,'rows','stable') returns the rows of C in the same
%  order that they appear in A, C=A[IA],A=C[IC] ��C����A�а���IA��ȡ

kp1=kp1(:,1:2);     %ȡǰ��������������y,x
[kp1,i,~]=unique(kp1,'rows','stable');  %stable����ԭ��˳���ţ���kp1�з��ظ����У�i�������ʣ�������
ds1=ds1(i,:);       %��˳��ȡds�ڵ�ֵ
kp2=kp2(:,1:2);     %ȡǰ��������������
[kp2,i,~]=unique(kp2,'rows','stable');
ds2=ds2(i,:);

matches = match(ds1,ds2); %������ƥ�䣬matches 2��K�У�K��ƥ������������2�зֱ���ͼһ��ͼ����ƥ������������
%kp:K-by-4,(y, x, scale, orientation)
kp1=kp1';   %����ת��
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
%data_origǰ���д洢����ͼһƥ�����������꣬4,5�д洢����ͼ�������������꣬ÿһ��ͼһͼ����������ƥ������
data_orig = [ kp1(2,matches(1,:)) ; 
              kp1(1,matches(1,:));
              ones(1,size(matches,2)) ; 
              kp2(2,matches(2,:)) ; 
              kp2(1,matches(2,:)) ;
              ones(1,size(matches,2)) ];
%***************************
%****ȥ���ظ������㣡����****
%***************************
data_orig=data_orig';
data_orig=unique(data_orig,'rows','stable');%����unique��stable����˳�򲻱�
data_orig=data_orig';
%data_orig��ʾ������������ƥ��õ���������������

%normalise2dpts��������һ����������Ӧ�Ծ���
%���������ǽ�һ��ԭ����nά��������һ��n+1ά��������ʾ��
% �����������(8,4,2)��(4,2,1)��ʾ�Ķ��Ƕ�ά��(4,2)��
% ���������α��ʽ[X Y H]���Ϳ�������ά�ѿ������꣬��
% [X Y H]��  = [x y 1]�� ������̳�Ϊ����������

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
% M     = 500;  % Number of hypotheses�ٶ� for RANSAC.
% thr   = 0.1;  % RANSAC threshold.

[ ~,res,~,~ ] = multigsSampling(100,data_norm,M,10);  %M��ransac�������¼��res
con = sum(res<=thr);   %res��thr�ȽϺ�õ�����0,1��ɵľ���Ȼ������͵õ�con
[ ~, maxinx ] = max(con);
inliers = find(res(:,maxinx)<=thr);
%inliers�õ�����M*1�ľ��󣬸þ����¼���ڵ�����

%inliers=inliers(1:2:end);%ʹ�������ɢ����Ҫ̫�ܼ�


if size(img1,1) == size(img2,1)
    % Show results of RANSAC.
    fprintf('  Showing results of RANSAC...');tic;
%     figure;
%     imshow([img1 img2]);
%     hold on;
    %����ͼ�϶�IMG1������ƥ��õĵ㻮��Ȧ
    %plot(data_orig(1,:),data_orig(2,:),'ro','LineWidth',2);
    %����ͼ�϶�IMG2������ƥ��õĵ㻮��Ȧ
    %plot(data_orig(4,:)+size(img1,2),data_orig(5,:),'ro','LineWidth',2);
    for i=1:length(inliers)
        %����ͼ�϶�IMG1�ϵ�i��ƥ��õ�inlier��(�ڵ�)����Ȧ��g��green o:circle linewidth:�߿�2
     %   plot(data_orig(1,inliers(i)),data_orig(2,inliers(i)),'go','LineWidth',2);
        %����ͼ�϶�IMG2�ϵ�i��ƥ��õ�inlier��(�ڵ�)����Ȧ
     %   plot(data_orig(4,inliers(i))+size(img1,2),data_orig(5,inliers(i)),'go','LineWidth',2);
        %��ƥ���inlier������
        %plot([data_orig(1,inliers(i)) data_orig(4,inliers(i))+size(img1,2)],[data_orig(2,inliers(i)) data_orig(5,inliers(i))],'g-');
    end
   % title('Ransac''s results');
    fprintf('done (%fs)\n',toc);
end


%���п�������ƥ���������,��ȡͼһ�������ڵ�inliers����Ϊ������
keypoint_img1=data_orig(1:2,inliers);
keypoint_img2=data_orig(4:5,inliers);
keypoint_norm=data_norm(:,inliers);

%-------------------------------------------------------------
% Step2.Finding Best-fitting Global homography (H).   best-fitting�����ϣ������
%-------------------------------------------------------------
penalty = zeros(length(inliers),1); %��ʼ�������������penaltyֵ��zeros(m,n),����m��n��0����ÿ��ѡȡ��seed�㣬�õ��penaltyֵ��һ
flag=zeros(length(inliers),1); %��ʼ�������������flag��־ֵ����ʾ��û�б�ѡΪseed�㣬0Ϊû�б�ѡ��
%Hg_fitness_threshold=2
Hg_fitness_threshold=3;
sum_threshold=0.01;             %����֮��sum����ֵ

edge_cost=1;
Hg_stack=[];           %�洢�ҵ��ĵ�Ӧ�Ծ���
Hg_fitness_stack=[];   %�洢����������Ӧ�Ծ���
edge_cost_stack=[];    %�洢edge_cost
while edge_cost>0.05
    if mean(flag)==1   %���е����ѡȡ
        break;
    end    
    %���ѡȡһ��seed������
    rng('shuffle');                     %����ʱ�����seed��ÿ�ζ�����ͬ
    seednum=randi(length(inliers));     %seednumΪ����������������1-length����
    
    %һ�������㱻ѡ��Ϊ��Ч��seed��ʱ
    %��֮ǰӦ��û�б�ѡΪseed���
    %��������penaltyֵӦ�õ�������������penalty�ľ�ֵ��
    if flag(seednum)==1
       continue;
    else
       flag(seednum)=1;   %���ø�seednum�㱻�ѱ�ѡΪseed��
    end   
    if penalty(seednum)>mean(penalty)
        continue;
    end

    %����seed��������K��������%
    seedfpt=data_orig(4:5,inliers(seednum));    %����seed feature point�����꣬��img2Ϊ��׼
    index=knnsearch(keypoint_img2',seedfpt','K',length(inliers));   %index���������seednum�������������˳��
    %����������ļ��ϵ�����idx_featrues_group
    featrues_group_idx=index(1:9); %��ʼ�������㼯�ϣ���ֻ����seed�������ڵ�3����ɵ�4����ļ���
    
    %Ѱ����������
    is_Hg_exist=0;
    for i=10:length(index)
        %-------------
        %--����group--
        %-------------
        featrues_group_idx=[featrues_group_idx,index(i)]; %�����㼯����ÿ������һ��������
        %--------
        %--��Hg--  ��������֪����ô�������
        %--------
        [ h,A,D1,D2 ] = feval(fitfn,keypoint_norm(:,featrues_group_idx));%FEVAL Execute the specified function.fitfn = 'homography_fit';  
        Hg = T2\(reshape(h,3,3)*T1);    %Hg��ȫ�ֵ�Ӧ�Ծ���H��    %reshape���µ��������������������ά��
                                                                 %������ԭ����ĳ�ȡ�ǰ����г�ȡ�ģ���ԭ���鰴�г�ȡ����ȡ��Ԫ�����Ϊ��������У� 
        Hg_fitness=fitness(A,h);
        Hg_fitness;    %Hg����϶�
        %�ж�H�Ƿ����������ֵ
        if Hg_fitness>Hg_fitness_threshold 
            if i==10  
                is_Hg_exist=0;
                break;
            else
                %��������������ֵ������Hg��Ϊѡ�е�Hg����step3
                featrues_group_idx=featrues_group_idx(:,1:i-1);
                [ h,A,D1,D2 ] = feval(fitfn,keypoint_norm(:,featrues_group_idx));
                Hg = T2\(reshape(h,3,3)*T1);    %Hg��ȫ�ֵ�Ӧ�Ծ���H��    
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
    
    %������������ĵ��������penaltyֵ+1��penalty�����аѺ�featrues_group_idx�����м�һ
    penalty(featrues_group_idx)=penalty(featrues_group_idx)+1;
     
    %�鿴��ͼ�Ƿ�ת  ��������
    TL = Hg\[1;1;1];    %���Ͻ�
    %TL(1)ָ�����е�һ��Ԫ�أ������������к��У�����TL(1��3)ָ��һ�е����е�Ԫ��
    %�����TL���й�һ����TL��3x1����һ������2x1
    TL = round([ TL(1)/TL(3) ; TL(2)/TL(3) ]);
    BL = Hg\[1;size(img2,1);1];     %���½�
    BL = round([ BL(1)/BL(3) ; BL(2)/BL(3) ]);
    TR = Hg\[size(img2,2);1;1];     %���Ͻ�
    TR = round([ TR(1)/TR(3) ; TR(2)/TR(3) ]);
    BR = Hg\[size(img2,2);size(img2,1);1];      %���½�
    BR = round([ BR(1)/BR(3) ; BR(2)/BR(3) ]);
    %{
    if TL(1)>TR(1) || TL(2)>BL(2) || BR(1)<BL(1) || BR(2)<TR(2)%�����Hg��img2������ת�������
        continue;
    end
    %}
    % Canvas size.(����)
    cw = max([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) - min([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) + 1; %������
    ch = max([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) - min([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) + 1; %����߶�

    % Offset for left image.��img1���Ͻǵ���canvas�����꣩
    off = [ 1 - min([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) + 1 ; 1 - min([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) + 1 ];
    % Warping source image with global homography 
    warped_img1 = uint8(zeros(ch,cw,3));
    warped_img1(off(2):(off(2)+size(img1,1)-1),off(1):(off(1)+size(img1,2)-1),:) = img1;    %��img1��ֵ�������ϵ�λ��
    warped_img2 = imagewarping(double(ch),double(cw),double(img2),Hg,double(off));
    warped_img2 = reshape(uint8(warped_img2),size(warped_img2,1),size(warped_img2,2)/3,3);  %��/3
    
%     imwrite(warped_img1,'warped_img1.png','png');
%     imwrite(warped_img2,'warped_img2.png','png');
%     
    
    
      
    
    %-------------------------------------------------------------
    % Step3.Evaluate the alignment quality of H
    %-------------------------------------------------------------
    %--------------------------
    %��H��������ǰ����Ԥɸѡ��ȥ��һЩ�����Ʊ任Hs���̫���H
    %--------------------------
    %--------
    %--��Hs-- Homography Screening   
    %--------
    %round ��������
    %�ĸ��ǵ�����   
    % C1=[1;1;1]; 
    % C2=[size(img2,2);1;1];
    % C3=[1;size(img2,1);1];
    % C4=[size(img2,2);size(img2,1);1];
    C1bar=Hg*C1;%����
    C1bar=round([ C1bar(1)/C1bar(3) ; C1bar(2)/C1bar(3) ]);
    C2bar=Hg*C2;%����
    C2bar=round([ C2bar(1)/C2bar(3) ; C2bar(2)/C2bar(3) ]);
    C3bar=Hg*C3;%����
    C3bar=round([ C3bar(1)/C3bar(3) ; C3bar(2)/C3bar(3) ]);
    C4bar=Hg*C4;%����
    C4bar=round([ C4bar(1)/C4bar(3) ; C4bar(2)/C4bar(3) ]);        
    save('corners.mat','C1bar','C2bar','C3bar','C4bar','img2','img1')
    %����Hg����ϵ����Ʊ任Hs
    %��������С���˷����||AX-b||^2=0
    %һ�ױ�Ҫ������A'AX-A'b=0
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
    %--��Hg��Hs����֮��sum--
    %-----------------------
    %�󾭹�Hs�任���Ľ�����
    C1s=Hs*[1;1;1];   %����
    C2s=Hs*[size(img1,2);1;1];    %����
    C3s=Hs*[1;size(img1,1);1];    %����
    C4s=Hs*[size(img1,2);size(img1,1);1];     %����   
    Sum=norm(C1bar-C1s)+norm(C2bar-C2s)+norm(C3bar-C3s)+norm(C4bar-C4s);
    %����֮�͹�һ��
    img1_size=size(img1,1)*size(img1,2);
    norm_sum=Sum/img1_size;
    %�ж�H�Ƿ�����Ԥ�����õ���ֵ
    if norm_sum>sum_threshold 
       continue;
    end
    
    %-------------------------------------
    %��Hg������������
    %-------------------------------------
    %BW2 = imfill(BW,'holes')
    %����ֵͼ���еĿն����� �磬 ��ɫ�ı������и���ɫ��ԲȦ�� �����ԲȦ�����򽫱���䡣
    %matlab�к���im2bwʹ����ֵ��threshold���任���ѻҶ�ͼ��grayscale image��ת���ɶ�ֵͼ��
    %һ����������ָֻ�д��ڣ�0�������ף�255��������ɫ��ͼ�� ��Ȼ�� Ҳ��������������������ɫ����ϡ�
    
    %�����ص�����mask
    w1 = imfill(im2bw(uint8(warped_img1), 0),'holes');
    w2 = imfill(im2bw(uint8(warped_img2), 0),'holes');
    
    mask = w1 & w2;            %��Ĥ���ص�����Ϊ1����������Ϊ0
    
    %�ص����߽�һ��Ϊ�����ֱ���source��sink����
    % getborder��ȡ��Ĥ�ı߽�
    border1=getborder(w1);
    border2=getborder(w2);
    border=border1&border2;
    border1=border1-border;   
    source=border2&mask;%��ͼ�ı�Ե���ص����Ĳ�������source
    sink=border1&mask;  %��ͼ�ı�Ե���ص����Ĳ�������sink
    %Ѱ���ص�����ľ��ΰ���rect
    %[i,j]=find(A) ���ؾ���A�з���Ԫ�����ڵ���i,��j
    [gy, gx] = find(mask); %Ѱ��mask��Ϊ��ĵ�
    miny = min(gy);
    minx = min(gx);
    maxy = max(gy);
    maxx = max(gx);
    %-------------------------------
    %���ص�����mask�ڵ�����������
    %-------------------------------
    %��ʼ������������
    kp_img=zeros(ch,cw);
    for i=1:length(index)
        %������������
        point1=keypoint_img1(:,index(i)); %img1�е�����������
        point2=keypoint_img2(:,index(i)); %img2����������
        %����point2����Hg�任��img1����ϵ�ϵĵ������
        point2=pt_after_H(point2,Hg);
        %��img1����ϵ��canvas����ϵ��X=x+off(1)-1;Y=y+off(2)-1; ��Ϊ����ôת��  canvas����
        %��img1��img2�ϵ�������ת��canvas���λ����Ϊ1
        kp_img(uint8(point1(2))+off(2)-1,uint8(point1(1))+off(1)-1)=1;
        kp_img(point2(2)+off(2)-1,point2(1)+off(1)-1)=1;
    end
    %ȥ���ص����������������
    kp_img=mask&kp_img;
    kp_img=kp_img(miny:maxy,minx:maxx);%��canvas��С��С�������С
    [kpy,kpx]=find(kp_img);%kpx,kpy������
    %�ص����������㼯��
    kp=[kpy,kpx];%n*2,�õ��ص���������������
    
    mask=mask(miny:maxy,minx:maxx);%mask��С�������С
    
    wd=maxx-minx+1;%����Ŀ�
    ht=maxy-miny+1;%����ĸ�

    %���㺯��fc�����Ŷȣ�ͼ
    all_one=ones(ht,wd);
    [i,j]=find(all_one);
    ps=[i j];
    fc_dim1=fc(ps,kp);%1*m
    fc_img=reshape(fc_dim1,ht,wd);
    fc_img=fc_img.*mask;%ֻ�����ص������ڵĵ��fc
    
  %  imwrite(fc_img,'fc_img.png','png');
    
    %�����ֵͼEd
    %canny����Ե���õ���Եͼ
    gray_warped_img1=rgb2gray(warped_img1);           %rgb2gray��һ�ֺ����������ǽ����ɫͼ��ת��Ϊ�Ҷ�ͼ��
    gray_warped_img2=rgb2gray(warped_img2);
    warped_img1_edge = edge(gray_warped_img1,'canny'); 
    warped_img2_edge = edge(gray_warped_img2,'canny'); %ʹ��canny�������ͼ���Ե���ڱ�Ե�ĵط�����ֵΪ1
    Ed=bitxor(warped_img1_edge,warped_img2_edge);      %��λ�����ͬΪ0����ͬΪ1��Ϊ�˵õ�
    Ed=Ed(miny:maxy,minx:maxx);
    
%     figure;
%     imshow(Ed);



    %����fc*Edͼ
    fced_img=fc_img.*Ed;
    %����(x,y)��(x+1,y)�Ĵ���ͼcost_x
    cost_x=fced_img+rowmove(fced_img,-1)+1e-6;       % x>0,����x x<0,����x
    infinity=ones(ht,wd)*1e10;                  %infinity�����
    cost_x=cost_x.*mask+rowmove(infinity,-1).*(~mask);  %�������ص������ۣ�fced��+�����ڷ��ص������ۣ������
    cost_x(:,wd)=0;
 
    %����(x,y)��(x,y+1)�Ĵ���ͼcost_y
    cost_y=fced_img+colmove(fced_img,-1)+1e-6;
    cost_y=cost_y.*mask+colmove(infinity,-1).*(~mask);
    cost_y(ht,:)=0;
    
    %------------------------------------
    %����graph����Ҫ������dc,nb
    %------------------------------------
    %������----------------------
    %����source��sink��ڵ������Ĵ���dc
    %������----------------------
    %rect�����Լ������ڷ��ص����ֵĵ����ص���source����
    %rect�Ҳ�ߵ�������sink����
    %��source�����ĵ�
    %{ 
    source=~mask;%�����ڷ��ص�����Ϊ1������Ϊ0   
    source(1:ht,1)=1;%����Ϊ1
    source(1:ht,wd)=0;%�Ҳ��Ϊ0
    source=source*1e10;   %���������
    %��sink�����ĵ�
    sink=zeros(ht,wd);
    sink(1:ht,wd)=1;  %�Ҳ���ϵĵ���Ϊ1
    %sink=sink&mask;   %ֻ�����Ҳ����λ���ص�����ĵ�
    sink=sink*1e10;   %���������
    %}
    
    %????
    source=source(miny:maxy,minx:maxx);
    sink=sink(miny:maxy,minx:maxx);
    source=source+~mask;%source������ڷ��ص���������
    source=source*1e10;   %��������󣬴��ۼ����ڵ��Ȩֵ���߽���������
    sink=sink*1e10;   %���������
    %����dc��datacost������
    source=source';  %����ת�úͣ�����'��Ϊ�˵õ�һ��Ȩ�أ����е�ǰ��˳���ǰ�ͼ������б�
    sink=sink';
    dc=[source(:)';sink(:)'];%2*(wd*ht)�ľ���  A=[3 4 2;1 5 3;4 7 1];A(:),���һ�У�314457231����ת�ñ����
     
    %������------------------------------------
    %���ø��ڵ�֮��Ĵ���nb   ������������sc��smoothcostƽ����
    %������------------------------------------
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
    nb=sparse(I,J,S,wd*ht,wd*ht);   %ϡ�����I��J�У�ֵΪS
    %ϡ����󣬾��ǿ���������ɢ�ģ�Ҳ����˵������������У��������Ԫ������Ԫ�ء�
    %----------------------
    %��graph-cut����С���ۣ�ͼ��
    %----------------------
    hinc = BK_Create(wd*ht,2*wd*ht);   %����Ŀ����
    BK_SetUnary(hinc,dc);              %dc�������õ�,���������
    BK_SetNeighbors(hinc,nb);          %nb�������õ�,dc,���ƽ����
    edge_cost=BK_Minimize(hinc);       %ʹ������С������
    
    
    L=BK_GetLabeling(hinc);            %�õ��ĸ��������ĸ���
    maskone=zeros(wd*ht,1);
    maskone=im2uint8(maskone)+1;
    B=L-maskone;    
    C=reshape(B,wd,ht);
    
    C=C';                   %source��ͼ 0,sink��ͼ 1,�����ص������mask,ԭmask���������ŵ����ڵĴ�С
    Cc=C.*255;
%     figure;
%     imshow(Cc);
    D=uint8(~C);
    
    %imwrite(C,'warped_img1.png','png');
    
    
    Hg_stack=[Hg_stack;Hg];
    edge_cost_stack=[edge_cost_stack;edge_cost];
    Hg_fitness_stack=[Hg_fitness_stack;Hg_fitness];
    %��ʾfeatrues group
    if size(img1,1) == size(img2,1)
        fprintf('  Showing results of features group');tic;
       % figure;
      %  imshow([img1 img2]);
       % hold on;   
        %����ͼ�϶�IMG1������ƥ��õĵ㻮��Ȧ
     %   plot(keypoint_img1(1,index),keypoint_img1(2,index),'ro','LineWidth',2);
        %����ͼ�϶�IMG2������ƥ��õĵ㻮��Ȧ
     %   plot(keypoint_img2(1,index)+size(img1,2),keypoint_img2(2,index),'ro','LineWidth',2);

        %����ͼ�϶�IMG1�ϵ�i��ƥ��õ�inlier��(�ڵ�)����Ȧ
     %   plot(keypoint_img1(1,featrues_group_idx),keypoint_img1(2,featrues_group_idx),'go','LineWidth',2);
        %����ͼ�϶�IMG2�ϵ�i��ƥ��õ�inlier��(�ڵ�)����Ȧ
      %  plot(keypoint_img2(1,featrues_group_idx)+size(img1,2),keypoint_img2(2,featrues_group_idx),'go','LineWidth',2);

      %  title('features group''s results');
        fprintf('done (%fs)\n',toc);
    end
end
[~,min_index]=min(edge_cost_stack);
Hg=Hg_stack(3*min_index-2:3*min_index,:);



%{
%��ʾfeatrues group
if size(img1,1) == size(img2,1)
    fprintf('  Showing results of features group');tic;
    figure;
    imshow([img1 img2]);
    hold on;    
    %{
    %����ͼ�϶�IMG1������ƥ��õĵ㻮��Ȧ
    plot(data_orig(1,:),data_orig(2,:),'ro','LineWidth',2);
    %����ͼ�϶�IMG2������ƥ��õĵ㻮��Ȧ
    plot(data_orig(4,:)+size(img1,2),data_orig(5,:),'ro','LineWidth',2);
    %}
    
    %����ͼ�϶�IMG1������ƥ��õĵ㻮��Ȧ
    plot(keypoint_img1(1,index),keypoint_img1(2,index),'ro','LineWidth',2);
    %����ͼ�϶�IMG2������ƥ��õĵ㻮��Ȧ
    plot(keypoint_img2(1,index)+size(img1,2),keypoint_img2(2,index),'ro','LineWidth',2);
   
    %����ͼ�϶�IMG1�ϵ�i��ƥ��õ�inlier��(�ڵ�)����Ȧ
    plot(keypoint_img1(1,featrues_group_idx),keypoint_img1(2,featrues_group_idx),'go','LineWidth',2);
    %����ͼ�϶�IMG2�ϵ�i��ƥ��õ�inlier��(�ڵ�)����Ȧ
    plot(keypoint_img2(1,featrues_group_idx)+size(img1,2),keypoint_img2(2,featrues_group_idx),'go','LineWidth',2);
   
    title('features group''s results');
    fprintf('done (%fs)\n',toc);
end
%}
TL = Hg\[1;1;1];    %���Ͻ�
TL = round([ TL(1)/TL(3) ; TL(2)/TL(3) ]);
BL = Hg\[1;size(img2,1);1];     %���½�
BL = round([ BL(1)/BL(3) ; BL(2)/BL(3) ]);
TR = Hg\[size(img2,2);1;1];     %���Ͻ�
TR = round([ TR(1)/TR(3) ; TR(2)/TR(3) ]);
BR = Hg\[size(img2,2);size(img2,1);1];      %���½�
BR = round([ BR(1)/BR(3) ; BR(2)/BR(3) ]);
cw = max([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) - min([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) + 1; %������
ch = max([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) - min([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) + 1; %����߶�
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


%��ʾwarped_img1��warped_img2
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


