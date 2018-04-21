function border=getborder(mask)
    %获取掩膜的边界
    ht=size(mask,1);
    wd=size(mask,2);
    border=zeros(ht,wd);
    for i=1:ht
        former_point=mask(i,1);%第i行第1个点
        for j=2:wd %从第i行第2个点开始
            point=mask(i,j);    %当前点
            if point==1 && former_point==0 %左侧边
                border(i,j)=1;  %该点为边界点
            end
            if point==0 && former_point==1 %右侧边
                border(i,j-1)=1;    %前一点为边界点
            end
            former_point=point;
        end
    end
    for j=1:wd
        former_point=mask(1,j);%第j列第1个点开始
        for i=2:ht %从第j列第2个点开始
            point=mask(i,j);    %当前点
            if point==1 && former_point==0 %上侧边
                border(i,j)=1;  %该点为边界点
            end
            if point==0 && former_point==1 %下侧边
                border(i-1,j)=1;    %前一点为边界点                
            end
            former_point=point;
        end
    end
end