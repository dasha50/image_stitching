function border=getborder(mask)
    %��ȡ��Ĥ�ı߽�
    ht=size(mask,1);
    wd=size(mask,2);
    border=zeros(ht,wd);
    for i=1:ht
        former_point=mask(i,1);%��i�е�1����
        for j=2:wd %�ӵ�i�е�2���㿪ʼ
            point=mask(i,j);    %��ǰ��
            if point==1 && former_point==0 %����
                border(i,j)=1;  %�õ�Ϊ�߽��
            end
            if point==0 && former_point==1 %�Ҳ��
                border(i,j-1)=1;    %ǰһ��Ϊ�߽��
            end
            former_point=point;
        end
    end
    for j=1:wd
        former_point=mask(1,j);%��j�е�1���㿪ʼ
        for i=2:ht %�ӵ�j�е�2���㿪ʼ
            point=mask(i,j);    %��ǰ��
            if point==1 && former_point==0 %�ϲ��
                border(i,j)=1;  %�õ�Ϊ�߽��
            end
            if point==0 && former_point==1 %�²��
                border(i-1,j)=1;    %ǰһ��Ϊ�߽��                
            end
            former_point=point;
        end
    end
end