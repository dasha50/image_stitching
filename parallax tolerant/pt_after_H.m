function position=pt_after_H(point,H)
    %��ͼ�ĵ�任����ͼ��
    point=[point;1];
    position = H\point;    %���Ͻ�
    position = round([ position(1)/position(3) ; position(2)/position(3) ]);
end