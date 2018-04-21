function position=pt_after_H(point,H)
    %右图的点变换到左图上
    point=[point;1];
    position = H\point;    %左上角
    position = round([ position(1)/position(3) ; position(2)/position(3) ]);
end