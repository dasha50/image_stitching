function B=colmove(A,y)
%y>0,обрфy
%y<0,иорфy
[m,n]=size(A);
if y>0
   C=A(1:m-y,:);
   D=zeros(y,n);
   B=cat(1,D,C);
else if y<0
        C=A(abs(y)+1:m,:);
        D=zeros(abs(y),n);
        B=cat(1,C,D);
    end
end
end