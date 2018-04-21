function fitness=fitness(A,h)
    %A:2N*9
    %h:9*1
    fit=A*h;
    fitness=norm(fit)^2;
end