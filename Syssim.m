function dr = syssim(t,R)


SystemStateVectorLength = length(R);

RR = zeros(SystemStateVectorLength,1);

h = 0;
for i = 1:n
    for h = 1:3 
        RR(((i-1)*6)+h) = R(((i-1)*6)+h+3);
    end
    h = 0;
    for h = 4:6
        RR(((i-1)*6)+h) = G*(sum(m)-m(i))*;
    end
    h = 0;
end


%for i = 1, we want to subtract rx1 from all the rx#




dr = dr';