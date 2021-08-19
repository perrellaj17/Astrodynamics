function rr = syssim(t,R, mu,bc, we)

rr = zeros(6,1);

rmag = sqrt(R(1)^2 + R(2)^2 + R(3)^2)^3;

for i = 1:6
    if i<=3
        rr(i) = R(i+3);
    else 
        rr(i) = -(R(i-3).*mu)/rmag;
    end
    
end







