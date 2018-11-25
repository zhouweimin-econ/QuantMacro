function DC = mydctmx(n)
DC=NaN(n);
for j=0:n-1
    DC(1,j+1)=1/sqrt(n);
    for i=1:n-1
        DC(i+1,j+1)=sqrt(2/(n))*cos(pi*(j+1/2)*i/((n)));
    end
end
end