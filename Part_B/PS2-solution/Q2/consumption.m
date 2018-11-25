function F = consumption(c,k,ap,K,emp,z)
sigma = 2;alfa=0.36;Gamma = 0.65;gamma = 3.5;d0g=0;d1g=1;d0b=0;d1b=1;

F = (c^(-sigma*emp*(1-alfa)*z*(K/((d0g+d1g*log(K))*z + (d0b+d1b*log(K))*(1-z)))^(alfa))/Gamma)^(-1/gamma)*(1-alfa)*z*(K/((d0g+d1g*log(K))*z + (d0b+d1b*log(K))*(1-z)))^(alfa)*emp ...
    - ap + alfa*z*(K/((d0g+d1g*log(K))*z + (d0b+d1b*log(K))*(1-z)))^(alfa-1).*k - c;

end