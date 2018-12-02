function F=tran(c_0)
global a_0 beta T R;
F=[-1 R 1]*([(beta*R)^.5 0 0;-1 R 1;0 0 1]^T)*[c_0 a_0 1]';