%%2008-02-19 Pulse-labeling HX of SNase for PPOE2 model:

clear

global k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 kch
global residue_case


fitPara_k=[209.83	23.796	7340.4	1.0953e-009	49.741	0.0018621	3322.9	3.4916	2446	924.22];

k1=fitPara_k(1);  k2=fitPara_k(2);  k3=fitPara_k(3);  k4=fitPara_k(4); 
k5=fitPara_k(5);  k6=fitPara_k(6);  k7=fitPara_k(7);  k8=fitPara_k(8); k9=fitPara_k(9);  k10=fitPara_k(10);
kch=0.01;

tp=0.05;
Hu=1-exp(-kch*tp);

i=1;

% for tf=[0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100]; %folding time
for tf=0.01:0.01:3  
    [t,y] = ode15s('snase_folding_ppoe2',tf,[1 0 0 0 0 0]);
    
    sizer=size(y);
    aft_tf=y(sizer(1),:);
    
    result(i,1)=tf;
    for residue_case=1:3
        [t,y] = ode15s('snase_folding_ppoe2_hx', tp, [aft_tf 0]);
        sizer=size(y);
        result(i,residue_case+1)=y(sizer(1),7)/Hu;  %calculate the Amplitude
        subplot(2,2,residue_case)
        semilogx(tf, y(sizer(1),7)/Hu, 'or')
        axis([0.9e-2 5 -0.05 1.05])
        xlabel('Folding time (s)')
        ylabel('Amplitude')
        hold on
    end
    i=i+1;
end

result
