%%2008-02-19 Pulse-labeling HX of SNase for IUP model:

clear

global k1 k2 k3 k4 k5 k6 k7 k8 kch
global residue_case


fitPara_k=[133.04	5.8526e-005	13.782	9.1504e-006	83.919	0.00010986	2.4617	0.0049486];

k1=fitPara_k(1);  k2=fitPara_k(2);  k3=fitPara_k(3);  k4=fitPara_k(4); 
k5=fitPara_k(5);  k6=fitPara_k(6);  k7=fitPara_k(7);  k8=fitPara_k(8);
kch=0.01;

tp=0.05;
Hu=1-exp(-kch*tp);

i=1;

% for tf=[0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100]; %folding time
for tf=0.01:0.01:3
    
    [t,y] = ode15s('snase_folding_iup',tf,[1 0 0 0]);
    
    sizer=size(y);
    aft_tf=y(sizer(1),:);
    
    result(i,1)=tf;
    for residue_case=1:4
        [t,y] = ode15s('snase_folding_iup_hx', tp, [aft_tf 0]);
        sizer=size(y);
        result(i,residue_case+1)=y(sizer(1),5)/Hu;  %calculate the Amplitude
        subplot(2,2,residue_case)
        semilogx(tf, y(sizer(1),5)/Hu, 'o')
        axis([0.9e-2 5 -0.05 1.05])
        xlabel('Folding time (s)')
        ylabel('Amplitude')
        hold on
    end
    i=i+1;
end

result
