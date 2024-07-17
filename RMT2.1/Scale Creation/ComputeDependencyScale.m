function [ S ] = ComputeDependencyScale( H, sigma )
N = size(H,1);
W=ceil(4*sigma);
%generate gaussian vector
acc = 0;
for i=1:2*W+1
    gaussFunc(i)= exp( - (i-W-1)^2 / (2*sigma^2));
    acc= acc+gaussFunc(i);
end
for i=1:2*W+1
    gaussFunc(i) = gaussFunc(i)/ acc;
end
Ho = H;
H = normalize_12132015(Ho);
HT = normalize_12132015(Ho');
%step2: calculte S matrix
i = floor(size(gaussFunc,2)/2)+1;
S = gaussFunc(1,i)*eye(N,N);
i=i+1;
j=2;
H1 = H;
H2 = HT;
while (j<W+1)
    %H1=normalize(H1);
    S=S+gaussFunc(1,i)*H1+gaussFunc(1,i)*H2;
    i=i+1;
    j=j+1;
    H1 = H1*H;
    H2 = H2*HT;
end





