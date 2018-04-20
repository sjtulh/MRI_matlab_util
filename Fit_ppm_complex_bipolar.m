% Projection onto Dipole Fields (PDF)
%   [p1, dp1, relres, p0]=Fit_ppm_complex(M)
%    
%   output
%   p1 - field map, may need further unwrapping
%   dp1 - a priori error estimate
%   relres - relative residual
%   p0 - initial phase
%
%   input
%   M - a multi-echo and could be a multi-channel dataset
%       echo needs to be the 4th dimension
%       channel needs to be the 5th dimension
%
%   When using the code, please cite 
%   T. Liu et al. MRM 2013;69(2):467-76
%   B. Kressler et al. IEEE TMI 2010;29(2):273-81
%   de Rochefort et al. MRM 2008;60(4):1003-1009
%
%   The coil combination method is similar to
%   MA. Bernstein et al. MRM 1994;32:330-334
%
%   Adapted from a linear fitting created by Ludovic de Rochefort
%   Modified by Tian Liu on 2011.06.01
%   Last modified by Alexey Dimov on 2016.05.12
% -----------------------------------------------------------
%   Fit separate phase offsets for even/odd echoes
%
% ---- By Zhe Liu, 2018/1/3 ---------------------------------


function [p1, dp1, relres, p0, p0_o2e]=Fit_ppm_complex_bipolar(M)

%Modification to handle one echo datasets - assuming zero phase at TE = 0;
%- AD, 05/12/2016
if size(M,4) == 1
    M = cat(4,abs(M),M);
end

if size(M,5)>1
% combine multiple coils together, assuming the coil is the fifth dimension
    M = sum(M.*conj( repmat(M(:,:,:,1,:),[1 1 1 size(M,4) 1])),5);  
    M = sqrt(abs(M)).*exp(1i*angle(M));
end


M= conj(M);
s0=size(M);
L_s0=length(s0);
nechos=size(M,L_s0);

M=reshape(M,[prod(s0(1:L_s0-1)),s0(L_s0)]);
s=size(M);

% Y=angle(M(:,1:min(3,nechos)));
Y_odd=angle(M(:,1:2:6));  % first 3 odd echoes
c=((Y_odd(:,2)-Y_odd(:,1)));
[m ind]=min([abs(c-2*pi),abs(c),abs(c+2*pi)],[],2);
c(ind==1)=c(ind==1)-2*pi;
c(ind==3)=c(ind==3)+2*pi;
for n=1:2
    cd=((Y_odd(:,n+1)-Y_odd(:,n)))-c;
    Y_odd(cd<-pi,(n+1):end)=Y_odd(cd<-pi,n+1:end)+2*pi;
    Y_odd(cd>pi,(n+1):end)=Y_odd(cd>pi,n+1:end)-2*pi;
end
Y_even=angle(M(:,2:2:6));  % first 3 even echoes
c=((Y_even(:,2)-Y_even(:,1)));
[m ind]=min([abs(c-2*pi),abs(c),abs(c+2*pi)],[],2);
c(ind==1)=c(ind==1)-2*pi;
c(ind==3)=c(ind==3)+2*pi;
for n=1:2
    cd=((Y_even(:,n+1)-Y_even(:,n)))-c;
    Y_even(cd<-pi,(n+1):end)=Y_even(cd<-pi,n+1:end)+2*pi;
    Y_even(cd>pi,(n+1):end)=Y_even(cd>pi,n+1:end)-2*pi;
end

% extra unknown, p0_o2e:    phase shift from odd to even echo
A = [1  0  0;
     0  1  1;
     1  0  2;];
     %0  1  3;
     %1  0  4;
     %0  1  5];
Y = reshape(permute(cat(3, Y_odd, Y_even), [1,3,2]), [prod(s0(1:L_s0-1)),6]); 
Y = Y(:, 1:size(A, 1));
ip = A(:,:)\Y(:,:)'; 
p0 = ip(1,:)';
p0_o2e = ip(2,:)';
p1 = ip(3,:)';

dp1 = p1;
tol = norm(p1(:))*1e-4;
iter = 0;
max_iter = 30;

% weigthed least square
% calculation of WA'*WA
v1=zeros(1,nechos); v1(1:2:end) = 1;
v2=zeros(1,nechos); v2(2:2:end) = 1;
v3=(0:(nechos-1));
% v3=floor((0:(nechos-1))/2);

% a11=sum(abs(M).^2.*(ones(s(1),1)*(v1.^2)),2);
% a12=sum(abs(M).^2.*(ones(s(1),1)*(v1.*v2)),2);
% a22=sum(abs(M).^2.*(ones(s(1),1)*(v2.^2)),2);
% % inversion
% d=a11.*a22-a12.^2;
% ai11=a22./d;
% ai12=-a12./d;
% ai22=a11./d;
a11=sum(abs(M).^2.*(ones(s(1),1)*(v1.^2)),2);
a12=sum(abs(M).^2.*(ones(s(1),1)*(v1.*v2)),2);
a13=sum(abs(M).^2.*(ones(s(1),1)*(v1.*v3)),2);
a22=sum(abs(M).^2.*(ones(s(1),1)*(v2.^2)),2);
a23=sum(abs(M).^2.*(ones(s(1),1)*(v2.*v3)),2);
a33=sum(abs(M).^2.*(ones(s(1),1)*(v3.^3)),2);
% co-factor
c11=(1)*(a22.*a33-a23.^2);
c12=(-1)*(a12.*a33-a23.*a13);
c13=(1)*(a12.*a23-a22.*a13);
c22=(1)*(a11.*a33-a13.^2);
c23=(-1)*(a11.*a23-a12.*a13);
c33=(1)*(a11.*a22-a12.^2);
% det
d=a11.*c11 + a12.*c12 + a13.*c13;
% inverse
ai11=c11./d;
ai12=c12./d;
ai13=c13./d;
ai22=c22./d;
ai23=c23./d;
ai33=c33./d;



while ((norm(dp1)>tol) &&(iter<max_iter))
    iter = iter+1;
    W = abs(M).*exp(1i*(p0*v1 + p0_o2e*v2 + p1*v3) );

    % projection
    pr1=sum(conj(1i*W).*(ones(s(1),1)*v1).*(M-W),2);
    pr2=sum(conj(1i*W).*(ones(s(1),1)*v2).*(M-W),2);
    pr3=sum(conj(1i*W).*(ones(s(1),1)*v3).*(M-W),2);

    dp0=real(ai11.*pr1+ai12.*pr2+ai13.*pr3);
    dp0_o2e=real(ai12.*pr1+ai22.*pr2+ai23.*pr3);
    dp1=real(ai13.*pr1+ai23.*pr2+ai33.*pr3);
    dp0(isnan(dp0))=0;
    dp0_o2e(isnan(dp0_o2e))=0;
    dp1(isnan(dp1))=0;
    
    %update
    p0 = p0+dp0;
    p0_o2e = p0_o2e+dp0_o2e;
    p1 = p1+dp1;
    

end

% error propagation
dp1=sqrt(ai33);
dp1(isnan(dp1)) = 0;
dp1(isinf(dp1)) = 0;

% relative residual
res = M - abs(M).*exp(1i*(p0*v1 + p0_o2e*v2 + p1*v3) );
relres = sum(abs(res).^2,2)./sum(abs(M).^2,2);
relres(isnan(relres)) = 0;

p1(p1>pi)=mod(p1(p1>pi)+pi,2*pi)-pi;
p1(p1<-pi)=mod(p1(p1<-pi)+pi,2*pi)-pi;

p0=reshape(p0,s0(1:L_s0-1));
p0_o2e=reshape(p0_o2e,s0(1:L_s0-1));
p1=reshape(p1,s0(1:L_s0-1));
dp1=reshape(dp1,s0(1:L_s0-1));
relres = reshape(relres,s0(1:L_s0-1));
    


