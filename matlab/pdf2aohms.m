function [ms,f,q]=pdf2aohms(xxs,yys,zz,q,dq)
% [ms,f,q]=pdf2aohms(xxs,yys,zz,q,dq)
% This function is to estimate the arbitrary order Hilbert marginal
% spectra using the statistical moments of joint pdf.
% Input
% xxs is the frequency (in log)
% yys is the modulation (in log)
% zz is the joint pdf
% Output
% ms is the corresponding arbitrary order Hilbert marginal spectra
% f is the frequency resolution

if nargin<3
    error('You should input at least three parameters!');
end

if nargin==3
    q=6;
    dq=0.5;
end
if nargin==4
    dq=0.5;
end

if isstruct(xxs)
    f=xxs.F;
else
    f=xxs(:,1);
end
if isstruct(yys)
    amp=abs(yys.A);
    damp=diff(yys.O);
else
    amp=yys(1,:);
    damp=diff([0 amp]);
    zz=zz';
end
Nf=length(f);
Na=length(amp);
qd=0:dq:q;
nord=length(qd);
ms=zeros(nord,Nf);

for i=1:nord
    ms(i,:)=sum(bsxfun(@times,zz',amp.^qd(i).*damp),2);
end
q=qd;