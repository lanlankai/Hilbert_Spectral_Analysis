function [xxs,yys,zz]=jpdfcm2(imf,ma,minw)
% [xxs,yys,zz]=jpdfcm2(imf,ma,minw)
% This function is to estimate the joint pdf  of instantaneous frequency
% \omega and modulation A
% Input
% imf is the imf from mhs or EMD
% ma is the range of the modulation
% Output
% xxs is the x grid in log
% yys is the y grid in log
% zz is the joint pdf in numbers
% 
% 
% To show the result:
%  surfc(xm,ym,zz,'Facecolor','interp','Edgecolor','none','facelighting','phong'); alpha(0.5)
%  contourf(xm,ym,zz,10)
% 
% See also jointpdf jpdf jpdfca

if nargin==1
    error('You should specify the range of the modulation.');
end
% check the input
if iscell(imf)
    Nt=length(imf);
    cflag=1;% the flag of the cell
    if ischar(imf{1})
        N=4;
    else
        N=1;
    end
else
    cflag=0;
    N=1;
    Nt=1;
end

if length(ma)==1
    ma=[-8 ma(1)];
end

if ma(2)>=10
    ma(2)=log10(ma(2));
end

switch cflag
    case 1  % cell from mhs or other codes
        if nargin==2
            NL=length(imf{N});
            minw=2/NL;
        end

        for i=N:Nt
            [xxs,yys,zz1]=jpdffc(imf{i},ma,minw);
            if i==N
                zz=zz1;
            else
                zz=zz+zz1;
            end
        end
        zz=zz/NL/(Nt-N+1);

    case 0 % imf from emd
        NL=length(imf);
        if nargin==2
            minw=2/NL;
        end

        
        
        [xxs,yys,zz]=jpdffc(imf,ma,minw);
        zz=zz/NL;
end

        









function [xxs,yys,zz]=jpdffc(imf,ma,minw)
% This function is to calculate the joint pdf of the instantaneous
% frequency and the instantaneous of amplitude
% Input
% imf is IMF mode from mhs or emd
% ma is the range of the amplitude
% Output
% xxs is the frequency resolution
% yys is the amplitude resolution
% zz is the joint pdf
% Written by Yongxiang Huang 29/06/2008
% 
% See also jpdf jointpdf jointpdfc

[instf,insta] = instantaneous(imf); % get the instantaneous frequency and amplitude
% [instf,insta]=NHHT(imf);


% amplitude
wmx=ma(2)-ma(1);
% wmn=10^-4;
NF=fix(wmx*10);
fa=linspace(ma(1),ma(2),NF);
FA.A=10.^fa;
dfa=fa(2)-fa(1);
fa=[fa-dfa/2 fa(end)+dfa/2];
fa=10.^fa;
FA.O=fa;




%*****
wmx=log10(0.5);
wmn=log10(minw);
NF=fix((wmx-wmn)*20);
ff=linspace(wmn,wmx,NF);
dff=ff(2)-ff(1);
FF.F=10.^ff;
ff=[ff-dff/2 ff(end)];
ff=10.^ff;
FF.O=ff;

zz=jointpdfc2(instf,insta,ff,fa); % get the joint pdf


% 
xxs=FF;
yys=FA;








% % % % % % % % % % % % % % % Normalized Hilbert transform% % % % % % % % % % % % % % % % % % % % % %
    function [Instf,Insta]= NHHT(imf)
        %         This function is to perform the so-called Normalized Hilbert
        %         Spectral Analysis

      
        Nimf=size(imf);
        Nimf(1)=Nimf(1)-1; % exclude the trend
        Nimf(2)=Nimf(2)-2;
        
        Instf=zeros(1,Nimf(1)*Nimf(2));
        Insta=zeros(1,Nimf(1)*Nimf(2));
        for i=1:Nimf
            [instf1,insta1]=NHSA(imf(i,:));
            Instf((i-1)*Nimf(2)+1:i*Nimf(2))=instf1;
            Insta((i-1)*Nimf(2)+1:i*Nimf(2))=insta1;
        end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % ?????
        function [Inst,FA]=NHSA(Ci)
            % This function is to perform the so-called Normalized Hilbert Spectral
            % Analysis
            %  see:  	ON INSTANTANEOUS FREQUENCY
            % NORDEN E. HUANG; ZHAOHUA WU; STEVEN R. LONG; KENNETH C. ARNOLD; XIANYAO CHEN; KARIN BLANK
            % Page 177 - 229
            % Input
            % Ci is one dimension of IMF
            % Output
            % Inst is the instantaneous frequency estimated by DQ method.
            % FA is the modulation
            %
            % Written by Yongxiang Huang
            % See also instantaneous,DQinstantaneous,DQ
            %
            
            
            
     
            
            Nc=length(Ci);% get the length of the IMF
            FA=ones(1,Nc); % initial the modulation
            
            t=1:Nc;
            [indmin,indmax,indzer]=extr(abs(Ci),t); % identify the extrema points
            
            if length(indmax)>4
                
                % perform the modulation and phase decomposition
                for i=1:10
                    [indmin,indmax,indzer]=extr(abs(Ci),t); % identify the extrema points
                    sl=spline(t(indmax),abs(Ci(indmax)),t); % estimate the envelope
                    Ci=Ci./sl;
                    FA=FA.*sl;
                end
            end
            FA=FA(2:Nc-1);
            xi= FA<0;
            FA(xi)=0;
            [Inst,insta]=instantaneous(Ci); % calculate the instantaneous frequency



% % % % % % % % % % Intantaneous Frequency% % % % % % % % % % % % % % %
            function [instf,insta] = instantaneous(data)
                % [instf,insta] = instantaneous(imf)
                % This function is to get the instantaneous frequency and amplitude
                % Input
                % imf is the imf mode from emd
                % Output
                % instf is the instantaneous frequency
                % insta is the instantaneous moudulation
                %
                % See also mrgnllog mrgnlrslog

                
                [knb,npt] = size(data);
                if knb>npt
                    data=data';
                    [knb,npt] = size(data);
                end
                
                data=hilbtm(data')'; % get the analysitical siginal by performing Hilbert transform
                
                if knb>1
                    f=zeros(knb-1,npt-2);
                    for i=1:knb-1 % exclude the residual
                        tmp=instfreq(data(i,:)')'; % Get the instantaneous frequency
                        f(i,:)=tmp;
                    end
                    
                    Nf=size(f);
                    a=abs(data);% get the amplitude
                    
                    instf=zeros(1,Nf(1)*Nf(2));
                    insta=zeros(1,Nf(1)*Nf(2));
                    
                    for i=1:knb-1 %exclude the residual
                        instf((i-1)*Nf(2)+1:i*Nf(2))=f(i,1:Nf(2));
                        insta((i-1)*Nf(2)+1:i*Nf(2))=a(i,2:Nf(2)+1);
                    end
                    
                else
                    tmp=instfreq(data(1,:)')'; % Get the instantaneous frequency
                    if length(tmp)<2
                        instf=0;
                        insta=0;
                    else
                        instf=tmp;
                        insta=abs(data(2:end-1));
                    end
                end
                
                
                
