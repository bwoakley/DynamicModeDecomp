close all;clear;  clc;
figure;set(gcf,'position',[100 100 800 700])
format long; 


kk=256;      %Grid size
N = 1050;   %N is largest Th to consider
dt = .002;   %time step
skip = 10;   %how many frames to skip. Plot ever t=skip*dt
            %skip 
            %= sav_flow_int from input.dat 
            %= time_step from res_to_dat.f
Pe = 10^3;   %Peclet number. I think Nu = 10^(-3) in input.dat is the diff const.
lambda = 1;   %rate
sigma = sqrt(1/(8*pi^2));   %standard deviation of the initial theta as in periodic.f
theta = zeros(kk,kk);     %initialize theta analytic


x=linspace(-.5,.5,kk+1);x=x(1:end-1);y=x;
%[Y X]=meshgrid(y,x)
[XX, YY] = meshgrid(x,y);

ss='Line';

v = VideoWriter('newfile.avi');
open(v)


for i = 1000 : 1 : N
     
    if i<10
    st=strcat('../Cases/',ss,'/Th0000',num2str(i));
    elseif i<100
    st=strcat('../Cases/',ss,'/Th000',num2str(i));
    elseif i<1000
    st=strcat('../Cases/',ss,'/Th00',num2str(i));
    elseif i<10000
    st=strcat('../Cases/',ss,'/Th0',num2str(i));
    else
    st=strcat('../Cases/',ss,'/Th',num2str(i));
    end


    fid=fopen(st,'rb');
    data=fread(fid,[1 1],'*float');
    data=fread(fid,[1 inf],'*double');
    
    
    fclose(fid);
    data=reshape(data,kk,kk,1);

    %If the data is close enough to 0, but is negative, then just round up to 0.
    if false
        %min(min(data))
        for j = 1:1:  size(data,1)
            for k = 1 : 1 :  size(data,2)
                if data(j,k) > -10^(-5) & data(j,k) < 0
                    data(j,k) = 0;
                end
            end
        end
        %min(min(data))
        %max(max(data))
    end

    %Activate to make movie of sim vs analytic
    if true
        %put theta sim in first plot
        subplot(1,2,1)
    
        pcolor(x,x,data(:,:,1)');shading interp;colorbar;
        set(gca,'fontsize',18)
        daspect([1,1,1])
        caxis([0 1])
        title('theta simulation')
        hold on;
        %C=contourc(x,x,data(:,:,1)',[.5 .5]);
        
    
         %What does below do?
        if false
            ll=size(C,2);
            li=0;
            ld=1;
            while ld<ll
                li=li+1;
                L(li)=C(2,ld);
                xc=C(1,ld+1:ld+L(li));
                yc=C(2,ld+1:ld+L(li));
                for ki=1:length(xc)-1
                    xkc(ki*2-1)=xc(ki);xkc(ki*2)=(xc(ki)+xc(ki+1))/2;
                    ykc(ki*2-1)=yc(ki);ykc(ki*2)=(yc(ki)+yc(ki+1))/2;
                end
                plot(xc,yc,'k');
                ld=ld+L(li)+1;
            end
        end
    
   
        t = 1*dt*skip*(i-1000);
    
        for j = 1 : size(data,1)
            for k = 1 : size(data,2)
        
                A = sigma^2+2*t*Pe^(-1)*exp(2*lambda*t);
                B = sigma^2*exp(2*lambda*t)+2*t*Pe^(-1);
    
                %theta(j,k) = sigma^2*exp(lambda*t-.5*(   2*XX(j,k)^2*t*exp(4*lambda*t)*Pe^(-1)...
                %    + 2*YY(j,k)^2*t*Pe^(-1) + (XX(j,k)^2 + YY(j,k)^2)*sigma^2*exp(2*lambda*t)  )...
                %    /(  ( sigma^2*exp(2*lambda*t)+2*t*Pe^(-1) )*( sigma^2 + 2*t*exp(2*lambda*t)*Pe^(-1) )   ) )...
                %    /( ( sigma^2*exp(2*lambda*t)+2*t*Pe^(-1) )^(.5)*( sigma^2 + 2*t*exp(2*lambda*t)*Pe^(-1) )^(.5) );
                
                %Or, more simply:
                %theta(j,k) = sigma^2*exp(lambda*t-.5*(1/(A*B))*( XX(j,k)^2*exp(2*lambda*t)*A ...
                %    +YY(j,k)^2*B ) )/(sqrt(A*B)); 
                
                %Somehow, I'm getting:
                theta(j,k) = sigma^2*exp(lambda*t-.5*(1/(A*B))*( XX(j,k)^2*exp(2*lambda*t)*A ...
                    +XX(j,k)*exp(2*lambda*t)*A -t*Pe^(-1)*(1/(2*sigma^2))*A ...
                    +YY(j,k)^2*B + YY(j,k)*B -t*Pe^(-1)*(1/(2*sigma^2))*exp(2*lambda*t)*B ) )/(sqrt(A*B)); 

                
            end
        end
    
        %min(min(theta))
        %max(max(theta))


        subplot(1,2,2)
        pcolor(x,x,theta);shading interp;colorbar;
    
        set(gca,'fontsize',18)
        daspect([1,1,1])
        caxis([0 1])
        title('theta analytic')
    
        %drawnow
        %Instead of drawnow, save the movie in F
    
        frame = getframe(gcf);
        writeVideo(v,frame);
    
    
        %writeVideo(v,data(:,:,1)')
        %writeVideo(v,C)
        
        hold off;
    end %This ends the if true/false to make movie

end

close(v)


%Just plot the initial condition.
if false

    for i = 1000 : 1 : 1000
             
        if i<10
        st=strcat('../Cases/',ss,'/Th0000',num2str(i));
        elseif i<100
        st=strcat('../Cases/',ss,'/Th000',num2str(i));
        elseif i<1000
        st=strcat('../Cases/',ss,'/Th00',num2str(i));
        elseif i<10000
        st=strcat('../Cases/',ss,'/Th0',num2str(i));
        else
        st=strcat('../Cases/',ss,'/Th',num2str(i));
        end
    
    
        fid=fopen(st,'rb');
        data=fread(fid,[1 1],'*float');
        data=fread(fid,[1 inf],'*double');
        
        
        fclose(fid);
        data=reshape(data,kk,kk,1);

        subplot(1,2,1)
    
        pcolor(x,x,data(:,:,1)');shading interp;colorbar;
        set(gca,'fontsize',18)
        daspect([1,1,1])
        caxis([0 1])
        title('theta simulation')
        hold on;
    
        for j = 1 : size(data,1)
            for k = 1 : size(data,2)
    
    
                theta(j,k) = exp(-1*(XX(j,k)^2 + YY(j,k)^2)/(2*sigma^2) );
            end
        end

        subplot(1,2,2)
        pcolor(x,x,theta);shading interp;colorbar;
    
        set(gca,'fontsize',18)
        daspect([1,1,1])
        caxis([0 1])
        title('theta analytic')

         

    end
    
    ERROR = zeros(kk,kk);     

    for j = 1 : size(data,1)
            for k = 1 : size(data,2)
            
                ERROR(j,k) = data(j,k) - theta(j,k);
    
            end
    end
    min(min(ERROR))
    max(max(ERROR))

end %End of if true/false



%Hi







