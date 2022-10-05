close all;clear;  clc;
figure;set(gcf,'position',[100 100 800 700])
format long; 


k=512;      %Grid size
N = 1010;   %N is largest Th to consider
dt = .002;   %time step
skip = 10;   %how many frames to skip. Plot ever t=skip*dt
            %skip 
            %= sav_flow_int from input.dat 
            %= time_step from res_to_dat.f
Pe = 10^3;   %Peclet number. I think Nu = 10^(-3) in input.dat is the diff const.
lambda = 1;   %rate
sigma = sqrt(2)/(2*pi);   %standard deviation of the initial theta as in periodic.f
theta = zeros(k,k);     %initialize theta analytic


x=linspace(-1,1,k+1);x=x(1:end-1);y=x;
%[Y X]=meshgrid(y,x)
[XX, YY] = meshgrid(x,y);

ss='Lint';

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
    data=reshape(data,k,k,1);

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

    % subplot(2,2,2)
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


            theta(j,k) = sigma^2*exp(lambda*t-.5*(   2*XX(k,j)^2*t*exp(4*lambda*t)*Pe^(-1)...
                + 2*YY(k,j)^2*t*Pe^(-1) + (XX(k,j)^2 + YY(k,j)^2)*sigma^2*exp(2*lambda*t)  )...
                /(  ( sigma^2*exp(2*lambda*t)+2*t*Pe^(-1) )*( sigma^2 + 2*t*exp(2*lambda*t)*Pe^(-1) )   ) )...
                /( ( sigma^2*exp(2*lambda*t)+2*t*Pe^(-1) )^(.5)*( sigma^2 + 2*t*exp(2*lambda*t)*Pe^(-1) )^(.5) );
               
        end
    end

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
     
end

close(v)

%vv = VideoWriter('newfile2.avi');
%open(vv)


%for i = 1000 : 1 : N

    

    %frame2 = getframe(gcf);
    %writeVideo(vv,frame2);
%end

%close(v)














