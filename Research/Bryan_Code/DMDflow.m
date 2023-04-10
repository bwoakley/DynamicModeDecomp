clear all;
close all;
clc;
format long;


make_movie = false;         %Set to true to plot movie of evals
if make_movie                   
    v = VideoWriter('temp.avi');
    open(v)
end

kk = 256;                    %Grid is kk by kk
rows = (kk^2)*2;             %Number of rows in the snapshots

r = 24;                       %Truncate to r singular values

N = 100;                      %How many state snapshots to use (the X1 and X2 will have N-1 columns)

pred = 1;                   %pred = number of time steps forward to predict.
if N+pred>500
    disp('Need more data')
end

no_Windows = 399;            %How many windows of length N to compute
shift = 0;                %Shift starting window. We will look from N = shift to N = shift+no_Windows


flowCase = 1;   %flowCase decides what flow to use. 
                % flowCase = 1 means 'turb'
                % flowCase = 2 means Linear flow (-x,y)
                % flowCase = 3 means Linear flow with time dependent amplitute (1-t/10). 



% Activate block to plot error for various r
 if false
    %Here I used:
    % pred = 90;                   
    % N = 10;
    flowCase = 1; %flowCase decides what flow to use. flowCase = 1 means 'turb'
    figure;
    rVec = [1,2,3,4,5,6,7,8,9]; %Max r to take is N-1
    
    for r = rVec
        
        [error, lambda, S] = DMDpred(pred,N,r,flowCase);
        plot(error)
        hold on;
        
    end
    
    hold off;
    legend('r=1','r=2','r=3','r=4','r=5','r=6','r=7','r=8','r=9')
    title('Error for different r truncations, given N = 10.')
    xlabel('Time steps after N')
    ylabel('Error')


    figure;
    plot(real(lambda),'*-')
    hold on; 
    plot(imag(lambda),'o-')
    hold off;

end

% Now fix r and vary N.
if false

    flowCase = 1; %flowCase decides what flow to use. flowCase = 1 means 'turb'
    r = 5;   %How many modes to use
    
    figure;
    nVec = [10,20,30,40,50]; %Max pred to take is 100-N
    
    for N = nVec
    
        pred = 100-N;
    
        [error, lambda, S] = DMDpred(pred,N,r,flowCase);
        plot(error)
        hold on;
        
    end
    
    hold off;
    legend('N=10','N=20','N=30','N=40','N=50')
    title('Error for different N, given r = 5.')
    xlabel('Number of time steps predicted')
    ylabel('Error')
end

% Plot the predicted flow
if false

    tt = 1; %time steps forward from N to observe
            %End time is T = tt*dt = 5*.02=.1
    upred = reshape( Xdmd(1:kk^2,tt) , kk,kk );
    vpred = reshape( Xdmd(kk^2+1:end,tt) , kk,kk );

    utrue = reshape( stateVecs(1:kk^2,N+tt) , kk,kk );
    vtrue = reshape( stateVecs(kk^2+1:end,N+tt) , kk,kk );


    figure1=figure('Position',[150 200 1200 300]);

    axes1 = axes('Parent',figure1,'Position',[.07 .1 .4 .8]);
    pcolor(real(upred)');shading interp;colorbar;daspect([1 1 1])%forecasting step
    title('Horizontal component of forecasted field at T=0.1')
    hold on;

    axes3 = axes('Parent',figure1,'Position',[.55 .1 .4 .8]);
    pcolor(real(upred)'-utrue');shading interp;colorbar;daspect([1 1 1])%forecasting step
    title('Error in horizontal component at T=0.1')
    hold off;


end

% Plot the singular values and DMD modes
if false

    figure2=figure('Position',[150 200 1400 300]);
    subplot(1,4,1)
    tempS = diag(S);
    semilogy([1:r],real(tempS(1:r)))
    %set(gca,'position',[0.07 .2 .14 .6])
    title('Singular values','fontsize',18)
    xlabel('Modal number');ylabel('Singular value')
    
    Ut=reshape(Phi(1:kk^2,1),kk,kk);
    Ut=Ut';
    subplot(1,4,2)
    pcolor(real(Ut));shading interp;daspect([1 1 1]);colorbar;hold on;
    xlabel('X');ylabel('Y')
    %set(gca,'position',[0.64 .1 .3 .8])
    title('DMD mode 2','fontsize',18)

    Ut=reshape(Phi(1:kk^2,2),kk,kk);
    Ut=Ut';
    subplot(1,4,3)
    pcolor(real(Ut));shading interp;daspect([1 1 1]);colorbar;hold on;
    xlabel('X');ylabel('Y')
    %set(gca,'position',[0.64 .1 .3 .8])
    title('DMD mode 3','fontsize',18)

    Ut=reshape(Phi(1:kk^2,3),kk,kk);
    Ut=Ut';
    subplot(1,4,4)
    pcolor(real(Ut));shading interp;daspect([1 1 1]);colorbar;hold on;
    xlabel('X');ylabel('Y')
    %set(gca,'position',[0.64 .1 .3 .8])
    title('DMD mode 4','fontsize',18)


    %  Plot DMD spectrum
    figure
    theta = (0:1:100)*2*pi/100;
    plot(cos(theta),sin(theta),'k--') % plot unit circle
    hold on, grid on
    scatter(real(lambda),imag(lambda),'ok')
    axis([-1.1 1.1 -1.1 1.1]);


end

%Plot DMD modes
    % U1W1 = zeros(kk, kk);
    % U1W2 = zeros(kk, kk);
    % s1 = zeros(1,no_Windows);
    % s2 = zeros(1,no_Windows);
    % s3 = zeros(1,no_Windows);
    % s4 = zeros(1,no_Windows);



%Now track top number of evals
top = 6;                %Number of top eval to track
topEval = zeros(top, no_Windows);
truncate = top;         %Number of top modes to keep when truncating Atilde

pVec = zeros(2,no_Windows);

for start_index = 1 : no_Windows

    i = start_index;

    [X1, X2, stateVecs] = DMDpullData( pred, N, flowCase, kk, start_index+shift );
    
    [Phi, lambda, b, Xdmd, S, Atilde] = DMD(X1,X2,pred,r);

    if make_movie
        %  Plot DMD spectrum
        figure;
        set(gcf,'position',[100 100 800 700])
        theta = (0:1:100)*2*pi/100;
        plot(cos(theta),sin(theta),'k--');      % plot unit circle
        hold on, grid on

%         scatter(real(lambda),imag(lambda),'ok');
        scale = 5;
        scale_2 = 10;
        sz = scale*(log(abs(b))+scale_2);
        for j = 1:size(sz)
            if sz(j)<0
                sz(j) = 1;
            end
        end
        
        scatter(real(lambda),imag(lambda),sz);

        axis([-1.1 1.1 -1.1 1.1]);
        title(['Window number ', num2str(start_index + shift)])
         
        frame = getframe(gcf);
        writeVideo(v,frame);
        
        close all;
    end

    if false   %Plot DMD modes
        %Atilde
        %[Z_r , D] = eig(Atilde);
    
        %Plot the e.vals.
        ss = lambda;
        s1(i) = real(ss(1));
        s2(i) = real(ss(2));
        s3(i) = real(ss(3));
        s4(i) = real(ss(4));
        
        %Plot the DMD modes
        figure;
        Ut=reshape(Phi(1:kk^2,1),kk,kk);
        Ut=Ut';
        if start_index == 1
            U1W1 = Ut;
        elseif start_index == 2
            U1W2 = Ut;
        end
            
        pcolor(real(Ut));shading interp;daspect([1 1 1]);colorbar;hold on;
        xlabel('X');ylabel('Y')
        set(gca,'position',[0.64 .1 .3 .8])
        title('DMD mode 1','fontsize',18)
    end


    if false   %Consider adjusting vectors in V,W so that they are always sampled from the right half of R^M. This should stop the sign switching of Atilde
        [Z_r , D] = eig(Atilde);
    
        Atilde;
    
        [W, S, V] = svd(X1, 'econ');
    
        for j = 1:10
            for k = 1:N-1
            headW(j,k) = W(j,k);
    
            end
        end
    
        headW
    
        S
    
        V
    
        W_r = W(:, 1:r); % truncate to rank-r
        S_r = S(1:r, 1:r);
        V_r = V(:, 1:r);
    end


    %Find top modes

    bAbs = abs(b);


    [sort_bAbs, idx] = sort(bAbs,'descend');

    %Plot the b's. Is there an exponential drop in the coeff b?
    if true
        R =length(b);
        x = linspace(1,R,R);
        y = log(sort_bAbs);
    %     figure;
    %     i
        p = polyfit(x,y,1);
        
        pVec(1,i) = p(1);
        pVec(2,i) = p(2);
    
    %     f = polyval(p,x);
    %     plot(x,y,'o',x,f,'-') 
    end

    for modeNumber = 1:top
        topEval(modeNumber,i) = lambda(idx(modeNumber));
    end

    truncIdx = idx(1:truncate);

%        sort_bAbs

%     b(truncIdx)
%     truncbAbs = bAbs(truncIdx)

    if truncate<r
        aa = bAbs(idx(truncate));
        bb = bAbs(idx(truncate+1));
        if abs(aa - bb) < aa/100
            disp('Warning: truncAtilde may be complex; the truncation contains an eigenmode but not its conjugate pair.')
        end
    end

    %Now truncate to these top modes by truncating lambda and Phi
    trunclambda = lambda(truncIdx);
    truncPhi = Phi(:, truncIdx);

    % Reconstructing DMD
    Nf=pred;                                        %number of time steps forward 
    x1=X2(:,end);                                   %Take the last snap shot from X2
    bbb=truncPhi\x1;                  %notice that bbb is slightly different from b            %Obtain initial condition from last snapshot
    temporal=zeros(truncate,Nf);
    for ii=1:Nf
        temporal(:,ii)=bbb.*trunclambda.^ii;%Raise by eigenvalue to power of timestep
    end
    truncXdmd=truncPhi*temporal;




    %Compare predictions to DNS

    % Activate block to Plot error. truncDMD vs DNS vs DMD
    if false
        
        error = zeros(1,pred);
        
        for iii = 1:pred
            tempX = Xdmd(:,iii) - stateVecs(:,N+iii) ;
            error(iii) = sumabs(real(tempX));
%             errorC(iii) = sumabs(imag(tempX));
        
        end
    
        uMax = max(abs(stateVecs(:,1)));
        
        error = error/(rows*uMax);     %Normalize it
        %error = error(2:end);   %The first entry is the last entry of X2, so let's drop it.
                   
        figure;
        plot(error)
%         hold on;
%         plot(errorC)
        title('DMD vs DNS')
%         hold off;



        error2 = zeros(1,pred);
        
        for iii = 1:pred
            tempX = truncXdmd(:,iii) - stateVecs(:,N+iii) ;
            error2(iii) = sumabs(real(tempX));
%             error2C(iii) = sumabs(imag(tempX));
        
        end
    
        uMax = max(abs(stateVecs(:,1)));
        
        error2 = error2/(rows*uMax);     %Normalize it
        %error = error(2:end);   %The first entry is the last entry of X2, so let's drop it.
                   
        figure;
        plot(error2)
%         hold on;
%         plot(error2C)
        title('truncDMD vs DNS')
%         hold off;





        error3 = zeros(1,pred);
        
        for iii = 1:pred
            tempX = truncXdmd(:,iii) - Xdmd(:,iii);
            error3(iii) = sumabs(real(tempX));
%             error3C(iii) = sumabs(imag(tempX));
        
        end
    
        uMax = max(abs(stateVecs(:,1)));
        
        error3 = error3/(rows*uMax);     %Normalize it
        %error = error(2:end);   %The first entry is the last entry of X2, so let's drop it.
                   
        figure;
        plot(error3)
%         hold on;
%         plot(error3C)
        title('truncDMD vs DMD')
%         hold off;


        figure;
        plot(error,'-')
        hold on;
        plot(error2,'--')
        hold off;
        legend('DMD','truncDMD')


    %     figure;
    %     plot(real(lambda),'*-')
    %     hold on; 
    %     plot(imag(lambda),'o-')
    %     hold off;
    %     title('Real and imaginary parts of the eigenvalues $\Lambda$ of $\widetilde{A}$','interpreter','latex')
    %     legend('Real part','Imaginary part')
    %     
    %     figure;
    %     plot(diag(S))
    %     title('Singular values of $X_1$','interpreter','latex')
    
    end



    %Now to improve the algorithm by learning on the evals
    
    


end


%Plot change in DMD mode 1
if false
    pcolor(real(U1 - U2));shading interp;daspect([1 1 1]);colorbar;hold on;
    xlabel('X');ylabel('Y')
    title('Difference of U1 and U2','fontsize',18)
    
    %Plot evals over windows
    tt = linspace(1,N);
    figure;
    scatter(tt,s1,'ok')
    
    figure;
    plot(s2)
    
    figure;
    plot(s3)
    
    figure;
    plot(s4)
end


%Plot top evals
plot_eval = false;    %Set true to plot top evals
if plot_eval
    
%     linS = {'-','--','-','--','-','--','-','--','-','--','-','--','-','--','-','--','-','--','-','--'};
    linS = {'-','--','-','--','-','--'};
    tim = linspace(1,no_Windows,no_Windows);
    
    
    figure;
    absTopEval = abs(topEval);
    
    for modeNumber = 1:top
        
        plot( tim, absTopEval(modeNumber,:),'linestyle', linS(modeNumber) )
        hold on;
    
    end
    axis([1 no_Windows .95 1.1])
    hold off;
    title('Modulus of top evals','fontsize',18)
%     legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20')
    legend('1','2','3','4','5','6')
    


    figure;
    angleTopEval = abs(angle(topEval));
    
    for modeNumber = 1:top
        
        plot( tim, angleTopEval(modeNumber,:),'linestyle', linS(modeNumber) )
        hold on;
    
    end
    hold off;
    title('Angle of top evals','fontsize',18)
%     legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20')
    legend('1','2','3','4','5','6')
    
end

figure;
plot(pVec(1,:))
title('Exp rate of b')

figure;
plot(pVec(2,:))
title('y int of log of b')

