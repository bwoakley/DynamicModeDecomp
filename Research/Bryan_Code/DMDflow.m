clear all;
close all;
clc;
format long;

kk = 256;                    %Grid is kk by kk
rows = (kk^2)*2;             %Number of rows in the snapshots

r = 30;                       %Truncate to r singular values

N = 40;                      %How many state snapshots to use (the X1 and X2 will have N-1 columns)
pred = 10;                   %pred = number of time steps forward to predict.
if N+pred>100
    disp('Need more data')
end

flowCase = 1;   %flowCase decides what flow to use. 
                % flowCase = 1 means 'turb'
                % flowCase = 2 means Linear flow (-x,y)
                % flowCase = 3 means Linear flow with time dependent amplitute (1-t/10). 

[ Xpred, Xdmd, stateVecs, Phi, lambda, S ] = DMDpred( pred, N, r, flowCase, kk );


%% Activate block to plot error for various r
%Here I used:
% pred = 90;                   
% N = 10; 
if false
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


%% Activate block to:
% Plot error and
% plot singular values and eigenvalues
%Consider :
%N = 10
%r = 5, 6, or 9
if true
    
    %Now analize the error at each time step
    error = zeros(1,pred);
    %errorC = zeros(1,pred);
    
    for i = 1:pred+1
    
        temp = Xpred(:,i) - stateVecs(:,N-1+i) ;
        error(i) = sumabs(real(temp));
        %errorC(i) = sumabs(imag(temp));
    
    end

    uMax = max(abs(stateVecs(:,1)));
    
    error = error/(rows*uMax);     %Normalize it
    error = error(2:end);   %The first entry is the last entry of X2, so let's drop it.
    
    % figure;
    % plot(error,'o-')
    % hold on;
    % plot(errorC,'*-')
    % hold off;
        
    %Now analize the error at each time step, for Xdmd prediction
    error2 = zeros(1,pred);
    % errorC = zeros(1,pred);
    
    for i = 1:pred
    
        temp = Xdmd(:,i) - stateVecs(:,N+i) ; %Xdmd(1) ~ stateVec(N+1)
        error2(i) = sumabs(real(temp));
    %     errorC(i) = sumabs(imag(temp));
    
    end
    
    error2 = error2/(rows*uMax);     %Normalize it
       
    % figure;
    % plot(error2,'o-')
    % hold on;
    % plot(errorC,'*-')
    % hold off;

    figure;
    plot(error)
    title('error')

    figure;
    plot(error2)
    title('error2')

    figure;
    plot(real(lambda),'*-')
    hold on; 
    plot(imag(lambda),'o-')
    hold off;
    title('Real and imaginary parts of the eigenvalues $\Lambda$ of $\widetilde{A}$','interpreter','latex')
    legend('Real part','Imaginary part')
    
    figure;
    plot(diag(S))
    title('Singular values of $X_1$','interpreter','latex')

end


%% Now fix r and vary N.
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


%% Plot the predicted flow
if true

    tt = 5; %time steps forward from N to observe
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

%% Plot the singular values and DMD modes
if true

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

end


                



