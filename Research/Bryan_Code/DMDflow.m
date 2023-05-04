clear all;
close all;
clc;
format long;
tic

kk = 256;                    %Grid is kk by kk

r = 24;                       %Truncate to r singular values

N = 100;                      %How many state snapshots to use (the X1 and X2 will have N-1 columns)

pred = 10;                   %pred = number of time steps forward to predict.
if N+pred>500
    disp('Need more data')
end

no_Windows = 1;            %How many windows of length N to compute
shift = 0;                %Shift starting window. We will look from N = shift to N = shift+no_Windows


flowCase = 4;   %flowCase decides what flow to use. 
                % flowCase = 1 means 'turb'
                % flowCase = 2 means Linear flow (-x,y)
                % flowCase = 3 means Linear flow with time dependent amplitute (1-t/10). 
                % flowCase = 4 means sin(r) e^(-kr) transported down the x-axis 

if flowCase == 4        % predicting a scalar
    isScalar = true;
    rows = kk^2;        %Number of rows in the snapshots
else                    % predicting a flow
    isScalar = false;
    rows = (kk^2)*2;    %Number of rows in the snapshots
end



make_movie = false;         %Set to true to plot movie of evals
if make_movie                   
    v = VideoWriter('temp.avi');
    open(v)
end



make_movie_prediction = true;         %Set to true to plot movie of predictions from DMD vs DNS
if make_movie_prediction                   
    v = VideoWriter('temp.avi');
    open(v)

    if isScalar == false
        disp('Warning: make_movie_prediction is only coded for scalar fields right now')
    end
end


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
top = 8;                %Number of top eval to track. Usually want top = even
topEval = zeros(top, no_Windows);
truncate = top;         %Number of top modes to keep when truncating Atilde

pVec = zeros(2,no_Windows);
AtildeVec = zeros(r^2,no_Windows);
truncAtildeVec = zeros(r^2,no_Windows);

plotBs = false; %Plot the b's. Is there an exponential drop in the coeff b?
MatNorm = false; %Compute matrix norm?

if true    %Iterate DMD over multiple windows.
    for start_index = 1 : no_Windows
    
        i = start_index;
    
        [X1, X2, stateVecs] = DMDpullData( pred, N, flowCase, kk, start_index+shift );
        

         if false % Centering data by subtracting mean
            
            if i == 1
                disp('*Centering data by subtracting mean*')
                
            end

            oneVec = ones(N-1, 1);

            mu1 = X1*oneVec/(N-1);
            X1 = X1 - mu1*oneVec' ;

            mu2 = X2*oneVec/(N-1);
            X2 = X2 - mu2*oneVec' ;
        end


        [Phi, lambda, b, Xdmd, S, Atilde] = DMD(X1,X2,pred,r);
        
        AtildeVec(:,i) = reshape(Atilde,[],1);
    
    %     Atilde
    
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
    
        if false %Turn this on to trucate Atilde to the top modes
    
            b;
    
            %Find top modes
            bAbs = abs(b);
       
            [sort_bAbs, idx] = sort(bAbs,'descend');
        
            if plotBs
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
        
            
            %We can also truncate Atilde itself by restricting to those top freq. Perhaps this truncAtilde will sign switch less...
            [Z_r , D] = eig(Atilde);
    %         Z_r*D/Z_r - Atilde
    %         [Ztilde, Dtilde, Ztilde2] = svd(Atilde, 'econ')
    %         Phi'*Phi
            truncZ = Z_r(:,truncIdx);
            tempD = diag(D);
            truncD = diag( tempD(truncIdx));
            truncAtilde = truncZ*truncD*pinv(truncZ);
            
            imagTruncAtilde = sum(sum(imag(truncAtilde)));
            if imagTruncAtilde > 1e-10
                disp('truncAtilde seems to be complex, taking only the real part.')
            end
            truncAtilde = real(truncAtilde);
    
            truncAtildeVec(:,i) = reshape(truncAtilde,[],1);
    
    
    %         E = sum(sum(abs( Atilde-truncAtilde )));
    
    %         [ truncZtilde , truncDtilde ] = eig(truncAtilde);
            
    %         truncZ - truncZtilde
    %         truncD-truncDtilde
    
    
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
        
        if MatNorm    %Compare Frob norm of diff: abs vs raw
            if i > 1
                
                errs(i-1)=sum(sum(abs(Atilde-prevAtilde)));
                errs2(i-1)=sum(sum(abs((abs(Atilde)-abs(prevAtilde)))));
               
            end
            prevAtilde = Atilde;
        end
        
        if false   %For fixed window, use the fixed DMD modes and see how the flow snapshots' coeff evolve over that window.
    
            windowBvec = zeros(r,N);
    
            for j = 1:N
    
                x1 = stateVecs(:,j);
                windowB = Phi\x1;
                windowBvec(:,j) = windowB; 
                
            end
    
            figure;
            %oneVec = ones(1,N);
    
            for j = 1:r
    
                plot(abs(windowBvec(j,:)))       
                hold on;
    
                %expFit = abs(lambda(j))*oneVec;
                expFit = zeros(1,N);
                for m = 1:N
    
                    expFit(m) = abs(lambda(j))^(m-1)*abs(windowBvec(j,1)) ;
    
                end
                plot(expFit,'--')
                hold on;
                %expFit(1)
    
                AvgErrorBs = sum(abs( abs(windowBvec(j,:)) - expFit))/N;
                %disp
                message = ['Average error of coeff b is ', num2str(AvgErrorBs),  ' for r = ', num2str(j) ];
                disp(message)
            end
            hold off;
    
            title('Coeff b along the window of length N=100')
            
             %windowBvec(1:10,1)
    
    
            figure;     
            for j = 1:r
    
                plot(log10(abs(windowBvec(j,:))))       
                hold on;
    
                %expFit = abs(lambda(j))*oneVec;
                expFit = zeros(1,N);
                for m = 1:N
    
                    expFit(m) = abs(lambda(j))^(m-1)*abs(windowBvec(j,1)) ;
    
                end
                plot(log10(expFit),'--')
                hold on;
                %expFit(1)
    
%                 sum(abs( abs(windowBvec(j,:)) - expFit))
    
            end
%              ylim([-15 6])

            hold off;
    
            title('Log of coeff b along the window of length N=100')
    
    
    
        end
    


        if true         % If we know the future coefficients b, how well can we use the dynamic modes to approximate DNS?

            windowBvec = zeros(r,N+pred);
            
            projErrorVec = zeros(1,N+pred);

            for j = 1:N+pred
    
                x1 = stateVecs(:,j);
                windowB = Phi\x1;
                windowBvec(:,j) = windowB; 
                
                exactB_Xdmd = Phi*windowB;

                projErrorVec(j) = sum(abs( exactB_Xdmd  - x1  ));
            
                %complexError = sum(abs(imag(exactB_Xdmd)))

                if j > N
                    if make_movie_prediction
    
                        figure;
                        set(gcf,'position',[100 100 800 700])
                        
                        theta = reshape(abs(exactB_Xdmd),kk,kk);
    
                        pcolor(theta); shading interp; drawnow;
        
                        frame = getframe(gcf);
                        writeVideo(v,frame);
                        
                        close all;
    
                    end
                end


            end
    
            %Take average error over spatial domain
            uMax = max(abs(stateVecs(:,N)));
            normalizeConst = rows*uMax ;

            projErrorVec = projErrorVec/normalizeConst; %Normalize the error   


            if false      %Plot error and log of error

                figure;
                plot(projErrorVec)
                title('Error of the projection onto dynamic modes')
    
                figure;
                plot(log10(projErrorVec),'-o','MarkerSize',4)
                title('Log of the error of the projection onto dynamic modes')

            end



            %Now compare to DMD
            error = zeros(1,pred+1);
                    
            for iii = 1:pred+1
                tempX = Xdmd(:,iii) - stateVecs(:,N+iii-1) ;
                error(iii) = sum(sum(abs(abs(tempX))));
            %             errorC(iii) = sumabs(imag(tempX));
            
            end
        
           
            error = error/normalizeConst;     %Normalize the error

            if false     % Plot error
                figure;
    
                plot(error)
                hold on;
    
                endProjErrorVec = projErrorVec(N:end);
                plot(endProjErrorVec)
    
                title('Error of DMD and exactB DMD ')
                legend('DMD', 'exactB DMD')
    
                hold off;
    
                % Plot log of error
                figure;
    
                plot(log10(error))
                hold on;
    
                plot(log10(endProjErrorVec))
    
                title('Log of error of DMD and exactB DMD ')
                legend('DMD', 'exactB DMD')
    
                hold off;

            end


            if false    % How much of an improvement is exactB DMD over DMD? Plot the difference:
                
                figure;
                plot( abs(error - endProjErrorVec) )
                title('DMD error minus exactB DMD error')
    
                figure;
                plot( log10(abs(error - endProjErrorVec) ) )
                title('Log of DMD error minus exactB DMD error')
    
                figure;
                plot( error./endProjErrorVec )
                title('Ratio of DMD to exactB DMD')

            end


            


        end


    end
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
    %axis([1 no_Windows .95 1.1])
    hold off;
    title('Modulus of top evals','fontsize',18)
%     legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20')
%     legend('1','2','3','4','5','6')
    legend('1','2')
 


    figure;
    angleTopEval = abs(angle(topEval));
    
    for modeNumber = 1:top
        
        plot( tim, angleTopEval(modeNumber,:),'linestyle', linS(modeNumber) )
        hold on;
    
    end
    hold off;
    title('Angle of top evals','fontsize',18)
%     legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20')
%     legend('1','2','3','4','5','6')
    legend('1','2')



    %Decompose top eval pair using fft
    if true
        Eval1 = absTopEval(1,:);

        L = no_Windows;
        t = (1:L)/L;        % Let's rescale time to be just one time unit...

        S = Eval1;

%         figure;
%         plot(t,S)
        
        Y = fft(S);

        rescaledY = Y/L;


        Yabs = abs(rescaledY); %Just looking at modulus. Which frequencies have most energy?

        last = ceil(L/2); %Last positive freq
           
        Ypos = Yabs(1:last);
        
        [sort_Ypos, idx_Ypos] = sort(Ypos,'descend');   %idx_Ypos contains the indices's with most energy
        top_freq = idx_Ypos - 1;
                
        
        
        shiftY = fftshift(Y);
        
        if mod(L,2) == 0
            middle = ceil(L/2) + 1;
        else 
            middle = ceil(L/2);
        end
        
        
        
        
        
        kkkk = 6;     %Use top kk frequencies to approx Eval1. Note kk should be < L/2
        if kkkk > last
            disp('kkkk too large')
        end
        
        Sapprox = zeros(1,L);
        
        for j = 1:kkkk
        
            freq = top_freq(j)
            Sapprox = Sapprox +  shiftY(middle+freq)*exp((2*pi*1i)*freq*t) + shiftY(middle-freq)*exp((-2*pi*1i)*freq*t);
          
        
        end

        Sapprox = Sapprox/(2*L);

        tempSapprox = Sapprox;     %For some reason, Sapprox is shifted one to the left. Let's shift it back...
        Sapprox(1) = S(1);
        for j = 2:L
            Sapprox(j) = tempSapprox(j-1);
        end
                
        % Sapprox
        
        figure; 
        plot(t, S)
        hold on;
        plot(t, Sapprox,'--')
        hold off;
        legend('orig','approx')


    end

    
end


%Plot the b's. Is there an exponential drop in the coeff b?
if plotBs

    figure;
    plot(pVec(1,:))
    title('Exp rate of b')
    
    figure;
    plot(pVec(2,:))
    title('y int of log of b')

end

%Plot MatNorm of diff
if MatNorm
    figure;
    plot(errs)
    xlabel('Time');title('$||\tilde{A}_{i+1}-\tilde{A}_i||$','Interpreter', 'latex');ylabel('error');
    set(gca,'fontsize',18)
    
    
    figure;
    plot(errs2)
    xlabel('Time');title('$|| |\tilde{A}_{i+1}| - |\tilde{A}_i| ||$','Interpreter', 'latex');ylabel('error');
    set(gca,'fontsize',18)

%     figure;
%     plot(errs)
%     hold on;
%     plot(errs2)
%     hold off;
%     legend('$||\tilde{A}_{i+1}-\tilde{A}_i||$','$|| |\tilde{A}_{i+1}| - |\tilde{A}_i| ||$','Interpreter', 'latex')
    
%     figure;
%     plot( abs( errs-errs2 ) )
%     title('errs-errs2')
end

if false %Plot the entries of Atilde
    figure;
    for j = 1:r^2
        
        plot( AtildeVec(j,:)  )     %There are some entries that are approx 1, others seem to be clumped around 0
    
        %plot( truncAtildeVec(j,:)  )     %There are some entries that are approx 1, others seem to be clumped around 0

        %plot( AtildeVec(j,:) - AtildeVec(j,1) )
        hold on;
    
    end
    hold off;
end

if false  %Save the Atilde to a csv file
    AtildeVecTran = AtildeVec';
    
    writematrix(AtildeVecTran)
end

if false    %Plot the error of DMD vs DNS vs futureAtilde+oldFlow




    % DMD
    start_index = 1;
    [X1, X2, stateVecs] = DMDpullData( pred, N, flowCase, kk, start_index+shift );
    [Phi, lambda, b, Xdmd, S, Atilde] = DMD(X1,X2,pred,r);

    %Analize error:
    error = zeros(1,pred+1);
            
    for iii = 1:pred+1
        tempX = Xdmd(:,iii) - stateVecs(:,N+iii-1) ;
        error(iii) = sum(sum(abs(abs(tempX))));
    %             errorC(iii) = sumabs(imag(tempX));
    
    end

    uMax = max(abs(stateVecs(:,N)));
    
    error = error/(rows*uMax);     %Normalize it
    %error = error(2:end);   %The first entry is the last entry of X2, so let's drop it.
               
    figure;
    plot(error)
    %         hold on;
    %         plot(errorC)
    title('DMD vs DNS')
    %      hold off;










    % Now compute DMD with alternative algorithms:
    % 'current' refers to using future/preditions of Atilde, with old modes Phi,
    %       giving 'improvedDMD'. However, this does not seem to be an improvement.
    if false

    
       % This is hard coded to N=100, r=24
       if r ~= 24 || N ~= 100
           disp('Warning: should set N=100, r=24')
       end
    
    
        %    currentAtilde = Atilde;  %Take the Atilde from the last loop as the currentAtilde.
    
       AtildeVec24Tran = importdata('AtildeVec24.csv');
       AtildeVec24 = AtildeVec24Tran';
    
        %Check that currentAtilde = AtildeVec24(:,shift+1)
        %     AtildeTemp = AtildeVec24(:,shift+1);
        % 
        %     AtildeTemp2 = reshape(AtildeTemp,r,r);
        % 
        %     errorAtilde = currentAtilde - AtildeTemp2;
        %     sum(sum(abs(errorAtilde)))
    
    
        %Instead of recomputing Phi, just use the old DMD modes Phi, but only update lambda
            %Now I need the most recent V_r and X2. We copy the following from DMD.m
            %             [W, S, V] = svd(X1, 'econ');
            %         
            %             adjustV = true;
            %             if adjustV
            %                 sizeX1 = size(X1);
            %                 for j = 1:sizeX1(2)    %adjusting vectors in V,W so that they are always sampled from the right half of R^M. This should stop the sign switching of Atilde
            %                     if V(1,j) < 0
            %                         V(:,j) = -1*V(:,j);
            %                         W(:,j) = -1*W(:,j);
            %                     end
            %                 end
            %             else
            %                 disp('consider adjusting V')
            %             end
            %                 
            %             ss=diag(S);
            %             ind=find(ss>ss(1)*1e-10);
            %         
            %             %  Compute DMD (Phi are eigenvectors)
            %             max_r = length(ind);  % truncate at r modes, which is singular value >1e-10 max
            %             r_orig = r;
            %             r = min([r_orig, size(W,2), max_r]);
            %             if r < r_orig
            %                 disp(['Singular value(s) < ss(1)*1e-10 detected, truncating from r=', num2str(r_orig), ' to r=', num2str(r)])
            %             end
            %             
            %             W_r = W(:, 1:r); % truncate to rank-r
            %             S_r = S(1:r, 1:r);
            %             V_r = V(:, 1:r);
            %End of copying from DMD.m
    
    
        origLambda = lambda;
        origAtilde = Atilde;
        origX1 = X1;
        origX2 = X2;
        origPhi = Phi;
        

        checkXdmd = zeros(rows,pred+1);
        improvedXdmd = zeros(rows,pred+1);
    
       
       
    
        for j = 1:pred
    
            j;
    
            AtildeTemp = AtildeVec24(:,shift+j);
            currentAtilde = reshape(AtildeTemp,r,r) ;
            
    
            if true
                alpha = 1;  % Using the future Atilde seems to be an overcorrection, and performs worse than DMD. So, instead use a weighted average...? 
                    % Let alpha be the weight coeff. 
                    % alpha = 1 uses only future Atilde, 
                    % alpha = 0 uses only origAtilde
    
                currentAtilde = alpha*currentAtilde + (1-alpha)*origAtilde;
    
            end
            
    
    
            [currentZ_r , currentD] = eig(currentAtilde);
        
            
    
            %Phi = X2 * V_r / S_r * Z_r; % DMD modes        %Maybe: instead of recomputing Phi, just use the old DMD modes Phi, but only update lambda
    
            currentLambda = diag(currentD);           
    
    
            
            
            if j == 1
                temp1 = currentLambda;
            else
                temp2 = currentLambda;
                max(abs(temp1-temp2));
    
                diffLam = sum(sum(abs(temp1-temp2)));
                avgDiffLam = diffLam/(r^2);
                avgDiffLam;
                
                temp1 = currentLambda;
            end
            %j
            %lambda
    
    
    
            if j == 1                           
    
    
                x1=origX2(:,end);                        %Take the last snap shot from X2
                b=origPhi\x1;                          %Obtain initial condition from last snapshot
                improvedXdmd(:,1)=origPhi*b ;              %j==1 should be just retrieving the final snapshot. The error with DNS should be approx 0. This is a check.
                improvedXdmd(:,2)=origPhi*(b.*currentLambda) ;      %j==2 uses AtildeVec24(:,shift+1) to predict
                
                checkXdmd(:,1)=origPhi*b ; 
                checkXdmd(:,2)=origPhi*(b.*origLambda) ;
    
                prevB = b.*currentLambda;
                checkPrevB = b.*origLambda;
            else
                %x1=improvedXdmd(:,j);                  %Use the approximation at j...
                %b=origPhi\x1;                            %together with the approx DMD modes from AtildeVec24(:,shift+j). 
                                                     %     Note that these are not the true DMD modes because they use old flow snapshots in V_r, S_r, and X2 
                %improvedXdmd(:,j+1)=origPhi*(b.*currentLambda) ;    %Use AtildeVec24(:,shift+j) to predict the next flow snapshot
                
    
    
                %Just evolve the prev coeff using the new lambda
                improvedXdmd(:,j+1)=origPhi*(prevB.*currentLambda) ;    %Use AtildeVec24(:,shift+j) to predict the next flow snapshot
                
                checkXdmd(:,j+1)=origPhi*(checkPrevB.*origLambda) ;   
    
                prevB=prevB.*currentLambda;
                checkPrevB = checkPrevB.*origLambda;
    
            end
            
            
    
            diffB = sum(sum(abs(prevB-checkPrevB)));
            percentDiffB = max( abs(prevB-checkPrevB)./abs(prevB)  );
    
            
        end
    
        %sum(sum(abs(Xdmd-checkXdmd)))
    
        %Now plot the error of impDMD vs DNS
        error2 = zeros(1,pred+1);
                
        for iii = 1:pred+1
            tempX = improvedXdmd(:,iii) - stateVecs(:,N+iii-1) ;
            error2(iii) = sum(sum(abs(abs(tempX))));
            %             errorC(iii) = sumabs(imag(tempX));
        
        end
    
        uMax = max(abs(stateVecs(:,N)));
        
        error2 = error2/(rows*uMax);     %Normalize it
        %error = error(2:end);   %The first entry is the last entry of X2, so let's drop it.
                   
        figure;
        plot(error2)
        %         hold on;
        %         plot(errorC)
        title('improvedDMD vs DNS')
        %         hold off;
    end




    % Let's try an alternative algorithm 'feedback':
    % 'feedback' refers to plugging our predictions of the flow snapshot back into feedbackX1 and feedbackX2,
    %       using that to generate feedbackAtilde and then a new predition.
    if true
        feedbackXdmdPred = zeros(rows,pred+1);
        for j = 1:pred
        
            start_index = j;
            
            if j == 1
    
                [feedbackX1, feedbackX2, feedbackStateVecs] = DMDpullData( 1, N, flowCase, kk, start_index+shift );
                [feedbackPhi, feedbackLambda, feedbackB, feedbackXdmd, feedbackS, feedbackAtilde] = DMD(feedbackX1,feedbackX2,1,r);
                
                feedbackXdmdPred(:,1:2) = feedbackXdmd(:,1:2);
    
            else
    
                feedbackX1 = feedbackX2;
                feedbackX2 = [feedbackX2(:,2:end), feedbackXdmd(:,2) ] ;
    
                [feedbackPhi, feedbackLambda, feedbackB, feedbackXdmd, feedbackS, feedbackAtilde] = DMD(feedbackX1,feedbackX2,1,r);
            
                feedbackXdmdPred(:,j+1) = feedbackXdmd(:,2);
    
    
            end
    
        end
        %Now plot the error of fbDMD vs DNS
        error3 = zeros(1,pred+1);
                
        for iii = 1:pred+1
            tempX = feedbackXdmdPred(:,iii) - stateVecs(:,N+iii-1) ;
            error3(iii) = sum(sum(abs(abs(tempX))));
            %             errorC(iii) = sumabs(imag(tempX));
        
        end
    
        uMax = max(abs(stateVecs(:,N)));
        
        error3 = error3/(rows*uMax);     %Normalize it
        %error = error(2:end);   %The first entry is the last entry of X2, so let's drop it.
                   
        figure;
        plot(error3)
        %         hold on;
        %         plot(errorC)
        title('feedbackDMD vs DNS')
        %         hold off;


        figure;
        plot(error)
        hold on;
        plot(error3)
        title('Error')
        legend('DMD', 'feedbackDMD')
        hold off;

        ratioDMDtoFeedback = error./error3;
        figure;
        plot(ratioDMDtoFeedback)
        title('Ratio of DMD to FeedbackDMD error')

    end












    % Let's try an alternative algorithm 'feedbackPlusFuture' or 'fbpf'
    % 'feedback' refers to plugging our predictions of the flow snapshot back into feedbackX1 and feedbackX2,
    %       using that to generate feedbackPhi.
    % 'PlusFuture' refers to using the future/predicted Atilde instead of feedbackAtilde
    if false
    
        % This is hard coded to N=100, r=24
        if r ~= 24 || N ~= 100
            disp('Warning: should set N=100, r=24')
        end
    
        AtildeVec24Tran = importdata('AtildeVec24.csv');
        AtildeVec24 = AtildeVec24Tran';
    
        fbpfXdmdPred = zeros(rows,pred+1);
        for j = 1:pred
        
            start_index = j;
            AtildeTemp = AtildeVec24(:,shift+j);
            currentAtilde = reshape(AtildeTemp,r,r) ;
            
            if j == 1
    
                [feedbackX1, feedbackX2, feedbackStateVecs] = DMDpullData( 1, N, flowCase, kk, start_index+shift );
                [feedbackPhi, feedbackLambda, feedbackB, feedbackXdmd, feedbackS, feedbackAtilde] = DMD(feedbackX1,feedbackX2,1,r);
                
                fbpfXdmdPred(:,1:2) = feedbackXdmd(:,1:2);
    
            else
    
                feedbackX1 = feedbackX2;
                feedbackX2 = [feedbackX2(:,2:end), feedbackXdmd(:,2) ] ;
    
                X1 = feedbackX1;
                X2 = feedbackX2;
                Atilde = currentAtilde;
    
                % DMD
                [W, S, V] = svd(X1, 'econ');
                
                adjustV = true;
                if adjustV
                    sizeX1 = size(X1);
                    for j = 1:sizeX1(2)    %adjusting vectors in V,W so that they are always sampled from the right half of R^M. This should stop the sign switching of Atilde
                        if V(1,j) < 0
                            V(:,j) = -1*V(:,j);
                            W(:,j) = -1*W(:,j);
                        end
                    end
                else
                    disp('consider adjusting V')
                end
    
                
                
                ss=diag(S);
                ind=find(ss>ss(1)*1e-10);
                
                %  Compute DMD (Phi are eigenvectors)
    %             max_r = length(ind);  % truncate at r modes, which is singular value >1e-10 max
    %             r_orig = r;
    %             r = min([r_orig, size(W,2), max_r]);
    %             if r < r_orig
    %                 disp(['Singular value(s) < ss(1)*1e-10 detected, truncating from r=', num2str(r_orig), ' to r=', num2str(r)])
    %             end
                
                W_r = W(:, 1:r); % truncate to rank-r
                S_r = S(1:r, 1:r);
                V_r = V(:, 1:r);
                
                Atilde = W_r' * X2 * V_r / S_r; % low -rank dynamics
                
                [Z_r , D] = eig(Atilde);
                Phi = X2 * V_r / S_r * Z_r; % DMD modes
                lambda = diag(D); % discrete -time eigenvalues
                                    
               
                % Reconstructing DMD
                Nf=1;%number of time steps forward 
                x1=X2(:,end);%Take the last snap shot from X2
                b=Phi\x1;%Obtain initial condition from last snapshot
                temporal=zeros(r,Nf+1);
                for i=1:Nf+1
                    temporal(:,i)=b.*lambda.^(i-1);%Raise by eigenvalue to power of timestep
                end
                Xdmd=Phi*temporal;
    
    
    
            
                fbpfXdmdPred(:,j+1) = Xdmd(:,2);
    
    
                feedbackXdmd = Xdmd;
    
            end
    
        end
        %Now plot the error of fbDMD vs DNS
        error4 = zeros(1,pred+1);
                
        for iii = 1:pred+1
            tempX = fbpfXdmdPred(:,iii) - stateVecs(:,N+iii-1) ;
            error4(iii) = sum(sum(abs(abs(tempX))));
            %             errorC(iii) = sumabs(imag(tempX));
        
        end
    
        uMax = max(abs(stateVecs(:,N)));
        
        error4 = error4/(rows*uMax);     %Normalize it
        %error = error(2:end);   %The first entry is the last entry of X2, so let's drop it.
                   
        figure;
        plot(error4)
        %         hold on;
        %         plot(errorC)
        title('feedbackPlusFutureDMD vs DNS')
        %         hold off;
    end

end



if make_movie || make_movie_prediction
    close(v)
end



toc