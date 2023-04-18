clear all;
close all;
clc;
format long;

AtildeVecTran = importdata('AtildeVec24.csv');

r=24;

if true %Plot the entries of Atilde
    figure;
    for j = 1:r^2
        
       % plot( AtildeVec(j,:)  )     %There are some entries that are approx 1, others seem to be clumped around 0
        
%         plot( AtildeVecTran(50:150,j)  )     %There are some entries that are approx 1, others seem to be clumped around 0
        plot( AtildeVecTran(50:end ,j)  )     %There are some entries that are approx 1, others seem to be clumped around 0

        %plot( truncAtildeVec(j,:)  )     %There are some entries that are approx 1, others seem to be clumped around 0

        %plot( AtildeVec(j,:) - AtildeVec(j,1) )
%         hold on;
        title('j=',j)
    
        pause;
    end
    hold off;
end


% figure;
% plot( AtildeVecTran(:,5)  )     %There are some entries that are approx 1, others seem to be clumped around 0
