% test effect of broken tracks on reversal analysis via intertime ecdf
clear
close all

% generate some poisson events
N = 40;
R = poissrnd(0.1,N,1000);
figure, hold on
for row = 1:N
    iR{row,:} = diff(find(R(row,:)));
    ecdf(iR{row,:},'Bounds','on','function','survivor')
end

iRcat = horzcat(iR{:});
[f, x] = ecdf(iRcat,'function','survivor'); % this figure checks that concatenation of multiple interreversal times doesn't affect the analysis
stairs(x,f,'g-','LineWidth',2)
set(gca,'YScale','log')
title('individual poisson sequences and concatenated interevent times')
%% check impact of right-censored data
figure, hold on
ecdf(iRcat,'Bounds','off','function','survivor')
set(gca,'YScale','log')
legend('original data')

iRnans = iRcat;
cutoff = 30;
iRnans(iRnans >= cutoff*(1 + rand(size(iRnans)))) = NaN; % randomly sets large values to NaN

figure, hold on
ecdf(iRcat,'Bounds','off','function','survivor')
ecdf(iRnans,'Bounds','off','function','survivor') % this should turn out skewed
set(gca,'YScale','log')
legend('original data','some missing data')

iRcens = iRcat;
censored = isnan(iRnans);
iRcens(censored) = cutoff + rand(size(iRcat(censored))).*(iRcat(censored) - cutoff); % randomly assigns a time in the right range for when censoring happened
figure, hold on
ecdf(iRcat,'Bounds','off','function','survivor')
ecdf(iRnans,'Bounds','off','function','survivor')
ecdf(iRcens,'Bounds','off','function','survivor','censoring',censored) % this should correct for the censoring
set(gca,'YScale','log')
legend('original data','some missing data','missing darta censored')