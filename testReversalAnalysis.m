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

%% check impact of right-censored data
figure, hold on
ecdf(iRcat,'Bounds','on','function','survivor')

iRnans = iRcat;
cutoff = 30;
iRnans(iRnans > cutoff*(1 + rand(size(iRnans)))) = NaN; % randomly sets large values to NaN

ecdf(iRnans,'Bounds','on','function','survivor') % this should turn out skewed

iRcens = iRcat;
censored = isnan(iRnans);
iRcens(censored) = cutoff + rand(size(iRcat(censored))).*(iRcat(censored) - cutoff); % randomly assigns a time in the right range for when censoring happened
ecdf(iRcens,'Bounds','on','function','survivor','censoring',censored) % this should correct for the censoring
set(gca,'YScale','log')