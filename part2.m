function part2(win,abp,feat,otimes,COVals,CVals)

MAP = feat(1+(win-1)/2:end-(win-1)/2,6);

delP = zeros(length(feat),1);
for i = 1:length(feat)
    delP(i) = abp(otimes(i+1))-abp(otimes(i));
end
delP = delP(1+(win-1)/2:end-(win-1)/2);

Tn = feat(1+(win-1)/2:end-(win-1)/2,7);

TPR=MAP./(COVals-CVals.*(delP./Tn));

timemin=(1:length(TPR))';
n = 1;
min = intmax;
while abs(timemin(n) - 720) < min
    min = abs(timemin(n) - 720);
    n = n+1;
end

figure;
plot(timemin(1:n),TPR(1:n));
xlabel('time (min)');
ylabel('TPR (mmHg)');
axis([0 720]);


end
