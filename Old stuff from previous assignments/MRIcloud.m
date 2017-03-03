%% Joe and the Joes - MRICloud (Left Cluster 3)

SPGnormvol = [7855 7549 9448 8297 5735];
SPGalzvol = [8085 6016 10561 13703 6184];
PoCGnormvol = [9364 10124 16581 9458 8734];
PoCGalzvol = [12510 7853 10487 11541 9020];
PrCGnormvol = [9949 9442 13010 11685 10200];
PrCGalzvol = [14235 11582 13352 12704 11605];
SMGnormvol = [8882 10906 14651 9113 8269];
SMGalzvol = [8856 8089 7897 9720 6357];
STGnormvol = [13039 11986 12182 12523 10356];
STGalzvol = [11673 11500 12921 12510 11029];
Insnormvol = [6671 5619 4720 5596 5703];
Insalzvol = [5198 5284 6328 6243 5291];


normalvol = {SPGnormvol PoCGnormvol PrCGnormvol SMGnormvol Insnormvol STGnormvol};
alzvol = {SPGalzvol PoCGalzvol PrCGalzvol SMGalzvol Insalzvol STGalzvol};
titles = {'SPG', 'PoCG', 'PrCG', 'SMG', 'Insula', 'STG'};

for i = 1:length(normalvol)
    meannorm = mean(normalvol{i});
    meanalz = mean(alzvol{i});
    stdnorm = std(normalvol{i});
    stdalz = std(alzvol{i});
    ratio = abs(meannorm - meanalz) / ((stdnorm^2+stdalz^2)/2)^0.5;
    reject = ttest(normalvol{i},alzvol{i});
    meannormstr = strcat('Normal mean:', num2str(meannorm), ' mm3');
    meanalzstr = strcat('Alzheimer mean:', num2str(meanalz), ' mm3');
    stdnormstr = strcat('Normal std:', num2str(stdnorm));
    stdalzstr = strcat('Alzheimer std:', num2str(stdalz));
    ratiostr = strcat('Ratio:', num2str(ratio));
    if reject == 1
        tteststr = strcat('The difference is significant.');
    else
        tteststr = strcat('The difference is not significant.');
    end

    figure;
    histogram(normalvol{i},'FaceColor','g')
    hold on;
    histogram(alzvol{i},'FaceColor','r')
    title(titles{i});
    legend('Normal','Alzheimers');
    xlabel('Volume (mm3)');
    ylabel('Count');
    x = get(gca,'xlim');
    y = get(gca,'ylim');
    text(0.56*(x(2)-x(1))+x(1),0.8*(y(2)-y(1))+y(1),meannormstr);
    text(0.56*(x(2)-x(1))+x(1),0.75*(y(2)-y(1))+y(1),meanalzstr);
    text(0.56*(x(2)-x(1))+x(1),0.7*(y(2)-y(1))+y(1),stdnormstr);
    text(0.56*(x(2)-x(1))+x(1),0.65*(y(2)-y(1))+y(1),stdalzstr);
    text(0.56*(x(2)-x(1))+x(1),0.6*(y(2)-y(1))+y(1),ratiostr);
    text(0.56*(x(2)-x(1))+x(1),0.55*(y(2)-y(1))+y(1),tteststr);
end