a=project4AvRon;
[max,maxid]=findpeaks(a);

posmax=[];
for i=1:length(maxid)
    if(max(i)>0)
        posmax =[posmax;maxid(i)];
    end
end


rec_lengths=[0,0];
for i = 1:(length(posmax)-1)
    if(posmax(i+1)- posmax(i) > 700)
        
        rec_lengths = [rec_lengths; posmax(i+1),(posmax(i))];
    end
end
%xy = size(rec_lengths);
%rec_lengths(1,2) = posmax(1);
rec_lengths(1,1) = posmax(1);
rec_lengths = [rec_lengths; 0,posmax(end)];
%rec_lengths(end,2) = posmax(end);

for i=1:length(rec_lengths)-1
rec_lengths(i,2)=rec_lengths(i+1,2);
end
rec_lengths=rec_lengths(1:end-1,:);

%plot(a)
%hold on

%plot(rec_lengths(:,1),50,'*');
%plot(rec_lengths(:,2),54,'*');
%plot(posmax,60,'*')

q=zeros(length(rec_lengths),1);
for i = 1:length(rec_lengths)
    q(i)=rec_lengths(i,2)-rec_lengths(i,1);
end

fig=figure;
[counts, bins] = hist(q,25);
plot(bins, counts); %# get a line plot of the histogram
title('Threshold = 0.01')
xlabel('Active time (sec * 10^-5)')
ylabel('frequency, bins=25')
hgexport(fig,'-clipboard')