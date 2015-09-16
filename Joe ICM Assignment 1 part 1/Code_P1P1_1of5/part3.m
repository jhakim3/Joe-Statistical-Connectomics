function part3(patnum,index,abp, onsettimes, feat) 

%This function uses the onsettimes and abp waveform to overlay features on
%a 125 Hz sampled ABP waveform.
%Inputs: a list of patient numbers, the position in that list, the abp
%waveform, the onsettimes of the individual peaks, and an array of the
%features (including systole, diastole etc.)

hrs=[10,11]; %hrs for title

%for k = [1,2,3,4]; %each patient

    tc = 125*3600; %seconds in hour, times 125 Hz sampling

    for j =[1,2]%each hour, 10 and 11
        figure;
        systime=feat(:,1);
        times=find(systime>(tc*hrs(j)));%tc*hrs(j) is 10 or 11 hours
        for i = times
            temp=systime(i);%first measured time that's above 10 hours
        end
        tstart=temp(1);%first heartbeat
        tend=temp(20);  %20th heartbeat

        t=transpose([tstart:tend]);
        plot(t,abp(tstart:tend));%ABP waveform
        hold on;

        %overlaying features on ABP waveform
        %plot(res([times(1):times(20)],1),res([times(1):times(20)],2),'*');%systole marker
        %plot(res([times(1):times(20)],3),res([times(1):times(20)],4),'*');%diastole marker
        a=onsettimes([times(1):times(20)]);
        plot(a,abp(a),'*');%onset of wave pulse time
        a=feat([times(1):times(20)],9);%end of systole, .3 sqrt(RR) method (A)
        plot(a,abp(a),'*');
        a=feat([times(1):times(20)],11);%end of systole, 1st min slope method (B)
        plot(a,abp(a),'*');
        
        title(sprintf('Patient %d hour %d',patnum(index), hrs(j)));
        legend('ABP','onset-time','end systole (A)','end systole (B)');
        xlabel('125 * sec')
        ylabel('ABP (mmHg)')

        hold off;
    end
%end