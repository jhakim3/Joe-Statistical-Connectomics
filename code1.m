%function void = code1(fin)


a='D:\Downloads D\s01502-2719-07-27-20-26_ABP.txt';
b='D:\Downloads D\s00708-2868-08-02-12-43_ABP.txt';
c='D:\Downloads D\s00214-3084-11-28-16-23_ABP.txt';
d='D:\Downloads D\s00020-2567-03-30-17-47_ABP.txt';

filelist=cellstr([a;b;c;d]);

for k = [1,2,3,4]
    fin=char(filelist(k));


    ABPgross=load(fin);
    abp=ABPgross(:,2);
    res=abpfeature(abp,wabp(abp));

    tc = 125*3600;

    tin=[tc*10,tc*11];

    for j =[1,2]
        figure;
        systime=res(:,1);
        times=find(systime>tin(j));
        for i = times
            temp=systime(i);
        end
        tstart=temp(1);
        tend=temp(20);  

        t=transpose([tstart:tend]);
        plot(t,abp(tstart:tend));
        hold on;

        plot(res([times(1):times(20)],1),res([times(1):times(20)],2),'*');
        plot(res([times(1):times(20)],3),res([times(1):times(20)],4),'*');
        a=res([times(1):times(20)],9);
        plot(a,abp(a),'*');
        a=res([times(1):times(20)],11);
        plot(a,abp(a),'*');
        hold off;
    end
end