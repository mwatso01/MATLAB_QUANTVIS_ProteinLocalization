clear all
close all
clc



counter = 1;

for n = 2:8
    
    % 
    filename = ['all_cell_',num2str(n),'_35kPA_Histogram.mat'];
    contained = exist(filename,'file');
    
    if contained == 2
        load(filename)
        r = .05:.1:.95;

        fitdata = fit(r',mean_bins','poly3');
        p1 = fitdata.p1;
        p2 = fitdata.p2;
        p3 = fitdata.p3;
        p4 = fitdata.p4-mean(mean_bins);
        
        
        figure(1)
        set(gca,'fontsize',24)
        plot(0:.01:1,fitdata(0:.01:1),'linewidth',2)
        hold on
        box off
        axis([0 1 0 1.2])
        daspect([1 1.2 1])
        xlabel('r*')
        ylabel('Relative Intensity')
        drawnow
        rValues = roots([p1 p2 p3 p4]);
        for j = 1:length(rValues)
            if rValues(j) >= 0.01
                rMean = rValues(j);
            end
        end
        r_mean(counter) = rMean;
        cell_number(counter) = n;
        counter = counter + 1;
    end
end

r_mean = r_mean';
cell_number = cell_number';

data_for_excel = [cell_number,r_mean];



