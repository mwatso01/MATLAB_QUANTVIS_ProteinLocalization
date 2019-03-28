clear all
close all
clc


% Define counters
counter1 = 1;
counter2 = 1;
counter3 = 1;
counter4 = 1;
counter5 = 1;
counter6 = 1;
counter7 = 1;
counter8 = 1;
counter9 = 1;
counter10 = 1;

for n = 2:8                 % MODIFY FOR NUMBER OF CELL OBJECTS
    
    filename = ['all_cell_',num2str(n),'_35kPA_protein.mat']; % Change if you want plot of nuclei
    contained = exist(filename,'file');
    

    
    if contained == 2
        load(filename)

        for i = 1:length(normal_r)
            if normal_r(i) >= 0 && normal_r(i) < 0.1
                bin1(counter1) = norm_intensity(i);
                counter1 = counter1 + 1;
            elseif normal_r(i) >= 0.1 && normal_r(i) < 0.2
                bin2(counter2) = norm_intensity(i);
                counter2 = counter2 + 1;
            elseif normal_r(i) >= 0.2 && normal_r(i) < 0.3
                bin3(counter3) = norm_intensity(i);
                counter3 = counter3 + 1;
            elseif normal_r(i) >= 0.3 && normal_r(i) < 0.4
                bin4(counter4) = norm_intensity(i);
                counter4 = counter4 + 1;
            elseif normal_r(i) >= 0.4 && normal_r(i) < 0.5
                bin5(counter5) = norm_intensity(i);
                counter5 = counter5 + 1;
            elseif normal_r(i) >= 0.5 && normal_r(i) < 0.6
                bin6(counter6) = norm_intensity(i);
                counter6 = counter6 + 1;
            elseif normal_r(i) >= 0.6 && normal_r(i) < 0.7
                bin7(counter7) = norm_intensity(i);
                counter7 = counter7 + 1;
            elseif normal_r(i) >= 0.7 && normal_r(i) < 0.8
                bin8(counter8) = norm_intensity(i);
                counter8 = counter8 + 1;
            elseif normal_r(i) >= 0.8 && normal_r(i) < 0.9
                bin9(counter9) = norm_intensity(i);
                counter9 = counter9 + 1;
            elseif normal_r(i) >= 0.9 && normal_r(i) < 1.0
                bin10(counter10) = norm_intensity(i);
                counter10 = counter10 + 1;
            end
        end
    end
end
bin1_mean = mean(bin1);
bin1_std = std(bin1);

bin2_mean = mean(bin2);
bin2_std = std(bin2);

bin3_mean = mean(bin3);
bin3_std = std(bin3);

bin4_mean = mean(bin4);
bin4_std = std(bin4);

bin5_mean = mean(bin5);
bin5_std = std(bin5);

bin6_mean = mean(bin6);
bin6_std = std(bin6);

bin7_mean = mean(bin7);
bin7_std = std(bin7);

bin8_mean = mean(bin8);
bin8_std = std(bin8);

bin9_mean = mean(bin9);
bin9_std = std(bin9);

bin10_mean = mean(bin10);
bin10_std = std(bin10);

mean_bins1 = [bin1_mean bin2_mean bin3_mean bin4_mean bin5_mean bin6_mean bin7_mean bin8_mean bin9_mean bin10_mean];
mean_bins = mean_bins1/max(mean_bins1);
std_bins = [bin1_std bin2_std bin3_std bin4_std bin5_std bin6_std bin7_std bin8_std bin9_std bin10_std]/max(mean_bins1);
counters = [counter1 counter2 counter3 counter4 counter5 counter6 counter7 counter8 counter9 counter10];
sem_bins = std_bins./sqrt(counters);




figure(1)

errorbar( mean_bins,std_bins,'.','color',[0 0 0],'linewidth',2)
set(gca,'fontsize',18)
hold on
box off
bar(mean_bins,'facecolor',[0 0 0] + 0.5,'linewidth',2)
xticks(1:10)
xticklabels({'0-0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-0.6',...
    '0.6-0.7','0.7-0.8','0.8-0.9','0.9-1.0'})
xtickangle(65)
xlabel('r*')
ylabel('Relative Intensity')
axis([0 11 0 2])
daspect([11-0 2-0 1])

figure(2)
errorbar( mean_bins,sem_bins,'.','color',[0 0 0],'linewidth',2)
set(gca,'fontsize',18)
hold on
box off
bar(mean_bins,'facecolor',[0 0 0] + 0.5,'linewidth',2)
xticks(1:10)
xticklabels({'0-0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-0.6',...
    '0.6-0.7','0.7-0.8','0.8-0.9','0.9-1.0'})
xtickangle(65)
xlabel('r*')
ylabel('Relative Intensity')
axis([0 11 0 1.5])
daspect([11-0 1.5-0 1])



