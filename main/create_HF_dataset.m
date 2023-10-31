function create_HF_dataset(rhod)

% Load in US HF data 
filename= 'data_fig1.csv';
file_iv = '../data/jk/factor_data.csv';    % External instrument series
dataset_name = '../data/jk/data.csv';    % External instrument series

raw = readtable(filename);

% Create dates
for i = 1:size(raw,1)
    days_to_eom(i,:) = [days(datetime(raw{i,1},raw{i,2}, eomday(raw{i,1},raw{i,2}))- datetime(raw{i,1},raw{i,2},raw{i,3}))  ];
    dates(i)=datetime(raw{i,1},raw{i,2},raw{i,3});
    ff4_hf_disc(i) = rhod^(days_to_eom(i,:)) *raw{i,4};
    
    sp500_hf_disc(i) = rhod^(days_to_eom(i,:)) *raw{i,5};
    if i >2 && (raw{i,2} == raw{i-1,2})
        ff4_hf_new(i) = ff4_hf_new(i-1)  +  ff4_hf_disc(i); 
        sp500_hf_new(i) = sp500_hf_new(i-1) +sp500_hf_disc(i);
    else
        ff4_hf_new(i) =  ff4_hf_disc(i);
        sp500_hf_new(i) = sp500_hf_disc(i);
    end
end
tab1 = timetable(raw{:,4}, ff4_hf_new', raw{:,5}, sp500_hf_new' ,'RowTimes',dates ,'VariableNames',{'ff4_hf_orig', ['ff4_hf_disc_w_rhod=' num2str(rhod)],'sp500_hf_orig', ['sp500_hf_disc_w_rhod=' num2str(rhod)]  }' );
tab2 = retime(tab1,'monthly','sum');

tab3 = table(year(tab2.Time), month(tab2.Time), tab2{:,2},tab2{:,4} ,'VariableNames',{'year' 'month' 'ff4_hf' 'sp500_hf'});
writetimetable(tab2, file_iv);    % External instrument series

%Update final dataset
dataset = readtable(dataset_name); 
dataset{128:end,3:4}     = [tab2{:,2},tab2{:,4}];
writetable(dataset,dataset_name)
% summary(tab)
%         plot(dates,ff4_hf_new,'-r');
%         hold on 
%         plot(dates,raw{:,4},'-b');
%         legend('Discounted FFR HF shocks','Original JK FFR HF shocks')
%         corr(ff4_hf_new', raw{:,4})
%         probplot('normal',ff4_hf_new)

