function filtdata = filter_data_astrodir2026(data)

data_len = length(data);
win_len = 60;

win_mean = zeros(1,data_len-win_len+1);
win_mean(1) = mean(data(1:win_len));
sum_sqr_diff = zeros(1,data_len-win_len+1);
for i = 1 : win_len
    sum_sqr_diff(1) = sum_sqr_diff(1) + (data(i)-win_mean(1))^2;
end

for i = 2 : data_len-win_len+1
    temp_mean = win_mean(i-1) - (data(i-1)-win_mean(i-1))/(win_len-1);
    sum_sqr_diff(i) = sum_sqr_diff(i-1) - (data(i-1) - win_mean(i-1))*(data(i-1) - temp_mean);
    win_mean(i) = temp_mean + (data(i+win_len-1)-temp_mean)/win_len;
    sum_sqr_diff(i) = sum_sqr_diff(i) + (data(i+win_len-1) - temp_mean)*(data(i+win_len-1) - win_mean(i));
end
win_var = sum_sqr_diff ./ (win_len-1);
filtdata = [zeros(1,ceil((win_len-1)/2)),win_var,zeros(1,floor((win_len-1)/2))];
end
