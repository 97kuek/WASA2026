%% File: readXFLR.m
% xflr5空力解析データの読込 (readXFLR.m)
Relist = {'0.010','0.020','0.030','0.040','0.050','0.100',...
          '0.150','0.200','0.250','0.300','0.350','0.400'};

for i = 1:12
    dir_name   = 'airfoil/';
    foil_name  = 'pelafoil';
    datafile_1 = '_T1_Re';
    datafile_2 = Relist{i};
    datafile_3 = '_M0.00_N9.0.txt';
    dataname   = strcat(dir_name, foil_name, datafile_1, datafile_2, datafile_3);

    fpr = fopen(dataname);      % ファイルオープン
if fpr == -1
    error('readXFLR:FileNotFound', 'ファイルが開けません: %s', dataname);
end

    while true
        buffer = fscanf(fpr,'%c',[1,1]);
        if buffer == '-'
            break;
        end
    end
    buffer = fgetl(fpr); %#ok<*NASGU>

    if i == 1
        data_mat = (fscanf(fpr,'%f',[10,100]))';
    else
        data_mat(:,:,i) = (fscanf(fpr,'%f',[10,100]))';
    end
    fclose(fpr);
end
