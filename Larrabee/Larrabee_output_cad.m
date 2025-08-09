%% File: Larrabee_output_cad.m
%--------------------------------
%   SolidWorks用にファイル出力 (Larrabee_output_cad.m)
%--------------------------------

% 翼型datファイルの読み込み
% 翼型はairfoil/の中に入れる
fp = fopen('airfoil/dae51.dat');
airfoil = fgetl(fp);
airfoil = (fscanf(fp,'%f',[2,200]))';
fclose(fp);

% 複数翼型に対応予定
fp = fopen('airfoil/geminism.dat');
airfoil_2 = fgetl(fp);
airfoil_2 = (fscanf(fp,'%f',[2,200]))';
fclose(fp);

% 行列計算からベクトル計算にするためにx,yに分ける
airfoil_x = airfoil(:,1);
airfoil_y = airfoil(:,2);

% 位置rによるairfoilのベクトル作り
airfoil_x = airfoil_x .* linspace(1,1,length(r));
airfoil_y = airfoil_y .* linspace(1,1,length(r));

% 前縁から空力中心の位置に平行移動
airfoil_x = airfoil_x + 0.25;

% 角度変化とコード長に合わせて大きさ変化
deformation_x = chord .* cos(beta);
deformation_y = chord .* sin(beta);

cross_section_x = zeros(length(airfoil), length(r));
cross_section_y = zeros(length(airfoil), length(r));
for i = 1:length(r)
    for j = 1:length(airfoil)
        cross_section_x(j,i) = airfoil_x(j,i) .* deformation_x(i);
        cross_section_y(j,i) = airfoil_y(j,i) .* deformation_y(i);
    end
end

cross_section = zeros(length(airfoil) * length(r), 3);
for i = 1:length(r)
    for j = 1:length(airfoil)
        idx = (i-1)*length(airfoil) + j;
        cross_section(idx,1) = cross_section_x(j,i);
        cross_section(idx,2) = cross_section_y(j,i);
        cross_section(idx,3) = r(i);
    end
end

cad_filename = 'CADfile.csv';
fid = fopen(cad_filename, 'wt');
for i = 1:size(cross_section,1)
    fprintf(fid,'%f,%f,%f\n',cross_section(i,1),cross_section(i,2),cross_section(i,3));
end
fclose(fid);