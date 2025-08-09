function core
%-------------------------------------------------------------------------------
% Name:         "Wing Profile Manager"
% Purpose:      リブマスター図面出力（DXF）
% Notes:        ファイル先頭にメイン関数、末尾にローカル関数を配置。
%               Octave 記法を MATLAB 互換に修正。浮動小数比較にトレランスを導入。
%               法線方向符号 pm のバグを修正（else 側を -1）。
% Author:       なぽー @Luxion009（原著） / MATLAB 移植・修正: You @97kuek_
% Created:      9/24/2016
%-------------------------------------------------------------------------------

    % ===== メイン処理 =====
    clc;

    %-------------------- ユーザー設定読み込み --------------------
    setting = readmatrix("setting.csv");
    profilenum      = setting(1,1);
    plk_gap         = setting(2,1);
    str_gap         = setting(3,1);
    str_thi         = setting(4,1);
    cap_gap         = setting(5,1);
    te_l            = setting(6,1);
    plkend_sign_u   = setting(7,1);
    plkend_sign_l   = setting(8,1);
    sparpos_sign    = setting(9,1);

    %--------------------------------------------------------------
    % 各翼リブ寸法データ読み込み
    disp("Rib data input");
    wingdata      = readmatrix("rib.csv");
    foil_num      = wingdata(:,1);
    foil_name1    = wingdata(:,2);  % ベース翼型 index
    foil_name2    = wingdata(:,3);  % ミックス相手 index
    foil_chord    = wingdata(:,4);
    foil_ang      = wingdata(:,5);  % [deg]
    foil_mix      = wingdata(:,6);  % 0..1
    foil_sparpos  = wingdata(:,7);  % 相対桁位置（0..1）
    foil_spardia  = wingdata(:,8)./2;  % DXF CIRCLE の半径（入力が直径想定のため /2）
    foil_plk_ux   = wingdata(:,9);
    foil_plk_lx   = wingdata(:,10);

    % 上面ストリンガー位置読み込み（列: リブごと）
    usraw = readmatrix("u_strdata.csv");
    usdata = usraw';
    if ~isempty(usdata)
        usdata(1,:) = [];   % 先頭行を削除（元コード踏襲）
    end

    % 下面ストリンガー位置
    lsraw = readmatrix("l_strdata.csv");
    lsdata = lsraw';
    if ~isempty(lsdata)
        lsdata(1,:) = [];
    end

    disp("Profile data input");
    % 翼型データ読み込み（WPM_i.dat: 2列 [x y]）
    foil_p = cell(profilenum,1);
    fs     = zeros(profilenum,1);
    for i = 1:profilenum
        fn = sprintf("WPM_%d.dat", i);
        fp = fopen(fn,'r');
        if fp<0
            error('ファイルが開けません: %s', fn);
        end
        fgetl(fp); % 先頭行を捨てる（元コード踏襲）
        A = fscanf(fp,'%f',[2,Inf])';
        fclose(fp);
        foil_p{i} = foil_rev(A);  % 翼型反転判断
    end

    disp("Normalize profile data");
    % 翼型データを(0,0)原点に
    foil_p_nor = cell(profilenum,1);
    for i = 1:profilenum
        foil_p_nor{i} = foil_norm(foil_p{i});
    end

    disp("Genelating mix profile data");
    % 各翼データ生成（ミックス）
    mix_foil = cell(numel(foil_num),1);
    for i = 1:numel(foil_num)
        orgnum1   = foil_name1(i,1);
        orgnum2   = foil_name2(i,1);
        orgfoil1  = foil_p_nor{orgnum1};
        orgfoil2  = foil_p_nor{orgnum2};
        mix_foil{i} = foil_gen(orgfoil1, orgfoil2, foil_mix(i,1));
    end

    disp("Foil data output");

    % ===== 出力ループ =====
    tol = 1e-9;  % 浮動小数比較用

    for fnum = 1:numel(foil_num)
        fprintf('%d,', fnum);

        % 入力パラメータ
        chord   = foil_chord(fnum,1);
        sparpos = foil_sparpos(fnum,1);
        spardia = foil_spardia(fnum,1);
        alpha   = foil_ang(fnum,1);    % [deg]
        plk_ux  = foil_plk_ux(fnum,1);
        plk_lx  = foil_plk_lx(fnum,1);

        % ストリンガー位置（列 fnum）取り出し + 末端ゼロ除去
        str_ux = [];
        str_lx = [];
        if ~isempty(usdata)
            str_ux = usdata(:,fnum);
            str_ux = str_ux(~isnan(str_ux));
            zpos = find(str_ux==0,1,'first');
            if ~isempty(zpos)
                str_ux = str_ux(1:zpos-1);
            end
        end
        if ~isempty(lsdata)
            str_lx = lsdata(:,fnum);
            str_lx = str_lx(~isnan(str_lx));
            zpos = find(str_lx==0,1,'first');
            if ~isempty(zpos)
                str_lx = str_lx(1:zpos-1);
            end
        end

        % 出力ファイル
        fpw = fopen(sprintf("WPMout_%d.dxf", fnum),'wt');
        if fpw<0
            error('出力ファイルが開けません: WPMout_%d.dxf', fnum);
        end

        ed_foil = mix_foil{fnum};
        len = size(ed_foil,1);
        chordline = [0 0; 1 0];

        % 特徴点探索
        front_p = 0;  % 前縁 x 最小
        while true
            front_p = front_p + 1;
            if ed_foil(front_p,1) <= min(ed_foil(:,1))
                break;
            end
        end
        [~, u_middle_p] = max(ed_foil(:,2)); % 上面最大 y
        [~, l_middle_p] = min(ed_foil(:,2)); % 下面最小 y

        u_f_gap = front_p - u_middle_p - 1;
        f_l_gap = l_middle_p - front_p - 1;

        uf_buf = ceil(u_f_gap/10); % 無限傾き回避
        lf_buf = ceil(f_l_gap/10);

        sec_s  = u_middle_p + uf_buf;   % 第2セクション開始
        four_s = l_middle_p - lf_buf;   % 第4セクション開始

        % セクション用評価点（x または y 軸）
        xx_1 = linspace(ed_foil(sec_s,1), 1, 200);              % 第1: 上面 x→y
        xx_2 = linspace(ed_foil(sec_s,2), ed_foil(front_p,2),200); % 第2: 上面 逆関数 y→x
        xx_3 = linspace(ed_foil(front_p,2), ed_foil(four_s,2),200);% 第3: 下面 逆関数 y→x
        xx_4 = linspace(ed_foil(four_s,1), 1, 200);              % 第4: 下面 x→y

        % pchip 係数（単調性が崩れると NaN になる可能性あり）
        sf1 = pchip(ed_foil(1:sec_s,1),           ed_foil(1:sec_s,2));      % y(x)
        sf2 = pchip(ed_foil(sec_s:front_p,2),     ed_foil(sec_s:front_p,1));% x(y)
        sf3 = pchip(ed_foil(front_p:four_s,2),    ed_foil(front_p:four_s,1));
        sf4 = pchip(ed_foil(four_s:len,1),        ed_foil(four_s:len,2));

        % 値の評価
        sv1 = ppval(sf1, xx_1);
        sv2 = ppval(sf2, xx_2(1,2:end));
        sv3 = ppval(sf3, xx_3(1,2:end));
        sv4 = ppval(sf4, xx_4(1,2:end));

        % セクションの反転/整形
        xx_1_ed = fliplr(xx_1);   sv1_ed = fliplr(sv1);
        xx_2_ed = xx_2(1,2:end);  % y
        xx_3_ed = xx_3(1,2:end);  % y
        xx_4_ed = xx_4(1,2:end);  % x
        sv2_ed = sv2;  % 逆関数で得た x(y) 値をそのまま使用

        profx = [xx_1_ed, sv2_ed, sv3,   xx_4_ed]';
        profy = [sv1_ed,  xx_2_ed, xx_3_ed, sv4]';

        prof1 = [xx_1_ed', sv1_ed'];    % 上面（x,y）第1
        prof2 = [sv2_ed',  xx_2_ed'];   % 上面（x,y）第2 （x←y）
        prof3 = [sv3',     xx_3_ed'];   % 下面（x,y）第3 （x←y）
        prof4 = [xx_4_ed', sv4'];       % 下面（x,y）第4

        prof  = [profx profy];
        profu = [prof1; prof2]; % 上面
        profl = [prof3; prof4]; % 下面

        % ストリンガー/プランク分布（種）
        sp_ux = zeros(300,1);
        sp_lx = zeros(300,1);
        for i = 1:300
            sp_ux(i,1) = ((i-1)/299)^2;     % 上面寄り密
            sp_lx(i,1) = ((i-1)/299)^1.5;   % 下面寄り密
        end

        % ストリンガー座標 & プランク座標 統合
        sp_ux = [sp_ux; str_ux; plk_ux];
        sp_lx = [sp_lx; str_lx; plk_lx];
        sp_ux = sort(sp_ux);
        sp_lx = sort(sp_lx);

        % y 座標補間
        sp_u = [sp_ux, zeros(numel(sp_ux),1)];
        sp_l = [sp_lx, zeros(numel(sp_lx),1)];
        for i = 1:size(sp_u,1)
            sp_u(i,2) = interp1(profu(:,1), profu(:,2), sp_u(i,1), 'linear','extrap');
        end
        for i = 1:size(sp_l,1)
            sp_l(i,2) = interp1(profl(:,1), profl(:,2), sp_l(i,1), 'linear','extrap');
        end

        % 桁位置変更：キャンバーライン→桁中心点
        cambpoint = (0.001:0.001:0.999)';
        cambline  = zeros(numel(cambpoint),2);
        for i = 1:numel(cambpoint)
            yu = spline(sp_u(:,1), sp_u(:,2), cambpoint(i));
            yl = spline(sp_l(:,1), sp_l(:,2), cambpoint(i));
            cambline(i,1) = cambpoint(i);
            cambline(i,2) = abs(yu - yl);  % 厚み
        end
        spar_camb   = spline(cambline(:,1), cambline(:,2), sparpos);
        spar_camb_h = spar_camb/2;
        sparmidp    = [sparpos, spline(sp_l(:,1), sp_l(:,2), sparpos) + spar_camb_h];

        % 桁基準へ並進
        sp_ur     = [sp_u(:,1)-sparmidp(1), sp_u(:,2)-sparmidp(2)];
        sp_lr     = [sp_l(:,1)-sparmidp(1), sp_l(:,2)-sparmidp(2)];
        profr     = [prof(:,1)-sparmidp(1),  prof(:,2)-sparmidp(2)];
        profur    = [profu(:,1)-sparmidp(1), profu(:,2)-sparmidp(2)];
        proflr    = [profl(:,1)-sparmidp(1), profl(:,2)-sparmidp(2)];
        chordliner = [chordline(:,1)-sparmidp(1), chordline(:,2)-sparmidp(2)];

        % 絶対寸法へスケール
        profr      = profr      .* chord;
        profur     = profur     .* chord;
        proflr     = proflr     .* chord;
        chordliner = chordliner .* chord;
        sp_ur      = sp_ur      .* chord;
        sp_lr      = sp_lr      .* chord;

        % 後縁線（指定長）
        te_x0 = sp_lr(end,1) - 1;  % 初期値
        te_x  = te_x0;
        try
            te_x = fsolve(@(x) endform(x, sp_lr, te_l), te_x0, optimoptions('fsolve','Display','off'));
        catch
            % fsolve なし/失敗時は単純な後縁方向で代替
            % ここでは x=sp_lr(end,1)-te_l とする
            te_x = sp_lr(end,1) - te_l;
        end
        te_y   = interp1(sp_lr(:,1), sp_lr(:,2), te_x, 'linear','extrap');
        te_ang = atan2((sp_lr(end,2)-te_y),(sp_lr(end,1)-te_x)) + 0.5*pi;
        te_p   = [te_x, te_y; te_x + 10*cos(te_ang), te_y + 10*sin(te_ang)];

        % プランク/キャップ/ストリンガー座標生成
        plank_pu = [];
        plank_pl = [];
        cap_pu   = [];
        cap_pl   = [];
        str_pu   = {}; su = 1;
        str_pl   = {}; sl = 1;
        str_peu  = []; % プランク端用（上）
        str_pel  = []; % プランク端用（下）

        % 上面プランク・キャップ
        pu = 0; cu = 0; su_i = 1;
        for i = 1:size(sp_u,1)-1
            if sp_u(i,1) < plk_ux - tol
                if i == 1
                    pu = pu+1; plank_pu(pu,:) = plank_ex(sp_lr(1,:),    sp_ur(i+1,:), sp_ur(i,:),  plk_gap);
                else
                    pu = pu+1; plank_pu(pu,:) = plank_ex(sp_ur(i-1,:),  sp_ur(i+1,:), sp_ur(i,:),  plk_gap);
                end
            elseif abs(sp_u(i,1)-plk_ux) <= tol
                if i == 1
                    pu = pu+1; plank_pu(pu,:) = plank_ex(sp_ur(i,:),    sp_ur(i+1,:), sp_ur(i,:),  plk_gap);
                    cu = cu+1;   cap_pu(cu,:) = plank_ex(sp_ur(i,:),    sp_ur(i+1,:), sp_ur(i,:),  cap_gap);
                else
                    pu = pu+1; plank_pu(pu,:) = plank_ex(sp_ur(i-1,:),  sp_ur(i+1,:), sp_ur(i,:),  plk_gap);
                    cu = cu+1;   cap_pu(cu,:) = plank_ex(sp_ur(i-1,:),  sp_ur(i+1,:), sp_ur(i,:),  cap_gap);
                end
            else
                if i == 1
                    cu = cu+1;   cap_pu(cu,:) = plank_ex(sp_ur(i,:),    sp_ur(i+1,:), sp_ur(i,:),  cap_gap);
                else
                    cu = cu+1;   cap_pu(cu,:) = plank_ex(sp_ur(i-1,:),  sp_ur(i+1,:), sp_ur(i,:),  cap_gap);
                end
            end

            % ストリンガー
            if su_i <= numel(str_ux) && abs(sp_u(i,1) - str_ux(su_i)) <= tol
                str_pu{su_i} = str_ex(sp_ur(i-1,:), sp_ur(i+1,:), sp_ur(i,:), plk_gap, str_gap, str_thi);
                su_i = su_i + 1;
                if su_i > numel(str_ux)
                    su_i = numel(str_ux);
                end
            elseif abs(sp_u(i,1) - plk_ux) <= tol
                if plkend_sign_u == 1
                    str_peu = plankend_ex(sp_ur(i-1,:), sp_ur(i+1,:), sp_ur(i,:), plk_gap, str_gap, str_thi*2);
                end
            end
        end
        % 上面プランクの最後端にキャップ先端を追加（指定時）
        if plkend_sign_u ~= 1 && ~isempty(cap_pu)
            plank_pu(end+1,:) = cap_pu(1,:);
        end

        % 下面プランク・キャップ
        pl = 0; cl = 0; sl_i = 1;
        for i = 1:size(sp_l,1)-1
            if sp_l(i,1) < plk_lx - tol
                if i == 1
                    pl = pl+1; plank_pl(pl,:) = plank_ex(sp_lr(i+1,:), sp_ur(1,:),   sp_lr(i,:),  plk_gap);
                else
                    pl = pl+1; plank_pl(pl,:) = plank_ex(sp_lr(i+1,:), sp_lr(i-1,:),  sp_lr(i,:), plk_gap);
                end
            elseif abs(sp_l(i,1)-plk_lx) <= tol
                if i == 1
                    pl = pl+1; plank_pl(pl,:) = plank_ex(sp_lr(i+1,:), sp_lr(i,:),   sp_lr(i,:),  plk_gap);
                    cl = cl+1;   cap_pl(cl,:) = plank_ex(sp_lr(i+1,:), sp_lr(i,:),   sp_lr(i,:),  cap_gap);
                else
                    pl = pl+1; plank_pl(pl,:) = plank_ex(sp_lr(i+1,:), sp_lr(i-1,:),  sp_lr(i,:), plk_gap);
                    cl = cl+1;   cap_pl(cl,:) = plank_ex(sp_lr(i+1,:), sp_lr(i-1,:),  sp_lr(i,:), cap_gap);
                end
            else
                if i == 1
                    cl = cl+1;   cap_pl(cl,:) = plank_ex(sp_lr(i+1,:), sp_lr(i,:),   sp_lr(i,:),  cap_gap);
                else
                    cl = cl+1;   cap_pl(cl,:) = plank_ex(sp_lr(i+1,:), sp_lr(i-1,:),  sp_lr(i,:), cap_gap);
                end
            end

            % ストリンガー
            if sl_i <= numel(str_lx) && abs(sp_l(i,1) - str_lx(sl_i)) <= tol
                str_pl{sl_i} = str_ex(sp_lr(i+1,:), sp_lr(i-1,:), sp_lr(i,:), plk_gap, str_gap, str_thi);
                sl_i = sl_i + 1;
                if sl_i > numel(str_lx)
                    sl_i = numel(str_lx);
                end
            elseif abs(sp_l(i,1) - plk_lx) <= tol
                if plkend_sign_l == 1
                    str_pel = plankend_ex(sp_lr(i+1,:), sp_lr(i-1,:), sp_lr(i,:), plk_gap, str_gap, str_thi*2);
                end
            end
        end
        % 下面プランクの最後端にキャップ先端を追加（指定時）
        if plkend_sign_l ~= 1 && ~isempty(cap_pl)
            plank_pl(end+1,:) = cap_pl(1,:);
        end

        % --- DXF 出力 ---
        % 図面オフセット
        cntx = 400; cnty = 200; cnty1 = cnty;

        % プランク/キャップ厚による桁中心補正
        if sparpos_sign == 1
            if plk_ux > sparpos
                cnty1 = cnty1 - (plk_gap/2);
            else
                cnty1 = cnty1 - (cap_gap/2);
            end
            if plk_lx > sparpos
                cnty1 = cnty1 + (plk_gap/2);
            else
                cnty1 = cnty1 + (cap_gap/2);
            end
        end

        % 桁中心クロス座標
        sparcrossx = [cntx, 0;    cntx, 841];
        sparcrossy = [0,   cnty1; 1189, cnty1];

        % 迎角線（図面確認用）
        angline = [0,   -400*tan(alpha*pi/180) + cnty1; 
                   1100,  700*tan(alpha*pi/180) + cnty1];

        % DXF セクション開始
        fprintf(fpw,'  0\n'); fprintf(fpw,'SECTION\n');
        fprintf(fpw,'  2\n'); fprintf(fpw,'ENTITIES\n');

        % 桁穴
        circle_ex(fpw, cntx, cnty1, spardia);

        % 外形
        polyline_ex(fpw);
        for i = 1:size(profr,1)
            vertex_ex(fpw, profr(i,1)+cntx, profr(i,2)+cnty);
        end
        seqend_ex(fpw);

        % 上面プランク
        polyline_ex(fpw);
        for i = 1:size(plank_pu,1)
            vertex_ex(fpw, plank_pu(i,1)+cntx, plank_pu(i,2)+cnty);
        end
        seqend_ex(fpw);

        % 下面プランク
        polyline_ex(fpw);
        for i = 1:size(plank_pl,1)
            vertex_ex(fpw, plank_pl(i,1)+cntx, plank_pl(i,2)+cnty);
        end
        seqend_ex(fpw);

        % 上面キャップ
        polyline_ex(fpw);
        for i = 1:size(cap_pu,1)
            vertex_ex(fpw, cap_pu(i,1)+cntx, cap_pu(i,2)+cnty);
        end
        seqend_ex(fpw);

        % 下面キャップ
        polyline_ex(fpw);
        for i = 1:size(cap_pl,1)
            vertex_ex(fpw, cap_pl(i,1)+cntx, cap_pl(i,2)+cnty);
        end
        seqend_ex(fpw);

        % 上面ストリンガー
        for j = 1:numel(str_pu)
            if ~isempty(str_pu{j})
                polyline_ex(fpw);
                for i = 1:size(str_pu{j},1)
                    vertex_ex(fpw, str_pu{j}(i,1)+cntx, str_pu{j}(i,2)+cnty);
                end
                seqend_ex(fpw);
            end
        end
        % 上面プランク端（オプション）
        if plkend_sign_u == 1 && ~isempty(str_peu)
            polyline_ex(fpw);
            for i = 1:size(str_peu,1)
                vertex_ex(fpw, str_peu(i,1)+cntx, str_peu(i,2)+cnty);
            end
            seqend_ex(fpw);
        end

        % 下面ストリンガー
        for j = 1:numel(str_pl)
            if ~isempty(str_pl{j})
                polyline_ex(fpw);
                for i = 1:size(str_pl{j},1)
                    vertex_ex(fpw, str_pl{j}(i,1)+cntx, str_pl{j}(i,2)+cnty);
                end
                seqend_ex(fpw);
            end
        end
        % 下面プランク端（オプション）
        if plkend_sign_l == 1 && ~isempty(str_pel)
            polyline_ex(fpw);
            for i = 1:size(str_pel,1)
                vertex_ex(fpw, str_pel(i,1)+cntx, str_pel(i,2)+cnty);
            end
            seqend_ex(fpw);
        end

        % 後縁線
        polyline_ex(fpw);
        for i = 1:size(te_p,1)
            vertex_ex(fpw, te_p(i,1)+cntx, te_p(i,2)+cnty);
        end
        seqend_ex(fpw);

        % コード基準線
        polyline_ex(fpw);
        for i = 1:size(chordliner,1)
            vertex_ex(fpw, chordliner(i,1)+cntx, chordliner(i,2)+cnty);
        end
        seqend_ex(fpw);

        % 桁中心クロス（縦）
        polyline_ex(fpw);
        for i = 1:size(sparcrossx,1)
            vertex_ex(fpw, sparcrossx(i,1), sparcrossx(i,2));
        end
        seqend_ex(fpw);
        % 桁中心クロス（横）
        polyline_ex(fpw);
        for i = 1:size(sparcrossy,1)
            vertex_ex(fpw, sparcrossy(i,1), sparcrossy(i,2));
        end
        seqend_ex(fpw);

        % 迎角線
        polyline_ex(fpw);
        for i = 1:size(angline,1)
            vertex_ex(fpw, angline(i,1), angline(i,2));
        end
        seqend_ex(fpw);

        % DXF セクション終端
        fprintf(fpw,'  0\n'); fprintf(fpw,'ENDSEC\n');
        fprintf(fpw,'  0\n'); fprintf(fpw,'EOF\n');
        fclose(fpw);

        % クリア（MATLAB は関数スコープなので省略可だが、可読性のため）
        clear chord sparpos spardia alpha plk_ux plk_lx;
        clear ed_foil front_p u_middle_p l_middle_p u_f_gap f_l_gap uf_buf lf_buf;
        clear sec_s four_s xx_1 xx_2 xx_3 xx_4 sf1 sf2 sf3 sf4 sv1 sv2 sv3 sv4;
        clear profx profy prof prof1 prof2 prof3 prof4 profu profl;
        clear sp_ux sp_lx sp_u sp_l cambline spar_camb spar_camb_h sparmidp;
        clear sp_ur sp_lr profr profur proflr chordline chordliner te_x te_y te_ang te_p;
        clear plank_pu plank_pl cap_pu cap_pl str_pu str_pl str_peu str_pel;
        clear sparcrossx sparcrossy angline cntx cnty cnty1;
    end

    disp(" ");
    disp("Output complete");
end

% ============================================================================
% ここから下はローカル関数群
% ============================================================================
function orgfoil = foil_rev(orgfoil)
    % 翼型反転判断：前縁（x 最小）を境に体積を比較して上下を揃える
    fs = 0;
    while true
        fs = fs + 1;
        if orgfoil(fs,1) <= min(orgfoil(:,1)), break; end
    end
    vol1 = sum(orgfoil(1:fs,1));
    vol2 = sum(orgfoil(fs:end,1));
    if vol2 > vol1
        orgfoil = flipud(orgfoil);
    end
end

function ed_orgfoil = foil_norm(orgfoil)
    % 翼型を前縁が (0,0) へ来るよう平行移動
    fs = 0;
    while true
        fs = fs + 1;
        if orgfoil(fs,1) <= min(orgfoil(:,1)), break; end
    end
    ed_orgfoil = orgfoil;
    if orgfoil(fs,1) ~= 0
        ed_orgfoil(:,1) = orgfoil(:,1) - orgfoil(fs,1);
    end
    if orgfoil(fs,2) ~= 0
        ed_orgfoil(:,2) = orgfoil(:,2) - orgfoil(fs,2);
    end
end

function foildata = foil_gen(orgfoil1, orgfoil2, mixture)
    % 2枚の翼型を同一分割 x に補間し、mixture で線形合成（元コード互換: 179点）
    orgfoil{1} = orgfoil1;
    orgfoil{2} = orgfoil2;

    % x 分割（後縁→前縁→後縁。中央の重複を除く）
    x_foil = zeros(180,1);
    for j = 1:90
        x_foil(90+j,1) = ((j-1)/89)^2; % 0..1
    end
    x_foil(1:90,1) = flipud(x_foil(91:180));
    x_foil(90,:)   = [];  % 179 点

    for i = 1:2
        % 前縁探索（x 最小の位置）
        [~, fs] = min(orgfoil{i}(:,1));

        y_i = zeros(numel(x_foil),1);
        % 上面: 1..90 を y(x) で補間
        y_i(1:90)   = interp1(orgfoil{i}(1:fs,1),   orgfoil{i}(1:fs,2),   x_foil(1:90),   'linear','extrap');
        % 下面: 90..179 を y(x) で補間（前縁以降のデータ）
        y_i(90:end) = interp1(orgfoil{i}(fs:end,1), orgfoil{i}(fs:end,2), x_foil(90:end), 'linear','extrap');

        if i==1
            y1 = y_i; %#ok<NASGU>
        else
            y2 = y_i; %#ok<NASGU>
        end
    end

    % 合成
    foildata = [x_foil, mixture.*y1 + (1 - mixture).*y2];
end

function P = plank_ex(fp, rp, cp, gap)
    % プランク/キャップの基準点（外形から法線方向へ gap オフセット）
    % fp: 前点, rp: 後点, cp: 現在点
    ang = atan2((fp(1,2)-rp(1,2)), (fp(1,1)-rp(1,1))) + 0.5*pi;
    if ang < pi
        pm = 1;
    else
        pm = -1; % 元コードのバグ修正
    end
    P(1,1) = cp(1,1) + pm*gap*cos(ang);
    P(1,2) = cp(1,2) + pm*gap*sin(ang);
end

function S0 = str_ex(fp, rp, cp, plk_dep, str_dep, str_thi)
    % ストリンガー四点（底2点・上2点）
    ang = atan2((fp(1,2)-rp(1,2)), (fp(1,1)-rp(1,1))) + 0.5*pi;
    if ang < pi
        pm = 1;
    else
        pm = -1;
    end
    gap = plk_dep + str_dep;

    P0(1,1) = cp(1,1) + pm*gap*cos(ang);
    P0(1,2) = cp(1,2) + pm*gap*sin(ang);

    ang_t = ang - 0.5*pi;   % 厚み方向
    ang_l = ang - pi;       % 長手方向（底の左右）

    % 底点 2
    S0(2,1) = P0(1,1) + pm*(str_thi/2)*cos(ang - 0.5*pi + 0.5*pi);
    S0(2,2) = P0(1,2) + pm*(str_thi/2)*sin(ang - 0.5*pi + 0.5*pi);
    S0(3,1) = P0(1,1) + pm*(str_thi/2)*cos(ang + pi);
    S0(3,2) = P0(1,2) + pm*(str_thi/2)*sin(ang + pi);

    % 上点 2
    S0(1,1) = S0(2,1) + pm*(str_dep)*cos(ang_t);
    S0(1,2) = S0(2,2) + pm*(str_dep)*sin(ang_t);
    S0(4,1) = S0(3,1) + pm*(str_dep)*cos(ang_t);
    S0(4,2) = S0(3,2) + pm*(str_dep)*sin(ang_t);
end

function S0 = plankend_ex(fp, rp, cp, plk_dep, str_dep, str_thi)
    % プランク端ストリンガー
    ang = atan2((fp(1,2)-rp(1,2)), (fp(1,1)-rp(1,1))) + 0.5*pi;
    if ang < pi
        pm = 1;
    else
        pm = -1;
    end
    gap = plk_dep + str_dep;

    P0(1,1) = cp(1,1) + pm*gap*cos(ang);
    P0(1,2) = cp(1,2) + pm*gap*sin(ang);

    % 底 2 点
    S0(2,1) = P0(1,1) + pm*(str_thi/2)*cos(ang - 0.5*pi + 0.5*pi);
    S0(2,2) = P0(1,2) + pm*(str_thi/2)*sin(ang - 0.5*pi + 0.5*pi);
    S0(3,1) = P0(1,1) + pm*(str_thi/2)*cos(ang + pi);
    S0(3,2) = P0(1,2) + pm*(str_thi/2)*sin(ang + pi);

    ang_t = ang - 0.5*pi;

    if S0(3,1) > S0(2,1)
        S0(3,1) = P0(1,1);
        S0(3,2) = P0(1,2);
        % 上 2 点
        S0(1,1) = S0(2,1) + pm*(str_dep)*cos(ang_t);
        S0(1,2) = S0(2,2) + pm*(str_dep)*sin(ang_t);
        S0(4,1) = S0(3,1) + pm*(str_dep + plk_dep)*cos(ang_t);
        S0(4,2) = S0(3,2) + pm*(str_dep + plk_dep)*sin(ang_t);
    else
        S0(2,1) = P0(1,1);
        S0(2,2) = P0(1,2);
        S0(1,1) = S0(2,1) + pm*(str_dep + plk_dep)*cos(ang_t);
        S0(1,2) = S0(2,2) + pm*(str_dep + plk_dep)*sin(ang_t);
        S0(4,1) = S0(3,1) + pm*(str_dep)*cos(ang_t);
        S0(4,2) = S0(3,2) + pm*(str_dep)*sin(ang_t);
    end
end

function val = endform(x, foil, TE_length)
    % 後縁確定用: 点 (x, y(x)) から現後縁点までの距離 = TE_length
    y = interp1(foil(:,1), foil(:,2), x, 'linear','extrap');
    r = sqrt((x - foil(end,1))^2 + (y - foil(end,2))^2);
    val = r - TE_length;
end

% --- DXF 出力補助 ---
function polyline_ex(fpw)
    fprintf(fpw,'  0\n'); fprintf(fpw,'POLYLINE\n');
    fprintf(fpw,'  8\n'); fprintf(fpw,'0\n');
    fprintf(fpw,'  6\n'); fprintf(fpw,'CONTINUOUS\n');
    fprintf(fpw,' 62\n'); fprintf(fpw,'7\n');
    fprintf(fpw,' 66\n'); fprintf(fpw,'1\n');
    fprintf(fpw,' 10\n'); fprintf(fpw,'0\n');
    fprintf(fpw,' 20\n'); fprintf(fpw,'0\n');
    fprintf(fpw,' 30\n'); fprintf(fpw,'0\n');
    fprintf(fpw,' 70\n'); fprintf(fpw,'128\n');
    fprintf(fpw,' 40\n'); fprintf(fpw,'1.0\n');
    fprintf(fpw,' 41\n'); fprintf(fpw,'1.0\n');
end

function seqend_ex(fpw)
    fprintf(fpw,'  0\n'); fprintf(fpw,'SEQEND\n');
    fprintf(fpw,'  8\n'); fprintf(fpw,'0\n');
    fprintf(fpw,'  6\n'); fprintf(fpw,'CONTINUOUS\n');
    fprintf(fpw,' 62\n'); fprintf(fpw,'7\n');
end

function vertex_ex(fpw, xx, yy)
    fprintf(fpw,'  0\n'); fprintf(fpw,'VERTEX\n');
    fprintf(fpw,'  8\n'); fprintf(fpw,'0\n');
    fprintf(fpw,'  6\n'); fprintf(fpw,'CONTINUOUS\n');
    fprintf(fpw,' 62\n'); fprintf(fpw,'7\n');
    fprintf(fpw,' 10\n'); fprintf(fpw,'%f\n', xx);
    fprintf(fpw,' 20\n'); fprintf(fpw,'%f\n', yy);
    fprintf(fpw,' 30\n'); fprintf(fpw,'0\n');
    fprintf(fpw,' 40\n'); fprintf(fpw,'0\n');
    fprintf(fpw,' 41\n'); fprintf(fpw,'0\n');
end

function circle_ex(fpw, xx, yy, rad)
    fprintf(fpw,'  0\n'); fprintf(fpw,'CIRCLE\n');
    fprintf(fpw,'  8\n'); fprintf(fpw,'0\n');
    fprintf(fpw,'  6\n'); fprintf(fpw,'CONTINUOUS\n');
    fprintf(fpw,' 62\n'); fprintf(fpw,'7\n');
    fprintf(fpw,' 10\n'); fprintf(fpw,'%f\n', xx);
    fprintf(fpw,' 20\n'); fprintf(fpw,'%f\n', yy);
    fprintf(fpw,' 30\n'); fprintf(fpw,'0\n');
    fprintf(fpw,' 40\n'); fprintf(fpw,'%f\n', rad);
end
