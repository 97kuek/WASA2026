function core
%-------------------------------------------------------------------------------
% Name:         "Wing Profile Manager"
% Purpose:      ���u�}�X�^�[�}�ʏo�́iDXF�j
% Notes:        �t�@�C���擪�Ƀ��C���֐��A�����Ƀ��[�J���֐���z�u�B
%               Octave �L�@�� MATLAB �݊��ɏC���B����������r�Ƀg�������X�𓱓��B
%               �@���������� pm �̃o�O���C���ielse ���� -1�j�B
% Author:       �Ȃہ[ @Luxion009�i�����j / MATLAB �ڐA�E�C��: You @97kuek_
% Created:      9/24/2016
%-------------------------------------------------------------------------------

    % ===== ���C������ =====
    clc;

    %-------------------- ���[�U�[�ݒ�ǂݍ��� --------------------
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
    % �e�����u���@�f�[�^�ǂݍ���
    disp("Rib data input");
    wingdata      = readmatrix("rib.csv");
    foil_num      = wingdata(:,1);
    foil_name1    = wingdata(:,2);  % �x�[�X���^ index
    foil_name2    = wingdata(:,3);  % �~�b�N�X���� index
    foil_chord    = wingdata(:,4);
    foil_ang      = wingdata(:,5);  % [deg]
    foil_mix      = wingdata(:,6);  % 0..1
    foil_sparpos  = wingdata(:,7);  % ���Ό��ʒu�i0..1�j
    foil_spardia  = wingdata(:,8)./2;  % DXF CIRCLE �̔��a�i���͂����a�z��̂��� /2�j
    foil_plk_ux   = wingdata(:,9);
    foil_plk_lx   = wingdata(:,10);

    % ��ʃX�g�����K�[�ʒu�ǂݍ��݁i��: ���u���Ɓj
    usraw = readmatrix("u_strdata.csv");
    usdata = usraw';
    if ~isempty(usdata)
        usdata(1,:) = [];   % �擪�s���폜�i���R�[�h���P�j
    end

    % ���ʃX�g�����K�[�ʒu
    lsraw = readmatrix("l_strdata.csv");
    lsdata = lsraw';
    if ~isempty(lsdata)
        lsdata(1,:) = [];
    end

    disp("Profile data input");
    % ���^�f�[�^�ǂݍ��݁iWPM_i.dat: 2�� [x y]�j
    foil_p = cell(profilenum,1);
    fs     = zeros(profilenum,1);
    for i = 1:profilenum
        fn = sprintf("WPM_%d.dat", i);
        fp = fopen(fn,'r');
        if fp<0
            error('�t�@�C�����J���܂���: %s', fn);
        end
        fgetl(fp); % �擪�s���̂Ă�i���R�[�h���P�j
        A = fscanf(fp,'%f',[2,Inf])';
        fclose(fp);
        foil_p{i} = foil_rev(A);  % ���^���]���f
    end

    disp("Normalize profile data");
    % ���^�f�[�^��(0,0)���_��
    foil_p_nor = cell(profilenum,1);
    for i = 1:profilenum
        foil_p_nor{i} = foil_norm(foil_p{i});
    end

    disp("Genelating mix profile data");
    % �e���f�[�^�����i�~�b�N�X�j
    mix_foil = cell(numel(foil_num),1);
    for i = 1:numel(foil_num)
        orgnum1   = foil_name1(i,1);
        orgnum2   = foil_name2(i,1);
        orgfoil1  = foil_p_nor{orgnum1};
        orgfoil2  = foil_p_nor{orgnum2};
        mix_foil{i} = foil_gen(orgfoil1, orgfoil2, foil_mix(i,1));
    end

    disp("Foil data output");

    % ===== �o�̓��[�v =====
    tol = 1e-9;  % ����������r�p

    for fnum = 1:numel(foil_num)
        fprintf('%d,', fnum);

        % ���̓p�����[�^
        chord   = foil_chord(fnum,1);
        sparpos = foil_sparpos(fnum,1);
        spardia = foil_spardia(fnum,1);
        alpha   = foil_ang(fnum,1);    % [deg]
        plk_ux  = foil_plk_ux(fnum,1);
        plk_lx  = foil_plk_lx(fnum,1);

        % �X�g�����K�[�ʒu�i�� fnum�j���o�� + ���[�[������
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

        % �o�̓t�@�C��
        fpw = fopen(sprintf("WPMout_%d.dxf", fnum),'wt');
        if fpw<0
            error('�o�̓t�@�C�����J���܂���: WPMout_%d.dxf', fnum);
        end

        ed_foil = mix_foil{fnum};
        len = size(ed_foil,1);
        chordline = [0 0; 1 0];

        % �����_�T��
        front_p = 0;  % �O�� x �ŏ�
        while true
            front_p = front_p + 1;
            if ed_foil(front_p,1) <= min(ed_foil(:,1))
                break;
            end
        end
        [~, u_middle_p] = max(ed_foil(:,2)); % ��ʍő� y
        [~, l_middle_p] = min(ed_foil(:,2)); % ���ʍŏ� y

        u_f_gap = front_p - u_middle_p - 1;
        f_l_gap = l_middle_p - front_p - 1;

        uf_buf = ceil(u_f_gap/10); % �����X�����
        lf_buf = ceil(f_l_gap/10);

        sec_s  = u_middle_p + uf_buf;   % ��2�Z�N�V�����J�n
        four_s = l_middle_p - lf_buf;   % ��4�Z�N�V�����J�n

        % �Z�N�V�����p�]���_�ix �܂��� y ���j
        xx_1 = linspace(ed_foil(sec_s,1), 1, 200);              % ��1: ��� x��y
        xx_2 = linspace(ed_foil(sec_s,2), ed_foil(front_p,2),200); % ��2: ��� �t�֐� y��x
        xx_3 = linspace(ed_foil(front_p,2), ed_foil(four_s,2),200);% ��3: ���� �t�֐� y��x
        xx_4 = linspace(ed_foil(four_s,1), 1, 200);              % ��4: ���� x��y

        % pchip �W���i�P������������ NaN �ɂȂ�\������j
        sf1 = pchip(ed_foil(1:sec_s,1),           ed_foil(1:sec_s,2));      % y(x)
        sf2 = pchip(ed_foil(sec_s:front_p,2),     ed_foil(sec_s:front_p,1));% x(y)
        sf3 = pchip(ed_foil(front_p:four_s,2),    ed_foil(front_p:four_s,1));
        sf4 = pchip(ed_foil(four_s:len,1),        ed_foil(four_s:len,2));

        % �l�̕]��
        sv1 = ppval(sf1, xx_1);
        sv2 = ppval(sf2, xx_2(1,2:end));
        sv3 = ppval(sf3, xx_3(1,2:end));
        sv4 = ppval(sf4, xx_4(1,2:end));

        % �Z�N�V�����̔��]/���`
        xx_1_ed = fliplr(xx_1);   sv1_ed = fliplr(sv1);
        xx_2_ed = xx_2(1,2:end);  % y
        xx_3_ed = xx_3(1,2:end);  % y
        xx_4_ed = xx_4(1,2:end);  % x
        sv2_ed = sv2;  % �t�֐��œ��� x(y) �l�����̂܂܎g�p

        profx = [xx_1_ed, sv2_ed, sv3,   xx_4_ed]';
        profy = [sv1_ed,  xx_2_ed, xx_3_ed, sv4]';

        prof1 = [xx_1_ed', sv1_ed'];    % ��ʁix,y�j��1
        prof2 = [sv2_ed',  xx_2_ed'];   % ��ʁix,y�j��2 �ix��y�j
        prof3 = [sv3',     xx_3_ed'];   % ���ʁix,y�j��3 �ix��y�j
        prof4 = [xx_4_ed', sv4'];       % ���ʁix,y�j��4

        prof  = [profx profy];
        profu = [prof1; prof2]; % ���
        profl = [prof3; prof4]; % ����

        % �X�g�����K�[/�v�����N���z�i��j
        sp_ux = zeros(300,1);
        sp_lx = zeros(300,1);
        for i = 1:300
            sp_ux(i,1) = ((i-1)/299)^2;     % ��ʊ�薧
            sp_lx(i,1) = ((i-1)/299)^1.5;   % ���ʊ�薧
        end

        % �X�g�����K�[���W & �v�����N���W ����
        sp_ux = [sp_ux; str_ux; plk_ux];
        sp_lx = [sp_lx; str_lx; plk_lx];
        sp_ux = sort(sp_ux);
        sp_lx = sort(sp_lx);

        % y ���W���
        sp_u = [sp_ux, zeros(numel(sp_ux),1)];
        sp_l = [sp_lx, zeros(numel(sp_lx),1)];
        for i = 1:size(sp_u,1)
            sp_u(i,2) = interp1(profu(:,1), profu(:,2), sp_u(i,1), 'linear','extrap');
        end
        for i = 1:size(sp_l,1)
            sp_l(i,2) = interp1(profl(:,1), profl(:,2), sp_l(i,1), 'linear','extrap');
        end

        % ���ʒu�ύX�F�L�����o�[���C���������S�_
        cambpoint = (0.001:0.001:0.999)';
        cambline  = zeros(numel(cambpoint),2);
        for i = 1:numel(cambpoint)
            yu = spline(sp_u(:,1), sp_u(:,2), cambpoint(i));
            yl = spline(sp_l(:,1), sp_l(:,2), cambpoint(i));
            cambline(i,1) = cambpoint(i);
            cambline(i,2) = abs(yu - yl);  % ����
        end
        spar_camb   = spline(cambline(:,1), cambline(:,2), sparpos);
        spar_camb_h = spar_camb/2;
        sparmidp    = [sparpos, spline(sp_l(:,1), sp_l(:,2), sparpos) + spar_camb_h];

        % ����֕��i
        sp_ur     = [sp_u(:,1)-sparmidp(1), sp_u(:,2)-sparmidp(2)];
        sp_lr     = [sp_l(:,1)-sparmidp(1), sp_l(:,2)-sparmidp(2)];
        profr     = [prof(:,1)-sparmidp(1),  prof(:,2)-sparmidp(2)];
        profur    = [profu(:,1)-sparmidp(1), profu(:,2)-sparmidp(2)];
        proflr    = [profl(:,1)-sparmidp(1), profl(:,2)-sparmidp(2)];
        chordliner = [chordline(:,1)-sparmidp(1), chordline(:,2)-sparmidp(2)];

        % ��ΐ��@�փX�P�[��
        profr      = profr      .* chord;
        profur     = profur     .* chord;
        proflr     = proflr     .* chord;
        chordliner = chordliner .* chord;
        sp_ur      = sp_ur      .* chord;
        sp_lr      = sp_lr      .* chord;

        % �㉏���i�w�蒷�j
        te_x0 = sp_lr(end,1) - 1;  % �����l
        te_x  = te_x0;
        try
            te_x = fsolve(@(x) endform(x, sp_lr, te_l), te_x0, optimoptions('fsolve','Display','off'));
        catch
            % fsolve �Ȃ�/���s���͒P���Ȍ㉏�����ő��
            % �����ł� x=sp_lr(end,1)-te_l �Ƃ���
            te_x = sp_lr(end,1) - te_l;
        end
        te_y   = interp1(sp_lr(:,1), sp_lr(:,2), te_x, 'linear','extrap');
        te_ang = atan2((sp_lr(end,2)-te_y),(sp_lr(end,1)-te_x)) + 0.5*pi;
        te_p   = [te_x, te_y; te_x + 10*cos(te_ang), te_y + 10*sin(te_ang)];

        % �v�����N/�L���b�v/�X�g�����K�[���W����
        plank_pu = [];
        plank_pl = [];
        cap_pu   = [];
        cap_pl   = [];
        str_pu   = {}; su = 1;
        str_pl   = {}; sl = 1;
        str_peu  = []; % �v�����N�[�p�i��j
        str_pel  = []; % �v�����N�[�p�i���j

        % ��ʃv�����N�E�L���b�v
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

            % �X�g�����K�[
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
        % ��ʃv�����N�̍Ō�[�ɃL���b�v��[��ǉ��i�w�莞�j
        if plkend_sign_u ~= 1 && ~isempty(cap_pu)
            plank_pu(end+1,:) = cap_pu(1,:);
        end

        % ���ʃv�����N�E�L���b�v
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

            % �X�g�����K�[
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
        % ���ʃv�����N�̍Ō�[�ɃL���b�v��[��ǉ��i�w�莞�j
        if plkend_sign_l ~= 1 && ~isempty(cap_pl)
            plank_pl(end+1,:) = cap_pl(1,:);
        end

        % --- DXF �o�� ---
        % �}�ʃI�t�Z�b�g
        cntx = 400; cnty = 200; cnty1 = cnty;

        % �v�����N/�L���b�v���ɂ�錅���S�␳
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

        % �����S�N���X���W
        sparcrossx = [cntx, 0;    cntx, 841];
        sparcrossy = [0,   cnty1; 1189, cnty1];

        % �}�p���i�}�ʊm�F�p�j
        angline = [0,   -400*tan(alpha*pi/180) + cnty1; 
                   1100,  700*tan(alpha*pi/180) + cnty1];

        % DXF �Z�N�V�����J�n
        fprintf(fpw,'  0\n'); fprintf(fpw,'SECTION\n');
        fprintf(fpw,'  2\n'); fprintf(fpw,'ENTITIES\n');

        % ����
        circle_ex(fpw, cntx, cnty1, spardia);

        % �O�`
        polyline_ex(fpw);
        for i = 1:size(profr,1)
            vertex_ex(fpw, profr(i,1)+cntx, profr(i,2)+cnty);
        end
        seqend_ex(fpw);

        % ��ʃv�����N
        polyline_ex(fpw);
        for i = 1:size(plank_pu,1)
            vertex_ex(fpw, plank_pu(i,1)+cntx, plank_pu(i,2)+cnty);
        end
        seqend_ex(fpw);

        % ���ʃv�����N
        polyline_ex(fpw);
        for i = 1:size(plank_pl,1)
            vertex_ex(fpw, plank_pl(i,1)+cntx, plank_pl(i,2)+cnty);
        end
        seqend_ex(fpw);

        % ��ʃL���b�v
        polyline_ex(fpw);
        for i = 1:size(cap_pu,1)
            vertex_ex(fpw, cap_pu(i,1)+cntx, cap_pu(i,2)+cnty);
        end
        seqend_ex(fpw);

        % ���ʃL���b�v
        polyline_ex(fpw);
        for i = 1:size(cap_pl,1)
            vertex_ex(fpw, cap_pl(i,1)+cntx, cap_pl(i,2)+cnty);
        end
        seqend_ex(fpw);

        % ��ʃX�g�����K�[
        for j = 1:numel(str_pu)
            if ~isempty(str_pu{j})
                polyline_ex(fpw);
                for i = 1:size(str_pu{j},1)
                    vertex_ex(fpw, str_pu{j}(i,1)+cntx, str_pu{j}(i,2)+cnty);
                end
                seqend_ex(fpw);
            end
        end
        % ��ʃv�����N�[�i�I�v�V�����j
        if plkend_sign_u == 1 && ~isempty(str_peu)
            polyline_ex(fpw);
            for i = 1:size(str_peu,1)
                vertex_ex(fpw, str_peu(i,1)+cntx, str_peu(i,2)+cnty);
            end
            seqend_ex(fpw);
        end

        % ���ʃX�g�����K�[
        for j = 1:numel(str_pl)
            if ~isempty(str_pl{j})
                polyline_ex(fpw);
                for i = 1:size(str_pl{j},1)
                    vertex_ex(fpw, str_pl{j}(i,1)+cntx, str_pl{j}(i,2)+cnty);
                end
                seqend_ex(fpw);
            end
        end
        % ���ʃv�����N�[�i�I�v�V�����j
        if plkend_sign_l == 1 && ~isempty(str_pel)
            polyline_ex(fpw);
            for i = 1:size(str_pel,1)
                vertex_ex(fpw, str_pel(i,1)+cntx, str_pel(i,2)+cnty);
            end
            seqend_ex(fpw);
        end

        % �㉏��
        polyline_ex(fpw);
        for i = 1:size(te_p,1)
            vertex_ex(fpw, te_p(i,1)+cntx, te_p(i,2)+cnty);
        end
        seqend_ex(fpw);

        % �R�[�h���
        polyline_ex(fpw);
        for i = 1:size(chordliner,1)
            vertex_ex(fpw, chordliner(i,1)+cntx, chordliner(i,2)+cnty);
        end
        seqend_ex(fpw);

        % �����S�N���X�i�c�j
        polyline_ex(fpw);
        for i = 1:size(sparcrossx,1)
            vertex_ex(fpw, sparcrossx(i,1), sparcrossx(i,2));
        end
        seqend_ex(fpw);
        % �����S�N���X�i���j
        polyline_ex(fpw);
        for i = 1:size(sparcrossy,1)
            vertex_ex(fpw, sparcrossy(i,1), sparcrossy(i,2));
        end
        seqend_ex(fpw);

        % �}�p��
        polyline_ex(fpw);
        for i = 1:size(angline,1)
            vertex_ex(fpw, angline(i,1), angline(i,2));
        end
        seqend_ex(fpw);

        % DXF �Z�N�V�����I�[
        fprintf(fpw,'  0\n'); fprintf(fpw,'ENDSEC\n');
        fprintf(fpw,'  0\n'); fprintf(fpw,'EOF\n');
        fclose(fpw);

        % �N���A�iMATLAB �͊֐��X�R�[�v�Ȃ̂ŏȗ������A�ǐ��̂��߁j
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
% �������牺�̓��[�J���֐��Q
% ============================================================================
function orgfoil = foil_rev(orgfoil)
    % ���^���]���f�F�O���ix �ŏ��j�����ɑ̐ς��r���ď㉺�𑵂���
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
    % ���^��O���� (0,0) �֗���悤���s�ړ�
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
    % 2���̗��^�𓯈ꕪ�� x �ɕ�Ԃ��Amixture �Ő��`�����i���R�[�h�݊�: 179�_�j
    orgfoil{1} = orgfoil1;
    orgfoil{2} = orgfoil2;

    % x �����i�㉏���O�����㉏�B�����̏d���������j
    x_foil = zeros(180,1);
    for j = 1:90
        x_foil(90+j,1) = ((j-1)/89)^2; % 0..1
    end
    x_foil(1:90,1) = flipud(x_foil(91:180));
    x_foil(90,:)   = [];  % 179 �_

    for i = 1:2
        % �O���T���ix �ŏ��̈ʒu�j
        [~, fs] = min(orgfoil{i}(:,1));

        y_i = zeros(numel(x_foil),1);
        % ���: 1..90 �� y(x) �ŕ��
        y_i(1:90)   = interp1(orgfoil{i}(1:fs,1),   orgfoil{i}(1:fs,2),   x_foil(1:90),   'linear','extrap');
        % ����: 90..179 �� y(x) �ŕ�ԁi�O���ȍ~�̃f�[�^�j
        y_i(90:end) = interp1(orgfoil{i}(fs:end,1), orgfoil{i}(fs:end,2), x_foil(90:end), 'linear','extrap');

        if i==1
            y1 = y_i; %#ok<NASGU>
        else
            y2 = y_i; %#ok<NASGU>
        end
    end

    % ����
    foildata = [x_foil, mixture.*y1 + (1 - mixture).*y2];
end

function P = plank_ex(fp, rp, cp, gap)
    % �v�����N/�L���b�v�̊�_�i�O�`����@�������� gap �I�t�Z�b�g�j
    % fp: �O�_, rp: ��_, cp: ���ݓ_
    ang = atan2((fp(1,2)-rp(1,2)), (fp(1,1)-rp(1,1))) + 0.5*pi;
    if ang < pi
        pm = 1;
    else
        pm = -1; % ���R�[�h�̃o�O�C��
    end
    P(1,1) = cp(1,1) + pm*gap*cos(ang);
    P(1,2) = cp(1,2) + pm*gap*sin(ang);
end

function S0 = str_ex(fp, rp, cp, plk_dep, str_dep, str_thi)
    % �X�g�����K�[�l�_�i��2�_�E��2�_�j
    ang = atan2((fp(1,2)-rp(1,2)), (fp(1,1)-rp(1,1))) + 0.5*pi;
    if ang < pi
        pm = 1;
    else
        pm = -1;
    end
    gap = plk_dep + str_dep;

    P0(1,1) = cp(1,1) + pm*gap*cos(ang);
    P0(1,2) = cp(1,2) + pm*gap*sin(ang);

    ang_t = ang - 0.5*pi;   % ���ݕ���
    ang_l = ang - pi;       % ��������i��̍��E�j

    % ��_ 2
    S0(2,1) = P0(1,1) + pm*(str_thi/2)*cos(ang - 0.5*pi + 0.5*pi);
    S0(2,2) = P0(1,2) + pm*(str_thi/2)*sin(ang - 0.5*pi + 0.5*pi);
    S0(3,1) = P0(1,1) + pm*(str_thi/2)*cos(ang + pi);
    S0(3,2) = P0(1,2) + pm*(str_thi/2)*sin(ang + pi);

    % ��_ 2
    S0(1,1) = S0(2,1) + pm*(str_dep)*cos(ang_t);
    S0(1,2) = S0(2,2) + pm*(str_dep)*sin(ang_t);
    S0(4,1) = S0(3,1) + pm*(str_dep)*cos(ang_t);
    S0(4,2) = S0(3,2) + pm*(str_dep)*sin(ang_t);
end

function S0 = plankend_ex(fp, rp, cp, plk_dep, str_dep, str_thi)
    % �v�����N�[�X�g�����K�[
    ang = atan2((fp(1,2)-rp(1,2)), (fp(1,1)-rp(1,1))) + 0.5*pi;
    if ang < pi
        pm = 1;
    else
        pm = -1;
    end
    gap = plk_dep + str_dep;

    P0(1,1) = cp(1,1) + pm*gap*cos(ang);
    P0(1,2) = cp(1,2) + pm*gap*sin(ang);

    % �� 2 �_
    S0(2,1) = P0(1,1) + pm*(str_thi/2)*cos(ang - 0.5*pi + 0.5*pi);
    S0(2,2) = P0(1,2) + pm*(str_thi/2)*sin(ang - 0.5*pi + 0.5*pi);
    S0(3,1) = P0(1,1) + pm*(str_thi/2)*cos(ang + pi);
    S0(3,2) = P0(1,2) + pm*(str_thi/2)*sin(ang + pi);

    ang_t = ang - 0.5*pi;

    if S0(3,1) > S0(2,1)
        S0(3,1) = P0(1,1);
        S0(3,2) = P0(1,2);
        % �� 2 �_
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
    % �㉏�m��p: �_ (x, y(x)) ���猻�㉏�_�܂ł̋��� = TE_length
    y = interp1(foil(:,1), foil(:,2), x, 'linear','extrap');
    r = sqrt((x - foil(end,1))^2 + (y - foil(end,2))^2);
    val = r - TE_length;
end

% --- DXF �o�͕⏕ ---
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
