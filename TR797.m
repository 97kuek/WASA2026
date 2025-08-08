clear all                  % ワークスペースの変数をクリア  
tic;                       % 処理時間計測開始  

% 設計パラメータ  
Lift    = 1000;            % 揚力 [N]  
Uinf    = 7.21;            % 自由流速 [m/s]  
rho     = 1.111;           % 空気密度 [kg/m^3]  
span    = 34;              % 翼スパン全幅 [m]  
le      = span / 2;        % 半スパン長 [m]  
N       = 100;             % 半スパン分割数  
beta    = 0.90;            % 翼根曲げモーメント係数  
delta_S = linspace(le/N/2, le/N/2, N);  % 各パネル半幅 [m]  

% ---- Geometric Condition ----  
% yz平面上: y軸垂直, z軸水平方向  
y   = linspace(delta_S(1), le - delta_S(N), N);  % パネル中心y座標  
z   = zeros(1, N);                               % パネル中心z座標  
fai = zeros(1, N);                               % 局所ダイhedral角 [rad]  

% ---- Geometric variable number ----  
% 相互誘導係数計算用の幾何パラメータを事前計算  
for i = 1:N
    for j = 1:N
        % 回転・平行移動後の相対座標  
        ydash(i,j)   =  (y(i)-y(j))*cos(fai(j)) + (z(i)-z(j))*sin(fai(j));
        zdash(i,j)   = -(y(i)-y(j))*sin(fai(j)) + (z(i)-z(j))*cos(fai(j));
        y2dash(i,j)  =  (y(i)+y(j))*cos(fai(j)) - (z(i)-z(j))*sin(fai(j));
        z2dash(i,j)  =  (y(i)+y(j))*sin(fai(j)) + (z(i)-z(j))*cos(fai(j));
        % 各端点までの距離二乗  
        R2puls(i,j)     = (ydash(i,j)-delta_S(j))^2   + zdash(i,j)^2;
        R2minus(i,j)    = (ydash(i,j)+delta_S(j))^2   + zdash(i,j)^2;
        Rdash2puls(i,j) = (y2dash(i,j)+delta_S(j))^2 + z2dash(i,j)^2;
        Rdash2minus(i,j)= (y2dash(i,j)-delta_S(j))^2 + z2dash(i,j)^2;
    end
end

% ---- Influence coefficient Q の計算 ----  
% パネル j が i に与える誘導効果 Q(i,j)  
for i = 1:N
    for j = 1:N
        % 中間項を分割  
        term1 = ( ydash(i,j)-delta_S(j) )/R2puls(i,j) ...
              - ( ydash(i,j)+delta_S(j) )/R2minus(i,j);
        term2 = ( zdash(i,j)/R2puls(i,j) ) ...
              - ( zdash(i,j)/R2minus(i,j) );
        term3 = ( y2dash(i,j)-delta_S(j) )/Rdash2minus(i,j) ...
              - ( y2dash(i,j)+delta_S(j) )/Rdash2puls(i,j);
        term4 = ( z2dash(i,j)/Rdash2minus(i,j) ) ...
              - ( z2dash(i,j)/Rdash2puls(i,j) );
        % Q の定義式  
        Q(i,j) = -1/(2*pi) * ( ...
                  term1 * cos( fai(i)-fai(j) ) ...
                + term2 * sin( fai(i)-fai(j) ) ...
                + term3 * cos( fai(i)+fai(j) ) ...
                + term4 * sin( fai(i)+fai(j) ) ...
               );
    end
end

% ---- Normalization (無次元化) ----  
delta_sigma     = delta_S ./ le;    % 無次元パネル幅  
eta             = y        ./ le;   % 無次元 y 座標  
etadash         = ydash    ./ le;  
eta2dash        = y2dash   ./ le;  
zeta            = z        ./ le;  
zetadash        = zdash    ./ le;  
zeta2dash       = z2dash   ./ le;  
gamma2puls      = R2puls     /le^2;  
gamma2minus     = R2minus    /le^2;  
gammadash2puls  = Rdash2puls /le^2;  
gammadash2minus = Rdash2minus/le^2;  

% q(i,j): 無次元誘導係数  
for i = 1:N
    for j = 1:N
        t1 = ( etadash(i,j)-delta_sigma(j) )/gamma2puls(i,j) ...
           - ( etadash(i,j)+delta_sigma(j) )/gamma2minus(i,j);
        t2 = ( zetadash(i,j)/gamma2puls(i,j) ) ...
           - ( zetadash(i,j)/gamma2minus(i,j) );
        t3 = ( eta2dash(i,j)-delta_sigma(j) )/gammadash2minus(i,j) ...
           - ( eta2dash(i,j)+delta_sigma(j) )/gammadash2puls(i,j);
        t4 = ( zeta2dash(i,j)/gammadash2minus(i,j) ) ...
           - ( zeta2dash(i,j)/gammadash2puls(i,j) );
        q(i,j) = -1/(2*pi) * ( ...
                  t1 * cos( fai(i)-fai(j) ) ...
                + t2 * sin( fai(i)-fai(j) ) ...
                + t3 * cos( fai(i)+fai(j) ) ...
                + t4 * sin( fai(i)+fai(j) ) ...
               );
    end
end

% ---- Elliptic loading aerodynamic force ----  
BendingMoment_elpl = 2/3/pi * le * Lift;                     % 曲げモーメント (楕円分布)  
InducedDrag_elpl   = Lift^2/(2*pi*rho*Uinf^2*le^2);         % 誘導抵抗 (楕円分布)  
Vn_elpl            = Lift/(2*pi*rho*Uinf * le^2);           % 垂直誘導速度 (楕円分布)  

% ---- Creating the optimization equation ----  
c = 2 * cos(fai) .* delta_sigma;  
b = 3*pi/2 * ( eta.*cos(fai) + zeta.*sin(fai) ) .* delta_sigma;  
for i = 1:N
    for j = 1:N
        A(i,j) = pi * q(i,j) * delta_sigma(i);
    end
end

% ---- solve optimization problem ----  
AAA     = A + A.';                    % A の対称化  
ccc     = -c.';                       % 制約ベクトル c  
bbb     = -b.';                       % 制約ベクトル b  
LeftMat = [ AAA,   ccc, bbb;          % ブロック行列 (左辺)  
            ccc.',  0 ,   0;  
            bbb.',  0 ,   0 ];  
RightMat= [ zeros(N,1); -1; -beta ];  % 右辺ベクトル  
SolveMat= LeftMat \ RightMat;         % 連立方程式を解く  
g       = SolveMat(1:N);              % 局所循環係数 (無次元)  
mu      = SolveMat(N+1:N+2);          % Lagrange 乗数  

% ---- After Solve ----  
efficiency_inverse = g.' * A * g;                   % 逆効率  
efficiency         = 1 / efficiency_inverse;       % スパン効率 η  
Gamma              = g * Lift/(2*le*rho*Uinf);     % 実循環分布  
InducedDrag        = efficiency_inverse * InducedDrag_elpl;  % 総誘導抵抗  

% 局所揚力  
Local_Lift = 4 * rho * Uinf * (Gamma.') .* cos(fai);  
% 局所垂直誘導速度  
Vn = zeros(1, N);  
for i = 1:N
    for j = 1:N
        Vn(i) = Vn(i) + Q(i,j) * Gamma(j);
    end
end  
% 局所誘導抵抗  
Local_InducedDrag = rho * (Gamma.') .* Vn;  

% ---- Aerodynamic force when Elliptical circulation distribution ----  
Lift0_elpl             = 2*Lift/pi/le;  
Lift_elpl              = 4*Lift0_elpl * sqrt(1 - (y/le).^2);  
Gamma0_elpl            = Lift0_elpl/(rho*Uinf*cos(fai(1)));  
Gamma_elpl             = Gamma0_elpl * sqrt(1 - (y/le).^2);  
Local_InducedDrag_elpl = 2 * rho * Gamma_elpl * Vn_elpl;  

% ---- Bending Moment ----  
Local_BendingMoment       = zeros(1, N);  
Local_BendingMoment_elpl  = zeros(1, N);  
for i = 1:N
    tmp1 = 0; tmp2 = 0;
    for j = i:N
        tmp1 = tmp1 + Local_Lift(j)*(y(j)-y(i));
        tmp2 = tmp2 + Lift_elpl(j)*(y(j)-y(i));
    end
    Local_BendingMoment(i)      = tmp1;
    Local_BendingMoment_elpl(i) = tmp2;
end

% プロット・結果表示
figure('Units','pixels','Position',[100,100,1200,800]);
sgtitle('TR-797による解析結果');  % 全体タイトル

% 1. 循環分布
subplot(3,2,1);
plot(y, Gamma, y, Gamma_elpl);
xlabel('位置 [m]');
ylabel('循環');
title('循環分布');
legend(sprintf('\\beta = %.2f',beta), '楕円循環分布');
xlim([0, max(y)*1.05]);
grid on;

% 2. 揚力分布
subplot(3,2,2);
plot(y, Local_Lift, y, Lift_elpl);
xlabel('位置 [m]');
ylabel('揚力 [N]');
title('揚力分布');
legend(sprintf('\\beta = %.2f',beta), '楕円循環分布');
xlim([0, max(y)*1.05]);
grid on;

% 3. 曲げモーメント
subplot(3,2,3);
plot(y, Local_BendingMoment, y, Local_BendingMoment_elpl);
xlabel('位置 [m]');
ylabel(' [Nm]');
title('曲げモーメント');
legend(sprintf('\\beta = %.2f',beta), '楕円循環分布');
xlim([0, max(y)*1.05]);
grid on;

% 4. 誘起速度
subplot(3,2,4);
plot(y, linspace(Vn_elpl,Vn_elpl,N), y, Vn);
xlabel('位置 [m]');
ylabel('誘起速度 [m/s]');
title('誘起速度');
legend(sprintf('\\beta = %.2f',beta), '楕円循環分布');
xlim([0, max(y)*1.05]);
grid on;

% 5. 誘導抵抗
subplot(3,2,5);
plot(y, Local_InducedDrag, y, Local_InducedDrag_elpl);
xlabel('位置 [m]');
ylabel('誘導抵抗 [N]');
title('誘導抵抗');
legend(sprintf('\\beta = %.2f',beta), '楕円循環分布');
xlim([0, max(y)*1.05]);
grid on;

% Leave subplot(3,2,6) empty or use it for additional info:
subplot(3,2,6);
axis off;

% Save the combined figure
print('AllPlots.jpg','-djpeg','-r300');

% 入出力をコマンドウィンドウに表示
disp('---- Input ----');
fprintf('揚力: %.1f [N]\n', Lift);
fprintf('機速: %.2f [m/s]\n', Uinf);
fprintf('スパン長: %.2f [m]\n', span);
fprintf('分割数: %d\n', N);
fprintf('翼根曲げモーメント係数: %.3f\n', beta);

disp('---- Output ----');
fprintf('誘導抵抗: %.4f [N]\n', InducedDrag);
fprintf('翼効率: %.4f\n', efficiency);

disp('---- End ----');
toc;  % 処理時間計測終了  
