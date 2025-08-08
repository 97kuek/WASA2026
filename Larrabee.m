% Larrabee法によるプロペラ設計

clear; close all; clc;

%% パラメータ入力
B     = 2;                % ブレード枚数
R     = 1.5;              % プロペラ半径 [m]
rpm   = 140;              % 回転数 [rpm]
Omega = rpm/60*2*pi;      % 角速度 [rad/s]
V     = 7.21;             % 流速 [m/s]
Cl    = 0.60;             % 揚力係数
alpha = deg2rad(1.5);     % 迎え角 [rad]
Cd    = 0.015;            % 抗力係数
rho   = 1.111;            % 空気密度 [kg/m^3]
T     = 25;               % トルク [N·m]

% 非次元化パラメータ
Tc  = 2*T/(rho*V^2*pi*R^2); % 非次元化トルク
lam = V/(Omega*R);          % 非次元化流速

%% Goldstein 関数 G(r)
G = @(r) (2/pi) .* acos( exp( - (B/2*sqrt(lam^2+1)/lam) .* (1 - r./R) ) ) ...
       .* ( (Omega.*r./V).^2 ) ./ ( 1 + (Omega.*r./V).^2 );

%% I1, I2 計算用被積分関数（非次元半径 xi ∈ [0,1]）
i1 = @(xi) G(xi*R) ...
           .* (1 - Cd./(Cl.*(Omega.*(xi*R)./V))) ...
           .* xi;

i2 = @(xi) G(xi*R) ...
           .* (1 - Cd./(Cl.*(Omega.*(xi*R)./V))) ...
           .* xi ...
           ./( 1 + (Omega.*(xi*R)./V).^2 );

%% 積分による I1, I2 の評価
I1 = 4 * integral(i1, 0, 1);
I2 = 2 * integral(i2, 0, 1);

%% 誘起速度パラメータ ζ の計算
zeta       = (I1/(2*I2)) * (1 - sqrt(1 - 4*I2*Tc/(I1^2)));
v_induced  = zeta * V;  % （参考に計算）

%% 翼弦長分布 c(r)
c = @(r) (R/Cl) * zeta * 4*pi/B * lam ...
         .* G(r) ./ sqrt( 1 + (Omega.*r./V).^2 );

%% 迎え角 φ(r) と取り付け角 β(r)
phi  = @(r) atan2( lam*(1 + 0.5*zeta), r./R );
beta = @(r) phi(r) + alpha;

%% 効率 η の計算
epsilon     = atan(Cd/Cl);
eta_element = @(r) tan(phi(r)) ./ tan(phi(r) + epsilon) ...
                   * ( 1/(1 + 0.5*zeta) );

Tc_xi = @(r) 4*zeta .* G(r) .* (1 - Cd./(Cl.*(Omega.*r./V))) .* (r./R) ...
             - 2*zeta^2 .* G(r) .* (1 - Cd./(Cl.*(Omega.*r./V))) .* (r./R) ...
               ./ (1 + (Omega.*r./V).^2);

f_eta = @(r) eta_element(r) .* Tc_xi(r);

eta     = integral(f_eta, 0, R) / (Tc * R);
P_brake = T * V / eta;

%% 結果表示
fprintf('プロペラ効率 η = %.4f\n', eta);
fprintf('Brake power P_brake = %.4f W\n', P_brake);

%% プロットと画像保存
rs = linspace(0, R, 100);

% スクリプトの保存先フォルダを取得
[dirPath, ~, ~] = fileparts( mfilename('fullpath') );

% Goldstein factor
figure;
plot(rs, G(rs), 'LineWidth', 1.5);
xlabel('r [m]', 'FontSize', 15);
ylabel('Goldstein factor', 'FontSize', 15);
grid on;
saveas(gcf, fullfile(dirPath, 'Goldstein.png'));

% 弦長分布
figure;
plot(rs, c(rs), 'LineWidth', 1.5);
xlabel('r [m]', 'FontSize', 15);
ylabel('Chord [m]', 'FontSize', 15);
grid on;
saveas(gcf, fullfile(dirPath, 'chord.png'));

% 迎え角 φ と取り付け角 β
figure;
plot(rs, rad2deg(phi(rs)), 'DisplayName', '\phi', 'LineWidth', 1.5); hold on;
plot(rs, rad2deg(beta(rs)), 'DisplayName', '\beta', 'LineWidth', 1.5);
xlabel('r [m]', 'FontSize', 15);
ylabel('Pitch angle [deg]', 'FontSize', 15);
legend('Location', 'best');
grid on;
saveas(gcf, fullfile(dirPath, 'pitch.png'));
