%% DivergenceAnalysis.m

clear; clc; close all;

%% 1. CSV読み込みとパラメータ準備
% オリジナルの列ヘッダーを保持する
opts = detectImportOptions('wing.csv');
opts.VariableNamingRule = 'preserve';
T = readtable('wing.csv', opts);

% 列データ取得
xold    = T.("span") / 1000;   % span [mm] → [m]
GIp_old = T.GIp;               % GIp
he_old  = T.("T.C.");          % T.C.
c_old   = T.("c") / 1000;      % c [mm] → [m]

% 補間関数
fGIp = @(xq) interp1(xold, GIp_old, xq, 'linear', 'extrap');
fhe  = @(xq) interp1(xold, he_old,  xq, 'linear', 'extrap');
fc   = @(xq) interp1(xold, c_old,   xq, 'linear', 'extrap');

%% 2. 区間分割
N = 150;                
L = max(xold);          
x = linspace(0, L, N);
dx = x(2) - x(1);

GIp = fGIp(x);
he  = fhe(x);
c   = fc(x);

%% 3. 要素行列定義
NN   = dx * [2/6, 1/6; 1/6, 2/6];
NxNx = (1/dx) * [ 1, -1; -1,  1];

%% 4. ダイバージェンス速度範囲設定
Utry = 20;
Udiv = 100;
U    = linspace(0, Utry, Udiv);
rho  = 1.2;
diag_list = zeros(1, Udiv);

%% 5. 行列組立＆固有値計算
for k = 1:Udiv
    u = U(k);
    Kelastic = zeros(N, N);
    Kaero    = zeros(N, N);
    for i = 1:(N-1)
        Kelastic(i:i+1, i:i+1) = Kelastic(i:i+1, i:i+1) ...
            - GIp(i) * NxNx;
        aero_coeff = 0.5 * rho * u^2 * c(i)^2 * (he(i) - 0.25) * 2 * pi;
        Kaero(i:i+1, i:i+1) = Kaero(i:i+1, i:i+1) ...
            + aero_coeff * NN;
    end
    K_full = Kelastic + Kaero;
    K_reduced = K_full(2:end, 2:end);
    ev = eig(K_reduced);
    diag_list(k) = max(real(ev));
end

%% 6. データ保存
save('eigen_divergence.mat', 'U', 'diag_list');

%% 7. プロット＆PDF出力設定
hFig = figure('Units','inches','Position',[1 1 6 4]);
plot(U, diag_list, 'LineWidth', 1.5);
hold on;
yline(0, '--', 'LineWidth', 1);
xline(18.8, '--', 'LineWidth', 1);
hold off;
xlabel('Air speed (m/s)');
ylabel('Maximum real part of eigenvalues');
title('Divergence Eigenvalue vs. Air Speed');
grid on;

% PDF出力時にサイズ崩れしない設定
set(hFig, 'PaperPositionMode', 'auto');
print(hFig, 'Eigenvalues_of_divergence', '-dpdf', '-bestfit');
