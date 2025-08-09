%% File: Re_lookup.m
function data = Re_lookup(Re, Cl, ask, data_mat)
% Robust Re–Cl lookup with duplicate handling and safer interpolation
% ask = 1 -> Cd, ask = 2 -> alpha (deg)

Relist = [10000 20000 30000 40000 50000 100000 150000 200000 250000 300000 350000 400000];
Cllist = linspace(0, 1.5, 100);

switch ask
    case 1
        Zcol = 3; % Cd
    case 2
        Zcol = 1; % alpha (deg)
    otherwise
        error('ask must be 1 (Cd) or 2 (alpha).');
end

M = numel(Relist);
Zgrid_all = zeros(M, numel(Cllist));

for i = 1:M
    Cl_i = data_mat(:, 2, i);
    Z_i  = data_mat(:, Zcol, i);

    % 有効データのみ
    mask = isfinite(Cl_i) & isfinite(Z_i);
    Cl_i = Cl_i(mask);
    Z_i  = Z_i(mask);

    if isempty(Cl_i)
        Zi_on_grid = zeros(1, numel(Cllist));
    else
        % Cl を昇順ソート
        [Cl_i, ord] = sort(Cl_i(:));
        Z_i = Z_i(ord);

        % ほぼ同一の Cl を統合（平均）
        [Cl_u, ~, ic] = unique(round(Cl_i, 6), 'stable');
        Z_u = accumarray(ic, Z_i, [], @mean);

        if numel(Cl_u) == 1
            Zi_on_grid = repmat(Z_u, 1, numel(Cllist));
        else
            % 外挿は避け、クエリ範囲をクランプ
            Cq = min(max(Cllist, Cl_u(1)), Cl_u(end));
            % 'pchip' は形状保存で外れ値に強い
            Zi_on_grid = interp1(Cl_u, Z_u, Cq, 'pchip');
        end
    end

    Zgrid_all(i, :) = Zi_on_grid;
end

% 2次元補間で Re × Cl の格子へ
Zgrid = Zgrid_all'; % size = [numel(Cllist) x M]
data = interp2(Relist, Cllist, Zgrid, Re, Cl, 'linear');
end