%% File: Re_lookup.m
function data = Re_lookup(Re, Cl, ask, data_mat)
%-------------------------------------------------
% Cd = f(Re,Cl) または alpha = f(Re,Cl) を求めるための関数
% ask = 1 の時は Cd, ask = 2 の時は alpha
%-------------------------------------------------
Relist = [10000 20000 30000 40000 50000 100000 150000 200000 250000 300000 350000 400000];
Cllist = linspace(0, 1.5, 100);

if ask == 1
    for i = 1:12
        Cd_mat_new(i,:) = interp1(data_mat(:,2,i), data_mat(:,3,i), Cllist, 'cubic', 0);
    end
    data = interp2(Relist, Cllist, Cd_mat_new', Re, Cl);
elseif ask == 2
    for i = 1:12
        alpha_mat_new(i,:) = interp1(data_mat(:,2,i), data_mat(:,1,i), Cllist, 'cubic', 0);
    end
    data = interp2(Relist, Cllist, alpha_mat_new', Re, Cl);
else
    error('ask must be 1 (Cd) or 2 (alpha).');
end
end