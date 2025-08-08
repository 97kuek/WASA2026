%% File: Larrabee_output_fig.m
%--------------------------------
%   Larrabeeグラフ描写と画像書き出し (Larrabee_output_fig.m)
%--------------------------------
figure(1);
subplot(2,1,1);
plot(r, chord * 1000);
xlabel('r [m]');
ylabel('Chord [mm]');
xlim([0 R]);
grid on;

subplot(2,1,2);
plot(r, beta_deg);
xlabel('r [m]');
ylabel('\beta [deg]');
xlim([0 R]);
grid on;

print('-dpng', '-r100', 'result/chord_pitch.png');

figure(2);
plot(r, dTdr, r, dQdr);
ylabel('Thrust [N] , Torque [Nm]');
xlim([0 R]);
legend('Thrust[N]', 'Torque[Nm]');
grid on;
print('-dpng', '-r100', 'result/thrust.png');

figure(3);
plot(r, Re);
ylabel('Reynolds number');
xlim([0 R]);
grid on;
print('-dpng', '-r100', 'result/Re.png');