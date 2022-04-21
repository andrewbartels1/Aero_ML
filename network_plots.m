% % Manually import the log data from the training sessions % % 


FigHandle_01 = figure('Position', [100, 150, 350, 290]);
plot(1:901, [simpleBiLSTM_RSME, ValidationRSME_DeepLSTM(1:901,1), ValRSME_verysimpleLSTM(1:901)], 'LineWidth', 2.5)
legend('simple BiLSTM', 'Deep LSTM', 'simple LSTM')
xlabel('Iteration')
ylabel('RSME')
title('RMSE LSTM network comparison')
saveas(FigHandle_01, 'network RSME comparison plot', 'jpeg')


FigHandle_02 = figure('Position', [100, 150, 350, 290]);
plot(1:901, [simpleBiLSTM_ValLoss, ValidationLoss_DeepLSTM(1:901,1), ValLoss_verysimpleLSTM(1:901)], 'LineWidth', 2.5)
legend('simple BiLSTM', 'Deep LSTM', 'simple LSTM')
xlabel('Iteration')
ylabel('Validation Loss')
axis([0 910 0 200000])
title('Validation Loss LSTM network')
saveas(FigHandle_02, 'network Validation Loss comparison plot', 'jpeg')




