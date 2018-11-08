load('./parnative.mat')

errorbar(parnaiv(4,:),parnaiv(3,:),parnaiv(1,:)-parnaiv(3,:),parnaiv(2,:)-parnaiv(3,:),'ok','MarkerFaceColor','black')
xlabel('Parental Protease Stability')
ylabel('Naïve Protease Stability')