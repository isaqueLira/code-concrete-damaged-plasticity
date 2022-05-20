function [fcm,fc0,ftm,Eci,E0,Gf,Gch,gc,ecel,ecm] = eqInitialCompression(fck,leq)
fcm=fck+8; % Reistencia � Compress�o M�xima
fc0=0.4*fcm;% Resistencia � Compress�o do Limite Linear
ftm=0.3016*fck^(2/3); % Resistencia � Tra��o M�xima
Eci=10000*((fcm)^(1/3)); % M�dulo de Elasticidade Tangente
E0=round((10000*((fcm)^(1/3)))*(0.8+(0.2*(fcm/88))),2); % M�dulo de Elasticidade N�o Danificado
Gf=(0.073*(fcm^0.18)); % Energia de Fratura
Gch=((fcm/ftm)^2)*Gf; % Energia de Esmagamento
gc=round(Gch/leq,3); % Energia dissipada pelo dano na compress�o por unidade de volume
ecel=fc0/E0; % Deforma��o el�stica n�o danificada m�xima na Compress�o
ecm=(2*fcm)/E0; % Deforma��o total referente a Tens�o de Pico
end
