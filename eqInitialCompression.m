function [fcm,fc0,ftm,Eci,E0,Gf,Gch,gc,ecel,ecm] = eqInitialCompression(fck,leq)
fcm=fck+8; % Reistencia à Compressão Máxima
fc0=0.4*fcm;% Resistencia à Compressão do Limite Linear
ftm=0.3016*fck^(2/3); % Resistencia à Tração Máxima
Eci=10000*((fcm)^(1/3)); % Módulo de Elasticidade Tangente
E0=round((10000*((fcm)^(1/3)))*(0.8+(0.2*(fcm/88))),2); % Módulo de Elasticidade Não Danificado
Gf=(0.073*(fcm^0.18)); % Energia de Fratura
Gch=((fcm/ftm)^2)*Gf; % Energia de Esmagamento
gc=round(Gch/leq,3); % Energia dissipada pelo dano na compressão por unidade de volume
ecel=fc0/E0; % Deformação elástica não danificada máxima na Compressão
ecm=(2*fcm)/E0; % Deformação total referente a Tensão de Pico
end
