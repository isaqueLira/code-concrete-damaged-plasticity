function [RC,RT,DC,DT,arq] = results(fck,leq,b)
[Phi,e,fb0_fc0,Kc,Viscosidade] = plasticity();
[sigmac2,ech2,dc2] = secondCompressionBlock(fck,leq,b);
[sigmat2,etck,dt] = secondTensionBlock(fck,leq);

arq=fopen('concrete.txt','wt'); % Abertura do Arquivo de Texto
RC=[sigmac2; ech2]; % Resposta Tens�oxDeforma��o Inelastica na Compress�o
RT=[sigmat2; etck]; % Resposta Tens�oxDeforma��o Inelastica na Tra��o
DC=[dc2; ech2]; % Resposta DanoxDeforma��o Inelastica na Compress�o
DT=[dt; etck]; % Resposta DanoxDeforma��o Inelastica na Compress�o
E0=round((10000*((fck+8)^(1/3)))*(0.8+(0.2*((fck+8)/88))),2); % M�dulo de Elasticidade N�o Danificado

fprintf(arq,'*Material, name=CONCRETO_C26\n'); % Chama o input Material e o Nome do Material no Abaqus (Deve ser igual ao nome utilizado no pr�prio programa

fprintf(arq,'*Density\n'); % Chama a Densidade
fprintf(arq,'2.4e-09\n'); % Entra com a Densidade    

fprintf(arq,'*Elastic\n'); % Chama os parametros de Elasticidade
fprintf(arq,'%s, ',E0); % Entra com o M�dulo de Elasticidade
fprintf(arq,'0.2\n'); % Entra com o Coeficiente de Poisson

fprintf(arq,'*Concrete Damaged Plasticity\n'); % Chama os parametros de Plasticidade
fprintf(arq,'%s, ',Phi); % Entra com o Angulo de Dilata��o
fprintf(arq,'%s, ',e); % Entra com a Excentricidade
fprintf(arq,'%s, ',fb0_fc0); % Entra com o fb0/fc0
fprintf(arq,'%s, ',Kc); % Entra com Kc 
fprintf(arq,'%s\n',Viscosidade); % Entra com a viscosidade

fprintf(arq,'*Concrete Compression Hardening\n'); % Chama os parametros do Comportamento do Concreto na Compress�o 
fprintf(arq,'%4f, %4f\n',RC); % Entra com os dados de Tens�o e Deforma��o Inelastica na Compress�o

fprintf(arq,'*Concrete Tension Stiffening\n'); % Chama os parametros do Comportamento do Concreto na Tra��o
fprintf(arq,'%4f, %4f\n',RT); % Entra com os dados de Tens�o e Deforma��o Inelastica na Tra��o

fprintf(arq,'*Concrete Compression Damage\n'); % Chama os parametros do Comportamento do Dano do Concreto na Compress�o
fprintf(arq,'%4f, %4f\n',DC); % Entra com os dados de Dano e Deforma��o Inelastica na Compress�o

fprintf(arq,'*Concrete Tension Damage\n'); % Chama os parametros do Comportamento do Dano do Concreto na Tra��o
fprintf(arq,'%4f, %4f\n',DT); % Entra com os dados de Dano e Deforma��o Inelastica na Tra��o
end