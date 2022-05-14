%% ABERTURA
clc % limpeza do prompt de comando
clear all % limpeza das vari�veis armazenadas
close all % fecha todas figuras
fprintf('Seja Bem-vindo \n');
fprintf('Matlab \n');
fprintf('Program��o dos Par�metros de Dano do Concreto para o Concrete Damaged Plasticity Model (CDP) \n');
fprintf('Baseado em: Behnam et al. (2018), Alfarah et al. (2017), Genikomsou e Polak (2015), FIB (2010), Birtel e Mark (2006), Kr�tzig e P�lling (2004), Hordijk (1992)\n\n');

%% PLASTICIDADE

fprintf('        -PARAMETROS DE PLASTICIDADE- \n\n');
Phi=input('Entre com o  Angulo de Dilata��o em graus: ');
Varepsilon=input('Entre com o valor da Excentricidade da superficie pl�stica (e): ');
fb0_fc0=input('Entre com a rela��o entre a resistencia a compress�o do concreto biaxial e unixial (fb0/fc0): ');
Kc=input('Entre com a raz�o entre a tens�o desviadora em compress�o e tra��o uniaxial (Kc): ');
Viscosidade=input('Entre com a Viscosidade (u): ');

%% ENTRADA DE DADOS

fprintf('\n\n        -PARAMETROS DE DANO- \n\n');
fck=input('Entre com o valor do fck do concreto em MPa: ');
leq=input('Entre com o comprimento equivalente da malha de MEF em mm: ');
b=input('Entre com a rela��o entre a deforma��o plastica de compress�o e deforma��o de esmagamento: ');

%% EQUA��ES INICIAIS

fcm=fck+8; % Reistencia � Compress�o M�xima
fc0=0.4*fcm;% Resistencia � Compress�o do Limite Linear
ftm=0.3016*fck^(2/3); % Resistencia � Tra��o M�xima
ft0=ftm; % Resistencia � Tra��o do Limite Linear
Eci=10000*((fcm)^(1/3)); % M�dulo de Elasticidade Tangente
E0=round((10000*((fcm)^(1/3)))*(0.8+(0.2*(fcm/88))),2); % M�dulo de Elasticidade N�o Danificado
Poisson=0.2; % Coeficiente de Poisson do Concreto
Gf=(0.073*(fcm^0.18)); % Energia de Fratura
Gch=((fcm/ftm)^2)*Gf; % Energia de Esmagamento
gt=round(Gf/leq,3); % Energia dissipada pelo dano na tra��o por unidade de volume
gc=round(Gch/leq,3); % Energia dissipada pelo dano na compress�o por unidade de volume
wcc c=5.14*(Gf/ftm); % Abertura Cr�tica de Fissura (Behnam et al. (2018))

%% INTERATIVIDADE NA COMPRESS�O %%%

ecel=fc0/E0; % Deforma��o el�stica n�o danificada m�xima na Compress�o
ecm=(2*fcm)/E0; % Deforma��o total referente a Tens�o de Pico
yc=((pi^2)*fcm*ecm)/(2*(((Gch/leq)-(0.5*fcm)*((ecm*(1-b))+(b*(fcm/E0))))^2)); % Eq. 28 de Alfarah et al. (2017) 
gc2=((ecm-ecel)*fcm)/1.2; % Kr�tzig e P�lling (2004) 
gc3=gc-gc2; % Kr�tzig e P�lling (2004)
ecu=((yc*ecm)+2*ecm*(tan(gc3*sqrt(yc/(2*fcm*ecm))))*(sqrt(yc/(2*fcm*ecm))))/yc; % Alfarah et al. (2017)
ec=linspace(ecel,ecu,300); % Deforma��o total de compress�o
sigmac=zeros(length(ec),1)'; % Tens�o de compress�o
ech=zeros(length(ec),1)'; % Deforma��o inelastica de compress�o
ecpl=zeros(length(ec),1)'; % Deforma��o pl�stica de compress�o
dc=zeros(length(ec),1)'; % Dano de compress�o
b1=zeros(length(ec),1)'; % Raz�o deforma��o pl�stica de compress�o e deforma��o inelastica de compress�o
ac=((2*(fcm/fc0))-1)+(2*sqrt(((fcm/fc0)^2)-(fcm/fc0))); % Eq. 24 de Alfarah et al. (2017)
bc=((1.97*(fck+8))/Gch)*leq; % Eq. 18 de Alfarah et. al (2017)

for i=1:length(ec)
    if (ec(i)>=ecel) && (ec(i)<=ecm)
           sigmac(i)=round((((Eci*(ec(i)/fcm)-((ec(i)/ecm)^2)))/(1+(((Eci*(ecm/fcm))-2)*(ec(i)/ecm))))*fcm,3); % Tens�o de compress�o, fase 2
           ech(i)=ec(i)-(sigmac(i)/E0); 
           dc(i)=round(1-((1/(2+ac))*(((2*(1+ac))*exp(-bc*ech(i)))-(ac*exp((-2*bc)*ech(i))))),3); % Eq. 19 de Alfarah et. al (2017)
    elseif ec(i)>=ecm
        sigmac(i)=round((((2+(yc*fcm*ecm))/(2*fcm))-(ec(i)*yc)+((ec(i)^2)*(yc/(2*ecm))))^-1,3); % Tens�o de compress�o, fase 3
        ech(i)=ec(i)-(sigmac(i)/E0);
        dc(i)=round(1-((1/(2+ac))*(((2*(1+ac))*exp(-bc*ech(i)))-(ac*exp((-2*bc)*ech(i))))),3); % Eq. 19 de Alfarah et. al (2017)
    end
    if dc(i)<=0.994 % Valor estimado para garantir deforma��es plasticas sempre positivas
        ecpl(i)=ech(i)-((sigmac(i)/E0)*(dc(i)/(1-dc(i))));
    else
        dc(i)=0.994; % Valor estimado para garantir deforma��es plasticas sempre positivas
        ecpl(i)=ech(i)-((sigmac(i)/E0)*(dc(i)/(1-dc(i))));
    end
    b1(i)=ecpl(i)/ech(i);
end

b1_mean=mean(b1,'omitnan');

while abs(b1_mean-b)>=0.00001
    b=b1_mean;
    yc=((pi^2)*fcm*ecm)/(2*(((Gch/leq)-(0.5*fcm)*((ecm*(1-b))+(b*(fcm/E0))))^2)); % Eq. 28 de Alfarah et al. (2017) 
    ecu=((yc*ecm)+2*ecm*(tan(gc3*sqrt(yc/(2*fcm*ecm))))*(sqrt(yc/(2*fcm*ecm))))/yc; % Alfarah et al. (2017)
    ec=linspace(ecel,ecu,100); % Deforma��o total de compress�o
    sigmac=zeros(length(ec),1)'; % Tens�o de compress�o
    ech=zeros(length(ec),1)'; % Deforma��o inelastica de compress�o
    ecpl=zeros(length(ec),1)'; % Deforma��o pl�stica de compress�o
    dc=zeros(length(ec),1)'; % Dano de compress�o
    b1=zeros(length(ec),1)'; % Raz�o deforma��o pl�stica de compress�o e deforma��o inelastica de compress�o
for i=1:length(ec)
    if (ec(i)>=ecel) && (ec(i)<=ecm)
           sigmac(i)=round((((Eci*(ec(i)/fcm)-((ec(i)/ecm)^2)))/(1+(((Eci*(ecm/fcm))-2)*(ec(i)/ecm))))*fcm,3); % Tens�o de compress�o, fase 2
           ech(i)=ec(i)-(sigmac(i)/E0);
           dc(i)=round(1-((1/(2+ac))*(((2*(1+ac))*exp(-bc*ech(i)))-(ac*exp((-2*bc)*ech(i))))),3); % Eq. 19 de Alfarah et. al (2017)
    elseif ec(i)>=ecm
        sigmac(i)=round((((2+(yc*fcm*ecm))/(2*fcm))-(ec(i)*yc)+((ec(i)^2)*(yc/(2*ecm))))^-1,3); % Tens�o de compress�o, fase 3
        ech(i)=ec(i)-(sigmac(i)/E0);
        dc(i)=round(1-((1/(2+ac))*(((2*(1+ac))*exp(-bc*ech(i)))-(ac*exp((-2*bc)*ech(i))))),3); % Eq. 19 de Alfarah et. al (2017)
    end
    if dc(i)<=0.994 % Valor estimado para garantir deforma��es plasticas sempre positivas
        ecpl(i)=ech(i)-((sigmac(i)/E0)*(dc(i)/(1-dc(i))));
    else
        dc(i)=0.994; % Valor estimado para garantir deforma��es plasticas sempre positivas
        ecpl(i)=ech(i)-((sigmac(i)/E0)*(dc(i)/(1-dc(i))));
    end
    b1(i)=ecpl(i)/ech(i);
end
b1_mean=mean(b1,'omitnan');
end

%% BLOCO 1 - COMPRESS�O

b=b1_mean;
yc=((pi^2)*fcm*ecm)/(2*(((Gch/leq)-(0.5*fcm)*((ecm*(1-b))+(b*(fcm/E0))))^2)); % Eq. 28 de Alfarah et al. (2017) 
ecu=((yc*ecm)+2*ecm*(tan(gc3*sqrt(yc/(2*fcm*ecm))))*(sqrt(yc/(2*fcm*ecm))))/yc; % Alfarah et al. (2017)
ec1=linspace(0,ecu,300); % Deforma��o total de compress�o
sigmac1=zeros(length(ec1),1)'; % Tens�o de compress�o
for i=1:length(ec1)
    if ec1(i)<=ecel
        sigmac1(i)=round(ec1(i)*E0,2);
    elseif ec1(i)>ecel && ec1(i)<=ecm
        sigmac1(i)=round((((Eci*(ec1(i)/fcm)-((ec1(i)/ecm)^2)))/(1+(((Eci*(ecm/fcm))-2)*(ec1(i)/ecm))))*fcm,2);
    else
        sigmac1(i)=round((((2+(yc*fcm*ecm))/(2*fcm))-(ec1(i)*yc)+((ec1(i)^2)*(yc/(2*ecm))))^-1,2);
    end
end

%% BLOCO 2 - COMPRESS�O

ec2=linspace(ecm,ecu,50); % Deforma��o total de compress�o do trecho n�o-linear
sigmac2=zeros(length(ec2),1)'; % Tens�o de compress�o do trecho n�o-linear
ech2=zeros(length(ec2),1)'; % Deforma��o inelastica de compress�o
ecpl2=zeros(length(ec2),1)'; % Deforma��o pl�stica de compress�o
dc2=zeros(length(ec2),1)'; % Dano de compress�o
b2=zeros(length(ec2),1)'; % Raz�o deforma��o pl�stica de compress�o e deforma��o inelastica de compress�o 
for i=1:length(ec2)
    if ec2(i)==ecm
        sigmac2(i)=(((2+(yc*fcm*ecm))/(2*fcm))-(ec2(i)*yc)+((ec2(i)^2)*(yc/(2*ecm))))^-1;
        ech2(i)=0;
        dc2(i)=0;
        ecpl2(i)=0;
    else
        sigmac2(i)=(((2+(yc*fcm*ecm))/(2*fcm))-(ec2(i)*yc)+((ec2(i)^2)*(yc/(2*ecm))))^-1;
        ech2(i)=ec2(i)-(sigmac2(i)/E0);
        dc2(i)=1-((1/(2+ac))*(((2*(1+ac))*exp(-bc*ech2(i)))-(ac*exp((-2*bc)*ech2(i))))); % Eq. 19 de Alfarah et. al (2017)
        if dc2(i)<=0.994
            ecpl2(i)=ech2(i)-((sigmac2(i)/E0)*(dc2(i)/(1-dc2(i))));
            b2(i)=ecpl2(i)/ech2(i);
        else
            dc2(i)=0.994;
            ecpl2(i)=ech2(i)-((sigmac2(i)/E0)*(dc2(i)/(1-dc2(i))));
            b2(i)=ecpl2(i)/ech2(i);
        end
    end
end

%% BLOCO 1 - TRA��O

etm=ftm/E0; % Deforma��o referente a tens�o m�xima 
etu=etm+(wc/leq); % Deforma��o limite do amolecimento (Genikomsou e Polak (2015))
et1=linspace(0,etu,100); % Deforma��o de Tra��o
sigmat1=zeros(length(et1),1)'; % Tens�o de Tra��o
at=((2*(ftm/ftm))-1)+(2*sqrt(((ftm/ftm)^2)-(ftm/ftm))); % Eq. 25 de Alfarah et. al (2017)
bt=((0.453*(fck^(2/3)))/Gf)*leq; % Eq. 18 de Alfarah et. al (2017)
for i=1:length(et1) % Loop para a elabora��o das curvas tens�oxdeforma��o e danoxdeforma��o
    if et1(i)<=etm % Trecho Linear
        sigmat1(i)=et1(i)*E0; % Alfarah et al. (2017)
    else % Trecho N�o-Linear
        sigmat1(i)=(((1+(3*((et1(i)-etm)/(etu-etm)))^3)*(exp(-6.93*((et1(i)-etm)/(etu-etm)))))-(((et1(i)-etm)/(etu-etm))*(1+(3^3))*(exp(-6.93))))*ftm; % Adaptada de Hordijk (1992)
    end
end

%% BLOCO 2 - TRA��O

et2=linspace(etm,etu,50); % Deforma��o de Tra��o
sigmat2=zeros(length(et2),1)'; % Tens�o de Tra��o
etck=zeros(length(et2),1)'; % Deforma��o Inelastica de Tra��o
etpl=zeros(length(et2),1)'; % Deforma��o Pl�stica de Tra��o
dt=zeros(length(et2),1)'; % Dano de Tra��o
b3=zeros(length(et2),1)';
for i=1:length(et2) % Loop para a elabora��o das curvas tens�oxdeforma��o e danoxdeforma��o
    sigmat2(i)=(((1+(3*((et2(i)-etm)/(etu-etm)))^3)*(exp(-6.93*((et2(i)-etm)/(etu-etm)))))-(((et2(i)-etm)/(etu-etm))*(1+(3^3))*(exp(-6.93))))*ftm; % Adaptada de Hordijk (1992)
    etck(i)=et2(i)-(sigmat2(i)/E0); % Adaptada de Birtel e Mark (2006) e Alfarah et al. (2017)
    dt(i)=1-((1/(2+at))*(((2*(1+at))*exp(-bt*etck(i)))-(at*exp((-2*bt)*etck(i)))));
    etpl(i)=etck(i)-((sigmat2(i)/E0)*(dt(i)/(1-dt(i))));
    b3(i)=etpl(i)/etck(i);
end

%% Gr�ficos

figure(1); % Gr�ficos de tens�oxdeforma��o na Compress�o e Tra��o
subplot(2,1,1)
plot(ec1,sigmac1,'linewidth',2)
xlabel('Deforma��o de Compress�o','FontSize',14)
ylabel('Tens�o de Compress�o [MPa]','FontSize',14)
title('Gr�fico Tens�o x Deforma��o de Compress�o do Concreto','FontSize',16)

subplot(2,1,2)
plot(et1,sigmat1,'linewidth',2)
xlabel('Deforma��o de Tra��o','FontSize',14)
ylabel('Tens�o de Tra��o[MPa]','FontSize',14)
title('Gr�fico Tens�o x Deforma��o de Tra��o do Concreto','FontSize',16)

figure(2);
subplot(2,1,1)
plot(ech2,sigmac2,'linewidth',2)
xlabel('Deforma��o Inelastica de Compress�o','FontSize',14)
ylabel('Tens�o de Compress�o [MPa]','FontSize',14)
title('Gr�fico Tens�o x Deforma��o Inelastica de Compress�o do Concreto','FontSize',16)

subplot(2,1,2)
plot(etck,sigmat2,'linewidth',2)
xlabel('Deforma��o Inelastica de Tra��o','FontSize',14)
ylabel('Tens�o de Tra��o[MPa]','FontSize',14)
title('Gr�fico Tens�o x Deforma��o Inelastica de Tra��o do Concreto','FontSize',16)

figure(3);
subplot(2,1,1)
plot(ech2,dc2,'linewidth',2)
xlabel('Deforma��o Inelastica de Compress�o','FontSize',14)
ylabel('Dano em Compress�o','FontSize',14)
title('Gr�fico Dano x Deforma��o Inelastica de Compress�o do Concreto','FontSize',16)
axis([0 ecu 0 1.1])

subplot(2,1,2)
plot(etck,dt,'linewidth',2)
xlabel('Deforma��o Inelastica de Tra��o','FontSize',14)
ylabel('Dano em Tra��o','FontSize',14)
title('Gr�fico Dano x Deforma��o Inelastica de Tra��o do Concreto','FontSize',16)
axis([0 etu 0 1.1])

figure(4);
subplot(2,1,1)
plot(ech2,b2,'linewidth',2)
xlabel('Deforma��o Inelastica de Compress�o','FontSize',14)
ylabel('Dano em Compress�o','FontSize',14)
title('Gr�fico Dano x Deforma��o Inelastica de Compress�o do Concreto','FontSize',16)
axis([0 ecu 0 1])

subplot(2,1,2)
plot(etck,b3,'linewidth',2)
xlabel('Deforma��o Inelastica de Compress�o','FontSize',14)
ylabel('Dano em Compress�o','FontSize',14)
title('Gr�fico Dano x Deforma��o Inelastica de Compress�o do Concreto','FontSize',16)
axis([0 etu 0 1])

%% RESULTADOS

RC=[sigmac2; ech2]; % Resposta Tens�oxDeforma��o Inelastica na Compress�o
RT=[sigmat2; etck]; % Resposta Tens�oxDeforma��o Inelastica na Tra��o
DC=[dc2; ech2]; % Resposta DanoxDeforma��o Inelastica na Compress�o
DT=[dt; etck]; % Resposta DanoxDeforma��o Inelastica na Compress�o

arq=fopen('CDP.txt','wt'); % Abertura do Arquivo de Texto

fprintf(arq,strcat('*Material, name=CONCRETO_C26\n')); % Chama o input Material e o Nome do Material no Abaqus (Deve ser igual ao nome utilizado no pr�prio programa

fprintf(arq,'*Density\n'); % Chama a Densidade
fprintf(arq,'2.4e-09\n'); % Entra com a Densidade    

fprintf(arq,'*Elastic\n'); % Chama os parametros de Elasticidade
fprintf(arq,'%s, ',E0); % Entra com o M�dulo de Elasticidade
fprintf(arq,'0.2\n'); % Entra com o Coeficiente de Poisson

fprintf(arq,'*Concrete Damaged Plasticity\n'); % Chama os parametros de Plasticidade
fprintf(arq,'%s, ',Phi); % Entra com o Angulo de Dilata��o
fprintf(arq,'%s, ',Varepsilon); % Entra com a Excentricidade
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

%% FIM
fprintf('Fim da Programa��o! \n\n');
fprintf('\n\nResultados no Arquivo de Texto Criado na mesma Pasta destino do Arquivo .m. \n\n');
fprintf('Resultados em Gr�ficos nas Janelas que foram Aberta. \n\n');