function [wc,etm,etu,at,bt] = eqInitialTension(fck,leq)
[ftm,E0,Gf] = eqInitialCompression(fck,leq);

wc=5.14*(Gf/ftm); % Abertura Cr�tica de Fissura (Behnam et al. (2018))
etm=ftm/E0; % Deforma��o referente a tens�o m�xima 
etu=etm+(wc/leq); % Deforma��o limite do amolecimento (Genikomsou e Polak (2015))
at=((2*(ftm/ftm))-1)+(2*sqrt(((ftm/ftm)^2)-(ftm/ftm))); % Eq. 25 de Alfarah et. al (2017)
bt=((0.453*(fck^(2/3)))/Gf)*leq; % Eq. 18 de Alfarah et. al (2017)
end