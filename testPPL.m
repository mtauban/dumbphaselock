function  phi2  = testPPL( chi,samples,framesize, fluctuation)

% Le message est transmis en étant codé sur des symboles à 6 bits
symbols=[ [ 0 0 1 1 0 1 ] ;
          [ 0 0 1 1 1 0 ] ; 
          [ 0 1 0 0 1 1 ] ] ;
% Le message se termine par un bouchon, assemblage unique de 2 symboles sur
% 6 bits qui n'est pas réalisable par assemblage des symboles de données.
endsymbol = [ 1 0 1 1 0 0 1 1 1 0 0 0 ] ;

% Construction de TXdata, les données à transmettre
TXdata=[] ;      
for i=1:framesize ; 
   TXdata = cat(1,TXdata,symbols(floor(rand()*3)+1,:)') ; 
end
TXdata = cat(1,TXdata,endsymbol') ; 

% On ajuste le paramètre de taille pour correspondre au codage 6bits +
% bouchon à 12bits/
framesize=framesize*6+12 ; 

% On décide de noyer le signal dans du bruit
iniTX = cat(1,TXdata ) ;
TXdata = cat(1,round(rand(framesize, 1)), TXdata, rand(framesize, 1)) ; 
initime=-framesize*rand() ; 
TXtime=initime ;
timeend=size(TXdata,1) ;
RXdata=zeros(framesize,1) ; 

% Parametres du Phase Lock
 % On code notre position par rapport au déroulement du cycle par un entier
 % compris entre 0 et SIZE=INCREMENT*SAMPLES. Chaque échantillon nous fait
 % avancer la position de INCREMENT. Quand on remarque un retard ou une
 % avance, on corrige cette position en ajoutant ou retirant la moitié de
 % ce qui nous manque. Comme on a un signal périodique, quand on arrive à
 % SIZE/2, alors que nous étions en retard, on devient en avance par
 % rapport à l'autre période. Voilà pourquoi on travaille entre 0:SIZE/2 et
 % SIZE/2:SIZE.
PL_increment = 20 ;
PL_size = PL_increment * samples ; 
PL_trans= PL_size/2 ; 
PL_ramp=0 ; 
PL_value = 0 ; 
PL_count = 0 ; 
RXconf=0 ; 

% Boucle globale
while (TXtime<timeend) ; 
    
    % Dans un premier temps nous devons déterminer ce que l'éméteur à
    % transmis
    if (TXtime<0) TXtime=0 ; 
    end
    valueindex=floor(TXtime) ; 
    value=TXdata(valueindex+1) ;
    
    % Voilà la valeur émise, maintenant on peut travailler en simulant le
    % récepteur.
    RXget=value ; 
    
    
    % A chaque échantillon, on avance notre position 
    PL_ramp=PL_ramp+PL_increment ;
    % Compteur pour faire la moyenne sur les samples derniers bits.
    PL_count = PL_count + RXget  ; 
  
    % Si on remarque une transition, on vient corriger notre retard ou
    % notre avance en se dirigeant dans la bonne direction mais en ne
    % faisant que la moitié du chemin. J'ai librement adapté ce que j'ai pu
    % trouver sur le net où généralement c'est un incrément fixe qui est
    % choisi. 
    if (PL_value ~= RXget) 
     %   disp('tac') 
        PL_value= RXget ; 
        if (PL_ramp>PL_trans) 
            PL_ramp = PL_ramp + (PL_size-PL_ramp)/2 ; 
        else
            PL_ramp = PL_ramp/2 ;
            
        end 
    end
    
    % Nous avons fini un cycle selon notre indicateur, nous réalisons donc
    % la moyenne sur les échantillons du bit et nous enregistrons le bit.
    if (PL_ramp >= PL_size) 
       % disp ('toc') ;
        p=PL_count/samples ;
        PL_count=0 ; 
        if (p>chi) RXconf=1 ; 
        else RXconf=0 ; 
        end
        PL_ramp = PL_ramp - PL_size ; 
        
        % Ici j'enregistre toujours en dernière position et je shift mes
        % bits vers la gauche. J'utilise de fonctions complexes de Matlab
        % mais on pourra le réaliser facilement avec les opérateurs
        % binaires en C++ quand on travaillera sur des types adaptés.
        RXdata = circshift(RXdata,-1) ;
        RXdata(end) = RXconf ; 
    end
    
    % Ici je regarde si je trouve le symbole bouchon. Si je le trouve,
    % j'arrête tout. 
    if (RXdata(end-11:end)==endsymbol(:))
   %     disp('Found Pattern') ; %
        break
    end
    
    % Enfin on avance le temps en considérant une fluctuation sur la
    % fréquence d'horloge. 
    TXtime = TXtime + ((1/samples)*(1+normrnd(0,fluctuation)))  ;
end

% Pour quantifier le succès, nous faisons la retirons de 1 la sommes des différences au carré
% que nous normalisons par le nombre de bit. (Phi^2) varie entre 1 (accord
% parfait) et 2 (désaccord parfait). 
phi2=1-(norm(iniTX-RXdata)^2)/(framesize) ;
end

