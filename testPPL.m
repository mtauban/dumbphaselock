function  phi2  = testPPL( chi,samples,framesize, fluctuation)

% Le message est transmis en �tant cod� sur des symboles � 6 bits
symbols=[ [ 0 0 1 1 0 1 ] ;
          [ 0 0 1 1 1 0 ] ; 
          [ 0 1 0 0 1 1 ] ] ;
% Le message se termine par un bouchon, assemblage unique de 2 symboles sur
% 6 bits qui n'est pas r�alisable par assemblage des symboles de donn�es.
endsymbol = [ 1 0 1 1 0 0 1 1 1 0 0 0 ] ;

% Construction de TXdata, les donn�es � transmettre
TXdata=[] ;      
for i=1:framesize ; 
   TXdata = cat(1,TXdata,symbols(floor(rand()*3)+1,:)') ; 
end
TXdata = cat(1,TXdata,endsymbol') ; 

% On ajuste le param�tre de taille pour correspondre au codage 6bits +
% bouchon � 12bits/
framesize=framesize*6+12 ; 

% On d�cide de noyer le signal dans du bruit
iniTX = cat(1,TXdata ) ;
TXdata = cat(1,round(rand(framesize, 1)), TXdata, rand(framesize, 1)) ; 
initime=-framesize*rand() ; 
TXtime=initime ;
timeend=size(TXdata,1) ;
RXdata=zeros(framesize,1) ; 

% Parametres du Phase Lock
 % On code notre position par rapport au d�roulement du cycle par un entier
 % compris entre 0 et SIZE=INCREMENT*SAMPLES. Chaque �chantillon nous fait
 % avancer la position de INCREMENT. Quand on remarque un retard ou une
 % avance, on corrige cette position en ajoutant ou retirant la moiti� de
 % ce qui nous manque. Comme on a un signal p�riodique, quand on arrive �
 % SIZE/2, alors que nous �tions en retard, on devient en avance par
 % rapport � l'autre p�riode. Voil� pourquoi on travaille entre 0:SIZE/2 et
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
    
    % Dans un premier temps nous devons d�terminer ce que l'�m�teur �
    % transmis
    if (TXtime<0) TXtime=0 ; 
    end
    valueindex=floor(TXtime) ; 
    value=TXdata(valueindex+1) ;
    
    % Voil� la valeur �mise, maintenant on peut travailler en simulant le
    % r�cepteur.
    RXget=value ; 
    
    
    % A chaque �chantillon, on avance notre position 
    PL_ramp=PL_ramp+PL_increment ;
    % Compteur pour faire la moyenne sur les samples derniers bits.
    PL_count = PL_count + RXget  ; 
  
    % Si on remarque une transition, on vient corriger notre retard ou
    % notre avance en se dirigeant dans la bonne direction mais en ne
    % faisant que la moiti� du chemin. J'ai librement adapt� ce que j'ai pu
    % trouver sur le net o� g�n�ralement c'est un incr�ment fixe qui est
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
    
    % Nous avons fini un cycle selon notre indicateur, nous r�alisons donc
    % la moyenne sur les �chantillons du bit et nous enregistrons le bit.
    if (PL_ramp >= PL_size) 
       % disp ('toc') ;
        p=PL_count/samples ;
        PL_count=0 ; 
        if (p>chi) RXconf=1 ; 
        else RXconf=0 ; 
        end
        PL_ramp = PL_ramp - PL_size ; 
        
        % Ici j'enregistre toujours en derni�re position et je shift mes
        % bits vers la gauche. J'utilise de fonctions complexes de Matlab
        % mais on pourra le r�aliser facilement avec les op�rateurs
        % binaires en C++ quand on travaillera sur des types adapt�s.
        RXdata = circshift(RXdata,-1) ;
        RXdata(end) = RXconf ; 
    end
    
    % Ici je regarde si je trouve le symbole bouchon. Si je le trouve,
    % j'arr�te tout. 
    if (RXdata(end-11:end)==endsymbol(:))
   %     disp('Found Pattern') ; %
        break
    end
    
    % Enfin on avance le temps en consid�rant une fluctuation sur la
    % fr�quence d'horloge. 
    TXtime = TXtime + ((1/samples)*(1+normrnd(0,fluctuation)))  ;
end

% Pour quantifier le succ�s, nous faisons la retirons de 1 la sommes des diff�rences au carr�
% que nous normalisons par le nombre de bit. (Phi^2) varie entre 1 (accord
% parfait) et 2 (d�saccord parfait). 
phi2=1-(norm(iniTX-RXdata)^2)/(framesize) ;
end

