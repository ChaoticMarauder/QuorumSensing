%parameters
global betacomA betacomB betacomC betacomD betacomE, %basal synthesis rates
global maxcomA maxcomB maxcomC maxcomD maxcomE maxcomX, %maximal synthesis rate
global CSP, %extracellular concentration of CSP
global hmcomA hmcomB hmcomC hmcomD hmcomE hmcomX, %concentration of comE-P for half maximal synthesis
global deltaA deltaB deltaAB deltaC deltaD deltaE deltaX deltaDP deltaEP deltaCSP, %decay rates
global Vcell Vmedia, %volume of single cell, volume of the media
global lambda, %phosphorylation rate
global rho, %dephosphorylation rate
global deltaDCSP, %dissociation rate of comD and CSP complex
global bindDCSP, %complex formation rate of comD and CSP
global exportpreCSP, %export rate of pre-CSP
global hmpCSP, %comC concentration required for half-maximal export rate of pre-CSP
global kABf kABb, %reaction rates for forward and backward reaction of comAB formation

%variables
global XcomA XcomB XcomC XcomD XcomE XcomX XcomAB, %intracellular cconcentrations to keep track of
global XcomEP XcomDP XcomDCSP XcomEEP, %variables that dictate the synthesis of others

%modelling involving comE synthesis
dXcomEdt=betacomE+maxcomE*(XcomEP/(XcomEP+Vcell*hmcomE))-deltaE*XcomE-(lambda/Vcell)*XcomDP*XcomE+rho*XcomEP;
dXcomEPdt=(lambda/Vcell)*XcomDP*XcomE-rho*XcomEP-deltaEP*XcomEP;
%dXcomEEPdt=;

%modelling involving comAB synthesis
dXcomAdt=betacomA+maxcomA*(XcomEP/(XcomEP+Vcell*hmcomA))-kABf*XcomA*XcomB+kABb*XcomAB-deltaA*XcomA;
dXcomBdt=betacomB+maxcomB*(XcomEP/(XcomEP+Vcell*hmcomB))-kABf*XcomA*XcomB+kABb*XcomAB-deltaB*XcomB;
dXcomABdt=kABf*XcomA*XcomB-kABb*XcomAB-deltaAB*XcomAB;

%modelling involving comC synthesis
dXcomCdt=betacomC+maxcomC*(XcomEP/(XcomEP+Vcell*hmcomC))-deltaC*XcomC-exportpreCSP*(XcomC/XcomC+Vcell*hmpCSP);
%dXCSPdt=                             -deltaCSP*XCSP; once released out in
%the open XCSP becomes equal to to CSP because we have only a single cell

%modelling involving comD synthesis
dXcomDdt=betacomD+deltaDCSP*XcomDP-bindDCSP*XcomD*CSP-deltaD*XcomD+maxcomD*(XcomEP/(XcomEP+Vcell*hmcomD));
%dXcomDDdt=;
dXcomDPdt=bindDCSP*XcomD*CSP-deltaDCSP*XcomDP-deltaDP*XcomDP;

%modelling late competence factor
dXcomXdt=maxcomX*(XcomEP/(XcomEP+Vcell*hmcomX))-deltaX*XcomX;%betacomX is equivalent to zero since non-competent bacteria don't have comX mRNA

%modelling for a possible repressor DPR for the negative feedback