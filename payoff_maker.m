% this script is supposed to be used with the output file pydel of the
% python scripts which run various delay models across all lysogenic
% fractions and use pydel as the output file to store the payoff values. 
clear;
load fmoi_payoff
pydel = fmoi_payoff;
F1=1:99;
F2=F1;


PFM=zeros(length(F1));  
% phage_PFM=PFM;

for f2=F2
    for f1=F1
        PFM(f2==F2,f1==F1)=pydel(pydel(:,1)==f1&pydel(:,2)==f2,3);
        %         phage_PFM(f1==F1,f2==F2)=(phage_pydel(phage_pydel(:,1)==f1&phage_pydel(:,2)==f2,3));
    end
end
[v,i] =max(min(PFM,[],1))
[X,Y]=meshgrid(F1,F2);
figure;
surf(X,Y,PFM); %this interchange is done to ensure that the surface looks up the right way.
% xlabel('f1');ylabel('f2');
% 
% xlabel('f1');ylabel('f2');
