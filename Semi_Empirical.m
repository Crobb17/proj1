function [varagout] = ptable(varagin)

Ptable = readtable("PeriodicTable.csv");

%so it produces 2 graphs
if nargin == 0
atomicNumber = Ptable.atomicNumber;
massData = round(Ptable.atomicWeight);
[massCalc, ~, ~] = massformula(massData, atomicNumber);


figure(1)
plot(atomicNumber, massData*931.5, "rx--")
hold on

plot(atomicNumber, massCalc, "b+-")
title("Atomic Number vs. Mass")
xlabel("Atomic Number")
ylabel("Mass")
legend(["Mass on Periodic Table", "Mass from Semi-Empirical"], "location", "Northwest")
hold off

[~, ~, bepn] = massformula(massData, atomicNumber);
figure(2)
plot(atomicNumber, bepn, "k")
title("Binding Energy per Nucleon vs..Atomic Number")
xlabel("Atomic Number")
ylabel("Binding Energy per Nucleon")

elseif nargin == 1
    atomicNumber_2 = varagin{1};
%end

massCalc_2 = massformula(round(Ptable.AtomicWeight (atomicNumber_2)), atomicNumber_2);
atomInfo = Ptable(atomicNumber_2, :);

varagout{1} = massCalc_2;
varagout{2} = atomInfo;

elseif nargin == 2
    %
    input_1 = varagin{1} ;
    input_2 = varagin{2} ;
    
    if (class(input_1)=='double') & (class(input_1) == class(input_2))
        [isotope_mass, ~, isotope_bepn] = massoformula(input_2, input_1);
        
     if isotope_bepn >= 6.5
         sprintf("The inputted isotope is stable with a mass of %0.1d MeV/c^2", isotope_mass)
     else
         sprintf("The inputted isotope is not stable with a mass of %0.1d MeV/c^2", isotope_mass)
     end
     
     varagout{1} = isotope_mass ;
     
    elseif (class(input_1) == 'double' & (class(input_2) == 'string')
        SAX
        [~, ~, ~] = massformula (Ptable.AtomicWeight (input_1), input_1)
        columnNames = string(Ptable.Properties.VariableNames) ;
        
        if sum (contains(columnNames, input_2)) == 1
            columnTitle = find(contains(columnNames, input_2) == 1;
            specificValue = Ptable(input1.colomnTitle);
        end
    else
        disp("invalid number of inputs")
    end
end

function [mass, be, bepn] = massformula(A,Z)
%calculates the mass of an atmoic nuclei
%output units are MeV/c^2

aP = 11.2;
aA = 23.2;
mP = 938.3;
aV = 15.8;
mN = 939.57;
aS = 18.3;
aC = 0.714;
format bank
bE = Av.*A - as.*A.^(2/3) - ac.*Z.*(Z-1)./.^(1/3)-aA.*(A-2.*Z).^(2)./A;
if mod (A,2) == 0
    if mod (Z,2) == 0
        bE = bE + aP./(A.^(1/2));
    else
        bE = bE - ap./(A.^(1/2))
    end
end

bepn = be./A;
mass = Z*mP+(A-Z)*mN-bE;
end
%I worked on this project with max and Roman

