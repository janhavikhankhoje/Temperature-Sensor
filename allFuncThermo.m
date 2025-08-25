% SA: Ambient Sensor Value (ADC)
% SC: Ambient and Body Sensors Combined Value (ADC)
% TB: Body Temperature (°C)
[SC,SA,minRB,maxRB] = calcBounds();
[RB,SCgrid,SAgrid] = calcRB(SC,SA,minRB,maxRB);
[Temp] = TempPerRow(RB,SA);
[TB] = linearizeTemp(Temp,0.1);
[check] = checkTB(TB,Temp);
[SASCTB] = data(SA,SAgrid,SCgrid,TB);

function [TSCSA] = data(SA,SAgrid,SCgrid,TB)
    TSCSA = cell(1,length(SA));
    for i = 1:length(SA)
        [~,col,~] = find(SAgrid(i,:));
        sub = zeros(length(col),3);
        sub(:,1) = SAgrid(i,col);
        sub(:,2) = SCgrid(i,col);
        sub(:,3) = TB(i,col);
        TSCSA(i) = {sub};
    end
end
function [check] = checkTB(TB,Temp)
    [row,~] = size(TB);
    check = ones(1,row);
    for i = 1:row
        [~,col,ogtemp] = find(Temp(i,:));
        newtemp = TB(i,col);
        pos = find(newtemp > ogtemp);
        check(i) = length(pos);
    end
end
function [TB] = linearizeTemp(Temp,increment)
    TB = Temp;
    [row,~] = size(Temp);
    for i = 1:row
        [~,col,Temploop] = find(Temp(i,:));
        sumdiff = 0;
        [Temploop] = flip(Temploop,2); %order least to greatest
        lastswitch = 1;
        for j = 2:length(Temploop)
            % track sum of differences
            diff = Temploop(j) - Temploop(j-1);
            sumdiff = sumdiff+diff;
            % Temp is max for each SV. (since R is min).
            if (sumdiff >= increment) % is increment closer to j or j-1?
                Temploop(lastswitch:j) = shiftT(increment,j,Temploop,lastswitch);
                lastswitch = j;
                sumdiff = 0;
            else
            end
        end
        TB(i,col) = flip(Temploop,2);
    end
end

function [linearT] = shiftT(increment,j,Temploop,lastswitch)
    linearT = Temploop(lastswitch:j);
    [~,TatJ] = newTatJ(increment,j,Temploop,lastswitch);
    linearT(1) = TatJ-increment;
    slope = increment/(length(linearT)-1);
    x = 1:(length(linearT)-1);
    linearT(x+1) = linearT(1) + (x.*slope);
    
function [stepnum,TatJ] = newTatJ(increment,j,Temploop,lastswitch)
err1 = abs((( (Temploop(j-1)-Temploop(lastswitch)) -increment)/ increment )*100);
err2 = ( ((Temploop(j)-Temploop(lastswitch)) -increment) /increment )*100;
if (err2 > err1) % closer to j-1
    stepnum = round( Temploop(j-1)/increment );
    TatJ = (stepnum*increment);
elseif (err1 > err2) % closer to j
    stepnum = round( Temploop(j)/increment );
    TatJ = (stepnum*increment);
end
end

end
function [Temp] = TempPerRow(RB,SA)
% RB set --> abc value per row --> convert and replace RB with TB
Temp = RB;
for i = 1:length(SA) %loop through rows
    [~,col,RBloop] = find(RB(i,:));
    % min temp = last (biggest) RB. mid temp = median RB.
    [abc] = heartmat(308.15,311.65,315.15,RBloop(length(RBloop)),median(RBloop),RBloop(1));
    [Temploop] = Res2temp(abc(1),abc(2),abc(3),RBloop); % returns in Celsius
    Temp(i,col) = Temploop;
end
end

function [Temp] = Res2temp(a,b,c,R)
Temp = round( (((a + (b.*log(R)) + (c.*(log(R)).^3) ).^(-1)) - 273.15) ,5,'Significant');
end

function [RB,SC,SA] = calcRB(SC,SA,minRB,maxRB)
bit = 1023;
[SC,SA] = meshgrid(SC,SA);
VR = 3.3;
[row,col] = size(SA);
RB = zeros(row,col);
for i = 1:row
    SAl = SA(i,1);
    SCl = SC(1,:);
    RB(i,:) = round( ((SCl.*VR./bit).*((10000.*SAl.*VR./bit)/(5-(VR.*SAl./bit)))./(5-(SCl.*VR./bit))) ,8,'Significant');
end
out = find( (RB > maxRB)|(RB < minRB) ); 
SC(out) = 0;
SA(out) = 0;
RB(out) = 0;
end

function [abc] = heartmat(T1,T2,T3,R1,R2,R3)
eqn = [1 log(R1) (log(R1)).^3;
1 log(R2) (log(R2)).^3;
1 log(R3) (log(R3)).^3];
temps = [1/T1; 1/T2; 1/T3];
abc = linsolve(eqn,temps);
end
function [SC,SA,minRB,maxRB] = calcBounds()
% RA/RB = maxTemp:minTemp = minRofmaxT:maxRofminT
% RA = minR45T:maxR20T
% RB = minR42T:maxR35T
VR = 3.3;
bit = 1023;
minRA = 4090; %45C
[minSA] = SAfromR(minRA,VR,bit);
maxRA = 13188; % 20C
[maxSA] = SAfromR(maxRA,VR,bit);
maxRB = 6985; % 35C max
[minRB] = RTxformula(4,40,42,4998); % 40C min (42)
[maxSC] = SCfromRB(minSA,maxRB,VR,bit);
[minSC] = SCfromRB(maxSA,minRB,VR,bit);
SA = minSA:1:maxSA;
SC = minSC:1:maxSC;

function [RT] = RTxformula(ax,Tx,T,RTx)
RT = RTx .* (exp( (ax./100) .* ((Tx+273.15).^2) .*   ( (1./(T+273.15)) - (1./(Tx+273.15)) ) ));
end
function [SA] = SAfromR(RA,VR,bit)
SA = floor( ((5*RA)/(10000+RA))*(bit/VR) );
end
function [SC] = SCfromRB(SA,RB,VR,bit)
% RA = ( (10000*SA*VR/bit)/(5 - (VR*SA/bit)) );
SC = floor((5*RB)/( ((10000*SA*VR/bit)/(5 - (VR*SA/bit))) +RB)*(bit/VR));
end
end