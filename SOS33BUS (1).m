clear all;clc;

PDG2=500000; %Initializing for DG2 Size (in Watt)
Fit_Error=1000;
DG2_Location=30; % Location of DG2 according to Load Sensitiviy analysis

d=6; % Location of DG1 Based on ZBLF Study

%% Utility Scenario Data
P_Utility=0;
Q_Utility=0;

%% Load Factor Scenario
LF=1;

%% IEEE 33 Bus Network Data
Sbase=100e+6;  %...W....
Vbase=4245.66e+3;  %....v
Zbase=(Vbase^2)/Sbase;

%..data33bus=[LineNO FromNode ToNode R(ohm)  X(ohm)  PL_ToNode[kw] QL_ToNode[kvar] ]...

 IEEE33BusData=[1     1        2     0.0922   0.0477    100          60 
                2     2        3     0.493    0.2511    90           40
                3     3        4     0.366    0.1864    120          80
                4     4        5     0.3811   0.1941    60           30
                5     5        6     0.819    0.707     60           20
                6     6        7     0.1872   0.6188    200          100 
                7     7        8     0.7114   0.2351    200          100 
                8     8        9     1.03     0.74      60           20  
                9     9        10    1.044    0.74      60           20  
                10    10       11    0.1966   0.065     45           30  
                11    11       12    0.3744   0.1238    60           35  
                12    12       13    1.468    1.155     60           35  
                13    13       14    0.5416   0.7129    120          80  
                14    14       15    0.591    0.526     60           10  
                15    15       16    0.7463   0.545     60           20  
                16    16       17    1.289    1.721     60           20 
                17    17       18    0.732    0.574     90           40 
                18    2        19    0.164    0.1565    90           40  
                19    19       20    1.5042   1.3554    90           40  
                20    20       21    0.4095   0.4784    90           40
                21    21       22    0.7089   0.9373    90           40  
                22    3        23    0.4512   0.3083    90           50  
                23    23       24    0.898    0.7091    420          200 
                24    24       25    0.896    0.7011    420          200 
                25    6        26    0.203    0.1034    60           25   
                26    26       27    0.2842   0.1447    60           25   
                27    27       28    1.059    0.9337    60           20  
                28    28       29    0.8042   0.7006    120          70  
                29    29       30    0.5075   0.2585    200          600  
                30    30       31    0.9744   0.963     150          70   
                31    31       32    0.3105   0.3619    210          100 
                32    32       33    0.341    0.5302    60           40 
                ];
            
%% Calculating PLi & QLi for each Load Bus
PL=zeros(1,33); QL=zeros(1,33);
for ii=1:32
    PL(1,IEEE33BusData(ii,3))=IEEE33BusData(ii,6)*1000/(Sbase);
    QL(1,IEEE33BusData(ii,3))=IEEE33BusData(ii,7)*1000/(Sbase);
end

PL(1,1)=-P_Utility;
QL(1,1)=-Q_Utility;
PL=LF*PL;
QL=LF*QL;

a=PL(1,DG2_Location); % Saving PL , QL of DG2 Location for next purposes
b=QL(1,DG2_Location);

Total_Load=sum(PL)*Sbase;

phi=acos(PL(1,DG2_Location)/(PL(1,DG2_Location)^2+QL(1,DG2_Location)^2)); %Calculating phase angle of Load in DG2 Location

%% Calculating Admittance of Branches
Ybus=zeros(33,33);
for ii=1:33
    for k=1:33
        if ii==k
            for j=1:32
                if (IEEE33BusData(j,2)==ii) || (IEEE33BusData(j,3)==ii)
                    Ybus(ii,k)=Ybus(ii,k)+1/(IEEE33BusData(j,4)+1i*IEEE33BusData(j,5));
                end
            end
        else
            for j=1:32
                if (IEEE33BusData(j,2)==ii) && (IEEE33BusData(j,3)==k)
                    Ybus(ii,k)=-1/(IEEE33BusData(j,4)+1i*IEEE33BusData(j,5));
                elseif (IEEE33BusData(j,3)==ii) && (IEEE33BusData(j,2)==k)
                    Ybus(ii,k)=-1/(IEEE33BusData(j,4)+1i*IEEE33BusData(j,5));
                end
            end
        end
    end
end
Ybus=Zbase*Ybus;

Loss(1)=50000;
y=1;


%% Running SOS
while Fit_Error>0.001
 
y=y+1
PL(1,DG2_Location)=a;
QL(1,DG2_Location)=b;
PL(1,DG2_Location)=PL(1,DG2_Location)-PDG2/Sbase;
QL(1,DG2_Location)=QL(1,DG2_Location)-PDG2*tan(phi)/Sbase; % DG power factor and Load Power factor are the same

 

%% First Estimation of Voltages and Angles
V=zeros(1,33);delta=zeros(1,33);
for ii=1:33
    absV(ii)=1;
    delta(ii)=angle(V(ii));
end
for ii=1:33
    V(ii)=absV(ii)*(cos(delta(ii))+1i*sin(delta(ii)));
end

%% Newton Raphson Method
Error=10;

while Error > 0.05
    
%% Calculating Differentials
DPDd=zeros(32,32);DPDV=zeros(32,32);DQDd=zeros(32,32);DQDV=zeros(32,32);
iii=0;
for ii=1:33
    if ii~=d
        iii=iii+1;
    for j=2:33
        if ii==j
            for k=1:33
                if k~=ii
                    DPDd(iii,j-1)=DPDd(iii,j-1)+abs(V(ii))*abs(V(k))*abs(Ybus(ii,k))*sin(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                    DPDV(iii,j-1)=DPDV(iii,j-1)+abs(V(k))*abs(Ybus(ii,k))*cos(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                    DQDd(iii,j-1)=DQDd(iii,j-1)+abs(V(ii))*abs(V(k))*abs(Ybus(ii,k))*cos(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                    DQDV(iii,j-1)=DQDV(iii,j-1)-abs(V(k))*abs(Ybus(ii,k))*sin(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                else
                    DPDV(iii,j-1)=DPDV(iii,j-1)+2*abs(V(ii))*abs(Ybus(ii,ii))*cos(angle(Ybus(ii,ii)));
                    DQDV(iii,j-1)=DQDV(iii,j-1)-2*abs(V(ii))*abs(Ybus(ii,ii))*sin(angle(Ybus(ii,ii)));
                end
            end
        else
            DPDd(iii,j-1)=-abs(V(ii))*abs(V(j))*abs(Ybus(ii,j))*sin(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
            DPDV(iii,j-1)=abs(V(ii))*abs(Ybus(ii,j))*cos(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
            DQDd(iii,j-1)=-abs(V(ii))*abs(V(j))*abs(Ybus(ii,j))*cos(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
            DQDV(iii,j-1)=-abs(V(ii))*abs(Ybus(ii,j))*sin(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
        end
    end
    end
end


%% Calculating P Calculated and Q Calculated
Pcalc=zeros(1,32);Qcalc=zeros(1,32);
PLd=zeros(1,32);QLd=zeros(1,32);
iii=0;
for ii=1:33
    if ii~=d
        iii=iii+1;
        PLd(iii)=100*PL(ii);
        QLd(iii)=100*QL(ii);
    for k=1:33
        Pcalc(iii)=Pcalc(iii)+abs(V(ii))*abs(V(k))*abs(100*Ybus(ii,k))*cos(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
        Qcalc(iii)=Qcalc(iii)-abs(V(ii))*abs(V(k))*abs(100*Ybus(ii,k))*sin(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
    end
    end
end

 
%% Calculating State Variables in Newton Raphson
D=[DPDd DPDV;DQDd DQDV];
X=-inv([DPDd DPDV;DQDd DQDV])*[[Pcalc-PLd]';[Qcalc-QLd]'];
    
%% Calculating New Voltage and Delta
OldV=V;
for ii=1:32
    delta(ii+1)=delta(ii+1)+X(ii);
    absV(ii+1)=absV(ii+1)+X(2*ii);
end
for ii=1:33
    V(ii)=absV(ii)*(cos(delta(ii))+1i*sin(delta(ii)));
end
Error=sum(abs(V-OldV));
end

%% Calculating Ploss for each Branch
Ploss=zeros(1,32);
for ii=1:32
    %Ploss(1,ii)=((IEEE33BusData(ii,6)*1000/(Sbase))^2+(IEEE33BusData(ii,7)*1000/(Sbase))^2)*(IEEE33BusData(ii,4)/Zbase)/((abs(V(IEEE33BusData(ii,3))))^2);
    %Qloss(1,ii)=((IEEE33BusData(ii,6)*1000/(Sbase))^2+(IEEE33BusData(ii,7)*1000/(Sbase))^2)*(IEEE33BusData(ii,5)/Zbase)/((abs(V(IEEE33BusData(ii,3))))^2);
    Ploss(1,ii)=(IEEE33BusData(ii,4)/Zbase)*(abs((V(IEEE33BusData(ii,2))-V(IEEE33BusData(ii,3)))/(IEEE33BusData(ii,4)/Zbase+1i*IEEE33BusData(ii,5)/Zbase)))^2;
    Qloss(1,ii)=(IEEE33BusData(ii,5)/Zbase)*(abs((V(IEEE33BusData(ii,2))-V(IEEE33BusData(ii,3)))/(IEEE33BusData(ii,4)/Zbase+1i*IEEE33BusData(ii,5)/Zbase)))^2;
end

%% Calculating Total Active Loss
PLoss=sum(Ploss);

%% Mutalism Phase
PDG2J=rand*(Total_Load);
Mutual_Vector=(PDG2J+PDG2)/2;
BF1=randperm(2,1);
BF2=randperm(2,1);
PDG2new=PDG2+rand*(PDG2-Mutual_Vector*BF1);
PDG2Jnew=PDG2J+rand*(PDG2-Mutual_Vector*BF2);


% Running ZBLF again with PDG2new
PL(1,DG2_Location)=a;
QL(1,DG2_Location)=b;
PL(1,DG2_Location)=PL(1,DG2_Location)-PDG2new/Sbase;
QL(1,DG2_Location)=QL(1,DG2_Location)-PDG2new*tan(phi)/Sbase; % DG power factor and Load Power factor are the same

 

%% First Estimation of Voltages and Angles
V=zeros(1,33);delta=zeros(1,33);
for ii=1:33
    absV(ii)=1;
    delta(ii)=angle(V(ii));
end
for ii=1:33
    V(ii)=absV(ii)*(cos(delta(ii))+1i*sin(delta(ii)));
end

%% Newton Raphson Method
Error=10;
while Error > 0.05
%% Calculating Differentials
DPDd=zeros(32,32);DPDV=zeros(32,32);DQDd=zeros(32,32);DQDV=zeros(32,32);
iii=0;
for ii=1:33
    if ii~=d
        iii=iii+1;
    for j=2:33
        if ii==j
            for k=1:33
                if k~=ii
                    DPDd(iii,j-1)=DPDd(iii,j-1)+abs(V(ii))*abs(V(k))*abs(Ybus(ii,k))*sin(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                    DPDV(iii,j-1)=DPDV(iii,j-1)+abs(V(k))*abs(Ybus(ii,k))*cos(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                    DQDd(iii,j-1)=DQDd(iii,j-1)+abs(V(ii))*abs(V(k))*abs(Ybus(ii,k))*cos(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                    DQDV(iii,j-1)=DQDV(iii,j-1)-abs(V(k))*abs(Ybus(ii,k))*sin(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                else
                    DPDV(iii,j-1)=DPDV(iii,j-1)+2*abs(V(ii))*abs(Ybus(ii,ii))*cos(angle(Ybus(ii,ii)));
                    DQDV(iii,j-1)=DQDV(iii,j-1)-2*abs(V(ii))*abs(Ybus(ii,ii))*sin(angle(Ybus(ii,ii)));
                end
            end
        else
            DPDd(iii,j-1)=-abs(V(ii))*abs(V(j))*abs(Ybus(ii,j))*sin(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
            DPDV(iii,j-1)=abs(V(ii))*abs(Ybus(ii,j))*cos(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
            DQDd(iii,j-1)=-abs(V(ii))*abs(V(j))*abs(Ybus(ii,j))*cos(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
            DQDV(iii,j-1)=-abs(V(ii))*abs(Ybus(ii,j))*sin(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
        end
    end
    end
end


%% Calculating P Calculated and Q Calculated
Pcalc=zeros(1,32);Qcalc=zeros(1,32);
PLd=zeros(1,32);QLd=zeros(1,32);
iii=0;
for ii=1:33
    if ii~=d
        iii=iii+1;
        PLd(iii)=100*PL(ii);
        QLd(iii)=100*QL(ii);
    for k=1:33
        Pcalc(iii)=Pcalc(iii)+abs(V(ii))*abs(V(k))*abs(100*Ybus(ii,k))*cos(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
        Qcalc(iii)=Qcalc(iii)-abs(V(ii))*abs(V(k))*abs(100*Ybus(ii,k))*sin(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
    end
    end
end

 
%% Calculating State Variables in Newton Raphson
D=[DPDd DPDV;DQDd DQDV];
X=-inv([DPDd DPDV;DQDd DQDV])*[[Pcalc-PLd]';[Qcalc-QLd]'];
    
%% Calculating New Voltage and Delta
OldV=V;
for ii=1:32
    delta(ii+1)=delta(ii+1)+X(ii);
    absV(ii+1)=absV(ii+1)+X(2*ii);
end
for ii=1:33
    V(ii)=absV(ii)*(cos(delta(ii))+1i*sin(delta(ii)));
end
Error=sum(abs(V-OldV));
end

%% Calculating Ploss for each Branch
Ploss=zeros(1,32);
for ii=1:32
    %Ploss(1,ii)=((IEEE33BusData(ii,6)*1000/(Sbase))^2+(IEEE33BusData(ii,7)*1000/(Sbase))^2)*(IEEE33BusData(ii,4)/Zbase)/((abs(V(IEEE33BusData(ii,3))))^2);
    %Qloss(1,ii)=((IEEE33BusData(ii,6)*1000/(Sbase))^2+(IEEE33BusData(ii,7)*1000/(Sbase))^2)*(IEEE33BusData(ii,5)/Zbase)/((abs(V(IEEE33BusData(ii,3))))^2);
    Ploss(1,ii)=(IEEE33BusData(ii,4)/Zbase)*(abs((V(IEEE33BusData(ii,2))-V(IEEE33BusData(ii,3)))/(IEEE33BusData(ii,4)/Zbase+1i*IEEE33BusData(ii,5)/Zbase)))^2;
    Qloss(1,ii)=(IEEE33BusData(ii,5)/Zbase)*(abs((V(IEEE33BusData(ii,2))-V(IEEE33BusData(ii,3)))/(IEEE33BusData(ii,4)/Zbase+1i*IEEE33BusData(ii,5)/Zbase)))^2;
end

%% Calculating Total Active Loss
PLossnew=sum(Ploss);

%% Selecting the best answer between modified Xi and Xi
if PLossnew<PLoss
    PDG2=PDG2new;
    PLoss=PLossnew;
else
    PDG2=PDG2;
    PLoss=PLoss;
end

%% Commensalism Phase

PDG2J=rand*(Total_Load);
PDG2new=PDG2+(-rand+rand)*(PDG2-PDG2J);

% Running ZBLF again with PDG2new
PL(1,DG2_Location)=a;
QL(1,DG2_Location)=b;
PL(1,DG2_Location)=PL(1,DG2_Location)-PDG2new/Sbase;
QL(1,DG2_Location)=QL(1,DG2_Location)-PDG2new*tan(phi)/Sbase; % DG power factor and Load Power factor are the same

 

%% First Estimation of Voltages and Angles
V=zeros(1,33);delta=zeros(1,33);
for ii=1:33
    absV(ii)=1;
    delta(ii)=angle(V(ii));
end
for ii=1:33
    V(ii)=absV(ii)*(cos(delta(ii))+1i*sin(delta(ii)));
end

%% Newton Raphson Method
Error=10;
while Error > 0.05
%% Calculating Differentials
DPDd=zeros(32,32);DPDV=zeros(32,32);DQDd=zeros(32,32);DQDV=zeros(32,32);
iii=0;
for ii=1:33
    if ii~=d
        iii=iii+1;
    for j=2:33
        if ii==j
            for k=1:33
                if k~=ii
                    DPDd(iii,j-1)=DPDd(iii,j-1)+abs(V(ii))*abs(V(k))*abs(Ybus(ii,k))*sin(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                    DPDV(iii,j-1)=DPDV(iii,j-1)+abs(V(k))*abs(Ybus(ii,k))*cos(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                    DQDd(iii,j-1)=DQDd(iii,j-1)+abs(V(ii))*abs(V(k))*abs(Ybus(ii,k))*cos(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                    DQDV(iii,j-1)=DQDV(iii,j-1)-abs(V(k))*abs(Ybus(ii,k))*sin(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                else
                    DPDV(iii,j-1)=DPDV(iii,j-1)+2*abs(V(ii))*abs(Ybus(ii,ii))*cos(angle(Ybus(ii,ii)));
                    DQDV(iii,j-1)=DQDV(iii,j-1)-2*abs(V(ii))*abs(Ybus(ii,ii))*sin(angle(Ybus(ii,ii)));
                end
            end
        else
            DPDd(iii,j-1)=-abs(V(ii))*abs(V(j))*abs(Ybus(ii,j))*sin(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
            DPDV(iii,j-1)=abs(V(ii))*abs(Ybus(ii,j))*cos(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
            DQDd(iii,j-1)=-abs(V(ii))*abs(V(j))*abs(Ybus(ii,j))*cos(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
            DQDV(iii,j-1)=-abs(V(ii))*abs(Ybus(ii,j))*sin(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
        end
    end
    end
end


%% Calculating P Calculated and Q Calculated
Pcalc=zeros(1,32);Qcalc=zeros(1,32);
PLd=zeros(1,32);QLd=zeros(1,32);
iii=0;
for ii=1:33
    if ii~=d
        iii=iii+1;
        PLd(iii)=100*PL(ii);
        QLd(iii)=100*QL(ii);
    for k=1:33
        Pcalc(iii)=Pcalc(iii)+abs(V(ii))*abs(V(k))*abs(100*Ybus(ii,k))*cos(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
        Qcalc(iii)=Qcalc(iii)-abs(V(ii))*abs(V(k))*abs(100*Ybus(ii,k))*sin(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
    end
    end
end

 
%% Calculating State Variables in Newton Raphson
D=[DPDd DPDV;DQDd DQDV];
X=-inv([DPDd DPDV;DQDd DQDV])*[[Pcalc-PLd]';[Qcalc-QLd]'];
    
%% Calculating New Voltage and Delta
OldV=V;
for ii=1:32
    delta(ii+1)=delta(ii+1)+X(ii);
    absV(ii+1)=absV(ii+1)+X(2*ii);
end
for ii=1:33
    V(ii)=absV(ii)*(cos(delta(ii))+1i*sin(delta(ii)));
end
Error=sum(abs(V-OldV));
end

%% Calculating Ploss for each Branch
Ploss=zeros(1,32);
for ii=1:32
    %Ploss(1,ii)=((IEEE33BusData(ii,6)*1000/(Sbase))^2+(IEEE33BusData(ii,7)*1000/(Sbase))^2)*(IEEE33BusData(ii,4)/Zbase)/((abs(V(IEEE33BusData(ii,3))))^2);
    %Qloss(1,ii)=((IEEE33BusData(ii,6)*1000/(Sbase))^2+(IEEE33BusData(ii,7)*1000/(Sbase))^2)*(IEEE33BusData(ii,5)/Zbase)/((abs(V(IEEE33BusData(ii,3))))^2);
    Ploss(1,ii)=(IEEE33BusData(ii,4)/Zbase)*(abs((V(IEEE33BusData(ii,2))-V(IEEE33BusData(ii,3)))/(IEEE33BusData(ii,4)/Zbase+1i*IEEE33BusData(ii,5)/Zbase)))^2;
    Qloss(1,ii)=(IEEE33BusData(ii,5)/Zbase)*(abs((V(IEEE33BusData(ii,2))-V(IEEE33BusData(ii,3)))/(IEEE33BusData(ii,4)/Zbase+1i*IEEE33BusData(ii,5)/Zbase)))^2;
end

%% Calculating Total Active Loss
PLossnew=sum(Ploss);

%% Selecting the best answer between modified Xi and Xi
if PLossnew<PLoss
    PDG2=PDG2new;
    PLoss=PLossnew;
else
    PDG2=PDG2;
    PLoss=PLoss;
end

%% Parasitism Phase
PDG2J=rand*(Total_Load);
PDG2new=PDG2J;

% Running ZBLF again with PDG2new
PL(1,DG2_Location)=a;
QL(1,DG2_Location)=b;
PL(1,DG2_Location)=PL(1,DG2_Location)-PDG2new/Sbase;
QL(1,DG2_Location)=QL(1,DG2_Location)-PDG2new*tan(phi)/Sbase; % DG power factor and Load Power factor are the same

 

%% First Estimation of Voltages and Angles
V=zeros(1,33);delta=zeros(1,33);
for ii=1:33
    absV(ii)=1;
    delta(ii)=angle(V(ii));
end
for ii=1:33
    V(ii)=absV(ii)*(cos(delta(ii))+1i*sin(delta(ii)));
end

%% Newton Raphson Method
Error=10;
while Error > 0.05
%% Calculating Differentials
DPDd=zeros(32,32);DPDV=zeros(32,32);DQDd=zeros(32,32);DQDV=zeros(32,32);
iii=0;
for ii=1:33
    if ii~=d
        iii=iii+1;
    for j=2:33
        if ii==j
            for k=1:33
                if k~=ii
                    DPDd(iii,j-1)=DPDd(iii,j-1)+abs(V(ii))*abs(V(k))*abs(Ybus(ii,k))*sin(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                    DPDV(iii,j-1)=DPDV(iii,j-1)+abs(V(k))*abs(Ybus(ii,k))*cos(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                    DQDd(iii,j-1)=DQDd(iii,j-1)+abs(V(ii))*abs(V(k))*abs(Ybus(ii,k))*cos(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                    DQDV(iii,j-1)=DQDV(iii,j-1)-abs(V(k))*abs(Ybus(ii,k))*sin(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
                else
                    DPDV(iii,j-1)=DPDV(iii,j-1)+2*abs(V(ii))*abs(Ybus(ii,ii))*cos(angle(Ybus(ii,ii)));
                    DQDV(iii,j-1)=DQDV(iii,j-1)-2*abs(V(ii))*abs(Ybus(ii,ii))*sin(angle(Ybus(ii,ii)));
                end
            end
        else
            DPDd(iii,j-1)=-abs(V(ii))*abs(V(j))*abs(Ybus(ii,j))*sin(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
            DPDV(iii,j-1)=abs(V(ii))*abs(Ybus(ii,j))*cos(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
            DQDd(iii,j-1)=-abs(V(ii))*abs(V(j))*abs(Ybus(ii,j))*cos(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
            DQDV(iii,j-1)=-abs(V(ii))*abs(Ybus(ii,j))*sin(angle(Ybus(ii,j))-angle(V(ii))+angle(V(j)));
        end
    end
    end
end


%% Calculating P Calculated and Q Calculated
Pcalc=zeros(1,32);Qcalc=zeros(1,32);
PLd=zeros(1,32);QLd=zeros(1,32);
iii=0;
for ii=1:33
    if ii~=d
        iii=iii+1;
        PLd(iii)=100*PL(ii);
        QLd(iii)=100*QL(ii);
    for k=1:33
        Pcalc(iii)=Pcalc(iii)+abs(V(ii))*abs(V(k))*abs(100*Ybus(ii,k))*cos(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
        Qcalc(iii)=Qcalc(iii)-abs(V(ii))*abs(V(k))*abs(100*Ybus(ii,k))*sin(angle(Ybus(ii,k))-angle(V(ii))+angle(V(k)));
    end
    end
end

 
%% Calculating State Variables in Newton Raphson
D=[DPDd DPDV;DQDd DQDV];
X=-inv([DPDd DPDV;DQDd DQDV])*[[Pcalc-PLd]';[Qcalc-QLd]'];
    
%% Calculating New Voltage and Delta
OldV=V;
for ii=1:32
    delta(ii+1)=delta(ii+1)+X(ii);
    absV(ii+1)=absV(ii+1)+X(2*ii);
end
for ii=1:33
    V(ii)=absV(ii)*(cos(delta(ii))+1i*sin(delta(ii)));
end
Error=sum(abs(V-OldV));
end

%% Calculating Ploss for each Branch
Ploss=zeros(1,32);
for ii=1:32
    %Ploss(1,ii)=((IEEE33BusData(ii,6)*1000/(Sbase))^2+(IEEE33BusData(ii,7)*1000/(Sbase))^2)*(IEEE33BusData(ii,4)/Zbase)/((abs(V(IEEE33BusData(ii,3))))^2);
    %Qloss(1,ii)=((IEEE33BusData(ii,6)*1000/(Sbase))^2+(IEEE33BusData(ii,7)*1000/(Sbase))^2)*(IEEE33BusData(ii,5)/Zbase)/((abs(V(IEEE33BusData(ii,3))))^2);
    Ploss(1,ii)=(IEEE33BusData(ii,4)/Zbase)*(abs((V(IEEE33BusData(ii,2))-V(IEEE33BusData(ii,3)))/(IEEE33BusData(ii,4)/Zbase+1i*IEEE33BusData(ii,5)/Zbase)))^2;
    Qloss(1,ii)=(IEEE33BusData(ii,5)/Zbase)*(abs((V(IEEE33BusData(ii,2))-V(IEEE33BusData(ii,3)))/(IEEE33BusData(ii,4)/Zbase+1i*IEEE33BusData(ii,5)/Zbase)))^2;
end

%% Calculating Total Active Loss
PLossnew=sum(Ploss);

%% Selecting the best answer between modified Xi and Xi
if PLossnew<PLoss
    PDG2=PDG2new;
    PLoss=PLossnew;
else
    PDG2=PDG2;
    PLoss=PLoss;
end
Loss(y)=PLoss*Sbase;
Fit_Error=abs(Loss(y)-Loss(y-1));
end

%%

PDG1=(-P_Utility+PLoss*Sbase+Total_Load-PDG2)
PDG2

Total_Loss=PLoss*Sbase

