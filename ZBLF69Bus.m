clear all;clc;
%% Estimation of DG Location
for h=1:69
d=h;

%% Utility Scenario Data
P_Utility=0;
Q_Utility=0;

%% Load Factor Scenario
LF=1;


%% IEEE 69 Bus Network Data
Sbase=100e+6;  %...W....
Vbase=22800.66e+3;  %....v
Zbase=(Vbase^2)/Sbase;

%..data69bus=[LineNO FromNode ToNode R(ohm)  X(ohm)  PL_ToNode[kw] QL_ToNode[kvar] ]...

 IEEE33BusData=[1     1        2     0.0005   0.0012    0            0 
                2     2        3     0.0005   0.0012    0            0
                3     3        4     0.0015   0.0036    0            0
                4     4        5     0.0251   0.0294    0            0
                5     5        6     0.366    0.1864    2.6         2.2
                6     6        7     0.3811   0.1941    40.4         30 
                7     7        8     0.0922   0.047     75           54 
                8     8        9     0.0493   0.2707    28           19  
                9     9        10    0.819    0.74      60           20  
                10    10       11    0.1872   0.0619    145          104  
                11    11       12    0.7114   0.2351    145          104  
                12    12       13    1.03     0.34     8            5.5  
                13    13       14    1.044    0.345     8            5.5  
                14    14       15    1.058    0.3496    0            0  
                15    15       16    0.1966   0.065     45.5         30  
                16    16       17    0.3744   0.1238     60           35 
                17    17       18    0.0047    0.0016     60           35 
                18    18       19    0.3276    0.1083    0           0  
                19    19       20    0.2106    0.0696    1           0.6  
                20    20       21    0.3416    0.1129    114           81
                21    21       22    0.014     0.0046    5.3           3.5  
                22    22        23    0.1591   0.0526    0           0  
                23    23       24    0.3463    0.1145    28          20 
                24    24       25    0.7488    0.2475    0          0 
                25    25        26    0.3089    0.1021    14           10   
                26    26       27    0.1732   0.0572    14           10   
                27    3       28    0.0044    0.0108    26           18.6  
                28    28       29    0.064   0.1565    26          18.6  
                29    29       30    0.3978   0.1315    0          0  
                30    30       31    0.0702   0.0232     0          0   
                31    31       32    0.351   0.116    0          0 
                32    32       33    0.839    0.2816    14           10 
                33    33       34     1.708    0.5646     19.5          14
                34    34        35    1.474    0.4873     6           4
                35    4       36    0.0034   0.0084     0           0
                36    36       37    0.0851    0.2083     79           56.4
                37    37       38    0.2898    0.7091     384.7       274.5
                38    38       39    0.0822    0.2011     384.7        274.5
                39      8       40     0.0928   0.0473      40.5     28.3
                40      40      41      0.3319  0.1114      3.6     2.7
                41      9       42      0.174   0.0886      4.35    3.5
                42      42      43      0.203   0.1034      26.4    19
                43      43      44      0.2842  0.1447      24      17.2
                44      44      45      0.2813  0.1433      0       0
                45      45      46      1.59    0.5337      0       0
                46      46      47      0.7837  0.263       0       0
                47      47      48      0.3042  0.1006      100     72
                48      48      49      0.3861  0.1172      0       0
                49      49      50      0.5075  0.2585      1244    888
                50      50      51      0.0974  0.0496      32      23
                51      51      52      0.145   0.0738      0       0
                52      52      53      0.7105  0.3619      227     162
                53      53      54      1.041   0.5302      59      42
                54      11      55      0.2012  0.0611      18      13
                55      55      56      0.0047  0.0014      18      13
                56      12      57      0.7394  0.2444      28      20
                57      57      58      0.0047  0.0016      28      20
                58      3       59      0.0044  0.0108      26      18.55   
                59      59      60      0.064   0.1565      26      18.55
                60      60      61      0.1053  0.123       0       0
                61      61      62      0.0304   0.0355      24      17
                62      62      63      0.0018  0.0021      24      17
                63      63      64      0.7283  0.8509      1.2     1
                64      64      65      0.31    0.3623      0       0
                65      65      66      0.041   0.0478      6       4.3
                66      66      67      0.0092  0.0116      0       0
                67      67      68      0.1089  0.1373      39.22    26.3
                68      68      69      0.0009  0.0012      39.22   26.3];
            
%% Calculating PLi & QLi for each Load Bus
PL=zeros(1,69); QL=zeros(1,69);
for ii=1:68
    PL(1,IEEE33BusData(ii,3))=IEEE33BusData(ii,6)*1000/(Sbase);
    QL(1,IEEE33BusData(ii,3))=IEEE33BusData(ii,7)*1000/(Sbase);
end

PL(1,1)=-P_Utility;
QL(1,1)=-Q_Utility;
PL=LF*PL;
QL=LF*QL;

%PL=PL*Sbase;
%QL=QL*Sbase;
%% Calculating Admittance of Branches
Ybus=zeros(69,69);
for ii=1:69
    for k=1:69
        if ii==k
            for j=1:68
                if (IEEE33BusData(j,2)==ii) || (IEEE33BusData(j,3)==ii)
                    Ybus(ii,k)=Ybus(ii,k)+1/(IEEE33BusData(j,4)+1i*IEEE33BusData(j,5));
                end
            end
        else
            for j=1:68
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
 

%% First Estimation of Voltages and Angles
V=zeros(1,69);delta=zeros(1,69);
for ii=1:69
    absV(ii)=1;
    delta(ii)=angle(V(ii));
end
for ii=1:69
    V(ii)=absV(ii)*(cos(delta(ii))+1i*sin(delta(ii)));
end

%% Newton Raphson Method
Error=10;
while Error > 0.01
%% Calculating Differentials
DPDd=zeros(68,68);DPDV=zeros(68,68);DQDd=zeros(68,68);DQDV=zeros(68,68);
iii=0;
for ii=1:69
    if ii~=d
        iii=iii+1;
    for j=2:69
        if ii==j
            for k=1:69
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
Pcalc=zeros(1,68);Qcalc=zeros(1,68);
PLd=zeros(1,68);QLd=zeros(1,68);
iii=0;
for ii=1:69
    if ii~=d
        iii=iii+1;
        PLd(iii)=100*PL(ii);
        QLd(iii)=100*QL(ii);
    for k=1:69
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
for ii=1:68
    delta(ii+1)=delta(ii+1)+X(ii);
    absV(ii+1)=absV(ii+1)+X(2*ii);
end
for ii=1:69
    V(ii)=absV(ii)*(cos(delta(ii))+1i*sin(delta(ii)));
end
Error=sum(abs(V-OldV));
end

%% Calculating Ploss for each Branch
Ploss=zeros(1,68);
for ii=1:68
    %Ploss(1,ii)=((IEEE33BusData(ii,6)*1000/(Sbase))^2+(IEEE33BusData(ii,7)*1000/(Sbase))^2)*(IEEE33BusData(ii,4)/Zbase)/((abs(V(IEEE33BusData(ii,3))))^2);
    %Qloss(1,ii)=((IEEE33BusData(ii,6)*1000/(Sbase))^2+(IEEE33BusData(ii,7)*1000/(Sbase))^2)*(IEEE33BusData(ii,5)/Zbase)/((abs(V(IEEE33BusData(ii,3))))^2);
    Ploss(1,ii)=(IEEE33BusData(ii,4)/Zbase)*(abs((V(IEEE33BusData(ii,2))-V(IEEE33BusData(ii,3)))/(IEEE33BusData(ii,4)/Zbase+1i*IEEE33BusData(ii,5)/Zbase)))^2;
    Qloss(1,ii)=(IEEE33BusData(ii,5)/Zbase)*(abs((V(IEEE33BusData(ii,2))-V(IEEE33BusData(ii,3)))/(IEEE33BusData(ii,4)/Zbase+1i*IEEE33BusData(ii,5)/Zbase)))^2;
end

%% Calculating Total Active Loss
PLoss=sum(Ploss);
Z(h)=PLoss;
PDG=P_Utility-PLoss-sum(PL);
DG(h)=PDG;

end
for ii=1:69
    if Z(ii)==min(Z)
        DG_Location=ii;
        DG_Power=DG(ii)*Sbase;
        Total_Loss=min(Z)*Sbase;
    end
end
DG_Location
DG_Power
Total_Loss

