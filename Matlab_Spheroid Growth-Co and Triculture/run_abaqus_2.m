clear all
close all

% set tau and P0
mat_props(1) = -0.05; % TNF-A 0.05
mat_props(2) = 0.02; % TGB-B -.006 .03
mat_props(3) = 0.15; % IL-6 0.05  .175

%%%
% A = RAW4T1
% B = MC34T1
% C = MC3RAW4T1

% cycle through hydrogels
for i=1:3
    % cycle through groups
    for j=1:3
        
    % low stiffness (0.58 kPa)
    if i == 1 && j == 1
        mat_props(4) = 0.00011; % C10
        mat_props(5) = 4138; % D1
        
        mat_props(6) = 1; % PCRA  1.6;
        mat_props(7) = -0.27; % PCRB
        mat_props(8) = -0.71; % PCR6
        
        f = 'abaqus j=sph_axiM_RAW4T1_58 inp=sph_axi_RAW4T1.inp user=umat_NH_g21_PCR int';
        p = 'abaqus cae noGUI=Extract_coordA58.py';
        
    elseif i == 1 && j == 2
        mat_props(4) = 0.00011; % C10
        mat_props(5) = 4138; % D1
        
        mat_props(6) = -0.85; % PCRA
        mat_props(7) = 0.12; % PCRB
        mat_props(8) = -0.65; % PCR6
        
        f = 'abaqus j=sph_axiM_MC34T1_58 inp=sph_axi_MC34T1.inp user=umat_NH_g21_PCR int';
        p = 'abaqus cae noGUI=Extract_coordB58.py';
        
    elseif i == 1 && j == 3
        mat_props(4) = 0.00011; % C10
        mat_props(5) = 4138; % D1
        
        mat_props(6) = 1; % PCRA  3.84
        mat_props(7) = -0.61; % PCRB
        mat_props(8) = -0.46; % PCR6
        
        f = 'abaqus j=sph_axiM_MC3RAW4T1_58 inp=sph_axi_MC3RAW4T1.inp user=umat_NH_g21_PCR int';
        p = 'abaqus cae noGUI=Extract_coordC58.py';
        
    % medium stiffness (0.85 kPa)
    elseif i == 2 && j == 1
        mat_props(4) = 0.00016; % C10
        mat_props(5) = 2824; % D1
        
        mat_props(6) = 0.78; % PCRA
        mat_props(7) =-0.57; % PCRB
        mat_props(8) = 0.22; % PCR6
        
        f = 'abaqus j=sph_axiM_RAW4T1_85 inp=sph_axi_RAW4T1.inp user=umat_NH_g21_PCR int';
        p = 'abaqus cae noGUI=Extract_coordA85.py';
        
    elseif i == 2 && j == 2
        mat_props(4) = 0.00016; % C10
        mat_props(5) = 2824; % D1
        
        mat_props(6) = -0.62; % PCRA
        mat_props(7) = -0.60; % PCRB
        mat_props(8) = 1; % PCR6  1.12
        
        f = 'abaqus j=sph_axiM_MC34T1_85 inp=sph_axi_MC34T1.inp user=umat_NH_g21_PCR int';
        p = 'abaqus cae noGUI=Extract_coordB85.py';
        
    elseif i == 2 && j == 3
        mat_props(4) = 0.00016; % C10
        mat_props(5) = 2824; % D1
        
        mat_props(6) = 1; % PCRA  2.78
        mat_props(7) = -0.71; % PCRB
        mat_props(8) = 0.68; % PCR6
        
        f = 'abaqus j=sph_axiM_MC3RAW4T1_85 inp=sph_axi_MC3RAW4T1.inp user=umat_NH_g21_PCR int';
        p = 'abaqus cae noGUI=Extract_coordC85.py';
        
    % high stiffness (1.1 kPa)
    elseif i == 3 && j == 1
        mat_props(4) = 0.00021; % C10
        mat_props(5) = 2182; % D1
        
        mat_props(6) = 0.28; % PCRA
        mat_props(7) = -0.80; % PCRB
        mat_props(8) = 0.65; % PCR6

        f = 'abaqus j=sph_axiM_RAW4T1_11 inp=sph_axi_RAW4T1.inp user=umat_NH_g21_PCR int';
        p = 'abaqus cae noGUI=Extract_coordA11.py';
        
    elseif i == 3 && j == 2
        mat_props(4) = 0.00021; % C10
        mat_props(5) = 2182; % D1      
        
        mat_props(6) = -0.66; % PCRA
        mat_props(7) = -0.67; % PCRB
        mat_props(8) = -0.16; % PCR6
        
        f = 'abaqus j=sph_axiM_MC34T1_11 inp=sph_axi_MC34T1.inp user=umat_NH_g21_PCR int';
        p = 'abaqus cae noGUI=Extract_coordB11.py';
        
    elseif i == 3 && j == 3
        mat_props(4) = 0.00021; % C10
        mat_props(5) = 2182; % D1
        
        mat_props(6) = 0.50; % PCRA
        mat_props(7) = -0.52; % PCRB
        mat_props(8) = -0.18; % PCR6
        
        f = 'abaqus j=sph_axiM_MC3RAW4T1_11 inp=sph_axi_MC3RAW4T1.inp user=umat_NH_g21_PCR int';
        p = 'abaqus cae noGUI=Extract_coordC11.py';
        
    end

%% Edit the .inp with new material properties
edit_inp( mat_props, 'sph_param.inp' )

%% Run the Abaqus job
% Run the simulations
job1 = system(f);
pause(1)

% Warning if jobs don't work
if (job1==1)
    disp('****WARNING*****')
    disp('Problem for Material Parameter Set')
    c
end

%% Run the python script
system(p);
pause(1)

    end
end

%% Gather force displacement data
cA58 = csvread('coordsA58.txt'); cA85 = csvread('coordsA85.txt'); cA11 = csvread('coordsA11.txt');
cB58 = csvread('coordsB58.txt'); cB85 = csvread('coordsB85.txt'); cB11 = csvread('coordsB11.txt');
cC58 = csvread('coordsC58.txt'); cC85 = csvread('coordsC85.txt'); cC11 = csvread('coordsC11.txt');
% tA58 = csvread('timeA58.txt'); tA85 = csvread('timeA85.txt'); tA11 = csvread('timeA11.txt');
% tB58 = csvread('timeB58.txt'); tB85 = csvread('timeB85.txt'); tB11 = csvread('timeB11.txt');
% tC58 = csvread('timeC58.txt'); tC85 = csvread('timeC85.txt'); tC11 = csvread('timeC11.txt');

% Calculate volume
VA58 = (4/3).*pi.*cA58(end).^2; VA85 = (4/3).*pi.*cA85(end).^2; VA11 = (4/3).*pi.*cA11(end).^2;
VB58 = (4/3).*pi.*cB58(end).^2; VB85 = (4/3).*pi.*cB85(end).^2; VB11 = (4/3).*pi.*cB11(end).^2;
VC58 = (4/3).*pi.*cC58(end).^2; VC85 = (4/3).*pi.*cC85(end).^2; VC11 = (4/3).*pi.*cC11(end).^2;

% Calculate diameters
dA58 = 2*cA58(end) ; dA85 = 2*cA85(end) ; dA11 = 2*cA11(end) ;
dB58 = 2*cB58(end) ; dB85 = 2*cB85(end) ; dB11 = 2*cB11(end) ;
dC58 = 2*cC58(end) ; dC85 = 2*cC85(end) ; dC11 = 2*cC11(end) ;

%% figure
hfig1=figure(); set(hfig1,'color','w'); hold all
bar(1e3*[dA58 dB58 dC58]); ylim([0 70])

hfig1=figure(); set(hfig1,'color','w'); hold all
bar(1e3*[dA85 dB85 dC85]); ylim([0 70])

hfig1=figure(); set(hfig1,'color','w'); hold all
bar(1e3*[dA11 dB11 dC11]); ylim([0 70])

% for i=1:3
%     xt = t1{i};
%     xncell = ncell{i};
%     plot(xt, xncell, '-','linewidth',4); 
% end
% errorbar(exp_act_Ctr(:,1), exp_act_Ctr(:,2), exp_act_Ctr(:,4),'ksq','linewidth',2, 'markersize', 6, 'MarkerFaceColor','k')
% box on; haxsY=gca;
% set(haxsY,'tickdir','out','linewidth',2.5,'fontsize',22,...
%     'FontName','Arial','FontWeight','bold')
% pbaspect([1.25 1 1])


%% functions

function edit_inp(newVal, inpFile)
% EDIT_INP Edit an Abaqus *.inp file to update the material properties
%   EDIT_INP(NEWVAL, INPFILE) takes in a vector of new material properties
%   named NEWVAL and writes a new inp file which has the *PARAMETER at the
%   top of the file. This is later added to the "main" inp file using the
%   *INCLUDE option

fin = fopen(inpFile,'w');
fprintf(fin,'%s\n', '*parameter');
fprintf(fin,'%s\n', ['xA = ' num2str(newVal(1),'%E') ]);
fprintf(fin,'%s\n', ['xB = ' num2str(newVal(2),'%E') ]);
fprintf(fin,'%s\n', ['x6 = ' num2str(newVal(3),'%E') ]);
fprintf(fin,'%s\n', ['C10 = ' num2str(newVal(4),'%E') ]);
fprintf(fin,'%s\n', ['D1 = ' num2str(newVal(5),'%E') ]);
fprintf(fin,'%s\n', ['pcA = ' num2str(newVal(6),'%E') ]);
fprintf(fin,'%s\n', ['pcB = ' num2str(newVal(7),'%E') ]);
fprintf(fin,'%s\n', ['pc6 = ' num2str(newVal(8),'%E') ]);
fclose(fin);
end

