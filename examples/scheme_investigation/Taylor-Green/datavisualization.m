clear; close all 
LW=2.0;
timescale=32/pi;
nummethods=8;
alldata=cell(1,nummethods);
for i=1:nummethods
    alldata{1,i}=importdata(i);
end
linestylelist={"-o","-+","-*","-^","-x","-s","-*","-^"};
c_blue   = [71, 135, 224] / 255;
c_red    = [218, 62, 32] / 255;
c_green  = [82,147,47]/255;
c_yellow = [241,178,56]/255;
c_pink   = [149,20,83]/255;
c_purple = [138, 43, 226] / 255;
c_cyan   = [64, 224, 208] / 255;
c_orange = [255, 140, 0] / 255;
palette = {c_blue, c_red, c_green, c_yellow, c_pink, c_purple, c_cyan, c_orange};
%% Plot Kinetic Energy
% Total Kinetic Energy
figure
for i=1:nummethods
    if (i~=3 && i~=4 && i~=5)
        curr_data=alldata{i};
        idx = round(linspace(1, height(curr_data), 15)); 
        plot(curr_data.time*timescale,curr_data.KE/curr_data.KE(1),...
            linestylelist{i},...
            "LineWidth",LW,...
            'Color',palette{i},...
            'MarkerIndices',idx,...
            'MarkerSize', 10)
        hold on
    end 
end
hold off
xlabel("$tu/ \Delta x$","Interpreter","latex")
ylabel('KE/KE_0')
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',LW)% %% Post process weighted by its distance
% legend("S1","S2","S3","S4","S5","S6",'Location','Eastoutside')
xlim([1,100])
ylim([0.2,1.01])

% Liquid Kinetic Energy
figure
for i=1:nummethods
    if (i~=3 && i~=4 && i~=5)
        curr_data=alldata{i};
        idx = round(linspace(1, height(curr_data), 15)); 
        plot(curr_data.time*timescale,curr_data.KE_l/curr_data.KE(1),...
            linestylelist{i},...
            "LineWidth",LW,...
            'Color',palette{i},...
            'MarkerIndices',idx,...
            'MarkerSize', 10)
        hold on
    end
end
hold off
xlabel("$tu/ \Delta x$","Interpreter","latex")
ylabel('KE_l/KE_0')
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',LW)% %% Post process weighted by its distance
% legend({"1KE","1KE_l","1KE_g"})
xlim([1,100])
ylim([0.2,1.01])

% Gas Kinetic Energy
figure
for i=1:nummethods
    if (i~=3 && i~=4 && i~=5)
        curr_data=alldata{i};
        idx = round(linspace(1, height(curr_data), 15)); 
        plot(curr_data.time*timescale,curr_data.KE_g/curr_data.KE(1),...
            linestylelist{i},...
            "LineWidth",LW,...
            'Color',palette{i},...
            'MarkerIndices',idx,...
            'MarkerSize', 10)
        hold on
    end 
end
hold off
xlabel("$tu/ \Delta x$","Interpreter","latex")
ylabel('KE_g/KE_0')
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',LW)% %% Post process weighted by its distance
% legend({"1KE","1KE_l","1KE_g"})
xlim([1,100])
ylim([0.0005,0.06])



%% Plot U Momemtum
% Total U Momemtum
figure
for i=1:nummethods
    if (i~=3 && i~=4 && i~=5)
        curr_data=alldata{i};
        idx = round(linspace(1, height(curr_data), 15)); 
        plot(curr_data.time*timescale,curr_data.rhoU,...
            linestylelist{i},...
            "LineWidth",LW,...
            'Color',palette{i},...
            'MarkerIndices',idx,...
            'MarkerSize', 10)
        hold on
    end
end
hold off
xlabel("$tu/ \Delta x$","Interpreter","latex")
ylabel('$\rho U$',"Interpreter","latex")
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',LW)% %% Post process weighted by its distance
% legend("S1","S2","S3","S4","S5","S6",'Location','Eastoutside')
xlim([1,100])
% ylim([0.0005,1.01])

% Liquid Kinetic Energy
figure
for i=1:nummethods
    if (i~=3 && i~=4 && i~=5)
        curr_data=alldata{i};
        idx = round(linspace(1, height(curr_data), 15)); 
        plot(curr_data.time*timescale,curr_data.rhoU_l,...
            linestylelist{i},...
            "LineWidth",LW,...
            'Color',palette{i},...
            'MarkerIndices',idx,...
            'MarkerSize', 10)
        hold on
    end
end
hold off
xlabel("$tu/ \Delta x$","Interpreter","latex")
ylabel('$\rho_l U$',"Interpreter","latex")
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',LW)% %% Post process weighted by its distance
% legend({"1KE","1KE_l","1KE_g"})
xlim([1,100])
% ylim([0.0005,1.01])

% Gas Kinetic Energy
figure
for i=1:nummethods
    if (i~=3 && i~=4 && i~=5)
        curr_data=alldata{i};
        idx = round(linspace(1, height(curr_data), 15)); 
        plot(curr_data.time*timescale,curr_data.rhoU_g,...
            linestylelist{i},...
            "LineWidth",LW,...
            'Color',palette{i},...
            'MarkerIndices',idx,...
            'MarkerSize', 10)
        hold on
    end
end
hold off
ylabel('$\rho_g U$',"Interpreter","latex")
xlabel("$tu/ \Delta x$","Interpreter","latex")
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',LW)% %% Post process weighted by its distance
% legend({"1KE","1KE_l","1KE_g"})
xlim([1,100])
% ylim([0.0005,1.01])


%% Plot V Momemtum
% Total V Momemtum
figure
for i=1:nummethods
    if (i~=3 && i~=4 && i~=5)
        curr_data=alldata{i};
        idx = round(linspace(1, height(curr_data), 15)); 
        plot(curr_data.time*timescale,curr_data.rhoV,...
            linestylelist{i},...
            "LineWidth",LW,...
            'Color',palette{i},...
            'MarkerIndices',idx,...
            'MarkerSize', 10)
        hold on
    end
end
hold off
xlabel("$tu/ \Delta x$","Interpreter","latex")
ylabel('$\rho V$',"Interpreter","latex")
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',LW)% %% Post process weighted by its distance
% legend("S1","S2","S3","S4","S5","S6",'Location','Eastoutside')
xlim([1,100])
% ylim([0.0005,1.01])

% Liquid Kinetic Energy
figure
for i=1:nummethods
    if (i~=3 && i~=4 && i~=5)
        curr_data=alldata{i};
        idx = round(linspace(1, height(curr_data), 15)); 
        plot(curr_data.time*timescale,curr_data.rhoV_l,...
            linestylelist{i},...
            "LineWidth",LW,...
            'Color',palette{i},...
            'MarkerIndices',idx,...
            'MarkerSize', 10)
        hold on
    end
end
hold off
xlabel("$tu/ \Delta x$","Interpreter","latex")
ylabel('$\rho_l V$',"Interpreter","latex")
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',LW)% %% Post process weighted by its distance
% legend({"1KE","1KE_l","1KE_g"})
xlim([1,100])
% ylim([0.0005,1.01])

% Gas Kinetic Energy
figure
for i=1:nummethods
    if (i~=3 && i~=4 && i~=5)
        curr_data=alldata{i};
        idx = round(linspace(1, height(curr_data), 15)); 
        plot(curr_data.time*timescale,curr_data.rhoV_g,...
            linestylelist{i},...
            "LineWidth",LW,...
            'Color',palette{i},...
            'MarkerIndices',idx,...
            'MarkerSize', 10)
        hold on
    end
end
hold off
xlabel("$tu/ \Delta x$","Interpreter","latex")
ylabel('$\rho_g V$',"Interpreter","latex")
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',LW)% %% Post process weighted by its distance
% legend({"1KE","1KE_l","1KE_g"})
xlim([1,100])
% ylim([0.0005,1.01])


%% Plot W Momemtum
% Total U Momemtum
figure
for i=1:nummethods
    if (i~=3 && i~=4 && i~=5)
        curr_data=alldata{i};
        idx = round(linspace(1, height(curr_data), 15)); 
        plot(curr_data.time*timescale,curr_data.rhoW,...
            linestylelist{i},...
            "LineWidth",LW,...
            'Color',palette{i},...
            'MarkerIndices',idx,...
            'MarkerSize', 10)
        hold on
    end
end
hold off
xlabel("$tu/ \Delta x$","Interpreter","latex")
ylabel('$\rho W$',"Interpreter","latex")
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',LW)% %% Post process weighted by its distance
% legend("S1","S2","S3","S4","S5","S6",'Location','Eastoutside')
xlim([1,100])
% ylim([0.0005,1.01])

% Liquid Kinetic Energy
figure
for i=1:nummethods
    if (i~=3 && i~=4 && i~=5)
        curr_data=alldata{i};
        idx = round(linspace(1, height(curr_data), 15)); 
        plot(curr_data.time*timescale,curr_data.rhoW_l,...
            linestylelist{i},...
            "LineWidth",LW,...
            'Color',palette{i},...
            'MarkerIndices',idx,...
            'MarkerSize', 10)
        hold on
    end
end
hold off
xlabel("$tu/ \Delta x$","Interpreter","latex")
ylabel('$\rho_l W$',"Interpreter","latex")
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',LW)% %% Post process weighted by its distance
% legend({"1KE","1KE_l","1KE_g"})
xlim([1,100])
% ylim([0.0005,1.01])

% Gas Kinetic Energy
figure
for i=1:nummethods
    if (i~=3 && i~=4 && i~=5)
        curr_data=alldata{i};
        idx = round(linspace(1, height(curr_data), 15)); 
        plot(curr_data.time*timescale,curr_data.rhoW_g,...
            linestylelist{i},...
            "LineWidth",LW,...
            'Color',palette{i},...
            'MarkerIndices',idx,...
            'MarkerSize', 10)
        hold on
    end
end
hold off
xlabel("$tu/ \Delta x$","Interpreter","latex")
ylabel('$\rho_g W$',"Interpreter","latex")
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',LW)% %% Post process weighted by its distance
% legend({"1KE","1KE_l","1KE_g"})
xlim([1,100])
% ylim([0.0005,1.01])
%% Helper function
function [output] = importdata(id)
    opts = delimitedTextImportOptions("NumVariables", 17);
    
    % Specify range and delimiter
    opts.DataLines = [3, Inf];
    opts.Delimiter = " ";
    
    % Specify column names and types
    opts.VariableNames = ["time", "timestep", "KE", "KE_l", "KE_g", "rhoU", "rhoV", "rhoW", "rhoU_l", "rhoV_l", "rhoW_l", "rhoU_g", "rhoV_g", "rhoW_g", "VF", "rho", "Var17"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    opts.LeadingDelimitersRule = "ignore";
    
    % Specify variable properties
    opts = setvaropts(opts, "time", "TrimNonNumeric", true);
    opts = setvaropts(opts, "time", "ThousandsSeparator", ",");
    
    filename=strcat("/Users/zonghaozou/Repositories/NGA2/examples/scheme_investigation/Taylor-Green/result/",int2str(id),"/output_implicit.csv");
    % Import the data
    output = readtable(filename, opts);
   
    clear opts
end