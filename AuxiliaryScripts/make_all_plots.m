clear;
%%%%% CHANGE THESE TO MATCH YOUR NEEDS %%%%%
%outerDirs = ["NoRWA", "RWA"];
outerDirs = ["NoRWA"];
omega_0 = 0.360185;
p_z = 0.8097719;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Initialize data structures %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("Initializing data structures...\n")
numOuterDirs = length(outerDirs);
maxDirs = 0;
data = struct([]);
resData = struct([]);
posData = struct([]);
tempData = struct([]);
for i=1:numOuterDirs
    cd(outerDirs(i));
    numDirs = length(dir);
    if numDirs > maxDirs
        maxDirs = numDirs;
    end
    cd ..;
    data = [data, struct()];
    resData = [resData, resData()];
    posData = [posData, posData()];
    tempData = [tempData, tempData()];
end
dirStructs = struct([]);
for i = 1:length(outerDirs)
    data(i).data = 0;
    data(i).textdata = {};
    resData(i).data = [];
    posData(i).data = [];
    tempData(i).data = [];
    cd(outerDirs(i));
    innerDirs = dir;
    blankStructTemplate = struct();
    fn = fieldnames(innerDirs(1));
    numFields = length(fn);
    for j = 1:numFields
        blankStructTemplate.(fn{j}) = "";
    end
    numDirs = length(innerDirs);
    dirDiff = maxDirs - numDirs;
    fittedDirs = [innerDirs; repmat(blankStructTemplate,dirDiff,1)];
    dirStructs = [dirStructs, fittedDirs];
    cd ..;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:maxDirs
    if dirStructs(i,1).name(1) ~= '.' & dirStructs(i,1).name ~= "" & dirStructs(i,1).isdir == 1
        currentDir = dirStructs(i,1).name;
        matchingIndices = zeros(1,numOuterDirs-1);
        for j = 2:numOuterDirs
            for k = 1:maxDirs
                if matches(dirStructs(k,j).name, currentDir)
                    matchingIndices(j-1) = k;
                    break;
                end
            end
        end
        if ~ismember(0, matchingIndices)%%%%% Only make plot if the files are present everywhere %%%%%
            printStatement = strcat("Plotting ",currentDir, " ...\n");
            fprintf(printStatement);
            longestTraj = 0;
            for j = 1:numOuterDirs
                filePath = strcat(outerDirs(j),"/",currentDir,"/dataFile.txt");
                data(j) = importdata(filePath," ");
                resFilePath = strcat(outerDirs(j),"/",currentDir,"/resonanceVtime.txt");
                posFilePath = strcat(outerDirs(j),"/",currentDir,"/positionVtime.txt");
                tempFilePath = strcat(outerDirs(j),"/",currentDir,"/tempVtime.txt");
                if isfile(resFilePath)
                    plotRes = true;
                    resData(j).data = importdata(resFilePath," ");
                else
                    plotRes = false;
                end
                if isfile(posFilePath)
                    plotPos = true;
                    posData(j).data = importdata(posFilePath," ");
                    for k = 1:length(posData(j).data(:,2))
                        posData(j).data(k,2) = 10*(posData(j).data(k,2) - 0.5);
                    end
                else
                    plotPos = false;
                end
                if isfile(tempFilePath)
                    plotTemp = true;
                    tempData(j).data = importdata(tempFilePath," ");
                else
                    plotTemp = false;
                end
                %tmpData = importdata(filePath," ");
                %data(j) = tmpData;
                len = length(data(j).data);
                if len < longestTraj | longestTraj == 0
                    shortestIndex = j;
                    longestTraj = len;
                end
            end
            %%%%% Extraction/Calculation of data %%%%%
            tempCell = regexp(currentDir,"_","split");
            ECell = regexp(cell2mat(tempCell(1)),"E","split");
            EString = cast(cell2mat(ECell(2)),"char");
            E_0 = str2double(EString);
            wCell = regexp(cell2mat(tempCell(2)),"wRat","split");
            wString = cast(cell2mat(wCell(2)),"char");
            omega = str2double(wString)*omega_0;
            omega_r = 0.5*sqrt((omega_0-omega)^2+(E_0^2*p_z^2)/omega_0^2);
            inverse_omega_r_fs = pi/(41.341374575751*omega_r);
            period = num2str(inverse_omega_r_fs);
            RabiString = strcat("Rabi period = ",period," fs\n");
            %fprintf(RabiString);
            timeCell = data(shortestIndex).textdata(:,1);
            numTimeSteps = length(timeCell);
            time = zeros(numTimeSteps,1);
            popGS = zeros(numTimeSteps,numOuterDirs);
            popES = zeros(numTimeSteps,numOuterDirs);
            RabiGS = zeros(numTimeSteps,1);
            RabiES = zeros(numTimeSteps,1);
            for j = 1:numTimeSteps
                time(j) = str2double(cell2mat(timeCell(j)));
                RabiGS(j) = (((omega_0-omega)^2)/(4*omega_r^2))*sin(omega_r*time(j)*41.341374575751)^2 + cos(omega_r*time(j)*41.341374575751)^2;
                RabiES(j) = E_0^2*p_z^2/(4*omega_0^2*omega_r^2)*sin(omega_r*time(j)*41.341374575751)^2;
                %RabiGS(j) = 1 - RabiES(j);
                for k = 1:numOuterDirs
                    popGS(j,k) = str2double(cell2mat(data(k).textdata(j,6)));
                    popES(j,k) = str2double(cell2mat(data(k).textdata(j,7)));
                end
            end
            %%%%% Excited State Plots %%%%%
            f1 = figure('Visible','off','Position',[0 0 1618 1000]);
            hold on;
            for j = 1:numOuterDirs
                plot(time,popES(:,j),'LineWidth',8);
            end
            plot(time,RabiES,'LineWidth',4);
            set(gca,'FontSize',25);
            title('Excited State Population');
            xlabel('time (fs)');
            ylabel('Population');
            %legend("Numerical non-RWA Solution", "Numerical RWA Solution", 'Analytical Solution','Location','northwest');%%%%% Manually edit this line! %%%%%
            legend("Numerical non-RWA Solution", 'Analytical Solution', 'Location', 'northwest');%%%%% Manually edit this line! %%%%%
            figName = join(['plot_ES_pop_',currentDir,'.png']);
            saveas(f1, figName);
            %%%%% Ground State Plots %%%%%
            f2 = figure('Visible','off','Position',[0 0 1618 1000]);
            hold on;
            for j = 1:numOuterDirs
                plot(time,popGS(:,j),'LineWidth',8);
            end
            plot(time,RabiGS,'LineWidth',4);
            set(gca,'FontSize',25);
            title('Ground State Population');
            xlabel('time (fs)');
            ylabel('Population');
            %legend("Numerical non-RWA Solution", "Numerical RWA Solution", 'Analytical Solution','Location','southwest');%%%%% Manually edit this line! %%%%%
            legend("Numerical non-RWA Solution", 'Analytical Solution', 'Location', 'southwest');%%%%% Manually edit this line! %%%%%
            figName = join(['plot_GS_pop_',currentDir,'.png']);
            saveas(f2, figName);
            %%%%% Excited State Combo Plots %%%%%
            if plotRes
                f3 = figure('Visible','off','Position',[0 0 1618 1000]);
                hold on;
                for j = 1:numOuterDirs
                    plot(time,popES(:,j),'LineWidth',8);
                end
                plot(time,RabiES,'LineWidth',4);
                set(gca,'FontSize',25);
                title('Excited State Population and Resonant Frequency');
                xlabel('time (fs)');
                ylabel('Population');
                yyaxis right;
                for j = 1:numOuterDirs
                    plot(resData(j).data(:,1),resData(j).data(:,2),'LineWidth',4); 
                end
                ylabel('Resonant Frequency (a.u.)');
                %legend("Numerical non-RWA Solution", "Numerical RWA Solution", 'Analytical Solution','Location','northwest');%%%%% Manually edit this line! %%%%%
                legend("Numerical non-RWA Solution", 'Analytical Solution', 'Resonant Frequency','Location', 'northwest');%%%%% Manually edit this line! %%%%%
                figName = join(['plot_ES_vs_Res_',currentDir,'.png']);
                saveas(f3, figName);
                yyaxis left;
            end
            if plotPos
                f4 = figure('Visible','off','Position',[0 0 1618 1000]);
                hold on;
                for j = 1:numOuterDirs
                    plot(time,popES(:,j),'LineWidth',8);
                end
                plot(time,RabiES,'LineWidth',4);
                set(gca,'FontSize',25);
                title('Excited State Population and Internuclear Distance');
                xlabel('time (fs)');
                ylabel('Population');
                yyaxis right;
                for j = 1:numOuterDirs
                    plot(posData(j).data(:,1),posData(j).data(:,2),'LineWidth',4); 
                end
                ylabel('Internuclear Distance (\AA)', 'Interpreter', 'latex');
                %legend("Numerical non-RWA Solution", "Numerical RWA Solution", 'Analytical Solution','Location','northwest');%%%%% Manually edit this line! %%%%%
                legend("Numerical non-RWA Solution", 'Analytical Solution', 'Internuclear Distance','Location', 'northwest');%%%%% Manually edit this line! %%%%%
                figName = join(['plot_ES_vs_Dist_',currentDir,'.png']);
                saveas(f4, figName);
                yyaxis left;
            end
            if plotTemp
                f5 = figure('Visible','off','Position',[0 0 1618 1000]);
                hold on;
                for j = 1:numOuterDirs
                    plot(time,popES(:,j),'LineWidth',8);
                end
                plot(time,RabiES,'LineWidth',4);
                set(gca,'FontSize',25);
                title('Excited State Population and Temperature');
                xlabel('time (fs)');
                ylabel('Population');
                yyaxis right;
                for j = 1:numOuterDirs
                    plot(tempData(j).data(:,1),tempData(j).data(:,2),'LineWidth',4); 
                end
                ylabel('Temperature (K)');
                %legend("Numerical non-RWA Solution", "Numerical RWA Solution", 'Analytical Solution','Location','northwest');%%%%% Manually edit this line! %%%%%
                legend("Numerical non-RWA Solution", 'Analytical Solution', 'Temperature', 'Location', 'northwest');%%%%% Manually edit this line! %%%%%
                figName = join(['plot_ES_vs_Temp_',currentDir,'.png']);
                saveas(f5, figName);
                yyaxis left;
            end
        end
    end
end