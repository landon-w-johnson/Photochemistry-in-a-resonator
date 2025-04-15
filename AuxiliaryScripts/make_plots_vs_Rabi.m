clear;
%%%%% CHANGE THESE TO MATCH YOUR NEEDS %%%%%
outerDirs = ["NoRWA", "RWA"];
omega_0 = 0.36018531627654671;
p_z = 0.80977192915387741;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Initialize data structures %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("Initializing data structures...\n")
numOuterDirs = length(outerDirs);
maxDirs = 0;
data = struct([]);
for i=1:numOuterDirs
    cd(outerDirs(i));
    numDirs = length(dir);
    if numDirs > maxDirs
        maxDirs = numDirs;
    end
    cd ..;
    data = [data, struct()];
end
dirStructs = struct([]);
for i = 1:length(outerDirs)
    data(i).data = 0;
    data(i).textdata = {};
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
            wCell = regexp(cell2mat(tempCell(2)),"w","split");
            wString = cast(cell2mat(wCell(2)),"char");
            omega = str2double(wString);
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
            legend("Numerical non-RWA Solution", "Numerical RWA Solution", 'Analytical Solution','Location','northwest');%%%%% Manually edit this line! %%%%%
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
            legend("Numerical non-RWA Solution", "Numerical RWA Solution", 'Analytical Solution','Location','southwest');%%%%% Manually edit this line! %%%%%
            figName = join(['plot_GS_pop_',currentDir,'.png']);
            saveas(f2, figName);
        end
    end
end