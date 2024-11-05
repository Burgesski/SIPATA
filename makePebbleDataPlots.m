% Code to create pebble data plots used in Figure 3, Burgess et al.

% Note the path ../PebbleData/ will need editing/removing dependong on where the listed data files are located
pebbleFieldworkTable = readtable("../PebbleData/PebbleDataTables.xlsx", "sheet", "Fieldwork data");
pebblePubTable = readtable("../PebbleData/PebbleDataTables.xlsx", "sheet", "Published data");
pebbleConcentrationDistanceTable = readtable("../PebbleData/PebbleDataTables.xlsx", "sheet", "summary");

numberOfRegions = 4; % South, central, east and west, controls overall layout of the figure
locationsPerRegion = [2,5,1,2];
concentrationDataColumns = [1,5,9,13,17,21]; % Columns in PebbleDataTables.xlsx fieldwork data sheet with pebble concentration data

pubLocations = pebblePubTable.Properties.VariableNames; % get location names for published and field work data tables
fieldworkLocations = {"Hawksmoor Wood",	"Hulme Quarry", "Acresford Quarry", "Nottingham Castle", "Chester Canal", "Farndon Cliffs"};

allLocations = cat(2, pubLocations, fieldworkLocations); % Combine all locations into one cell array
pubPebbleConcentrations = table2array(pebblePubTable); % Extract published pebble concentrations
pubPebbleDim = size(pubPebbleConcentrations); % size of the published data matrix, dimension 1 is now of elements in longest column
fieldworkPebbleConcentrations = table2array(pebbleFieldworkTable(:, concentrationDataColumns));
fieldworkPebbleDim = size(fieldworkPebbleConcentrations); % size of the fieldwork data matrix, dimension 1 is now of elements in longest column

longestColumn = max(fieldworkPebbleDim(1), pubPebbleDim(1)); % Find the longest of two data tables in terms of column length
if pubPebbleDim(1) < longestColumn % Check of published column length is less and pad it to the same as field work column length if so
    pubPebbleConcentrations = [pubPebbleConcentrations; nan(longestColumn - pubPebbleDim(1), pubPebbleDim(2))];
end

% Both matrixes should have same column length now, padded with nan, so can merge into one matrix for plotting
allPebbleConcentrations = cat(2, pubPebbleConcentrations, fieldworkPebbleConcentrations);

% Create the figure element showing how mean varies with subsamples for each location, element E in figure 3, Burgess et al
figure
hold on
xlabel("Sub-sampling sample size");
ylabel("Mean");
grid on

for j = 1:size(allPebbleConcentrations, 2)

    onePebbleConcentrationSample = allPebbleConcentrations(:,j);
    onePebbleConcentrationSample = onePebbleConcentrationSample(~isnan(onePebbleConcentrationSample)); % Remove all the NaN

    runningMean = NaN(1,numel(onePebbleConcentrationSample));
    for k = 1:numel(onePebbleConcentrationSample)
        runningMean(k) = mean(onePebbleConcentrationSample(1:k));
    end
    plot(runningMean, "lineWidth",2, "DisplayName",allLocations{j})
   
end
set(gca, 'FontSize', 16)
legend("NumColumns",3)

% figure element A-D in figure 3, Burgess et al, pebble concentration histograms
figure("Position", [10 10 800 900])
tiledlayout(numberOfRegions, max(locationsPerRegion))
plotCount = 1;
tileCount = 1;
for j =1:numberOfRegions

    for k =1:max(locationsPerRegion)

        if k <= locationsPerRegion(j)
            nexttile(tileCount);
            oneLocationPebbleConcentrations = allPebbleConcentrations(:,plotCount); 
            oneLocationPebbleConcentrations = oneLocationPebbleConcentrations(~isnan(allPebbleConcentrations(:,plotCount)));
            binEdges = 0:20:100;
            histogram(oneLocationPebbleConcentrations, binEdges, "facecolor",[0.7961, 0.4039,0.3569]);
            ylim([0,10]);
            grid on
            xlabel("Pebble %")
            ylabel("Frequency")
            title(allLocations(plotCount));
            plotCount = plotCount + 1;
        end
        tileCount = tileCount + 1;
    end
end

% Create figure element F in figure 3, Burgess et al showing pebble distances from Armorican Massif source
pebbleConcentrationDistance = table2array(pebbleConcentrationDistanceTable(:,7:8));
distanceData = pebbleConcentrationDistance(:,1);
concentrationData = pebbleConcentrationDistance(:,2);
figure
colour = [0.596, 0.984, 0.596; 0.596, 0.984, 0.596; ... % Green Wessex Basin points
    0.855, 0.898, 0.9412; 0.855, 0.898, 0.9412; 0.855, 0.898, 0.9412; 0.855, 0.898, 0.9412; 0.855, 0.898, 0.9412; ... % Pale blue Worcester Basin points
    1.00, 1.00, 0.773;  0.980, 0.784, 0.596;  0.980, 0.784, 0.596;]; % one pale yellow point for East Midlands, two pale orange points for Cheshire Basin
scatter(distanceData, concentrationData, 150, colour, "filled", "d", "MarkerEdgeColor",[0,0,0]);
grid on
set(gca, 'FontSize', 14)
xlabel("Distance from Armorician Massif (km)")
ylabel("Pebble concentration")

