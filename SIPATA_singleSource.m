% SIPATA - Sediment Input Point And Transport Analysis
% Version 1.0 started 2022, this version completed 17.10.2024 Peter Burgess, University of Liverpool
% See Burgess et al. 2024 Journal of the Geological Society of London for explanation of model, methods and the specific use to generate this output

% Single source, one transport axes best fit model runs, Burgess et al. Figure 4
figure('Position',[10, 10, 800, 800]);
oneTransportAxis = 1; % So 1 for one transport axis or 0 for multiple independent axes of sediment transport
sedimentSourceCount = 1;
inputDataFileName = 'outcropPebbleDataAll.csv';
subplot(3,1, 1);
pebblePrediction(inputDataFileName, sedimentSourceCount, oneTransportAxis);
inputDataFileName = 'outcropPebbleDataNoDistal.csv';
subplot(3,1, 2);
pebblePrediction(inputDataFileName, sedimentSourceCount, oneTransportAxis);
inputDataFileName = 'outcropPebbleDataNoIntermediate.csv';
subplot(3,1, 3);
pebblePrediction(inputDataFileName, sedimentSourceCount, oneTransportAxis);

function pebblePrediction(inputDataFileName, sedimentSourceCount, oneTransportAxis)

    % Read outcrop data file
    outcropData = readtable(inputDataFileName);
    outcropDistance = outcropData{:,4}; % NB squicgly brackets convert from table to double
    outcropPebbleConcentration = outcropData{:,5}; % NB squicgly brackets convert from table to double

    % Define and initialise model parameters for multiple source at x=0,400 and 480km, two lithologies, but pebbles from three input points
    rng(42);
    maxDist = 700;
    
    MCIterations = 200000; % 200000 is a good compromise between compute time and quality of fit
    MCIterationsOneTenth = MCIterations / 10.0;
    randScaleFactorLithoStartX = 2.0;
    randScaleFactorLithoVol = 2.0;
    randScaleFactorLithoDeposProportion = 10.0;
    errorRecord = nan(1,MCIterations);

    % Define the baseline parameters for two-lithology single source at x=0
    lithoNames = ["Sand", "Pebbles"];
    % Starting values for best fit single source models
    lithoStartX = [0, 0];
    lithoVol = [1000, 500];
    lithoDeposProportion = [0.001, 0.007];

    fprintf("MC iteration loop progress:");
    tic
    
    for j = 1:MCIterations

        % Set a vector of random values to adjust the volume and deposition proportion for the sediment source
        
        % lithoStartXAdjustment = [1, 1, proximalSourcesStartX];
        % lithoVolAdjustment = [1, 1, rand(1,(sedimentSourceCount-1)*2) .* randScaleFactorLithoVol];
        % lithoDeposProportionAdjustment = [1, 1, rand(1,(sedimentSourceCount-1)*2) .* randScaleFactorLithoDeposProportion];
   
        lithoVolAdjustment = rand(1,2) .* randScaleFactorLithoVol;
        lithoDeposProportionAdjustment = rand(1,2) .* randScaleFactorLithoDeposProportion;
        
        
        % Multiply the standard values by the random adjustments to calculate a single iteration of the model
        adjustedLithoVol = lithoVol .* lithoVolAdjustment;
        adjustedLithoDeposProportion = lithoDeposProportion .* lithoDeposProportionAdjustment;

        [~, error] = calculateLithologyProportions(sedimentSourceCount, maxDist, lithoStartX, adjustedLithoVol, adjustedLithoDeposProportion, outcropDistance, outcropPebbleConcentration, oneTransportAxis);

        nonNaNErrorRecord = errorRecord(~isnan(errorRecord));
        if isempty(nonNaNErrorRecord) || error < min(nonNaNErrorRecord) % empety error record indicates first iteration
        
            bestFitLithoVol = lithoVol .* lithoVolAdjustment;
            bestFitLithoDeposProportion = lithoDeposProportion .* lithoDeposProportionAdjustment;
        end

        errorRecord(j) = error;
        
        if mod(j, MCIterationsOneTenth) == 0
            fprintf("#");
        end
    end
    
    fprintf("\n");
    toc

    fprintf("%d sediment sources\n", sedimentSourceCount);
    fprintf("Outcrop data from %s\n", inputDataFileName);
    fprintf("Minimum error %5.4f from best fit parameters:\n", min(errorRecord(errorRecord > 0.0)));
    fprintf("Sediment volume ");
    fprintf("%4.0f ", bestFitLithoVol);
    fprintf("\n");
    fprintf("Sediment deposition proportion");
    fprintf("%5.4f ",bestFitLithoDeposProportion);
    fprintf("\n\n");

    [deposLitho, error] = calculateLithologyProportions(sedimentSourceCount, maxDist, lithoStartX, bestFitLithoVol, bestFitLithoDeposProportion, outcropDistance, outcropPebbleConcentration, oneTransportAxis);
    plotResults(oneTransportAxis, sedimentSourceCount, maxDist, lithoNames, deposLitho, outcropDistance, outcropPebbleConcentration, error)
end

function [deposLitho, error] = calculateLithologyProportions(sedimentSourceCount, maxDist, lithoStartX, lithoVol, lithoDeposProportion, outcropDistance, outcropPebbleConcentration, oneTransportAxis)

    lithoN = sedimentSourceCount * 2; % calculate number of lithologies -  2 lithologies per source
    deposLitho = zeros(maxDist, lithoN); % Initialise to avoid slow size-change memory allocation in loop below 

    for j = 1:lithoN % Calculate volume chnage with x distance for 2 lithologies for each sediment source
        for x = 1:maxDist
            if x >= lithoStartX(j) % check the sediment entry position along the profile
                deposLitho(x,j) = lithoVol(j) * lithoDeposProportion(j); % Record deposition of the specificed per-x-increment proportion
                lithoVol(j) = lithoVol(j) - deposLitho(x,j); % Remove deposited volume from remaining sediment volume for lithology j
            end
        end
    end

    % Calculate proportions of sediment volumes at each x position for each lithology
    totalThickness = sum(deposLitho,2); % Sum thickness at each x location along the profile
    j = 1:lithoN;
    deposLitho(:,j) = deposLitho(:,j) ./ totalThickness;
    
    % Calculate the error assuming either one transport axis, so all sediment sources are mixed together, therefore calculate one pebble proportion 
    if oneTransportAxis
        
        totalPebblesProportion = sum(deposLitho(:,2:2:lithoN), 2);
        errors = zeros(1, numel(outcropDistance));
        for j = 1:numel(outcropDistance)
            xIndex = round(outcropDistance(j)); % Position of the outcrop data points on the x distance transect
            errors(j) = abs(totalPebblesProportion(xIndex) - outcropPebbleConcentration(j));
        end
        
    else % or sediment sources remain separate, so calculate error compared to source that fits closest to observed pebble proportion at each location
        deposLitho(deposLitho == 0) = NaN; % Remove zeros by setting zero values to NaN
        % Calculate the error for pebble proportions between data points and model
        errors = zeros(1, numel(outcropDistance));
        for j = 1:numel(outcropDistance)
            % extract the n sediment proportions for jth outcrop data point x position
            xIndex = round(outcropDistance(j)); % get the vector index for the outcrop x km position
            modelPebbleProportions = deposLitho(xIndex,2:2:lithoN); % Extract the modelled proportions at x position
            [~,index] = min(abs(outcropPebbleConcentration(j) - modelPebbleProportions)); % Finds the position of the closest matching value in the modelPebbleProportions vector
            errors(j) = abs(modelPebbleProportions(index) - outcropPebbleConcentration(j));
        end
    end
    
    error = sum(errors);
end


function plotResults(oneTransportAxis, sedimentSourceCount, maxDist, lithoName, deposLitho, outcropDistance, outcropPebbleConcentration, error)
    
    lithoN = sedimentSourceCount * 2;
    if oneTransportAxis
        totalPebblesProportion = sum(deposLitho(:,2:2:lithoN), 2);
    end
    deposLitho(deposLitho == 0) = NaN;
    
    for j = 1:lithoN
        switch j
            case 1
                colourVector = [218/255,165/255,32/255]; % goldenrod
                lineWidth = 2.0;
                lineStyle = ":";
            case 2
                colourVector = [199/255,14/255,14/255];
                lineWidth = 2.0;
                lineStyle = "-";
            otherwise
                colourVector = [0.9,  j / (lithoN * 2), 0.3]; % Yellow hard to see on plots so lithoN*2 to ensure lower green value for orange lines, not yellow
                lineWidth = 1.0;
                if rem(j,2)
                    lineStyle = ":";
                else
                    lineStyle = "-";
                end
        end
        line(1:maxDist,deposLitho(:,j), "color", colourVector, "LineWidth", lineWidth, "LineStyle", lineStyle, "DisplayName", lithoName(j));
    end
    hold on;
    
    if oneTransportAxis
        line(1:maxDist,totalPebblesProportion, "color", [146/255, 43/255, 142/255], "LineWidth", 2.0, "DisplayName", "Total pebbles");
    end

    % Plot the outcrop pebble proportion data points
    scatter(outcropDistance, outcropPebbleConcentration, "b", "filled", "DisplayName", "Outcrops");
    
    xlabel("Distance (km)");
    ylabel("Proportion of sediment volume")
    ylim([0,1]);
    grid on
    legend
    drawnow
end
