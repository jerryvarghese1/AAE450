clear all
clc

%% Data Entry
%Each subcategory of options is defined below as a matrix. Each row of the
%matrix represents an option. The cost FOA value is in the first column,
%the risk in the second column, and the science return value in the third
%column.

%Mission Design
missionDesign = [10   0  10;    %Uranus Orbit + Moon Observation
                 9.16 1  9.05;  %Uranus Orbit + Single Moon Flyby
                 5.38 10 2.73]; %Uranus Orbit + Multi Moon Flyby/Observation 

%Earth Data Return
earthDataReturn = [2.9    4.3   5.7;   %34m antenna full time
                   1.59   5.88  8.2;   %70m antenna
                   2.59   4.52  6.02;  %34m Antenna majority time with 70m for essentials
                   2.4    4.26  5.74;  %34m for cruise, 70m for Uranus science portion
                   0.0085 1.525 1.96;  %26m Antenna Full time
                   2.41   4.175 5.645];%Combination of three(26m NEO, 34m majority, 70m essential)
               
%Data Reception
dataReception = [5   3.7  3.2;  %4m (Cassini)
                 4.9 3.87 3.47; %5m parabolic (falcon fairing)
                 4.2 3.95 4.4;  %7m folding reflectarry (falcon fairing)
                 1.2 4.3  6.2]; %20m circular folding reflectarray 

%Data Frequency
dataFrequency = [0.472 14.448 23.619; %KA Band full time(DSN) 
                 0.847 8.019  13.209; %X band full time(DSN)
                 0.532 12.823 20.960; %KA band science, X band all else
                 0.612 11.074 18.096; %X band science, KA  band all else
                 1.574 4.359  7.496]; %S band full time
                 
%Number of Spacecraft
numberOfSC = [6.618 7.786 0.671;  %1 spacecraft
              3.950 4.001 1.150;  %2 spacecraft
              1.900 2.084 1.250]; %Multiple spacecraft (DART)
          
%Mission Duration (After reaching Uranus)
missionDuration = [11.000 14.830 0.750;   %Hours
                 11.000 13.900 2.500;   %Days
                 11.000 11.550 6.500;   %Months
                 10.000 10.000 10.000]; %Years

%Power System
powerSystem = [2.000 10.000 0.000;  %Solar Panels
               8.000 10.000 2.000;  %RTG
               7.000 6.0000 10.000; %Kilopower
               9.000 5.000 4.000];  %Dynamic RPS
         
count = 1;  %Initialize a counter variable to keep track of the total number of permutations

%% Gimme the Loops
for i = 1:size(missionDesign,1) %Seperate look that goes down the length of each category to permute every possible option
    for j = 1:size(earthDataReturn,1)
        for k = 1:size(dataReception,1)
            for l = 1:size(dataFrequency,1)
                for m = 1:size(numberOfSC,1)
                    for n = 1:size(missionDuration,1)
                        for o = 1:size(powerSystem,1)
                            
                            %Calculate every possible score for Cost
                            score(count,1) = missionDesign(i,1) + ...
                                         earthDataReturn(j,1) + ...
                                         dataReception(k,1) + ...
                                         dataFrequency(l,1) + ...
                                         numberOfSC(m,1) + ...
                                         missionDuration(n,1) + ...
                                         powerSystem(o,1);
                                        
                            %Calculate every possible score for Risk
                            score(count,2) = missionDesign(i,2) + ... 
                                         earthDataReturn(j,2) + ...
                                         dataReception(k,2) + ...
                                         dataFrequency(l,2) + ...
                                         numberOfSC(m,2) + ...
                                         missionDuration(n,2) + ...
                                         powerSystem(o,2);
                                        
                            %Calculate every possible score for Science Return
                            score(count,3) = missionDesign(i,3) + ...
                                         earthDataReturn(j,3) + ... 
                                         dataReception(k,3) + ...
                                         dataFrequency(l,3) + ...
                                         numberOfSC(m,3) + ...
                                         missionDuration(n,3) + ...
                                         powerSystem(o,3);                                       
                            
                            %Uncomment these for the pretty rainbow colored 3d plot
                            figure(1)
                            plot3(score(count,1),score(count,2),score(count,3),'.')
                            hold on
                            
                            tags(count,:) = [i j k l m n o]; %Create a matrix of 'tags' that record which option choices correspond to each score
                            count = count+1; %iterate the counter variable
                        
                        end
                    end
                end
            end
        end
    end
end


for i = 1:length(tags)
    if (tags(i,1) == 1 || tags(i,1) == 2) && tags(i,5) == 1
        score(i,:) = 0;
    end
end


%Create a list of aggregate scores for each configuration
aggScore = sum(score,2);

%Find the top 5 aggregate scores and list their index location
[best5, locs] = maxk(aggScore,5);

%Select the tag arrays corresponding to each of the top 5 aggregate scores
topfive = tags(locs,:);

%Formatted output of the top 5 configurations
printOut(topfive)

%% Plotting 

%Plot each of the three figures of merit against each other

%Risk vs Cost
subplot(1,3,1)
plot(score(:,2), score(:,1),'.')
title('Risk vs. Cost')
xlabel('Risk')
ylabel('Cost')
grid on

%Science Return vs Cost
subplot(1,3,2)
plot(score(:,3), score(:,1),'.')
title('Science Return vs. Cost')
xlabel('Science Return')
ylabel('Cost')
grid on

%Risk vs Science Return
subplot(1,3,3)
plot(score(:,2), score(:,3),'.')
title('Risk vs. Science Return')
xlabel('Risk')
ylabel('Science Return')
grid on

%% Output formatting

function printOut(tags)

%Strings for every option

%Mission Design
mdes = ["Uranus Orbit + Moon Observation"
        "Uranus Orbit + Single Moon Flyby"
        "Uranus Orbit + Multi Moon Flyby/Observation"];

%Earth Data Return
dret = ["34m antenna full time"
        "70m antenna"
        "34m Antenna majority time with 70m for essentials"
        "34m for cruise, 70m for Uranus science portion"
        "26m Antenna Full time"
        "Combination of three(26m NEO, 34m majority, 70m essential)"];
               
%Data Reception
drec = ["4m (Cassini)"
        "5m parabolic (falcon fairing)"
        "7m folding reflectarry (falcon fairing)"
        "20m circular folding reflectarray"];

%Data Frequency
dfreq = ["KA Band full time(DSN)"
         "X band full time(DSN)"
         "KA band science, X band all else"
         "X band science, KA  band all else"
         "S band full time"];
                 
%Number of Spacecraft
num = ["1 spacecraft"
       "2 spacecraft"
       "Multiple spacecraft (DART)"];
          
%Mission Duration
mdur = ["Hours"
        "Days"
        "Months"
        "Years"];

%Power System
psys = ["Solar Panels"
        "RTG"
        "Kilopower"
        "Dynamic RPS"
        "Fusion Reactor"];

fprintf('Top 5 Configurations:\n\n')

%Prints the top 5 options, selecting the string for each category using the
%previously created 'tags' arrays
for i = 1:5
fprintf(['%d. Mission Design:    %s\n',...  
          '   Earth Data Return: %s\n',... 
          '   Data Reception:    %s\n',... 
          '   Data Frequency:    %s\n',...
          '   Number of S/C:     %s\n',...
          '   Mission Duration:  %s\n',... 
          '   Power System:      %s\n\n'],i,...
            mdes(tags(i,1)),...
            dret(tags(i,2)),...
            drec(tags(i,3)),...
            dfreq(tags(i,4)),...
            num(tags(i,5)),...
            mdur(tags(i,6)),...
            psys(tags(i,7)))
end

end