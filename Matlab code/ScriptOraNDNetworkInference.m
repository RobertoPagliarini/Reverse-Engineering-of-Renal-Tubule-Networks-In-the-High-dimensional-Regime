%Clear the workspace
clear

%Creating folders
mkdir('RenalTubuleSegmentsData')

mkdir('RenalTubuleSegmentsNetworksOra')

%%%%%%%Generating Ntwork Models%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Segments = {'S1','S2','S3','DTL1','DTL2','DTL3','ATL','mTAL','cTAL','DCT','CNT','CCD','OMCD','IMCD'};

%Importing data from xlsx file
[num, txt, raw] = xlsread('KTEA_proteome_individual_sample.xlsx',1);

num(:,1) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumSamples = [3 3 3 3 3 3 3 4 3 3 3 4 3 3];

CumSum = cumsum(NumSamples);

GGNetResults = {};

GenesInModels = {};

for i = 1:length(NumSamples)
    
    if i == 1
        
        data = num(:,1:3);
        
    else
        
        data = num(:, CumSum(i-1)+1:CumSum(i));
        
    end
    
    dataclean = [];

    GeneInModel = {};
    
    for j = 1:length(data)
        
        if NumSamples(i) == 3
    
            if sum((data(j,:)) == zeros(1,3)) == 0
        
                dataclean(end+1,:) = data(j,:);
        
                GeneInModel{end+1,1} = txt{j+1,3};
            
            end
            
        else
            
            if sum((data(j,:)) == zeros(1,4)) == 0
        
                dataclean(end+1,:) = data(j,:);
        
                GeneInModel{end+1,1} = txt{j+1,3};
                
            end
            
        end
    
    end

    save(strcat(Segments{i},'.mat'),'data','dataclean',"GeneInModel");

    movefile(strcat(Segments{i},'.mat'),'RenalTubuleSegmentsData');
    
    [MatricesOraQuantile] = ORASchafferShrinkageNetworkInference(dataclean','Ora','Q',GeneInModel);

    [MatricesOraQuantileND] = ORASchafferShrinkageNetworkInference(dataclean','Ora','N',GeneInModel);

    save(strcat('AdjacencyMatrices_',Segments{i},'.mat'),'MatricesOraQuantile',"MatricesOraQuantileND");
    
    movefile(strcat('AdjacencyMatrices_',Segments{i},'.mat'),'RenalTubuleSegmentsNetworksOra');
    
end
