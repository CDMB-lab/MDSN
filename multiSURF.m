function [multiSURFScores,locations_multiSURF,factor]=multiSURF(DataFile)

%%
inFile=[DataFile,'.mat'];
factor=[];
%% download data
tmp=load(inFile);
pts=tmp.pts;
class=tmp.class;
IsFactor=0;
if isfield(tmp,'factor')
    factor=tmp.factor;
    IsFactor=1;
end
[SampleNum, SNPNum]=size(pts);
clear tmp;

%% Distance, the number of different SNPs between two samples
DistanceMatrix=zeros(SampleNum);
for i=1:SampleNum
    SampleI=repmat(pts(i,:),SampleNum,1);
    diff=SampleI~=pts;
    DistanceI=sum(transpose(diff));
    DistanceMatrix(i,:)=DistanceI;
    DistanceMatrix(:,i)=DistanceI';
end

[SortedDistance,SortedNeighbors]=sort(DistanceMatrix);
SortedDistance(1,:)=[];
SortedNeighbors(1,:)=[];
% clear DistanceMatrix SampleI diff DistanceI

%% Nearest neighbors
thresholds=mean(SortedDistance,1)-std(SortedDistance,1)/2;
thresholdsMatrix=repmat(thresholds,SampleNum-1,1);
NearestNeighborsTag=SortedDistance<=thresholdsMatrix;
NearestNeighborsSum=sum(NearestNeighborsTag);
NearestNeighbors(2,:)=SortedNeighbors(NearestNeighborsTag);

Tag=0;
for i=1:SampleNum
    for j=1:NearestNeighborsSum(i)
        Tag=Tag+1;
        NearestNeighbors(1,Tag)=i;
        
        if class(i)==class(NearestNeighbors(2,Tag))
            NearestNeighbors(3,Tag)=0; % hit
        else
            NearestNeighbors(3,Tag)=1; % miss
        end
    end
end

Hits=NearestNeighbors(1:2,NearestNeighbors(3,:)==0);
Misses=NearestNeighbors(1:2,NearestNeighbors(3,:)==1);
HitsNum=size(Hits,2);
MissesNum=size(Misses,2);
% clear thresholds thresholdsMatrix NearestNeighborsTag NearestNeighborsSum ...
%     NearestNeighbors Tag SortedDistance SortedNeighbors

%% multiSURF Scores 
multiSURFScores=zeros(SNPNum,2);
MissesScore=zeros(1,SNPNum);
HitsScore=zeros(1,SNPNum);
for i=1:SNPNum
    for j=1:MissesNum
        if pts(Misses(1,j),i)~=pts(Misses(2,j),i)
            MissesScore(1,i)=MissesScore(1,i)+1;
        end
    end
    for j=1:HitsNum
        if pts(Hits(1,j),i)~=pts(Hits(2,j),i)
             HitsScore(1,i)=HitsScore(1,i)+1;
        end
    end
end

%%%%%%%%%%%%%%%%%%
%% mean(M)-mean(H)
meanMisses=MissesScore/MissesNum;
meanHits=HitsScore/HitsNum;
scores=meanMisses-meanHits;
[multiSURFScores(:,2),multiSURFScores(:,1)]=sort(scores,'descend');
locations_multiSURF=[];
if IsFactor
    [~,locations_multiSURF]=ismember(factor,multiSURFScores(:,1));
end

Outfile=[DataFile,'_multiSURF.mat'];
save(Outfile,'factor','multiSURFScores','locations_multiSURF');

%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%
% %% abs(mean(M)-mean(H))
% meanMisses=MissesScore/MissesNum;
% meanHits=HitsScore/HitsNum;
% scores=abs(meanMisses-meanHits);
% [multiSURFScores(:,2),multiSURFScores(:,1)]=sort(scores,'descend');
% locations_multiSURF=[];
% if IsFactor
%     [~,locations_multiSURF]=ismember(factor,multiSURFScores(:,1));
% end
% %%%%%%%%%%%%%%%%%%
% 
% clear Hits Misses
% 
% %% STIR Scores
% STIRScores=zeros(SNPNum,3);
% StdMisses=((1-meanMisses).^2.*MissesScore+(0-meanMisses).^2.*(MissesNum-MissesScore))/MissesNum;
% StdHits=((1-meanHits).^2.*HitsScore+(0-meanHits).^2.*(HitsNum-HitsScore))/HitsNum;
% Sp=sqrt(((MissesNum-1)*StdMisses+(HitsNum-1)*StdHits)/(MissesNum+HitsNum-2));
% scores2=scores./(Sp*sqrt(1/MissesNum+1/HitsNum));
% [STIRScores(:,2),STIRScores(:,1)]=sort(scores2,'descend');
% locations_STIR=[];
% if IsFactor
%     [~,locations_STIR]=ismember(factor,STIRScores(:,1));
% end
% 
% %% Pvalue
% STIRScores(:,3)=1-tcdf(STIRScores(:,2),MissesNum+HitsNum-2);
% TopN=sum(STIRScores(:,3)<=Pvalue);
% Detection=[];
% if IsFactor
% %    Detect=double(locations_STIR<=TopN);
%     Detection=factor(locations_STIR<=TopN);
% end
% SNPs_STIR=STIRScores(1:TopN,1);
% 
% clear pts class i j meanHits meanMisses MissesNum MissesScore HitsNum HitsScore ...
%     scores scores2 Sp StdHits StdMisses IsFactor
% Outfile=[DataFile,'_Relief.mat'];
% save(Outfile,'SNPs_STIR','TopN','STIRScores','multiSURFScores','Detection','locations_STIR','locations_multiSURF','factor','multiSURFScoresNA','locations_multiSURFNA')
