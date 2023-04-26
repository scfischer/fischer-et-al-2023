(* ::Package:: *)

(*Supplementary material Sabine C. Fischer, Elena Corujo-Sim\[OAcute]n, Joaqu\[IAcute]n Lilao, Ernst H.K. Stelzer and Silvia Mu\[NTilde]oz-Descalzo, "The transition from local to global patterns governs the differentiation of mouse blastocysts"*)


BeginPackage["naf`",{"cns`"}];


label::usage="extract embryo label from file path"


convertFromMINS::usage="convert data obtained from a segmentation with MINS in the format for the prostprocessing pipeline"


postProcessingPipelinePositions::usage="postprocessing pipeline (Schmitz et al Scientific Reports 2017), adapted for segmentation results that only give the positions of nuclei, not the whole segment; 
needs cns`in CellNucleiSegmentation_v8 available from https://www.physikalischebiologie.de/downloads"


rescaleCoordSphere::usage="embryos are slightly squashed during imaging, coordinates of cell centroids are rescaled, assuming that the embryo is a ball"


rescaleCoordScale::usage="embryos are slightly squashed during imaging, coordinates of cell centroids are rescaled according to a given scaling parameter"


compareValNeighboursFun::usage="determines for a cell its level of expression and the levels of expression of the neighbours"


compareValPopNeighboursFun::usage="determines for a cell its level of expression and the levels of expression and the population type of the neighbours"


compareValPopCentroidNeighboursFun::usage="determines for a cell its level of expression and the levels of expression, the population type and centroids of the neighbours"


prepareDataMoran::usage="determine values and weights for calculation of Moran's I from nucleiFeatures and Delaunay Cell Graph for an embryo"


calcMoranI::usage="calculates Moran's index for a set of values and corresponding weights; weights can e.g. be given as entries of adjacency matrix of corresponding graph if 
values are vertex properties of a graph"


Begin["`Private`"]; (* Begin Private Context *) 


label[fName_,splitCharacter1_:"_",splitPosition1_:1,splitCharacter2_:"-",splitPosition2_:-1]:=StringSplit[StringSplit[FileBaseName@fName,splitCharacter1][[splitPosition1]],splitCharacter2][[splitPosition2]]


Options[convertFromMINS]={"XinMicrons"->0.395,"YinMicrons"->0.395,"ZinMicrons"->1.0,"NumberMeasurements"->18,"CellFate"->"TE/ICM","CellFateThresholdRemove"->0,
"CellFateReplacementRule"->(Rule["TE/ICM",Switch[Round[#],1,"ICM+Epi",2,"TE",3,"PrE"]]&),
"ImagingChannels"->{"CH1-Avg"->"Dapi-Avg","CH1-Sum"->"Dapi-Sum","CH2-Avg"->"Brightfield-Avg","CH2-Sum"->"Brightfield-Sum","CH3-Avg"->"Gata6-Avg","CH3-Sum"->"Gata6-Sum",
"CH4-Avg"->"ABC-Avg","CH4-Sum"->"ABC-Sum","CH5-Avg"->"Nanog-Avg","CH5-Sum"->"Nanog-Sum"}};


convertFromMINS[inputFileName_,outputFileName_,opts:OptionsPattern[]]:=Block[{
dataRaw,nucleiMeasures,nucleiMeasuresRelevant,nucleiMeasuresCorrected,numMeasures=OptionValue["NumberMeasurements"],cellFate=OptionValue["CellFate"],
removeThreshold=OptionValue["CellFateThresholdRemove"],cellFateReplRule=OptionValue["CellFateReplacementRule"],
replaceChannels=OptionValue["ImagingChannels"],xum=OptionValue["XinMicrons"],yum=OptionValue["YinMicrons"],zum=OptionValue["ZinMicrons"]},
	If[FileExtension[inputFileName]=="csv",
		dataRaw=Import[inputFileName],
		dataRaw=Import[inputFileName][[1]]
	];
		nucleiMeasures=MapThread[Rule,{StringDelete[#," "]&/@dataRaw[[1,1;;numMeasures]],#}]&/@dataRaw[[2;;-1,1;;numMeasures]];
		nucleiMeasuresRelevant=Select[nucleiMeasures,(cellFate/.#)>removeThreshold&]/.replaceChannels;		
		nucleiMeasuresCorrected={"NucleiFeatures"->((DeleteCases[#,Rule["",""]]&/@nucleiMeasuresRelevant)/.Rule["Size",s_]->Rule["Size",s*xum*yum*zum]/.Rule[cellFate,t_]->cellFateReplRule[t]/.{a__,PatternSequence["X"->x_,"Y"->y_,"Z"->z_],b__}->{a,"Centroid"->{xum*x,yum*y,zum*z},b})};
	Export[outputFileName,nucleiMeasuresCorrected]
]


rescaleCoordSphere[coord_]:=Block[{extension,radius},
Sow[extension=Table[Apply[Plus,{Minus@Min[#[[All,i]]],Max[#[[All,i]]]}],{i,1,3}]&@coord];(*extension of the embryo along the x,y, and z axis*)
Sow[radius=r/.Flatten@Solve[RegionMeasure[Ball[{0,0,0},r]]==RegionMeasure@Ellipsoid[{0,0,0},extension/2],r,Reals]];
(*get the radius of a sphere with the same volume as the ellipsoid defined by the extension*)
Table[{Rescale[coord[[i,1]],Through[{Min,Max}[coord[[All,1]]]],{#-radius,#+radius}&@Mean@Through[{Min,Max}[coord[[All,1]]]]],
Rescale[coord[[i,2]],Through[{Min,Max}[coord[[All,2]]]],{#-radius,#+radius}&@Mean@Through[{Min,Max}[coord[[All,2]]]]],
Rescale[coord[[i,3]],Through[{Min,Max}[coord[[All,3]]]],{#-radius,#+radius}&@Mean@Through[{Min,Max}[coord[[All,3]]]]]},{i,1,Length[coord]}]
]


rescaleCoordScale[coord_,scale_]:=Block[{extension,adjustment},
Sow[extension=Table[Apply[Plus,{Minus@Min[#[[All,i]]],Max[#[[All,i]]]}],{i,1,3}]&@coord];(*extension of the embryo along the x,y, and z axis*)
Sow[adjustment=(extension*scale/2)];
(*get the radius of a sphere with the same volume as the ellipsoid defined by the extension*)
Table[{Rescale[coord[[i,1]],Through[{Min,Max}[coord[[All,1]]]],{#-adjustment[[1]],#+adjustment[[1]]}&@Mean@Through[{Min,Max}[coord[[All,1]]]]],
Rescale[coord[[i,2]],Through[{Min,Max}[coord[[All,2]]]],{#-adjustment[[2]],#+adjustment[[2]]}&@Mean@Through[{Min,Max}[coord[[All,2]]]]],
Rescale[coord[[i,3]],Through[{Min,Max}[coord[[All,3]]]],{#-adjustment[[3]],#+adjustment[[3]]}&@Mean@Through[{Min,Max}[coord[[All,3]]]]]},{i,1,Length[coord]}]
]


Options[postProcessingPipelinePositions]={"Alpha"-> 30,"OutlierDistanceThreshold"-> 20,"EdgeDistanceThreshold"->30,"IncludeSubfolders"->True,"AlphaShape"->True};
postProcessingPipelinePositions[folder_String,pattern_String,opts:OptionsPattern[]]:=Block[
{
alphaShapeQ=OptionValue["AlphaShape"],\[Alpha]=OptionValue["Alpha"],outlierDist=OptionValue["OutlierDistanceThreshold"],edgeDistance=OptionValue["EdgeDistanceThreshold"],
nucleiFeatures,expressionNames,centroids,surf,nucleiFeaturesUpdated,c,files,featuresAll,minExpression,maxExpression,time,dcg,pcg,alphaShape,delCellGraph,
pcgFeatures={"PCGNeighborCount","PCGClusteringCoefficient","PCGMinNeighborDistance","PCGMaxNeighborDistance","PCGMeanNeighborDistance","PCGStandardDeviationNeighborDistance"},
dcgFeatures={"DCGNeighborCount","DCGClusteringCoefficient","DCGMinNeighborDistance","DCGMaxNeighborDistance","DCGMeanNeighborDistance","DCGStandardDeviationNeighborDistance"}
},
files=If[OptionValue["IncludeSubfolders"],FileNames[pattern,{folder},Infinity],FileNames[pattern,{folder}]];
featuresAll=("NucleiFeatures"/.Import[#])&/@files;
expressionNames=Flatten[StringCases[featuresAll[[1,1,All,1]],__~~"-Avg"]];
{minExpression,maxExpression}=Transpose[Through[{Min,Max}[Select[#,NumberQ[#]==True&]]]&/@Transpose[Flatten[expressionNames/.featuresAll,1]]];
Map[
(
Print["dataset: ",#];
time=Round@First@AbsoluteTiming[

nucleiFeatures="NucleiFeatures"/.Import[#];(*import segmentation data*)
expressionNames=Flatten[StringCases[nucleiFeatures[[1,All,1]],__~~"-Avg"]];


nucleiFeaturesUpdated=If[alphaShapeQ,
	(Sow[alphaShape=Check[cns`Private`alphaShapes[First@ConnectedComponents@cns`Private`delaunayCellGraph["Centroid"/.nucleiFeatures, outlierDist], \[Alpha]],"alpha shape failed"]];(*compute alpha shape from nuclei centroids*)
	surf=RegionBoundary@alphaShape; (*approximate surface as the region boundary of the alpha shape*)
	cns`Private`appendSurfaceDistance[nucleiFeatures, surf]),
	nucleiFeatures
];



(*generate proximity cell graph and delaunay cell graph*)
pcg = cns`Private`proximityCellGraph["Centroid"/.nucleiFeaturesUpdated,edgeDistance];
dcg = cns`Private`delaunayCellGraph["Centroid"/.nucleiFeaturesUpdated,edgeDistance];

(*append graph features*)
nucleiFeaturesUpdated=appendNormalisedExpressionPerExp[cns`Private`appendGraphFeatures[cns`Private`appendGraphFeatures[nucleiFeaturesUpdated,pcg,pcgFeatures],dcg,dcgFeatures],expressionNames,minExpression,maxExpression];
(*export new data to mx file*)
Export[FileNameJoin[{DirectoryName[#],FileBaseName[#]<>"_Final.mx"}],
{
"NucleiFeatures"->nucleiFeaturesUpdated,
If[alphaShapeQ,
"GlobalFeatures"->{"CellNumberStage"->Switch[Length[nucleiFeatures],x_/;x<32,3.0,x_/;x<65,3.5,x_/;x<=90,4.0,x_/;x>90,4.5],"ProximityCellGraph"->pcg,"DelaunayCellGraph"->dcg,"Surface"->surf,"SurfaceArea"->Area@surf,"Volume"->RegionMeasure[alphaShape],"Centroid"->RegionCentroid[surf],"MinDistanceSurface"->RegionDistance[surf,RegionCentroid[surf]]},
"GlobalFeatures"->{"CellNumberStage"->Switch[Length[nucleiFeatures],x_/;x<32,3.0,x_/;x<65,3.5,x_/;x<=90,4.0,x_/;x>90,4.5],"ProximityCellGraph"->pcg,"DelaunayCellGraph"->dcg}
]}];
];
Export[FileNameJoin[{DirectoryName[#],StringTake[FileBaseName@(#),First@First@StringPosition[FileBaseName@(#),"Features"]-1]<>"Summary_Final.xlsx"}],
Prepend[Prepend[(Append[Rule@@@Transpose@{Keys@Options@postProcessingPipelinePositions,OptionValue[#]&/@Keys@Options@postProcessingPipelinePositions},"PostProcessingTimeSeconds"->time])/.Rule->List,{}],{"Summary","",DateString[]}]];
)&,
files
];
];


computeNormalisedExpressionPerExp[objects_List,names_List,minima_List,maxima_List]:=Block[{expressionValues},
expressionValues=Transpose[names/.objects];
Transpose[Table[Rescale[#,{minima[[i]],maxima[[i]]}]&/@expressionValues[[i]],{i,1,Length[expressionValues]}]]
]


appendNormalisedExpressionPerExp[objects_List,names_List,minima_List,maxima_List]:=Block[{expressionValuesNorm,namesNorm},
expressionValuesNorm=computeNormalisedExpressionPerExp[objects,names,minima,maxima];
namesNorm=Map[#<>"Norm"&,names];
Table[objects[[cur]]~Join~MapThread[Rule,{namesNorm,expressionValuesNorm[[cur]]}],{cur,1,Length@objects,1}]
]


compareValNeighboursFun[nucleiFeatures_,cellGraphs_,exprCell_,exprNeigh_]:=Block[{neighbours},
Table[neighbours=Cases[({Union[#[[All,1,1]]],#[[All,1,2]]}&/@GatherBy[ArrayRules[AdjacencyMatrix[cellGraphs[[k,2]]]],#[[1,1]]&]),{{_Integer},{___Integer}}];
{cellGraphs[[k,1]],Table[{{"TE/ICM",exprCell}/.(Select[nucleiFeatures[[k]],("Centroid"/.#)==(VertexList[cellGraphs[[k,2]]][[neighbours[[i,1,1]]]])&][[1]]),Table[exprNeigh/.(Select[nucleiFeatures[[k]],("Centroid"/.#)==(VertexList[cellGraphs[[k,2]]][[j]])&][[1]]),{j,neighbours[[i,2]]}]},{i,1,Length[neighbours]}]},{k,1,Length[cellGraphs]}]
]


(*compareValPopNeighboursFun[nucleiFeatures_,cellGraphs_,exprCell_,exprNeigh_]:=Block[{neighbours},
(*Table[neighbours=Cases[({Union[#[[All,1,1]]],#[[All,1,2]]}&/@GatherBy[ArrayRules[AdjacencyMatrix[cellGraphs[[k,2]]]],#[[1,1]]&]),{{_Integer},{___Integer}}];
{cellGraphs[[k,1]],Table[{{"TE/ICM","Quadrant",exprCell}/.(Select[nucleiFeatures[[k]],("Centroid"/.#)==(VertexList[cellGraphs[[k,2]]][[neighbours[[i,1,1]]]])&][[1]]),Table[{exprNeigh,"TE/ICM","Quadrant"}/.(Select[nucleiFeatures[[k]],("Centroid"/.#)==(VertexList[cellGraphs[[k,2]]][[j]])&][[1]]),{j,neighbours[[i,2]]}]},{i,1,Length[neighbours]}]},{k,1,Length[cellGraphs]}]*)
]*)


compareValPopNeighboursFun[nucleiFeatures_,cellGraphs_,exprCell_,exprNeigh_,fateLabel_:"Quadrant"]:=Block[{neighbours},
Table[neighbours=Cases[({Union[#[[All,1,1]]],#[[All,1,2]]}&/@GatherBy[ArrayRules[AdjacencyMatrix[cellGraphs[[k,2]]]],#[[1,1]]&]),{{_Integer},{___Integer}}];
{cellGraphs[[k,1]],Table[{{"TE/ICM",fateLabel,"Centroid",exprCell}/.(Select[nucleiFeatures[[k]],("Centroid"/.#)==(VertexList[cellGraphs[[k,2]]][[neighbours[[i,1,1]]]])&][[1]]),Table[{exprNeigh,"TE/ICM",fateLabel,"Centroid"}/.(Select[nucleiFeatures[[k]],("Centroid"/.#)==(VertexList[cellGraphs[[k,2]]][[j]])&][[1]]),{j,neighbours[[i,2]]}]},{i,1,Length[neighbours]}]},{k,1,Length[cellGraphs]}]
]


compareValPopNeighboursOrganoidsFun[nucleiFeatures_,cellGraphs_,exprCell_,exprNeigh_]:=Block[{neighbours},
Table[neighbours=Cases[({Union[#[[All,1,1]]],#[[All,1,2]]}&/@GatherBy[ArrayRules[AdjacencyMatrix[cellGraphs[[k,2]]]],#[[1,1]]&]),{{_Integer},{___Integer}}];
{cellGraphs[[k,1]],Table[{{"Label","Quadrant",exprCell}/.(Select[nucleiFeatures[[k]],("Centroid"/.#)==(VertexList[cellGraphs[[k,2]]][[neighbours[[i,1,1]]]])&][[1]]),Table[{exprNeigh,"Quadrant"}/.(Select[nucleiFeatures[[k]],("Centroid"/.#)==(VertexList[cellGraphs[[k,2]]][[j]])&][[1]]),{j,neighbours[[i,2]]}]},{i,1,Length[neighbours]}]},{k,1,Length[cellGraphs]}]
]


(* ::Input::Initialization:: *)
prepareDataMoran[nucleiFeatures_,dcgGraph_,protein_,population_:"Quadrant",throphectIdentifier_:"TE/ICM"]:=Block[{nucleiFeaturesICM,icmGraph,weights,values},
nucleiFeaturesICM=Select[nucleiFeatures,(throphectIdentifier/.#)!="TE"&];
icmGraph=Subgraph[dcgGraph,"Centroid"/.nucleiFeaturesICM];
weights=AdjacencyMatrix[icmGraph];
values=If[StringContainsQ[#,protein],1,0]&/@(population/.nucleiFeaturesICM);
{values,weights}
]


(* ::Input::Initialization:: *)
calcMoranI[values_,weights_]:=Block[{diffToMean,numerator,denominator},
diffToMean=values-Mean[values];
numerator=Outer[Times,diffToMean,diffToMean]*weights;
denominator=diffToMean^2;
N@Length[values]/Total[weights,2]*Total[numerator,2]/Total[denominator,2]
]


End[];(* End Private Context *)


EndPackage[];
