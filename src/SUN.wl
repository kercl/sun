(* ::Package:: *)

BeginPackage["SUN`"];

	Irrep::usage = 
		"Prepare an irredcible representation identified by a Dynkin label.";
		
	IrrepDimension::usage = 
		"Compute the dimension of a representation.";
		
	IrrepRLBasisMatrices::usage = 
		"Compute the basis matrices of raising/lowering operators for the given irreducible representation.";

	IrrepRotMomentumBasisMatrices::usage =
		"Compute the basis matrices of rotation momentum operators for the given irreducible representation";

	IrrepDynkinLabel::usage = 
		"Get the Dynkin label of an irreducible representation.";
	
Begin["Private`"];
	sunLibCore = FindLibrary["sun_core"];

	Irrep/:MakeBoxes[
		Irrep[dynkin__Integer],form:(StandardForm|TraditionalForm)]:=Module[
			{above = {{BoxForm`SummaryItem[{
					   "\[GothicS]\[GothicU]("<>ToString[Length@List@dynkin + 1]<>") Irrep: ", 
					   "D("<>StringRiffle[List@dynkin,","]<>")"}]},
					  {BoxForm`SummaryItem[{"Dimension Irrep: ", IrrepDimension[List@dynkin]}]}},
			 below = {{BoxForm`SummaryItem[{"Dimension Lie Algebra: ", (Length[List@dynkin]+1)^2-1}]}},
			 icon = Graphics[({{PointSize[Medium],Point@#},Line@#}&@{{1,0},{1/2,2/Sqrt[3]},{-(1/2),2/Sqrt[3]},{-1,0},{-(1/2),-(2/Sqrt[3])},{1/2,-(2/Sqrt[3])},{1,0}})~Join~{
					  Line[{{-1,0},{1,0}}],Line[{{1/2,2/Sqrt[3]},{-(1/2),-(2/Sqrt[3])}}],Line[{{1/2,-(2/Sqrt[3])},{-(1/2),2/Sqrt[3]}}]},
					  AspectRatio->1, ImageSize -> Dynamic[{Automatic, 3.5 CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}]]

			},
			BoxForm`ArrangeSummaryBox[
				Irrep,(*head*)
				List[dynkin],(*interpretation*)
				icon,(*icon,use None if not needed*)
				above,(*always shown content*)
				below,(*expandable content*)
				form,
				"Interpretable"->True]];

	IrrepDynkinLabel[Irrep[dynkin__Integer]]:=List@dynkin;

	IrrepDimension[irrep_List] := Module[{
		topRow=Append[Reverse@Accumulate@Reverse@irrep, 0], 
		DimIrrep=LibraryFunctionLoad[sunLibCore,"ml_dimension_irrep",{{Integer,1}},Integer]
	},
		DimIrrep[topRow]
	];

	IrrepDimension[irrep_Irrep] := IrrepDimension[IrrepDynkinLabel[irrep]];

	IrrepBuildLoweringRootOperators[topRow_List,basis_,dim_]:=Module[{
		matrices={}, bID,
		RootLowering=LibraryFunctionLoad[sunLibCore,"ml_root_lowering",{Integer,Integer}, {Integer,_}]
	},
		matrices=Table[
			{Transpose[#[[3;;]]]->#[[1]]/#[[2]]}&@RootLowering[basis,l],
			{l,Length@topRow-1}];
		SparseArray[#,{dim,dim}]&/@matrices
	];
	
	LToPQ[l_]:=Module[{
		w=Floor[(Sqrt[8l+1]-1)/2], t},
		t = (w(w+1))/2;
		{l-t+1, w+1}
	];		

	IrrepRLBasisMatrices[irrep_Irrep] := Module[{
		topRow=Append[Reverse@Accumulate@Reverse@IrrepDynkinLabel[irrep], 0],
		dim=IrrepDimension[irrep],
		basis,cartan,lowering,
		loweringRoot,
		SUNInitBasis=LibraryFunctionLoad[sunLibCore, "ml_init_basis", {{Integer,1}, Integer}, "Void"],
		SUNCartan=LibraryFunctionLoad[sunLibCore, "ml_cartan_matrix", {Integer,Integer}, {Integer,_}]
	},
		basis=CreateManagedLibraryExpression["SUNGTBasis", SUNGTBasis];
		SUNInitBasis[topRow, ManagedLibraryExpressionID@basis];
		loweringRoot=IrrepBuildLoweringRootOperators[topRow,ManagedLibraryExpressionID@basis,dim];
		cartan=Table[SparseArray[
			Transpose[{Range@dim,Range@dim}]->SUNCartan[ManagedLibraryExpressionID@basis, l],
			{dim,dim}],
			{l,Length@topRow-1}];
		lowering=Table[Module[{chain=Range@@LToPQ[l-1], A},
			A=loweringRoot[[First@chain]];
			Do[A=B.A-A.B, {B, Rest@loweringRoot[[chain]]}];
			A],
		{l, (Length[topRow](Length[topRow]-1))/2}];
		Join[Flatten[Transpose@{lowering,ConjugateTranspose/@lowering}, 1], cartan]
	];
	
	IrrepRLBasisMatrices[dynkin__Integer] := IrrepRLBasisMatrices@Irrep[dynkin];

	IrrepRotMomentumBasisMatrices[irrep_Irrep] := Module[{
		Y=IrrepRLBasisMatrices[irrep],
		dynkin=IrrepDynkinLabel[irrep]
	},
		Join[Flatten[Table[{
			(Y[[i]]+Y[[i+1]])/2, (Y[[i]]-Y[[i+1]])/(2I)
		},
		{i, 1, Length[Y]-Length[dynkin], 2}],1], Y[[Length[Y]-Length[dynkin]+1;;]]]
	];
	
	IrrepRotMomentumBasisMatrices[dynkin__Integer]:=IrrepRotMomentumBasisMatrices@Irrep[dynkin];
	
End[];
EndPackage[];
