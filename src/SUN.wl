(* ::Package:: *)

BeginPackage["SUN`"];

	Irrep::usage = 
		"Irredcible representation identified by a Dynkin label.";
	
	IrrepBasis::usage =
		"Set of basis matrices of a particular irreducible representation.";
	
	IrrepDimension::usage = 
		"Compute the dimension of a representation.";
	
	LieAlgebraBasisMatrices::usage =
		"Generate set of basis matrices spanning the particular representation";
	
	IrrepDynkinLabel::usage = 
		"Get the Dynkin label of an irreducible representation.";
	
	SU3BasisMatrices::usage =
		"Get \[GothicS]\[GothicU](3) Irrep basis";
	
	SU2BasisMatrices::usage =
		"Get \[GothicS]\[GothicU](2) Irrep basis";
	
Begin["Private`"];
	sunLibCore = FindLibrary["sun_core"];

	Irrep/:MakeBoxes[
		Irrep[dynkin__Integer],form:(StandardForm|TraditionalForm)]:=Module[
			{above = {{BoxForm`SummaryItem[{
					   "\[GothicS]\[GothicU]("<>ToString[Length@List@dynkin + 1]<>") Irrep: ", 
					   "D("<>StringRiffle[List@dynkin,","]<>")"}]},
					  {BoxForm`SummaryItem[{"Dimension Irrep: ", IrrepDimensionDynkin[List@dynkin]}]}},
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

	IrrepBasis/:MakeBoxes[
		IrrepBasis[matrices_List, dynkin_List],form:(StandardForm|TraditionalForm)]:=Module[
			{above = {{BoxForm`SummaryItem[{
					   "\[GothicS]\[GothicU]("<>ToString[Length@dynkin + 1]<>") representation basis: ", 
					   "D("<>StringRiffle[dynkin,","]<>")"}]},
					  {BoxForm`SummaryItem[{"Dimension Irrep: ", IrrepDimensionDynkin[dynkin]}]}},
			 below = {},
			 icon = Graphics[({{PointSize[Medium],Point@#},Line@#}&@{{1,0},{1/2,2/Sqrt[3]},{-(1/2),2/Sqrt[3]},{-1,0},{-(1/2),-(2/Sqrt[3])},{1/2,-(2/Sqrt[3])},{1,0}})~Join~{
					  Line[{{-1,0},{1,0}}],Line[{{1/2,2/Sqrt[3]},{-(1/2),-(2/Sqrt[3])}}],Line[{{1/2,-(2/Sqrt[3])},{-(1/2),2/Sqrt[3]}}]},
					  AspectRatio->1, ImageSize -> Dynamic[{Automatic, 3.5 CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}]]

			},
			BoxForm`ArrangeSummaryBox[
				IrrepBasis,(*head*)
				matrices,(*interpretation*)
				icon,(*icon,use None if not needed*)
				above,(*always shown content*)
				below,(*expandable content*)
				form,
				"Interpretable"->True]];

	IrrepBasis[matrices_List, dynkin_List][k_] := matrices[[k]];
	Unprotect[Length];
	Length[IrrepBasis[matrices_List, dynkin_List]] := Length[matrices];
	Protect[Length];

	IrrepDynkinLabel[Irrep[dynkin__Integer]]:=List@dynkin;

	IrrepDimensionDynkin[irrep_List] := Module[{
		topRow=Append[Reverse@Accumulate@Reverse@irrep, 0], 
		DimIrrep=LibraryFunctionLoad[sunLibCore,"ml_dimension_irrep",{{Integer,1}},Integer]
	},
		DimIrrep[topRow]
	];

	IrrepDimension[irrep_Irrep] := IrrepDimensionDynkin[IrrepDynkinLabel[irrep]];

	IrrepBuildLoweringRootOperators[topRow_List, basis_, dim_]:=Module[{
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
		loweringRoot, res,
		SUNInitBasis=LibraryFunctionLoad[sunLibCore, "ml_init_basis", {{Integer,1}, Integer}, "Void"],
		SUNCartan=LibraryFunctionLoad[sunLibCore, "ml_cartan_matrix", {Integer,Integer}, {Integer,_}]
	},
		basis=CreateManagedLibraryExpression["SUNGTBasis", SUNGTBasis];
		SUNInitBasis[topRow, ManagedLibraryExpressionID@basis];
		loweringRoot=IrrepBuildLoweringRootOperators[topRow,ManagedLibraryExpressionID@basis,dim];
		cartan=Table[SparseArray[
			Transpose[{Range@dim,Range@dim}]->1/2 SUNCartan[ManagedLibraryExpressionID@basis, l],
			{dim,dim}],
			{l,Length@topRow-1}];
		lowering=Table[Module[{chain=Range@@LToPQ[l-1], A},
			A=loweringRoot[[First@chain]];
			Do[A=B.A-A.B, {B, Rest@loweringRoot[[chain]]}];
			A],
		{l, (Length[topRow](Length[topRow]-1))/2}];
		res = Join[Flatten[Transpose@{lowering,ConjugateTranspose/@lowering}, 1], cartan];
		IrrepBasis[res, IrrepDynkinLabel[irrep]]
	];

	IrrepAngMomBasisMatrices[irrep_Irrep] := Module[{
		Y=IrrepRLBasisMatrices[irrep],
		dynkin=IrrepDynkinLabel[irrep],
		res
	},
		res = Join[Flatten[Table[{
			(Y[i]+Y[i+1])/2, (Y[i]-Y[i+1])/(2I)
		},
		{i, 1, Length[Y]-Length[dynkin], 2}],1], Y[Length[Y]-Length[dynkin]+1;;]];
		IrrepBasis[res, dynkin]
	];

	Options[LieAlgebraBasisMatrices]={BasisType->"AngularMomentum"};
	LieAlgebraBasisMatrices[irrep_Irrep, OptionsPattern[]]:=Switch[OptionValue[BasisType],
		"AngularMomentum", IrrepAngMomBasisMatrices[irrep],
		"LoweringRaising", IrrepRLBasisMatrices[irrep]
	];

	Options[SU3BasisMatrices]={BasisType->"GellMann"};
	SU3BasisMatrices[p_, q_, OptionsPattern[]] := Module[{
		coords,
		X=If[OptionValue[BasisType]=="GellMann",
			 LieAlgebraBasisMatrices@Irrep[p, q],
			 LieAlgebraBasisMatrices[Irrep[p, q], BasisType->OptionValue[BasisType]]]
	},
		If[OptionValue[BasisType]=="GellMann",
			X = X[{1, 2, 7, 3, 4, 5, 6, 8}];
			X[[8]]=(Sqrt[Tr[X[[3]].X[[3]]]/Tr[#.#]]#)&@(X[[8]]-Tr[X[[3]].X[[8]]]/Tr[X[[8]].X[[8]]] X[[3]]);
			IrrepBasis[X, {p, q}],
			X]
	];
	
	Options[SU2BasisMatrices]={BasisType->"AngularMomentum"};
	SU2BasisMatrices[p_, OptionsPattern[]]:=LieAlgebraBasisMatrices[Irrep[p], BasisType->OptionValue[BasisType]]
End[];
EndPackage[];
