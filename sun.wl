(* ::Package:: *)

Needs["Combinatorica`"]


GTPatternToRAlignedMatrix[pattern_]:=Module[{
	matrix = ConstantArray[0,{Length[pattern], Length[pattern]}]
},
	For[i=1, i<=Length[pattern], i++,
		matrix[[i,i;;]] = pattern[[i]]
	];
	matrix
]


GTPatternToLAlignedMatrix[pattern_]:=Module[{
	matrix = ConstantArray[0,{Length[pattern], Length[pattern]}]
},
	For[i=1, i<=Length[pattern], i++,
		matrix[[i,;;Length[pattern]-i+1]] = pattern[[i]]
	];
	matrix
]


GTPatternRecoverFlattened[pattern_]:=Module[{
	N = (-1 + Sqrt[1 + 8 Length[pattern]]) / 2,
	indices = None
},
	If[\[Not](N \[Element] Integers), Return[{}]];
	indices = Prepend[Accumulate[Range[N,1,-1]], 0];
	Return[Table[pattern[[indices[[i]]+1;;indices[[i+1]]]], {i,1,N}]]
]


GTPatternMinimumFromTopRow[topRow_]:=
	Table[topRow[[i;;]], {i,Length[topRow]}]


GTPatternMaximumFromTopRow[topRow_]:=
	Table[topRow[[;;-i]], {i,Length[topRow]}]


IncrementArrayWithLimit[array_,minArray_,maxArray_]:=Module[{
	lowPoints = Position[Thread[maxArray>array],True]//Flatten,
	arrayRet = array
},
	If[Length[lowPoints] == 0, Return[{minArray,True}]];
	arrayRet[[ lowPoints[[1]] ]]++;
	arrayRet[[;;lowPoints[[1]]-1]] = arrayRet[[lowPoints[[1]]]];
	{arrayRet,False}
]


IncrementPatternTopFixed[pattern_]:=Module[{
	carry = True,
	colFromRight = 1,
	newCol = None,
	arrayRAligned = GTPatternToRAlignedMatrix[pattern]//Reverse,
	minArrayRAligned = GTPatternToRAlignedMatrix[GTPatternMinimumFromTopRow[ pattern[[1]] ]]//Reverse
},
	While[carry,
		If[colFromRight >= Length[pattern], Return[{}]];
		
		newCol = IncrementArrayWithLimit[
			arrayRAligned[[colFromRight;;Length[pattern]-1,-colFromRight]],
			minArrayRAligned[[colFromRight;;Length[pattern]-1,-colFromRight]],
			arrayRAligned[[colFromRight+1;;,-colFromRight-1]]
		];
		carry = newCol[[2]];
		arrayRAligned[[colFromRight;;Length[pattern]-1,-colFromRight]] = newCol[[1]];
		If[carry, arrayRAligned[[;;,-colFromRight]] = minArrayRAligned[[;;,-colFromRight]] ];
		colFromRight++;
	];
	Table[arrayRAligned[[i,-i;;]],{i,Length[pattern]} ]//Reverse
]


GTPatternBasisFromFixedTopRow[toprow_]:=Module[{
	currentPattern = GTPatternMinimumFromTopRow[toprow],
	basis = {}
},
	While[Length[currentPattern] > 0,
		AppendTo[basis, Flatten[currentPattern]];
		currentPattern = IncrementPatternTopFixed[currentPattern]
	];
	basis
]


BinaryRangeSearch[l_,k_]:=Module[{
	pos = BinarySearch[l, k],
	min = None,
	max = None
},
	If[\[Not](pos \[Element] Integers),Return[{}]];
	min = pos;
	max = pos;
	
	For[max=pos,max <= Length[l], max++,
			If[l[[max]] != k, Break[]]
		];
	For[min=pos,min >= 1, min--,
			If[l[[min]] != k, Break[]]
		];
	{min+1,max-1}
]


GTPatternIndexOf[basis_,pattern_]:=Module[{
	candidates = {1,Length[basis]},
	aFlattened = Flatten[pattern]
},
	If[\[Not](basis[[1,;;Length[pattern] ]] == pattern[[1]]), Return[None] ];
	For[i=Length[pattern]+1, i<=(Length[pattern](Length[pattern]+1))/2,i++,
		candidates = BinaryRangeSearch[ basis[[candidates[[1]];;candidates[[2]], i]], aFlattened[[i]] ] + candidates[[1]] - 1;
		If[Length[candidates] == 0, Return[None]];
	];
	If[candidates[[1]] == candidates[[2]], Return[candidates[[2]]], Return[None]];
]


BuildSUNRaisingOp[basis_]:=Module[{

},
	None
]


BuildSUNAlgebra[rep_]:=Module[{
	carrierSpaceBasis = GTPatternBasisFromFixedTopRow[Append[rep,0]],
	N = 0,
	X = {},
	H = {},
	\[Sigma] = {},
	diagonals = {}
},
	N = Length[rep]+1;
	\[Sigma] = Append[#,0]&/@(Total[#,{2}]&/@GTPatternRecoverFlattened/@carrierSpaceBasis);
	diagonals = \[Sigma][[;;,2;;N]] - (1/2) * (\[Sigma][[;;,1;;N-1]] + \[Sigma][[;;,3;;N+1]]) // Transpose;
	H = (SparseArray[
		Table[{i,i}->#[[i]],{i,1,Length[carrierSpaceBasis]}],
		{Length[carrierSpaceBasis],Length[carrierSpaceBasis]}]&)/@diagonals;
	H
]
