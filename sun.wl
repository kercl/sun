(* ::Package:: *)

GTPatternToRAlignedMatrix[pattern_]:=Module[{
	matrix = ConstantArray[0,{Length[pattern], Length[pattern]}],i
},
	For[i=1, i<=Length[pattern], i++,
		matrix[[i,i;;]] = pattern[[i]]
	];
	matrix
]


GTPatternToLAlignedMatrix[pattern_]:=Module[{
	matrix = ConstantArray[0,{Length[pattern], Length[pattern]}], i
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


bsearchMax[list_List, elem_]:=Module[{
	n0 = 1, n1 = Length[list],
	m
},
	While[n0 <= n1,
		m = Floor[(n0 + n1)/2];
		If[list[[m]] == elem, 
			While[list[[m]] == elem, m++]; 
			Return[m - 1]];
		If[list[[m]] < elem, n0 = m + 1, n1 = m - 1]
	];
	If[list[[m]] < elem, m, m - 1] 
];


bsearchMin[list_List, elem_]:=Module[{
	n0 = 1, n1 = Length[list], 
	m
},
	While[n0 <= n1,
		m = Floor[(n0 + n1)/2];
		If[list[[m]] == elem, 
			While[list[[m]] == elem, m--]; 
			Return[m + 1]];
		If[list[[m]] < elem, n0 = m + 1, n1 = m - 1]
	];
	If[list[[m]] > elem, m, If[m==Length[list],None,m + 1]]
]


bsearchRange[list_List, elem_]:=Module[{
	idxmin = bsearchMin[list, elem],
	idxmax = bsearchMax[list, elem]
},
	If[idxmin > Length[list] || idxmax > Length[list], Return[{}]];
	If[list[[idxmin]] != elem, Return[{}]];
	{idxmin, idxmax}
]


GTPatternIndexOfBsearch[basis_,pattern_]:=Module[{
	candidates = {1,Length[basis]},
	aFlattened = Flatten[pattern], i
},
	If[\[Not](basis[[1,;;Length[pattern] ]] == pattern[[1]]), Return[None] ];
	For[i=Length[pattern]+1, i<=(Length[pattern](Length[pattern]+1))/2,i++,
		candidates = bsearchRange[ basis[[candidates[[1]];;candidates[[2]], i]], aFlattened[[i]] ] + candidates[[1]] - 1;
		If[Length[candidates] == 0, Return[None]];
Print[candidates];
	];
	If[candidates[[1]] == candidates[[2]], Return[candidates[[2]]], Return[None]];
]


GTPatternIndexOf[basis_,pattern_]:=Module[{},
	GTPatternIndexOfFlattened[basis, Flatten[pattern]]
]


GTPatternIndexOfFlattened[basis_,pattern_]:=Module[{
	ret
},
	ret = Position[basis, pattern]//Flatten;
	If[Length[ret]==1, First[ret], None]
]


BuildSUNRaisingOp[basis_, l_]:=Module[{
	entries = {}, col, row, k, m, incIdx,
	N1,N2,D0
},
	For[col=1, col<=Length[basis], col++,
		For[k=1, k<=l, k++,
			m = basis[[col]];
			incIdx = IntegerPart[-l(l+1)/2] + k - 1;
			m[[incIdx]]++;
			row = GTPatternIndexOfFlattened[basis, m];
			If[row == None, Continue[]];
			m[[incIdx]]--;
			m = GTPatternRecoverFlattened[m]//GTPatternToLAlignedMatrix//Reverse;
			
			N1 = Times@@(m[[l+1,;;l+1]] - m[[l,k]] + k - Range[l+1]);
			N2 = If[l==1,1,Times@@(m[[l-1,;;l-1]] - m[[l,k]] + k - 1 - Range[l-1])];
			D0 = (m[[l,;;l]] - m[[l,k]] + k - Range[l])(m[[l,;;l]] - m[[l,k]] + k - 1 - Range[l]);
			D0 = (Times@@D0[[;;k-1]])(Times@@D0[[k+1;;]]);
			AppendTo[entries,{col,row}->-((N1 N2)/D0)];
		];
	];
	SparseArray[entries, {Length[basis],Length[basis]}]
]


BuildSUNRaisingOp[basis_, l_]:=Module[{
	entries = {}, col, row, k, m, incIdx,
	N1,N2,D0
},
	For[col=1, col<=Length[basis], col++,
		For[k=1, k<=l, k++,
			m = basis[[col]];
			incIdx = IntegerPart[-l(l+1)/2] + k - 1;
			m[[incIdx]]++;
			row = GTPatternIndexOfFlattened[basis, m];
			If[row == None, Continue[]];
			m[[incIdx]]--;
			m = GTPatternRecoverFlattened[m]//GTPatternToLAlignedMatrix//Reverse;
			
			N1 = Times@@(m[[l+1,;;l+1]] - m[[l,k]] + k - Range[l+1]);
			N2 = If[l==1,1,Times@@(m[[l-1,;;l-1]] - m[[l,k]] + k - 1 - Range[l-1])];
			D0 = (m[[l,;;l]] - m[[l,k]] + k - Range[l])(m[[l,;;l]] - m[[l,k]] + k - 1 - Range[l]);
			D0 = (Times@@D0[[;;k-1]])(Times@@D0[[k+1;;]]);
			AppendTo[entries,{col,row}->Sqrt[-((N1 N2)/D0)]//Simplify];
		];
	];
	SparseArray[entries, {Length[basis],Length[basis]}]
]


SUNLieAlgebraIrep[rep_]:=Module[{
	carrierSpaceBasis = GTPatternBasisFromFixedTopRow[Append[rep,0]],
	N, X, H, \[Sigma], i, diagonals, UMatrices
},
	N = Length[rep]+1;
	If[rep==ConstantArray[0,N-1],Return@ConstantArray[{{0}},N^2-1]];
	\[Sigma] = Append[#,0]&/@(Total[#,{2}]&/@GTPatternRecoverFlattened/@carrierSpaceBasis);
	diagonals = \[Sigma][[;;,2;;N]] - (1/2) * (\[Sigma][[;;,1;;N-1]] + \[Sigma][[;;,3;;N+1]]) // Transpose;
	H = (SparseArray[
		Table[{i,i}->-#[[i]],{i,1,Length[carrierSpaceBasis]}],
		{Length[carrierSpaceBasis],Length[carrierSpaceBasis]}]&)/@diagonals;
	(*H = If[normalize \[NotEqual] None, Orthogonalize[Normal/@H,normalize]//Simplify, None];*)
	(*H = SparseArray/@(Orthogonalize[Normal/@H, 2 Tr[#1.#2]&]//Simplify);*)
	UMatrices = Table[BuildSUNRaisingOp[carrierSpaceBasis, l], {l,N-1}];
	UMatrices = Join[UMatrices, #2.#1-#1.#2&@@#&/@Subsets[UMatrices,{2}] ];
	X = Join[1/2 (#+Transpose[#])&/@UMatrices, 1/(2I) (#-Transpose[#])&/@UMatrices];
	Join[X,H]
]


StandardIrrepNotation[n_]:=Reverse@Accumulate@Reverse@n


SU3[p_,q_]:=Module[{
	X=(SUNLieAlgebraIrep@StandardIrrepNotation@{p,q})[[{2,5,7,3,6,1,4,8}]]
},
	If[p!=0||q!=0,
		X[[8]]=(Sqrt[Tr[X[[3]].X[[3]]]/Tr[#.#]]#)&@(X[[8]]-Tr[X[[3]].X[[8]]]/Tr[X[[8]].X[[8]]] X[[3]])
	];
	Association[Table[i->X[[i]],{i,Length[X]}]]
]


SU2[p_]:=SUNLieAlgebraIrep@{p}
