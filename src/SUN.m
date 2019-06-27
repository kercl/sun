(* ::Package:: *)

BeginPackage["SUN`"]

  SUNIrrepDynkin::usage = 
	  "Generate a irredcible representation based on a Dynkin label."

  Begin["Private`"]

  SUNIrrepDynkin[dynkinLabel__] :=
    Module[{},
      SUNIrrep /: 
        MakeBoxes[obj : SUNIrrep[asc_?SUNIrrepQ], 
         form : (StandardForm | TraditionalForm)] := Module[{above, below},
         above = {(*example grid*)
           {BoxForm`SummaryItem[{"SU("<>ToString[Length[x]+1]<>") Irrep: ", 
                                 "D("<>StringRiffle[ToString/@{dynkinLabel}, ", "]<>")"}], SpanFromLeft},
            BoxForm`SummaryItem[{"Dimension: ", -1}]
           };
         below = {(*example column*)
           
           };
         BoxForm`ArrangeSummaryBox[
          SUNIrrep,(*head*)
          <|{"empty"->None}|>,(*interpretation*)
          None,(*icon use None if not needed*)
          (*above and below must be in a format suitable for Grid or Column*)
          above,(*always shown content*)
          below,(*expandable content*)
          form, 
          "Interpretable" -> True]];
       SUNIrrep
    ];

  End[]
EndPackage[]
