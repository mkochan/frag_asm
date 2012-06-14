(* frag_asm, A Simple Fragment Assembly Program, ver. 1.00
   Martin Kochan (mkochan@ksu.edu), 27 Nov 2005
*)

module Similarity = struct

  exception Negative_similarity
  
  (* Performs gapped Smith-Waterman alignment. "Gapped" means that strings being aligned
     may contain gaps, which are then considered as wildcards.
     Returns a 5-tuple, let it be (ascore,ai,aj,ash,afg), containing:
      ascore ... best local similarity score obtained
      ai, aj ... the resulting strings
      ash ...... the shift between them (aj respective to ai)
      afg ...... afg = (afgi, afgj) the pair of lists of fresh gaps appearing in ai, aj,
                 respectively.
     Note that in ai, aj we do not return only the overlapping, high-scoring, region of
     local similarity (like "classical" S-W would) but we also pre- and post-pend the
     low-scoring regions.
  *) 

  let gwater a b =
    let lena = String.length a in
    let lenb = String.length b in
    let m = Array.make_matrix (lena+1) (lenb+1) 0 in
    let s x y = if x = y then (if (x = '-' && y = '-') then 1 else 1) else -1 in
    let q = 3 in (* prize for a gap *)
    let rec outer best i =
      let rec inner best' j =
	if j > lenb then best'
	else
	  let p = s a.[i-1] b.[j-1] in      
	  let v = max 0 (max (m.(i-1).(j) - q) (max (m.(i).(j-1) - q) (m.(i-1).(j-1) + p))) in
	  m.(i).(j) <- v;	
	  if (v > (fun (x,y,z)->z) best') 
	  then (inner (i,j,v) (j+1)) 
	  else (inner best' (j+1))
      in
      if i > lena then best else
      let new_best = inner best 1 in
      outer new_best (i+1)
    in
    let (best_i, best_j, best_score) = outer (0,0,0) 1 in (*side effect: fills in the matrix m*)
    if best_score <= 10 then
      raise Negative_similarity
    else
    (* returns pair of starting indexes, the pair of strings, and a flag telling whether the
       starting coordinates are to be determined in this call (ie. the empty string pair was
       reached in the upper instance *)
      let rec backtrack i j =
	if i = 0 || j = 0 then ((i,j), ("",""), true, ([], []))
	else
	  if m.(i).(j) = m.(i-1).(j-1) + (s a.[i-1] b.[j-1]) then
	    let (r1,r2,r3,r4) = backtrack (i-1) (j-1) in
	    ((if r3=true then (i,j) else r1),
	     (fst r2 ^ Char.escaped a.[i-1], snd r2 ^ Char.escaped b.[j-1]),
	     false,
	     r4
	    )
	  else
	    if m.(i).(j) = m.(i-1).(j) - q then 
	      let (r1,r2,r3,r4) = backtrack (i-1) j in
	      ((if r3=true then (i,j) else r1),
	       (fst r2 ^ Char.escaped a.[i-1], snd r2 ^ "-"),
	       false,
	       (fst r4, j :: (snd r4))
	      )
	    else
	      if m.(i).(j) = m.(i).(j-1) - q then
		let (r1,r2,r3,r4) = backtrack i (j-1) in 
		((if r3=true then (i,j) else r1), 
		 (fst r2 ^ "-", snd r2 ^ Char.escaped b.[j-1]), 
		 false,
		 (i :: (fst r4), snd r4)
		)
	      else
		((i,j), ("", ""), true, ([], []))
      in
      let (starts, align, flag, fresh_gaps) = backtrack best_i best_j in
      let (start_i, start_j, end_i, end_j) =
	((fst starts)-1, (snd starts)-1, best_i-1, best_j-1) in
      (
       best_score,
       (String.sub a 0 start_i) ^ (fst align) ^
       (String.sub a (end_i+1) (max 0 (lena - (end_i+1)))),
       (String.sub b 0 start_j) ^ (snd align) ^
       (String.sub b (end_j+1) (max 0 (lenb - (end_j+1)))),
       start_i-start_j,
       (List.rev (fst fresh_gaps), List.rev (snd fresh_gaps))
      )

end
