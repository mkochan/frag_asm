(* frag_asm, A Simple Fragment Assembly Program, ver. 1.00
   Martin Kochan (mkochan@ksu.edu), 27 Nov 2005
   
   Assembles DNA sequences into contigs. Contigs are shown as layouts of mutually aligned   
   sequences. It is also possible to produce a file showing process of progressive merging
   of contigs.

   READ frag_asm.pdf BEFORE READING THIS FILE.
   HOW TO READ THIS FILE: FIRST TYPE DEFINITIONS, AND THEN FROM THE BOTTOM TO THE TOP.

*)
  
open Similarity;; (* contains gapped Smith-Waterman algorithm *) 
  
(************************************ Types definitions ****************************************) 

type layout = layout_row list      
and layout_row = {i: int;          (* i = identifier of read *)
                  shift: int;      (* shift = with respect to the reference row *)
		  chs: char list   (* chs = characters of the read (possibly with gaps added) *)
		 } 
and contig = int * layout          (* contig's id, it's layout *)
and edge = int * int               (* (i,j) 
				      edge in the similarity graph - this is acquired
				      during the proprocessing (using 14-mers) *)
and occur_of_14mer = int * int     (* (id, start)
				      describes occurence of a 14-mer in read with identifier id
				      and starting at position start *)
and band = int * int * int;;       (* (starti, startj, length)
				      describes region of 100% similarity between two sequences;
				      this region starts at starti in the first sequence,
				      at startj in the second, and has length length same in
				      both sequences *)

(********************************** Auxiliary functions ****************************************) 

let chars_of_string s =
  let len = String.length s in
  let rec f i acc = if i < 0 then acc else (f (i-1) (s.[i]::acc)) in f (len-1) [];;
    
let string_of_chars l = 
  let rec f i l s =
    match l with
       [] -> s
     | x::xs -> s.[i] <- x; f (i+1) xs s 
  in f 0 l (String.create (List.length l));;

(* For item of hash h keyed by k, which is a list, prepend value v to it (to the list). *)  
let prepend h k v =
  let rest = try Hashtbl.find h k with Not_found -> [] in
  Hashtbl.replace h k (v::rest);;
    
(* Remove repeating numbers in list l; p is predicate for determining equality  *)
let rec rem_stutters p l = match l with
   [] -> []
   | x::[] -> [x]
   | x::xx::xxs -> if p x xx then (rem_stutters p (xx::xxs)) else x::(rem_stutters p (xx::xxs));;
								       
(* Put the item of list l that satisfies predicat p into the head. Throw Failure(msg) if no
   such item is found.*)
let reshuffle p l =
  let (iteml, without_item) = 
    List.fold_right
	(fun x ans -> if (p x) then ([x], snd ans) else (fst ans, x::(snd ans))) l ([], []) in
  try 
    let item = List.hd iteml in item :: without_item
  with Failure(msg) -> failwith "reshuffle";;
      
(* Convert list of insert positions to absolute positions. "Insert positions" are indexes in a
   string where gaps, "-"s, should be put. "Absolute positions" designate indexes in the
   resulting text where those gaps should be placed. *)
let inserts_to_positions l =
  let rec f count done_as_far rest = match rest with
     [] -> done_as_far
   | x::xs -> f (count+1) ((x+count)::done_as_far) xs
  in List.rev (f 0 [] l);;

(******************************* Input & creation of singlets ***********************************) 

(* Return complementary strand (in list-of-chars format) *)
let comp_strand chs = List.map
    (function 'a' ->
      't' | 't' -> 'a' | 'c' -> 'g' | 'g' -> 'c' | _ -> raise (Invalid_argument "comp_dna") )
    (List.rev chs);;
    
(* Read a file and produce array of strings (ie. "reads") *)
let rec readfile filename = 
  let buff = Scanf.Scanning.from_file filename in
  let rec f i acc = 
    try
      let s = String.lowercase (Scanf.bscanf buff " %s" (fun x -> x)) in
      if s = "" then acc else f (i+1) (s::acc)
    with
	_ -> acc
  in Array.of_list (List.rev (f 0 []));; (* get reads from the file *)

(* Create a singlet, ie. a single-item conting, directly from input string array *)
let singlet i s = [ {i = i; shift = 0; chs = chars_of_string s} ];;

(* Create an array of singlets (these will be later gradually merged into contigs) *) 
let init_contig_array sa =
  Array.init (Array.length sa) (fun i -> (i, singlet i sa.(i)));;

(******************************************** Output *******************************************) 

(* contig -> string conversion *)
let string_of_contig (c : contig) =
  let (c_id, c_lay) = c in
  let find_min_shift = 
    List.fold_right (fun x ans -> if x.shift < ans then x.shift else ans) c_lay 10000 in
  "Contig #" ^ (string_of_int c_id) ^ ":\n" ^
  List.fold_right 
    (fun x ans -> 
      (String.make (x.shift-find_min_shift) ' ') ^ (string_of_chars x.chs) ^
      " (" ^ (string_of_int x.i) ^ ") \n" ^ ans
    ) c_lay "";;
    
let output_contigs out_stream (cs : contig list) =
  output_string
    out_stream (List.fold_right (fun x ans -> (string_of_contig x) ^ "\n" ^ ans) cs "");;
    
let output_progress progress_stream lay_i lay_j lay_res =
  output_string progress_stream (
  (string_of_contig lay_i) ^ "+++++ \n" ^
  (string_of_contig lay_j) ^ "====> \n" ^
  (string_of_contig lay_res) ^ "----------------\n\n"
 );;

(************** Preproccessing (heuristic), generation of edges in similarity graph ************) 

(* Produce hash of 14-mers from string array. The type of this hash is:
      (string, occur_of_14mer list) Hashtbl.t
   That means that upon asking for some 14-character string (a 14-mer), you will be given
   a list occurences of that 14-mer accross all strings from sa.
   Check the type occur_of_14mer to see how such occurence is represented.
*)
let get_14mers (sa:string array) =
  let totallen = Array.fold_right (fun x ans -> String.length x + ans) sa 0 in
  let h: (string, occur_of_14mer list) Hashtbl.t = Hashtbl.create totallen in
  Array.iteri 
    (fun n s ->
      let len = String.length s in
      if len <= 14 then ()
      else
	let rec loop i = 
	  if i > len-14 then () 
	  else (prepend h (String.sub s i 14) (n, i); loop (i+1)) 
	in
	loop 0;
    )
    sa; 
  h;;
    
(* Checks the hash of 14-mers to see if any 14-mers occur in multiple reads. That implies
   that these reads have places of high local similarity there. Report those pairs of reads
   (such pairs actually form edges in the "similarity graph").
   The output is yet another hash, it's keys being of type edge, and keyed values being lists
   of uncompressed bands. (See type band).
*)   
let edges_from_14mers (mers: (string, occur_of_14mer list) Hashtbl.t) (sa : string array) =
  let edges: (edge, band list) Hashtbl.t = Hashtbl.create (Array.length sa) in
  Hashtbl.iter 
    (fun key occs -> 
      let comp_by_fst = (fun (a,b) (c,d) -> compare a c) in
      let occs' = List.sort comp_by_fst occs in
      let occs'' = rem_stutters (fun a b -> comp_by_fst a b = 0) occs' in
      (fun folding_result -> ())
	(List.fold_right 
	   (fun (xn,xi) ans -> 
	     List.iter (fun (yn,yi) -> prepend edges (xn,yn) (xi,yi,14)) ans; (xn,xi)::ans
	   )
	   occs'' []
	)
    ) 
    mers;
  edges;;

(* Compress the bands in the hash edges.
   Bands that are shifted just by one character and otherwise overlapping, will get merged.
*)
let compress_edges edges: (edge, band list) Hashtbl.t =
  (* Compresses the 14-mers if there are some following each other; the output are bands that
     have length >= 14 chars. Constraint: l has to be ordered by xi1 where (xi1,xi2,xlen) is 
     item of l *)   
  let rec compress_bands l =
    match l with 
       [] -> []
     | x::[] -> [x]
     | x::xx::xxs -> 
	 let (xi1, xi2, xlen) = x in
	 let (xxi1, xxi2, xxlen) = xx in
	 if (xxi1 = (xi1+xlen-14)+1) then
	   compress_bands ((xi1, xi2, xlen+1)::xxs)
	 else
           x::(compress_bands (xx::xxs))
  in
  let edges' = Hashtbl.create (Hashtbl.length edges) in
  Hashtbl.iter 
    (fun key bs -> 
      Hashtbl.add edges' key (
      compress_bands (List.sort (fun (a,b,c) (d,e,f) -> compare a d) bs)
     )
    )
    edges;
  edges';;

(************************************** Mergining layouts **************************************)

(* Add gaps into layout-row r at positions in gl; positions in gl are given relative to
   the layout-common frame of reference. Attempts to add gaps before the row's beginning result
   in the row being shifted more to the right. *)
let add_gaps gs (r:layout_row) =
  let rec f1 gs r =
    match gs with
       [] -> ([], r)
     | gsx::gss -> 
         if gsx < r.shift then f1 gss {i = r.i; shift = r.shift+1; chs = r.chs} else (gs, r) 
  in let rec f2 i gs chs = match chs with
     [] -> []
   | ax::axs -> match gs with
	[] -> chs
      |	gx::gxs -> 
	  if i=gx then '-'::(f2 (i+1) gxs chs) else ax::(f2 (i+1) gs axs)
  in let (gs', r') = f1 gs r in
  {i = r'.i; shift = r'.shift; chs = (f2 r'.shift gs' r'.chs)};;
    
(* Shift the layout row r to the right by sh positions. *)
let shift_row sh r = {i = r.i; shift = r.shift + sh; chs = r.chs};;
    
(* Merges layouts - THIS IS *THE* MOST IMPORTANT FUNCTION
   Aligns the layouts li and lj according to their head rows (head elements of li and lj).
   It runs gapped Smith-Waterman produce pairwise alignment between the heads and changes the
   rest so that both layouts can be merged. See the file frag_asm.pdf to learn more.
*)
let merge_layouts (li : layout) (lj : layout) =
  let hi = List.hd li in
  let hj = List.hd lj in
  let (ascore,ai,aj,ash,afg') =
    Similarity.gwater (string_of_chars hi.chs) (string_of_chars hj.chs) in
  let (afgi, afgj) = (
    List.map ((+) hi.shift) (inserts_to_positions (fst afg')),
    List.map ((+) hj.shift) (inserts_to_positions (snd afg'))) in
  let li' = List.map (add_gaps afgi) li in
  let lj' = List.map (add_gaps afgj) lj in
  let lj'' = List.map (shift_row (hi.shift + ash - hj.shift)) lj' in
  (li' @ lj'' : layout);;

(************************************** Storing contigs *****************************************)
(*
   NOTE: The contigs are stored in a reference array, ca.
   It's structural invariant: When i-th read is present in a contig ci, then ca[i] = ci.
   Thus actually all positions in ca that correspond to ci's reads contain reference to ci.*)

(* Updates reference array (of contigs), ca, so that all contigs with id's equal to
   id_i or or id_j "become" contig (id_i, lay').
   That effectively merges the original contigs. *)
let update_ca id_i id_j (lay' : layout) (ca : contig array) =
  Array.iteri
    (fun ix c ->
      let c_id = fst c in if (c_id = id_i || c_id = id_j) then ca.(ix) <- (id_i, lay') else ()
    ) ca;
  ca;;
    
(* Find out individual (unique) contigs from the reference array of contigs.*)
let get_unique_contigs_list (ca : contig array) =
  let sorted_contigs = List.sort (fun cx cy -> compare (fst cx) (fst cy)) (Array.to_list ca) in
  rem_stutters (fun cx cy -> (fst cx) = (fst cy)) sorted_contigs;;

(**************************************** Central loop ******************************************) 

(* The fragment assembly.
   Take array of reads, list of proproccessed edges and a stream to scribble the partial results
   to. Return list of contigs. See the file frag_asm.tex to see what gets done.
*)   
let assemble sa (edges: (edge, band list) Hashtbl.t) progress_stream =
  let total_band_len bs = List.fold_right (fun x ans -> ((fun (a,b,c) -> c) x) + ans) bs 0 in
  let sorted_edges =
    List.map 
      fst
      (List.sort
	 (fun (a,b) (c,d) -> compare d b)
	 (Hashtbl.fold (fun key v ans -> (key, total_band_len v) :: ans) edges [])
      ) in
  let (ca : contig array) = init_contig_array sa in
  let ca' =
    List.fold_left 
      (fun ca e ->
        let (i,j) = e in
        let (id_i, lay_i_unshuffled) = ca.(i) in
        let (id_j, lay_j_unshuffled) = ca.(j) in
        if id_i = id_j then ca
        else 
          let lay_i = reshuffle (fun x -> x.i=i) lay_i_unshuffled in
          let lay_j = reshuffle (fun x -> x.i=j) lay_j_unshuffled in
          try
            let l' = merge_layouts lay_i lay_j in
            output_progress progress_stream (id_i,lay_i) (id_j,lay_j) (id_i,l'); 
            update_ca id_i id_j l' ca
          with
            Similarity.Negative_similarity -> ca
      )
      ca sorted_edges in
  get_unique_contigs_list ca';;

(************************************** The main() function *************************************) 
    
let main () =
  if  (Array.length Sys.argv) = 2 && Sys.argv.(1) = "--help"
  then
    print_string
      ("Fragment assembly demonstration program\n"
       ^ "Usage: frag_asm <in_file> <out_file> <progress_file>")
  else
    if (Array.length Sys.argv) <> 4
    then print_string "Incorrect number of parameters (type frag_asm --help).\n"
    else
      let in_filename = Sys.argv.(1) in
      let out_filename = Sys.argv.(2) in
      let progress_filename = Sys.argv.(3) in
      let sa = readfile in_filename in
      (*let sa_len = Array.length sa in*)
      let mers14 = get_14mers sa in
      let edges = compress_edges (edges_from_14mers mers14 sa) in
      let out_stream = open_out out_filename in
      let progress_stream = open_out progress_filename in
      output_string out_stream ("THERE ARE " ^ (string_of_int (Array.length sa)) ^ " READS.\n\n");
      output_contigs out_stream (assemble sa edges progress_stream); (* <--- THE ASSEMBLY! *)
      close_out out_stream;
      close_out progress_stream;
      exit 0;;
      
main ();;
