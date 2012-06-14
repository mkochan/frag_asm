all:
	ocamlc -o frag_asm similarity.ml frag_asm.ml

clean:
	rm -f *.cm*

