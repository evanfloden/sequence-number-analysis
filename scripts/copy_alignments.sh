for i in ../../../sequence_sets/*.fa; 
    do BASE=$(basename $i); arBASE=(${BASE//./ }); 
       cp ../../../../../results/alignments/${arBASE[0]}.${arBASE[1]}.${arBASE[2]}.dpa.1000.CLUSTALO_DPA.with.CLUSTALO.tree.aln ${arBASE[0]}.${arBASE[1]}.${arBASE[2]}.regressive.1000.CLUSTALO.with.CLUSTALO.tree.aln; 
    done
