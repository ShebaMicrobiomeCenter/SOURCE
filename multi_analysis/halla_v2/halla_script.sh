# single simple halla run
halla \
      -x x.txt \
      -y y.txt \
      -o halla_res \
      -m spearman \
      --fdr_alpha 0.25 &
	  

## running all halla
for f in *; 
do
  halla \
  -x $f'/x.txt' \
  -y $f'/y.txt' \
  --x_dataset_label ${f%_vs_*} \
  --y_dataset_label ${f##*_vs_} \
  -o $f'/halla_res_alpha0.25' \
  -m spearman \
  --fdr_alpha 0.25
done

# copy ony the hallograms to a new destination, for easier reading
for f in full_pairs/*/halla_res_alpha0.25/hallagram.pdf; do
    f2=${f%/halla_res*}
	f2=${f2##*full_pairs/}
	cp -v "$f" res_plots_q25/"${f2//\//_}"
done

mkdir 
# copy ony the hallograms to a new destination, for easier reading
for f in full_pairs_uMets/*/halla_res_alpha0.25/hallagram.pdf; do
    f2=${f%/halla_res*}
	f2=${f2##*full_pairs/}
	cp -v "$f" res_plots_q25_uMets/"${f2//\//_}"
done

