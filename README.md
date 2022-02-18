# tagmmdi_wd
Paper specific scripts.

## Integrating LOPIT and GO data

To reproduce the comparison between the MDI and TAGM models using the data compared by Breckels et al., 2015, please run the following lines from the terminal, setting ``save_dir`` and ``script_dir`` to appropriate values. The values given below should work if ``tagmmdi_wd`` is the working directory.

```{bash, validationStudy}
script_dir="./Scripts/";
save_dir="./CV_output/";
n_seeds=10;
R=15000;
thin=50;
burn=5000;

for ((s=1; s<=${n_seeds}; s++));
do
  Rscript ${script_dir}mdiTagmCVSingleLoop.R --datasets 'dunkley2006 dunkley2006goCC' --R ${R} --thin ${thin} --burn ${burn} --seed ${s} --test_size 0.70 --K 100 --n_chains 5 --save_dir ${save_dir} --categorical_column_threshold 5
  
  Rscript ${script_dir}mdiTagmCVSingleLoop.R --datasets 'groen2014r1 groen2014r1goCC' --R ${R} --thin ${thin} --burn ${burn} --seed ${s} --test_size 0.70 --K 100 --n_chains 5 --save_dir ${save_dir} --categorical_column_threshold 5
  
  Rscript ${script_dir}mdiTagmCVSingleLoop.R --datasets 'E14TG2aS1 E14TG2aS1goCC' --R ${R} --thin ${thin} --burn ${burn} --seed ${s} --test_size 0.70 --K 100 --n_chains 5 --save_dir ${save_dir} --categorical_column_threshold 5
  
  Rscript ${script_dir}mdiTagmCVSingleLoop.R --datasets 'HEK293T2011 HEK293T2011goCC' --R ${R} --thin ${thin} --burn ${burn} --seed ${s} --test_size 0.70 --K 100 --n_chains 5 --save_dir ${save_dir} --categorical_column_threshold 5
  
  Rscript ${script_dir}mdiTagmCVSingleLoop.R --datasets 'tan2009r1 tan2009r1goCC' --R ${R} --thin ${thin} --burn ${burn} --seed ${s} --test_size 0.70 --K 100 --n_chains 5 --save_dir ${save_dir} --categorical_column_threshold 5
done;
```
