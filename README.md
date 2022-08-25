# tagmmdi_wd
Paper specific scripts.

## Simulation study

Generate the data by calling:

```{bash generateSims}
Rscript ./Scripts/generateSimulationDataV3.R --save_dir "./Simulations/DataPhiModel/" --plot_data TRUE;
```

Analyse this data by calling

```{bash simStudy}
scenarios=("Gaussian" "MVT" "LogPoisson")
n_sims=100;

for scn in "${scenarios[@]}"
do
  for(( ii=1; ii<=${n_sims}; ii++))
  do
    # The case with K known
    Rscript ./Scripts/analyseSimulation.R --R 15000 --thin 100 --burn 7500 --index ${ii} --n_chains 10 --K 6 --n_clust_unsupervised 50 --scn ${scn} --save_dir ./Simulations/OutputPhiModel/ --data_dir ./Simulations/DataPhiModel/;
    
    # Using an overfitted mixture
    Rscript ./Scripts/analyseSimulation.R --R 15000 --thin 100 --burn 7500 --index ${ii} --n_chains 10 --K 60 --n_clust_unsupervised 50 --scn ${scn} --save_dir ./Simulations/OutputPhiModel/ --data_dir ./Simulations/DataPhiModel/;
  done
done
```

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

## Toxoplasma gondii

First process the input data

```{bash processTgondiiInputs}
Rscript tGondiiMDIInputs.R --K 125 --save_dir "./T_gondii/Output/" --data_dir "./T_gondii/Original_data/";
```

Then call MDI using

```{bash tGondiiCC}
# The long chains
n_chains=9;
for(( ii=1; ii<=${n_chains}; ii++))
do
  Rscript ./Scripts/callMDITGondiiModelOnly.R --data_file ./T_gondii/Data/TGondiiMDI_K_125_input.rds --R 45000 --thin 250 --seed ${ii} --K 125 --n_chains 1 --save_dir ./T_gondii/Output/;
done;

# Consensus clustering
n_chains=50;
for(( ii=1; ii<=${n_chains}; ii++))
do
  Rscript ./Scripts/callMDITGondiiModelOnly.R --data_file ./T_gondii/Data/TGondiiMDI_K_125_input.rds --R 12000 --thin 1000 --seed ${ii} --K 125 --n_chains 1 --save_dir ./T_gondii/ConsensusClustering/;
done;
```

Then to process the output of the consensus clustering

```{r processCC}
Rscript ./Scripts/processTGondiiCCoutput.R --data_dir ./T_gondii/Data/ --R 12000 --K 125 --save_dir ./T_gondii/ConsensusClustering/ --model_output_dir ./T_gondii/ConsensusClustering/
```
