# GBRseabirds

This project models and predicts the foraging niches of tropical seabirds across the globe and their associated foraging ranges. It performs analyses for the accompanying paper: **Refining seabird marine protected areas by predicting habitat inside foraging range - a case study from the global tropics**, authored by Mark Miller, Graham Hemson, Maria Beger, Bradley Congdon and others.

### core scripts

**tracking_compile.R** Compiles GPS tracking data contributed for 60 colonies across the globe.<br />
**tracking_preprocess.R** Cleans, standardises and prepares tracking data for analyses.<br />
**identify_foraging.R** Uses the Expectation Maximisation binary Clustering (EMbC) algorthm to assign foraging behaviour to tracking datapoints.<br />
**oceano_preprocess.R** Compiles oceanographic raster data, carries out kernel analyses on tracking dataset foraging datapoints and then extracts oceanogrpahic data for presence and pseudo-absence locations.<br />
**modelling_kernUD_tune.R** Preliminary Environmental Niche Modelling running random forests with varying hyperparameters setups to identify optimal hyperparameters to tune final models.<br />
**modelling_kernUD.R** Runs Environmental Niche Modelling with random forests using tuned hyperparameters. Trains colony-specific and multi-colony models then predicts to other colonies and performs model validation. Also predicts multi-colony models to each tracked colony to compare with known foraging areas (50%UDs) for refinement framework validation.<br />
**paper_code.R** # Code takes outputs from random forest ENMs and transferability analyses and:calculates summary foraging radius values; conducts global Framework validation analyses; conducts foraging radius refinement on the GBR; calculates summary tables for paper generates non-GIS paper figures.<br />
