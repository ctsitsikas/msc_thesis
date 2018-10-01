Scripts
PDF: function that computes Gumbel PDF 

allowances_with_waves: allowances computation including wave change

subtract_distributions: subtracts end century-historical waves

comparison_PDF_at_locations: plots GTSR and GESLA2 PDFs at selected locations

wave_data_gumbel_models: makes Gumbel fits to seperate models, computes the ensemble mean and makes gumbel fit to it

gtsr_gesla2_comparison: comparison between datasets in terms of gumbel parameters and PDFs

nearest_points_computation: computes nearest points between two datasets with a distance limit

slr_gtsr_waves_at_ipcc_locations: makes plot of Sea Level Change PDF with all allowances at selected locations

percent_of_events_freq_change: computes percentage of points affected by specific change in exceedance frequency

Europe_allow_plus_waves_checks: analyzes allowances including waves in Europe, meridional trend and contributions by mean+std of wave change

downloadend, downloadhist: computes and downloads significant wave height maxima projections from the CSIRO server for mid-century and end-century


Use

'gtsr_gesla2_comparison' gives plots comparing the Gumbel parameters of GTSR and GESLA 2 at their nearest points, like difference, ratio, rmse. Then, 'comparison_PDF_at_locations' plots the Gumbel pdfs at nine specific locations to compare.
Using 'nearest_points_computation' we compute the nearest points between GTSR and sea level change grid. Then, we use the output in the 'allow_normal_GTSR' script to compute the GTSR allowances. The nearest point computation for these specific datasets can also be done within the 'allow_normal_GTSR' script. The same process stands for reproducing GESLA2 allowances, changing the input to the gesla mat file.
To include waves, 'downloadend' and 'downloadhist' access the CSIRO server and download annual maxima of significant wave height by 8 models for two periods, one historical and one for the end of century. For this data, 'wave_data_gumbel_models' makes Gumbel fits for individual models as well as for ensemble mean. Using 'combine_subtract_distributions' the ensemble end century and historical distributions are subtracted and that gives a wave maxima change distribution (P2). Nearest points between GTSR and waves grid can be computed with 'nearest_points_computation'. Then the output is plugged in 'allowances_with_waves' and the new set of allowances is computed. 'Europe_allow_plus_waves_checks' computes and plots the contribution by wave change mean and variance.
For the plot of the nine selected locations with SLC and all allowances 'slr_gtsr_waves_at_ipcc_locations' is used. 
'plottingallowancescontinents' is used to make the meridional plots of allowances for every continent in both sets of allowances(changing the dataset input and titles). 'plottingallowances_waves' and 'plottingallowanceseurope' are used to plot global and european maps with allowances. In case of Europe, the difference between GTSR and GTSR+wave allowances is computed and plotted within the script.
'percent_of_events_freq_change' is used to compute the percentage of points affected by a change in exceedance frequency. The input set of allowances and frequency change can be adjusted.


Note
Scripts run using MATLAB, it is more convenient to run gradually section by section with Ctrl+Enter. Paths of data files have to be altered suitably. Allowances computation additionally requires compiling a Fortran code. See Aimee Slangen's folder.
csirolib is a matlab library with map plotting tools.

Christos Tsitsikas 
email: ctsitsikas@hotmail.com