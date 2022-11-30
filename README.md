# 1858_cloudy
Computing Cloudy line ratios 

### Instructions 

To run the cloudy grid, run 

```
python3 1858_run_cloudy.py
```

you'll need to make sure that the kwarg cloudy_path is set correctly in the routine that runs cloudy. I haven't checked this recently because my cloudy install is failing right now. 

You can then make the plots using plot_ratios.py or PlotRatios.ipynb. 

I've also included the line ratio data I get from Cloudy in files .lines. 

NOTE: currently I've done this for one density and BB temperature, but you could clearly adapt these routines to be more flexible and run larger grids.

If you want to run in parallel, I can provide a script to do this. 

Finally, convert_ip.py is a useful script for converting xi to U (different ionization parameters) which is used by the 1858_run_cloudy.py routine. 