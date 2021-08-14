# TDMTW
A full moment tensor solver

Focal mechanisms is computed using the full-waveform time-domain moment tensor (TDMT) technique (Dreger and Helmberger 1993), which implemented at the Geological Survey of Israel (Wetzler et al., 2019). Data are extracted by a magnitude dependent time windows starting 80 s before the event origin time, corrected for instrument response and the horizontal components are rotated to great circle path. Green’s functions are computed using the frequency-wavenumber integration code (FKPROG) of Saikia (1994) based on the Israel velocity model of Gitterman et al. (2002) with a 10 Hz sampling rate. This velocity model can be updated and recalculated using your own velocity model. See below for more details.

Green’s functions are band-pass filtered in the frequency band of 0.06–0.1 Hz for earthquake with magnitude 3.0≤MW≤4.0 (Fig. S4), and 0.5–1.0 Hz for MW≤2.9 in order to capture the high frequency seismic energy content of smaller magnitude earthquakes. This can be modified manually using -fmin and -fmax on the command line.

The best result is achieved through a grid-search on the depth and choosing the moment tensor solution and centroid depth for which the variance reduction (VR) is at maximum.  

Dreger, D.S., Helmberger, D. V., 1993. Determination of source parameters at regional distances with three- component sparse network data. J. Geophys. Res. 98, 8107–8125. https://doi.org/10.1029/93JB00023

Gitterman, Y., Pinsky, V., Shapira, A., Ergin, M., Kalafat, D., Gurbuz, G., Solomi, K., 2002. Improvement in detection, location, and identification of small events through joint data analysis by seismic stations in the Middle East/Eastern Mediterranean region 13. https://doi.org/DTRA01-00-C-0119

Saikia, C.K., 1994. Modified frequency‐wavenumber algorithm for regional seismograms using Filon’s quadrature: modelling of Lg waves in eastern North America. Geophys. J. Int. 118, 142–158. https://doi.org/10.1111/j.1365-246X.1994.tb04680.x

Wetzler, N., Shalev, E., Göbel, T., Amelung, F., Kurzon, I., Lyakhovsky, V., Brodsky, E.E., 2019. Earthquake swarms triggered by groundwater extraction near the Dead Sea Fault. Geophys. Res. Lett. 46, 2019GL083491. https://doi.org/10.1029/2019GL083491



# Run
 As a test case I chose one of the M4 along the main Dead Sea Transform - a left lateral strike-slip fault system.
 Use this line on your command line from your /MT folder:
 
 python3.6 tdmtw.py -ot 2021-06-15T23:08:54 -lat 30.099 -lon 35.178 -d 21 -m 4.1 -c GFZ
 
 -ot: Origine time
 
 -lat: Latitude of the earthquake
 
 -lon: Longitude of the earthquake
 
 -d: Depth
 
 -m: Magnitude
 
 -c: the FDSN server
 
 For more information about this event:
 https://earthquake.co.il/en/earthquake/feltInfo.php?ID=202106152307


![Main](https://user-images.githubusercontent.com/88764899/129444678-5f9478a5-4dad-4169-b254-eb3101704fe5.png)

# Solution
Two solutions are 
![Figure_2011](https://user-images.githubusercontent.com/88764899/129444688-54977b42-5c54-4845-a90a-ad0273cb503a.png)

# Green's Functions
You can make your own GF by running GF/Py/mk_Green_functions.py

Velocity model is used from GF/models/ .Use Gitt02.csv as a tamplate

Or you can download the Israel GF from my dropbox from https://www.dropbox.com/sh/i4v38eukn001sbj/AABj35JHfR-7MSk9GfOrD9Xma?dl=0

Documantation is still in progress... More shold be added in the near future

Please contact me directly if you have any questions: nadav.wetzler@gmail.com

Cheers,

Nadav

