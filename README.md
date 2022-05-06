# Thermal Structure of OTFs

This project aims to investigate the role of the rheologic-control thermal structure in the mid-ocean ridge-transform system in the estimation of gravity and seismicity. Here we consider the effects of rheology and hydrothermal circulation on the T field.

# Work flow

1. RMBA calculation and gravity-derived crustal thickness prediction

(1) Mantle bouger anormaly(MBA). We intergrate the multibeam high resolution bathymetry data with ETOPO 1 data to calculate MBA (FAA->BGA->MBA), following the same approach in Guo et al. (2022).

(2) Thermal correction of the mantle. Here, we use two well-known approach: age-dependent thermal modeling and passive  upwelling mantle flow modeling. The former uses the half-spaceing cooling process and applies an ideally 2D static age-based T field without the active, dynamic mantle flow. In the latter, a 3D nature of mantle flow is formed below the ridge-transform-ridge boundary. After the thermal field becomes steady-state, we extract the thermal structure under the mid-ocean ridge-transform system (MORTs) which is the same area as the bathymetry. In the 3D geodynamic flow model, we applied different rheologic settings (i.e., isoviscous, nonlinear viscosity, viscoplastic, viscoelastoplastic) and compare them to the reference age-dependent T structure.

(3) Residual mantle bouger anormaly (RMBA). We calculate the RMBA field by subtracting the thermal correction of lithospheric cooling from MBA.

(4) Gravity-derived crustal thickness. We invert the RMBA filed to the crustal thickness (or Moho depth) field for all interested MORTs, following the same approach in Chappell and Kusznir (2008).

2. Estimation of the maximum depth and temperature of the seismicity in the MORTs

(1) Compiling the distribution of earthquakes (magnitude >5?) for each OTFs.

(2) Extracting the transect of OTF and plotting the Depth-Distance-T maps with confined seismogenic zone (i.e., *C isotherms).

(3) Comparing T field with other observations such as heat flow, rock samples (?) in each OTFs to validate the estimated T field of seismicity from thermal models.

More...
