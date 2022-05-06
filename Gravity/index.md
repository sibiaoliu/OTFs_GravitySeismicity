# Results of 3d thermal models for different OTF systems with spreading rate spectrum.

Here, we run 15 sets of OTF systems in which the spreading rate differs from >10 cm/yr (fast) and 10-5 cm/yr (intermediate) to 5-2 cm/yr (slow) and <2 cm/yr (ultra-slow).

​	Fast: Clipperton, Gofar (south), Garrett, Discovery

​	Intermediate:  Chile Ridge, Vlamingh, Blanco, Pitman, SWIR-88S

​	Slow: Marie Celeste (east), Oceanographer, Atlantis, (*Romanche,) Kane, Marathon

​	Ultra-slow: Atlantis_II, Marion, Du Troit


Model simulation time depends on the main transform margin length and half-spreading rate.
For example, in the Atlantis case, Main Transform Margin length (MTL)=62 km, Half-spreading rate (v_half)=1.2 cm/yr, so the model ends at 6.46 Ma (=1.25 * MTL / v_half). 


Tested parameter - **rheology option**.

1. <u>Viscous</u>: only dislocation creep (disl.) applies to the model.

2. Viscous/<u>isoviscous</u>: constant viscosity (10^21 Pas) in the model.

3. <u>Visco-plastic</u>: disl. + plasticity (plas.).

4. <u>Viscous-elasto-plastic</u>: disl. + plas. + elasticity (elas.).

Previous gravity studies: 
Age-based approach for thermal correction
    Using Half-space cooling thermal model, which applies an ideally 2-D static age-based T field (Geodynamic book) without the active, dynamic mantle flow. Also, it does not account for cooling across the transform (3-D) - ( Georgen et al., 2001; Wang et al., 2011; Zhang et al., 2018)
    The use of crustal age rather than passive £ow modeling for the calculation of thermal correction, however, may lead to an underestimation of the thermal effects associated with large fracture zones such as the Andrew Bain FZ, because the e°ects of reduced upwelling and lateral asthenospheric £ow were not considered. -Geprgen et al., 2001

Passive upwelling mantle flow approach for thermal correction, considering 3D nature of flow around transform faults
    Using Isoviscous - ( Morgan and Forysth, 1988: <-- used this approach in Kuo & Forysth, 1988;Lin et al., 1990; Escartin & Lin, 1995; Gregg et al, 2007; Blackman et al., 2008)
To estimate the steady-state T field in the ridge-transform system:
    Using T-dependent Viscosity - (Behn et al., 2007)
    Using T-dependent vp -  (Behn et al., 2007; Gregg et al., 2009)

 
