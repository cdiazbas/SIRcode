# SIRcode
 SIR: Stokes Inversion based on Response functions
 
 SIR is a general-purpose code capable of dealing with gradients of the physical quantities with height. It admits one and two-component model atmospheres. It allows the recovery of the stratification of the temperature, the magnetic field vector, and the line of sight velocity through the atmosphere, and the micro- and macroturbulence velocities - which are assumed to be constant with depth. It is based on the response functions, which enter a Marquardt nonlinear least-squares algorithm in a natural way. Response functions are calculated at the same time as the full radiative transfer equation for polarized light is integrated, which determines values of many free parameters in a reasonable computation time. SIR demonstrates high stability, accuracy, and uniqueness of results, even when simulated observations present signal-to-noise ratios of the order of the lowest acceptable values in real observations.
 
For any question, send a mail to **brc@iac.es**.
 
ASCL Code Record: http://ascl.net/1212.008

Appears in: http://adsabs.harvard.edu/abs/1992ApJ...398..375R

-

### sirtools.py
python tools for SIR-files
    
    1.-  lambda_mA, stokesIQUV, [nL,posi,nN] = lperfil(filename)
    
    2.-  wperfil(filename, numberLine, lambda_mA, stokes)
    
    3.-  [tau, todoPlot] = lmodel8(filename, verbose=True)
    
    4.-  wmodel8(modelo, filename, verbose=False)
    
    5.-  mapa = readSIRMap(resultadoSir, magnitud)
    
    6.-  [height, width, nlambda] = shapeSIRMap(resultadoSir)
    
    7.-  mapa = readSIRProfileMap(resultadoSir, Nstoke)
    
    8.-  index = tauIndex(resultadoSir, logTau)
