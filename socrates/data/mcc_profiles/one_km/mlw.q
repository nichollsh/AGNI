netcdf mlw.q{                                                                   

dimensions:
    lat            =   1;
    lon            =   1;
    plev           =  33;


variables:
    float lat(lat);                                                           
             lat:units = "degree";                                            
             lat:title = "LATITUDE";                                         
    float lon(lon);                                                           
             lon:units = "degree";                                            
             lon:title = "LONGITUDE";                                        
    float plev(plev);                                                         
             plev:units = "Pa";                                               
             plev:title = "PRESSURE";                                        

    float q(plev,lon,lat);                                                     
             q:units = "None";                                                
             q:title = "MMR OF WATER VAPOUR";                                                

data:                                                                           
              lat =   .000000E+00;
              lon =   .000000E+00;
             plev =   .300000E-01,  .467000E+01,  .682000E+02,  .129000E+03,
                      .253000E+03,  .518000E+03,  .111000E+04,  .243000E+04,
                      .286000E+04,  .334000E+04,  .391000E+04,  .458000E+04,
                      .537000E+04,  .628000E+04,  .735000E+04,  .861000E+04,
                      .100700E+05,  .117800E+05,  .137800E+05,  .161000E+05,
                      .188200E+05,  .219900E+05,  .256800E+05,  .299200E+05,
                      .347300E+05,  .401600E+05,  .462700E+05,  .531300E+05,
                      .608100E+05,  .693800E+05,  .789700E+05,  .897300E+05,
                      .101800E+06;
                q =   .400000E-05,  .399943E-05,  .399821E-05,  .399770E-05,
                      .400000E-05,  .400051E-05,  .399888E-05,  .400000E-05,
                      .400086E-05,  .400739E-05,  .400757E-05,  .400216E-05,
                      .400460E-05,  .400197E-05,  .399832E-05,  .399856E-05,
                      .400000E-05,  .404233E-05,  .512239E-05,  .668740E-05,
                      .906969E-05,  .127002E-04,  .184185E-04,  .346395E-04,
                      .670241E-04,  .145600E-03,  .285757E-03,  .510053E-03,
                      .833132E-03,  .125677E-02,  .173578E-02,  .215146E-02,
                      .269024E-02;

}                                                                               