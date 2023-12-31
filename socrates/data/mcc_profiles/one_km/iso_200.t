netcdf iso_200.t{                                                               

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

    float t(plev,lon,lat);                                                     
             t:units = "K";                                                   
             t:title = "TEMPERATURE";                                                        

data:                                                                           
              lat =   .000000E+00;
              lon =   .000000E+00;
             plev =   .300000E-01,  .671000E+01,  .951000E+02,  .176000E+03,
                      .333000E+03,  .652000E+03,  .132000E+04,  .277000E+04,
                      .322000E+04,  .376000E+04,  .437000E+04,  .510000E+04,
                      .595000E+04,  .695000E+04,  .812000E+04,  .950000E+04,
                      .111000E+05,  .130000E+05,  .153000E+05,  .179000E+05,
                      .209000E+05,  .243000E+05,  .281000E+05,  .324000E+05,
                      .372000E+05,  .426000E+05,  .487000E+05,  .554000E+05,
                      .628000E+05,  .710000E+05,  .802000E+05,  .902000E+05,
                      .101300E+06;
                t =   .200000E+03,  .200000E+03,  .200000E+03,  .200000E+03,
                      .200000E+03,  .200000E+03,  .200000E+03,  .200000E+03,
                      .200000E+03,  .200000E+03,  .200000E+03,  .200000E+03,
                      .200000E+03,  .200000E+03,  .200000E+03,  .200000E+03,
                      .200000E+03,  .200000E+03,  .200000E+03,  .200000E+03,
                      .200000E+03,  .200000E+03,  .200000E+03,  .200000E+03,
                      .200000E+03,  .200000E+03,  .200000E+03,  .200000E+03,
                      .200000E+03,  .200000E+03,  .200000E+03,  .200000E+03,
                      .200000E+03;

}                                                                               
