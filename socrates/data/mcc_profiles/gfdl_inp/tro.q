netcdf tro.q{                                                                   

dimensions:
    lat            =   1;
    lon            =   1;
    plev           = 123;


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
             plev =   .100000E-03,  .100000E+00,  .116591E+00,  .135936E+00,
                      .158489E+00,  .184785E+00,  .215444E+00,  .251189E+00,
                      .292865E+00,  .341455E+00,  .398107E+00,  .464159E+00,
                      .541170E+00,  .630957E+00,  .735642E+00,  .857696E+00,
                      .100000E+01,  .116591E+01,  .135936E+01,  .158489E+01,
                      .184785E+01,  .215444E+01,  .251189E+01,  .292865E+01,
                      .341455E+01,  .398107E+01,  .464159E+01,  .541170E+01,
                      .630957E+01,  .735642E+01,  .857696E+01,  .100000E+02,
                      .116591E+02,  .135936E+02,  .158489E+02,  .184785E+02,
                      .215444E+02,  .251189E+02,  .292865E+02,  .341455E+02,
                      .398107E+02,  .464159E+02,  .541170E+02,  .630957E+02,
                      .735642E+02,  .857696E+02,  .100000E+03,  .116591E+03,
                      .135936E+03,  .158489E+03,  .184785E+03,  .215444E+03,
                      .251189E+03,  .292865E+03,  .341455E+03,  .398107E+03,
                      .464159E+03,  .541170E+03,  .630957E+03,  .735642E+03,
                      .857696E+03,  .100000E+04,  .116591E+04,  .135936E+04,
                      .158489E+04,  .184785E+04,  .215444E+04,  .251189E+04,
                      .292865E+04,  .341455E+04,  .398107E+04,  .464159E+04,
                      .541170E+04,  .630957E+04,  .735642E+04,  .857696E+04,
                      .100000E+05,  .120000E+05,  .140000E+05,  .160000E+05,
                      .180000E+05,  .200000E+05,  .220000E+05,  .240000E+05,
                      .260000E+05,  .280000E+05,  .300000E+05,  .320000E+05,
                      .340000E+05,  .360000E+05,  .380000E+05,  .400000E+05,
                      .420000E+05,  .440000E+05,  .460000E+05,  .480000E+05,
                      .500000E+05,  .520000E+05,  .540000E+05,  .560000E+05,
                      .580000E+05,  .600000E+05,  .620000E+05,  .640000E+05,
                      .660000E+05,  .680000E+05,  .700000E+05,  .720000E+05,
                      .740000E+05,  .760000E+05,  .780000E+05,  .800000E+05,
                      .820000E+05,  .840000E+05,  .860000E+05,  .880000E+05,
                      .900000E+05,  .920000E+05,  .940000E+05,  .960000E+05,
                      .980000E+05,  .100000E+06,  .101300E+06;
                q =   .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325000E-05,  .325000E-05,  .325000E-05,  .325000E-05,
                      .325011E-05,  .329790E-05,  .351826E-05,  .405698E-05,
                      .651015E-05,  .128093E-04,  .240200E-04,  .417677E-04,
                      .682389E-04,  .105819E-03,  .157008E-03,  .224327E-03,
                      .306015E-03,  .393705E-03,  .488854E-03,  .598790E-03,
                      .724830E-03,  .868115E-03,  .102977E-02,  .121089E-02,
                      .141255E-02,  .163582E-02,  .188171E-02,  .215125E-02,
                      .244539E-02,  .276507E-02,  .311118E-02,  .348453E-02,
                      .388574E-02,  .431514E-02,  .484042E-02,  .558224E-02,
                      .648043E-02,  .741370E-02,  .836521E-02,  .925203E-02,
                      .989498E-02,  .103772E-01,  .108894E-01,  .114386E-01,
                      .120261E-01,  .126539E-01,  .133237E-01,  .140378E-01,
                      .147985E-01,  .156047E-01,  .158728E-01;

}                                                                               
