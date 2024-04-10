netcdf gam3_nts.q{                                                          

dimensions:
    lat            =   3;
    lon            =   1;
    level          =  18;


variables:
    float lat(lat);                                                           
             lat:units = "degree";                                            
             lat:title = "LATITUDE";                                         
    float lon(lon);                                                           
             lon:units = "degree";                                            
             lon:title = "LONGITUDE";                                        
    int level(level);
             level:units = "None";
             level:title = "LEVEL";

    float plev(level,lon,lat);
             plev:units = "Pa";
             plev:title = "PRESSURE";
    float q(level,lon,lat);                                                     
             q:units = "none";                                                
             q:title = "MMR of  water vapour";                                               

data:                                                                           
              lat =  0.481890E+02, 0.000000E+00, -0.481890E+02;
              lon =  0.000000E+00;
            level =     1,     2,     3,     4,     5,     6,
                        7,     8,     9,    10,    11,    12,
                       13,    14,    15,    16,    17,    18;
             plev =  0.100000E+03, 0.100000E+03, 0.100000E+03, 0.300000E+03,
                     0.300000E+03, 0.300000E+03, 0.100000E+04, 0.100000E+04,
                     0.100000E+04, 0.300000E+04, 0.300000E+04, 0.300000E+04,
                     0.500000E+04, 0.500000E+04, 0.500000E+04, 0.700000E+04,
                     0.700000E+04, 0.700000E+04, 0.100000E+05, 0.905919E+04,
                     0.100000E+05, 0.150000E+05, 0.915024E+04, 0.150000E+05,
                     0.200000E+05, 0.100000E+05, 0.200000E+05, 0.202088E+05,
                     0.150000E+05, 0.203471E+05, 0.204119E+05, 0.200000E+05,
                     0.205516E+05, 0.250000E+05, 0.250000E+05, 0.250000E+05,
                     0.300000E+05, 0.300000E+05, 0.300000E+05, 0.400000E+05,
                     0.400000E+05, 0.400000E+05, 0.500000E+05, 0.500000E+05,
                     0.500000E+05, 0.700000E+05, 0.700000E+05, 0.700000E+05,
                     0.850000E+05, 0.850000E+05, 0.850000E+05, 0.944467E+05,
                     0.989738E+05, 0.983021E+05;
                q =  0.353000E-05, 0.345000E-05, 0.360000E-05, 0.341000E-05,
                     0.307000E-05, 0.341000E-05, 0.334000E-05, 0.287000E-05,
                     0.325000E-05, 0.315000E-05, 0.275000E-05, 0.308000E-05,
                     0.285000E-05, 0.251000E-05, 0.268000E-05, 0.263000E-05,
                     0.232000E-05, 0.237000E-05, 0.236000E-05, 0.237000E-05,
                     0.201000E-05, 0.318000E-05, 0.238000E-05, 0.264000E-05,
                     0.937000E-05, 0.239000E-05, 0.783000E-05, 0.104000E-04,
                     0.660000E-05, 0.932000E-05, 0.115000E-04, 0.220000E-04,
                     0.102000E-04, 0.324000E-04, 0.618000E-04, 0.271000E-04,
                     0.104000E-03, 0.281000E-03, 0.751000E-04, 0.374000E-03,
                     0.883000E-03, 0.278000E-03, 0.781000E-03, 0.178000E-02,
                     0.612000E-03, 0.214000E-02, 0.475000E-02, 0.160000E-02,
                     0.388000E-02, 0.904000E-02, 0.321000E-02, 0.552000E-02,
                     0.152000E-01, 0.584000E-02;

}