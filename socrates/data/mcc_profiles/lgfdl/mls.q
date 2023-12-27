netcdf mls.q{                                                                   

dimensions:
    lat            =   1;
    lon            =   1;
    plev           = 122;


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
             plev =   .500000E-01,  .108300E+00,  .126260E+00,  .147210E+00,
                      .171640E+00,  .200110E+00,  .233320E+00,  .272030E+00,
                      .317160E+00,  .369780E+00,  .431130E+00,  .502660E+00,
                      .586060E+00,  .683300E+00,  .796670E+00,  .928850E+00,
                      .108300E+01,  .126260E+01,  .147210E+01,  .171640E+01,
                      .200110E+01,  .233320E+01,  .272030E+01,  .317160E+01,
                      .369780E+01,  .431130E+01,  .502660E+01,  .586060E+01,
                      .683300E+01,  .796670E+01,  .928850E+01,  .108300E+02,
                      .126260E+02,  .147210E+02,  .171640E+02,  .200110E+02,
                      .233320E+02,  .272030E+02,  .317160E+02,  .369780E+02,
                      .431130E+02,  .502660E+02,  .586060E+02,  .683300E+02,
                      .796670E+02,  .928850E+02,  .108300E+03,  .126260E+03,
                      .147210E+03,  .171640E+03,  .200110E+03,  .233320E+03,
                      .272030E+03,  .317160E+03,  .369780E+03,  .431130E+03,
                      .502660E+03,  .586060E+03,  .683300E+03,  .796670E+03,
                      .928850E+03,  .108300E+04,  .126260E+04,  .147210E+04,
                      .171640E+04,  .200110E+04,  .233320E+04,  .272030E+04,
                      .317160E+04,  .369780E+04,  .431130E+04,  .502660E+04,
                      .586060E+04,  .683300E+04,  .796670E+04,  .928850E+04,
                      .110000E+05,  .130000E+05,  .150000E+05,  .170000E+05,
                      .190000E+05,  .210000E+05,  .230000E+05,  .250000E+05,
                      .270000E+05,  .290000E+05,  .310000E+05,  .330000E+05,
                      .350000E+05,  .370000E+05,  .390000E+05,  .410000E+05,
                      .430000E+05,  .450000E+05,  .470000E+05,  .490000E+05,
                      .510000E+05,  .530000E+05,  .550000E+05,  .570000E+05,
                      .590000E+05,  .610000E+05,  .630000E+05,  .650000E+05,
                      .670000E+05,  .690000E+05,  .710000E+05,  .730000E+05,
                      .750000E+05,  .770000E+05,  .790000E+05,  .810000E+05,
                      .830000E+05,  .850000E+05,  .870000E+05,  .890000E+05,
                      .910000E+05,  .930000E+05,  .950000E+05,  .970000E+05,
                      .990000E+05,  .100650E+06;
                q =   .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .400396E-05,  .400396E-05,
                      .400396E-05,  .400396E-05,  .404675E-05,  .465452E-05,
                      .959103E-05,  .208308E-04,  .406392E-04,  .726371E-04,
                      .120825E-03,  .173983E-03,  .218272E-03,  .269601E-03,
                      .328891E-03,  .396808E-03,  .474035E-03,  .561250E-03,
                      .654909E-03,  .746024E-03,  .846033E-03,  .956677E-03,
                      .107870E-02,  .121317E-02,  .136094E-02,  .154460E-02,
                      .180383E-02,  .208887E-02,  .239822E-02,  .273165E-02,
                      .308858E-02,  .346847E-02,  .387043E-02,  .429348E-02,
                      .473662E-02,  .519847E-02,  .567762E-02,  .617238E-02,
                      .668094E-02,  .720089E-02,  .772954E-02,  .826691E-02,
                      .881111E-02,  .936029E-02,  .991444E-02,  .104717E-01,
                      .110308E-01,  .114929E-01;

}                                                                               