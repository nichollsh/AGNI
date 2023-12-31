netcdf tro.q{                                                                   

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
                q =   .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325321E-05,  .325321E-05,  .325321E-05,  .325321E-05,
                      .325346E-05,  .334501E-05,  .368621E-05,  .441302E-05,
                      .850636E-05,  .169269E-04,  .308442E-04,  .523262E-04,
                      .836828E-04,  .127381E-03,  .185968E-03,  .261945E-03,
                      .349403E-03,  .437558E-03,  .539768E-03,  .657521E-03,
                      .791986E-03,  .944238E-03,  .111546E-02,  .130671E-02,
                      .151904E-02,  .175358E-02,  .201119E-02,  .229305E-02,
                      .259992E-02,  .293291E-02,  .329271E-02,  .368018E-02,
                      .409583E-02,  .453971E-02,  .514635E-02,  .602305E-02,
                      .694465E-02,  .789125E-02,  .884967E-02,  .966753E-02,
                      .101390E-01,  .106334E-01,  .111639E-01,  .117324E-01,
                      .123400E-01,  .129887E-01,  .136809E-01,  .144186E-01,
                      .152041E-01,  .158882E-01;

}                                                                               
