netcdf mls.t{                                                                   

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

    float t(plev,lon,lat);                                                     
             t:units = "K";                                                   
             t:title = "TEMPERATURE";                                                        

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
                t =   .210721E+03,  .210829E+03,  .210959E+03,  .211203E+03,
                      .211448E+03,  .211692E+03,  .211937E+03,  .212182E+03,
                      .212428E+03,  .212673E+03,  .212919E+03,  .213165E+03,
                      .213411E+03,  .213657E+03,  .213904E+03,  .214150E+03,
                      .214397E+03,  .214644E+03,  .214892E+03,  .215139E+03,
                      .215387E+03,  .215635E+03,  .215883E+03,  .216132E+03,
                      .216381E+03,  .216637E+03,  .216933E+03,  .217458E+03,
                      .218693E+03,  .220924E+03,  .223792E+03,  .226864E+03,
                      .230003E+03,  .233189E+03,  .236417E+03,  .239690E+03,
                      .243006E+03,  .246367E+03,  .249773E+03,  .253225E+03,
                      .256722E+03,  .260263E+03,  .263844E+03,  .267452E+03,
                      .271052E+03,  .274329E+03,  .275966E+03,  .275392E+03,
                      .273651E+03,  .271287E+03,  .268548E+03,  .265643E+03,
                      .262690E+03,  .259738E+03,  .256808E+03,  .253909E+03,
                      .251054E+03,  .248290E+03,  .245678E+03,  .243197E+03,
                      .240774E+03,  .238381E+03,  .236014E+03,  .233671E+03,
                      .231355E+03,  .229086E+03,  .226962E+03,  .225167E+03,
                      .223735E+03,  .222479E+03,  .221268E+03,  .220069E+03,
                      .218878E+03,  .217697E+03,  .216608E+03,  .215933E+03,
                      .215761E+03,  .215752E+03,  .215767E+03,  .216036E+03,
                      .217541E+03,  .220742E+03,  .224629E+03,  .228385E+03,
                      .231906E+03,  .235214E+03,  .238335E+03,  .241291E+03,
                      .244098E+03,  .246772E+03,  .249322E+03,  .251756E+03,
                      .254077E+03,  .256288E+03,  .258393E+03,  .260399E+03,
                      .262316E+03,  .264156E+03,  .265930E+03,  .267645E+03,
                      .269307E+03,  .270920E+03,  .272489E+03,  .274016E+03,
                      .275503E+03,  .276953E+03,  .278368E+03,  .279749E+03,
                      .281096E+03,  .282410E+03,  .283684E+03,  .284912E+03,
                      .286077E+03,  .287160E+03,  .288150E+03,  .289052E+03,
                      .289884E+03,  .290667E+03,  .291418E+03,  .292146E+03,
                      .292857E+03,  .293554E+03,  .294000E+03;

}                                                                               