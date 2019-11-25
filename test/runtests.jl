using msJ, Test
using Plots

function tests()
    @testset "Subset of tests"  begin
        info = MSJ.info("test.mzXML", verbose = true)
        @test info[1] == "parentFile: test.raw"                                        #1
        @test info[9] == "6 scans"                                                     #2
        @test info[10] == "MS1+"                                                       #3
        @test info[11] == "MS2+ 1255.5  CID(CE=18)"                                    #4
        @test info[12] == "MS3+ 902.33  PQD(CE=35)"                                    #5
        
        scans = MSJ.load("test.mzXML")
        @test eltype(scans)              == MSJ.MSscan                                 #6
        @test length(scans)              == 6                                          #7
        @test scans[1].num               == 1                                          #8
        @test scans[2].level             == 2                                          #9
        @test scans[3].polarity          == "+"                                        #10
        @test scans[2].activationMethod  == "CID"                                      #11
        @test scans[3].collisionEnergy   == 35.0                                       #12
        @test size(scans[1].int, 1)      == 22320                                      #13
        
        rt = MSJ.retention_time("test.mzXML")
        @test length(rt) == 6                                                          #14
        
        cr = MSJ.chromatogram("test.mzXML", method = MSJ.TIC() )
        @test length(cr.rt) == 6                                                       #15
        
        cr = MSJ.chromatogram("test.mzXML", method = MSJ.MZ([0, 500]))
        @test length(cr.rt) == 6                                                       #16
        
        cr = MSJ.chromatogram("test.mzXML", method = MSJ.∆MZ([1000, 1]))
        @test length(cr.rt) == 6                                                       #17
        
        cr = MSJ.chromatogram("test.mzXML", method = MSJ.BasePeak())
        @test length(cr.rt) == 6                                                       #18
        
        cr = MSJ.chromatogram("test.mzXML", MSJ.Polarity("+"), MSJ.Scan(2),MSJ.Precursor(1255.5), MSJ.Activation_Energy(18), MSJ.Activation_Method("CID"), MSJ.Level(2) )
        @test length(cr.rt) == 1                                                       #19

        rt = MSJ.retention_time(scans)
        @test length(rt) == 6                                                          #20
        
        cr = MSJ.chromatogram(scans, method = MSJ.TIC() )
        @test length(cr.rt) == 6                                                       #21
        
        cr = MSJ.chromatogram(scans, method = MSJ.MZ([0, 500]))
        @test length(cr.rt) == 6                                                       #22
        
        cr = MSJ.chromatogram(scans, method = MSJ.∆MZ([1000, 1]))
        @test length(cr.rt) == 6                                                       #23
        
        cr = MSJ.chromatogram(scans, method = MSJ.BasePeak())
        @test length(cr.rt) == 6                                                       #24
        
        cr = MSJ.chromatogram(scans, MSJ.Polarity("+"),MSJ.Scan(2),MSJ.Precursor(1255.5),MSJ.Activation_Energy(18),MSJ.Activation_Method("CID"),MSJ.Level(2) )
        @test (cr.rt, cr.ic) == ([0.7307], [9727.2])                                   #25
        
        cr = MSJ.chromatogram(scans, MSJ.Polarity(["+"]),MSJ.Scan([2,3]),MSJ.Precursor([1255.5, 902.33]),MSJ.Activation_Energy([18, 35]),MSJ.Activation_Method(["CID", "PQD"]),MSJ.Level([2, 3]) )
        @test (cr.rt, cr.ic) == ([0.7307, 2.1379], [9727.2, 11.3032])                  #26
        
        ms = MSJ.average("test.mzXML")
        @test length(ms.num) == 6                                                      #27
        
        ms = MSJ.average("test.mzXML", MSJ.Polarity("+"),MSJ.Scan(2),MSJ.Precursor(1255.5),MSJ.Activation_Energy(18),MSJ.Activation_Method("CID"),MSJ.RT(1),MSJ.IC([0, 1e4]))
        @test ms isa MSJ.MSscan                                                        #28
        @test ms.num == 2                                                              #29
                
        ms = MSJ.average(scans)
        @test length(ms.num) == 6                                                      #30
        
        ms = MSJ.average(scans, MSJ.Polarity("+"),MSJ.Scan(2),MSJ.Precursor(1255.5),MSJ.Activation_Energy(18),MSJ.Activation_Method("CID"),MSJ.RT(1),MSJ.IC([0, 1e4]))
        @test ms isa MSJ.MSscan                                                        #31
        @test ms.num == 2                                                              #32
        
        ms = MSJ.average(scans, MSJ.Polarity(["+"]),MSJ.Scan([2,3]),MSJ.Precursor([1255.5, 902.33]),MSJ.Activation_Energy([18, 35]),MSJ.Activation_Method(["CID", "PQD"]),MSJ.RT([1,2]),MSJ.IC([0, 1e4]))
        @test ms isa MSJ.MSscans                                                       #33
        @test ms.num == [2, 3]                                                         #34

        ms = MSJ.average("test.mzXML", MSJ.RT( [[1,2], [2,3]] ), stats = false )
        @test ms isa MSJ.MSscans                                                       #35
        @test ms.num == [2, 3, 4]                                                      #36

        ms = MSJ.average("test.mzXML", MSJ.Polarity(["+"]),MSJ.Scan([2,3]),MSJ.Precursor([1255.5, 902.33]),MSJ.Activation_Method(["CID", "PQD"]),MSJ.RT([1,2]),MSJ.IC([0, 1e4]))   #MSJ.Activation_Energy([18., 35.]),
        @test ms isa MSJ.MSscans                                                       #37
        @test ms.num == [2, 3]                                                         #38

        cr = MSJ.chromatogram("test.mzXML", MSJ.Polarity(["+"]),MSJ.Scan([2,3]),MSJ.Precursor([1255.5, 902.33]),MSJ.Activation_Method(["CID", "PQD"]),MSJ.Level([2, 3]) )   #MSJ.Activation_Energy([18.0, 35.0]),
        @test length(cr.rt) == 2                                                       #39

        ms = MSJ.average(scans, MSJ.RT( [[1,2], [2,3]] ), stats = false )
        @test ms isa MSJ.MSscans                                                       #40
        @test ms.num == [2, 3, 4]                                                      #41

        a = scans[1] / 2.
        @test a.tic == 2.540975e6                                                      #42

        a = scans[1] * 2.
        @test a.tic == 1.01639e7                                                       #43

        a = ms * 2.
        @test a.tic == 3.2120923354666666e6                                            #44

        a = 2. * scans[1]
        @test a.tic == 1.01639e7                                                       #45

        a = scans[1] * scans[2]
        @test a.tic == 4.943314404e10                                                  #46

        a = scans[2] * scans[1]
        @test a.tic == 4.943314404e10                                                  #47

        a = scans[1] - scans[2]
        @test a.tic == 5.0722228e6                                                     #48

        a = scans[2] - scans[1]
        @test a.tic == -5.0722228e6                                                    #49

        a = scans[1] - scans[4]
        @test a.tic == 273550.0                                                        #50

        b = ms - scans[1]
        @test b.num == [2,3,4]                                                         #51

        b = ms + scans[1]
        @test b.num == [2,3,4,1]                                                       #52

        b = scans[1] + ms 
        @test b.num == [1,2,3,4]                                                       #53

        b = scans[1] + scans[4] 
        @test b.num == [1,4]                                                           #54

        a = MSJ.avg(scans[1], scans[2])
        @test a.num == [1,2]                                                           #55

        a = MSJ.avg(scans[1], scans[4])
        @test a.num == [1,4]                                                           #56

        info = MSJ.info("test64.mzXML")
        @test info[2] == "MS1-"                                                        #57
        
        scans = MSJ.load("test64.mzXML")
        @test eltype(scans)              == MSJ.MSscan                                 #58

        info = MSJ.info("test.mzXMLM")
        @test info.msg == "File format not supported."                                 #59

        scans = MSJ.load("test.mzXMLL")
        @test info.msg == "File format not supported."                                 #60

        scans = MSJ.load("bad1.mzXML")
        @test scans.msg == "Not an mzXML file."                                        #61

        scans = MSJ.info("bad1.mzXML")
        @test scans.msg == "Not an mzXML file."                                        #62

        scans = MSJ.load("bad2.mzXML")
        @test scans[1].num == 0                                                        #63

        scans = MSJ.load("bad3.mzXML")
        @test scans[1].num == scans[2].num == 0                                        #64

        cr = MSJ.chromatogram("test.mzXML", method = MSJ.∆MZ([1, 2]))
        @test cr.msg == "Bad mz ± ∆mz values."                                         #65
        
        cr = MSJ.chromatogram(scans, method = MSJ.∆MZ([1, 2]))
        @test cr.msg == "Bad mz ± ∆mz values."                                         #66

        scans = MSJ.load("test.mzXML")
        @test MSJ.smooth(scans[1], method = MSJ.SG(7,15,0)) isa MSJ.MSscan             #67

        a = MSJ.avg(scans[1], scans[4])
        @test MSJ.smooth(a) isa MSJ.MSscans                                            #68

        a = scans[1] * scans[4]
        @test a.num == [1,4]                                                           #69

        a = (scans[2]+scans[3]) - scans[1]
        @test a.num == [2, 3]                                                          #70

        a = (scans[1] + scans[4]) - (scans[1] - scans[4])
        @test a.num == [1, 4]                                                          #71

        a = scans[1] + MSJ.avg(scans[2], scans[5])
        @test a.num == [1, 2, 5]                                                       #72

        a = MSJ.smooth(scans[1], method = MSJ.SG(5,9,0))
        @test a.num == 1                                                               #73

       a = MSJ.centroid(scans[1], method = MSJ.TBPD(:gauss, 4500., 0.2))               #74
       @test length(a.int) == 976

       @test typeof(plot(scans[1], method = :relative)) == Plots.Plot{Plots.GRBackend} #75
       @test typeof(plot(scans[1], method = :absolute)) == Plots.Plot{Plots.GRBackend} #76

       a = MSJ.avg(scans[2], scans[5])
       @test typeof(plot( a, method = :relative )) == Plots.Plot{Plots.GRBackend}      #77
       @test typeof(plot( a, method = :absolute )) == Plots.Plot{Plots.GRBackend}      #78

       cr = MSJ.chromatogram(scans)
       @test typeof(plot( cr, method = :relative )) == Plots.Plot{Plots.GRBackend}     #79
       @test typeof(plot( cr, method = :absolute )) == Plots.Plot{Plots.GRBackend}     #80

       a = MSJ.centroid(scans[1], method = MSJ.TBPD(:voigt, 4500., 0.2))
       @test length(a.int) == 961                                                      #81

       a = MSJ.centroid(scans[1], method = MSJ.TBPD(:lorentz, 4500., 0.2))
       @test length(a.int) == 969                                                      #82

       a = MSJ.centroid(scans[1], method = MSJ.TBPD(:other, 4500., 0.2))
       @test a.msg == "Unsupported peak profile. Use :gauss, :lorentz or :voigt."      #83

       a = MSJ.centroid(scans[1], method = MSJ.SNRA(1., 100))
       @test length(a.int) == 109                                                      #84

       s1 = MSJ.extract(scans, MSJ.Activation_Energy([18,35]))
       @test length(s1) == 4                                                           #85

       s1 = MSJ.extract("test.mzXML", MSJ.Activation_Energy(18))
       @test length(s1) == 2                                                           #86

       s1 = MSJ.extract(scans, MSJ.Scan(1))
       @test length(s1) == 1                                                           #87

       s1 = MSJ.extract("test.mzXML", MSJ.Scan(1))
       @test length(s1) == 1                                                           #88

       bs = MSJ.baseline_correction(scans, method = MSJ.TopHat(1))
       @test length(bs) == 6                                                           #89

       bs = MSJ.baseline_correction(scans[1], method = MSJ.TopHat(1))
       @test length(bs.int) == length(scans[1].int)                                    #90

       c = MSJ.centroid(scans, method = MSJ.TBPD(:gauss, 4500., 0.2)) ;
       bs = MSJ.baseline_correction( c, method = MSJ.LOESS(3))
       @test length(bs) == 6                                                           #91

       bs = MSJ.baseline_correction(c[1], method = MSJ.LOESS(3))
       @test length(bs.int) == length(c[1].int)                                        #92

       bs = MSJ.baseline_correction(scans, method = MSJ.IPSA(51,100))
       @test length(bs) == 6                                                           #93

       bs = MSJ.baseline_correction(scans[1], method = MSJ.IPSA(51,100))
       @test length(bs.int) == length(scans[1].int)                                    #94

       a = MSJ.smooth(scans, method = MSJ.SG(5,9,0))
       @test length(a) == 6                                                            #95

       c = MSJ.centroid(scans, method = MSJ.TBPD(:lorentz, 4500., 0.2)) ;
       d = MSJ.centroid(scans, method = MSJ.TBPD(:voigt, 4500., 0.2)) ;
       @test length(c) == length(d)                                                    #96

       a = MSJ.centroid(scans[3], method = MSJ.SNRA(1., 100))
       @test length(a.int) == 0                                                        #97

       bs = MSJ.baseline_correction(scans[1], method = MSJ.IPSA(50,100))
       @test length(bs.int) == length(scans[1].int)                                    #98

       cr = MSJ.chromatogram(scans, method = MSJ.BasePeak() )
       @test length(cr.rt) == 6                                                        #99

       f = MSJ.formula("CH3(13C)10H3Kr(NaH2)2")                                        #100
      @test f == Dict("Na" => 2,"Kr" => 1,"C" => 1,"13C" => 10,"H" => 10) 

       m = MSJ.masses("C254 H377 N65 O75 S6")                                          #101
      @test m == Dict("Monoisotopic" => 5729.60087099839, "Average" => 5733.55, "Nominal" => 5727.0) 

      I = MSJ.isotopic_distribution("CH4", 0.9999, charge = +1)                        #102
      @test I[2,1:end] == [16.03130012908, 0.9887541751052761, 1, 0, 4, 0]
      
      a = MSJ.simulate(I, 0.4, Npoints = 5)                                            #103
      @test a.int == [100.0, 36.34733624865424, 2.154581386492942, 1.147877462218771, 0.4130025252583078]

       m = MSJ.masses(f)                                                               #104
      @test m == Dict("Monoisotopic" => 282.0028349717, "Average" => 281.902086912, "Nominal" => 282.0) 

      I = MSJ.isotopic_distribution(f, 0.9999, charge = +1)                            #105
      @test I[2,1:end][1:2] == [282.0028349717, 0.5630635281692917]
      
      a = MSJ.simulate(I, 0.4, model=:lorentz, Npoints = 5)                            #106
      @test a.int == [1.0392077560077122, 17.523433547791814, 100.0, 1.8273069587773836, 0.41728492055754945]

    end
end
tests()
