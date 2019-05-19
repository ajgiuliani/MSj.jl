using msJ, Test

function tests()
    @testset "Subset of tests"  begin
        info = msJ.info("test.mzXML", verbose = true)
        @test info[1] == "parentFile: test.raw"                                        #1
        @test info[9] == "6 scans"                                                     #2
        @test info[10] == "MS1+"                                                       #3
        @test info[11] == "MS2+ 1255.5  CID(CE=18)"                                    #4
        @test info[12] == "MS3+ 902.33  PQD(CE=35)"                                    #5
        
        scans = msJ.load("test.mzXML")
        @test eltype(scans)              == msJ.MSscan                                 #6
        @test length(scans)              == 6                                          #7
        @test scans[1].num               == 1                                          #8
        @test scans[2].level             == 2                                          #9
        @test scans[3].polarity          == "+"                                        #10
        @test scans[2].activationMethod  == "CID"                                      #11
        @test scans[3].collisionEnergy   == 35.0                                       #12
        @test size(scans[1].int, 1)      == 22320                                      #13
        
        rt = msJ.retention_time("test.mzXML")
        @test length(rt) == 6                                                          #14
        
        rt, tic = msJ.chromatogram("test.mzXML", method = msJ.TIC() )
        @test length(rt) == 6                                                          #15
        
        rt, tic = msJ.chromatogram("test.mzXML", method = msJ.MZ([0, 500]))
        @test length(rt) == 6                                                          #16
        
        rt, tic = msJ.chromatogram("test.mzXML", method = msJ.∆MZ([1000, 1]))
        @test length(rt) == 6                                                          #17
        
        rt, tic = msJ.chromatogram("test.mzXML", method = msJ.BasePeak())
        @test length(rt) == 6                                                          #18
        
        xrt, xtic = msJ.chromatogram("test.mzXML", msJ.Polarity("+"), msJ.Scan(2),msJ.Precursor(1255.5), msJ.Activation_Energy(18), msJ.Activation_Method("CID"), msJ.Level(2) )
        @test length(xrt) == 1                                                         #19

        rt = msJ.retention_time(scans)
        @test length(rt) == 6                                                          #20
        
        rt, tic = msJ.chromatogram(scans, method = msJ.TIC() )
        @test length(rt) == 6                                                          #21
        
        rt, tic = msJ.chromatogram(scans, method = msJ.MZ([0, 500]))
        @test length(rt) == 6                                                          #22
        
        rt, tic = msJ.chromatogram(scans, method = msJ.∆MZ([1000, 1]))
        @test length(rt) == 6                                                          #23
        
        rt, tic = msJ.chromatogram(scans, method = msJ.BasePeak())
        @test length(rt) == 6                                                          #24
        
        xrt, xtic = msJ.chromatogram(scans, msJ.Polarity("+"),msJ.Scan(2),msJ.Precursor(1255.5),msJ.Activation_Energy(18),msJ.Activation_Method("CID"),msJ.Level(2) )
        @test (xrt, xtic) == ([0.7307], [9727.2])                                      #25
        
        xrt, xtic = msJ.chromatogram(scans, msJ.Polarity(["+"]),msJ.Scan([2,3]),msJ.Precursor([1255.5, 902.33]),msJ.Activation_Energy([18, 35]),msJ.Activation_Method(["CID", "PQD"]),msJ.Level([2, 3]) )
        @test (xrt, xtic) == ([0.7307, 2.1379], [9727.2, 11.3032])                     #26
        
        ms = msJ.msfilter("test.mzXML")
        @test length(ms.num) == 6                                                      #27
        
        ms = msJ.msfilter("test.mzXML", msJ.Polarity("+"),msJ.Scan(2),msJ.Precursor(1255.5),msJ.Activation_Energy(18),msJ.Activation_Method("CID"),msJ.RT(1),msJ.IC([0, 1e4]))
        @test ms isa msJ.MSscan                                                        #28
        @test ms.num == 2                                                              #29
                
        ms = msJ.msfilter(scans)
        @test length(ms.num) == 6                                                      #30
        
        ms = msJ.msfilter(scans, msJ.Polarity("+"),msJ.Scan(2),msJ.Precursor(1255.5),msJ.Activation_Energy(18),msJ.Activation_Method("CID"),msJ.RT(1),msJ.IC([0, 1e4]))
        @test ms isa msJ.MSscan                                                        #31
        @test ms.num == 2                                                              #32
        
        ms = msJ.msfilter(scans, msJ.Polarity(["+"]),msJ.Scan([2,3]),msJ.Precursor([1255.5, 902.33]),msJ.Activation_Energy([18, 35]),msJ.Activation_Method(["CID", "PQD"]),msJ.RT([1,2]),msJ.IC([0, 1e4]))
        @test ms isa msJ.MSscans                                                       #33
        @test ms.num == [2, 3]                                                         #34

        ms = msJ.msfilter("test.mzXML", msJ.RT( [[1,2], [2,3]] ), stats = false )
        @test ms isa msJ.MSscans                                                       #35
        @test ms.num == [2, 3, 4]                                                         #36

        ms = msJ.msfilter("test.mzXML", msJ.Polarity(["+"]),msJ.Scan([2,3]),msJ.Precursor([1255.5, 902.33]),msJ.Activation_Method(["CID", "PQD"]),msJ.RT([1,2]),msJ.IC([0, 1e4]))   #msJ.Activation_Energy([18., 35.]),
        @test ms isa msJ.MSscans                                                       #37
        @test ms.num == [2, 3]                                                         #38

        xrt, xtic = msJ.chromatogram("test.mzXML", msJ.Polarity(["+"]),msJ.Scan([2,3]),msJ.Precursor([1255.5, 902.33]),msJ.Activation_Method(["CID", "PQD"]),msJ.Level([2, 3]) )   #msJ.Activation_Energy([18.0, 35.0]),
        @test length(xrt) == 2                                                        #39

        
    end
end
tests()
