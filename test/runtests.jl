using msJ, Test

function tests()
    @testset "Subset of tests"  begin
        info = msJ.info("test.mzXML", verbose = true)
        @test info[1] == "parentFile: test.raw"
        @test info[9] == "6 scans"
        @test info[10] == "MS1+"
        @test info[11] == "MS2+ 1255.5  CID(CE=18)"
        @test info[12] == "MS3+ 902.33  PQD(CE=35)"
        
        scans = msJ.load("test.mzXML")
        @test eltype(scans)              == msJ.MSscan
        @test length(scans)              == 6
        @test scans[1].num               == 1
        @test scans[2].level             == 2
        @test scans[3].polarity          == "+"
        @test scans[2].activationMethod  == "CID"
        @test scans[3].collisionEnergy   == 35.0
        @test size(scans[1].int, 1)      == 22320
        
        rt = msJ.retention_time("test.mzXML")
        @test length(rt) == 6
        
        rt, tic = msJ.chromatogram("test.mzXML", method = msJ.TIC() )
        @test length(rt) == 6
        
        rt, tic = msJ.chromatogram("test.mzXML", method = msJ.MZ([0, 500]))
        @test length(rt) == 6
        
        rt, tic = msJ.chromatogram("test.mzXML", method = msJ.∆MZ([1000, 1]))
        @test length(rt) == 6
        
        rt, tic = msJ.chromatogram("test.mzXML", method = msJ.BasePeak())
        @test length(rt) == 6
        
        xrt, xtic = msJ.chromatogram("test.mzXML", msJ.Polarity("+"), msJ.Scan(2),msJ.Precursor(1255.5), msJ.Activation_Energy(18), msJ.Activation_Method("CID"), msJ.Level(2) )
        @test length(xrt) == 1
    end
end
tests()
