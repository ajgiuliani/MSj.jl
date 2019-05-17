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
        rt, tic = msJ.chromatogram("test.mzXML" )
        @test length(rt) == 6
        rt, tic = msJ.chromatogram("test.mzXML", method = msJ.MZ([0, 500]))
        @test length(rt) == 6
        rt, tic = msJ.chromatogram("test.mzXML", method = msJ.∆MZ([1000, 1]))
        @test length(rt) == 6
        rt, tic = msJ.chromatogram("test.mzXML", method = msJ.BasePeak())
        @test length(rt) == 6
        xrt4, xtic4 = msJ.chromatogram("test.mzXML", msJ.Precursor(1255.5))
        xrt5, xtic5 = msJ.chromatogram("test.mzXML", msJ.Activation_Energy(18))
        xrt6, xtic6 = msJ.chromatogram("test.mzXML", msJ.Activation_Method("CID"))
        xrt7, xtic7 = msJ.chromatogram("test.mzXML", msJ.Level(2));
        @test (xrt4 == xrt5 == xrt6 == xrt7) == true

        rt = msJ.retention_time(scans)
        @test length(rt) == 6
        rt, tic = msJ.chromatogram(scans)
        @test length(rt) == 6
        rt, tic = msJ.chromatogram(scans, method = msJ.MZ([0, 500]))
        @test length(rt) == 6
        rt, tic = msJ.chromatogram(scans, method = msJ.∆MZ([1000, 1]))
        @test length(rt) == 6
        rt, tic = msJ.chromatogram(scans, method = msJ.BasePeak())
        @test length(rt) == 6
        xrt4, xtic4 = msJ.chromatogram(scans, msJ.Precursor(1255.5))
        xrt5, xtic5 = msJ.chromatogram(scans, msJ.Activation_Energy(18))
        xrt6, xtic6 = msJ.chromatogram(scans, msJ.Activation_Method("CID"))
        xrt7, xtic7 = msJ.chromatogram(scans, msJ.Level(2));
        @test (xrt4 == xrt5 == xrt6 == xrt7) == true
     
        ms1 = msJ.msfilter("test.mzXML")
        ms2 = msJ.msfilter("test.mzXML", msJ.Precursor(1255.5))
        ms3 = msJ.msfilter("test.mzXML", msJ.Activation_Energy(18))
        ms4 = msJ.msfilter("test.mzXML", msJ.Activation_Method("CID"))
        ms5 = msJ.msfilter("test.mzXML", msJ.Level(2));
        @test (ms2.int == ms3.int == ms4.int == ms5.int) == true

        ms6 = msJ.msfilter(scans)
        @test (ms1.int == ms6.int) == true
        ms7 = msJ.msfilter(scans, msJ.Precursor(1255.5))
        ms8 = msJ.msfilter(scans, msJ.Activation_Energy(18))
        ms9 = msJ.msfilter(scans, msJ.Activation_Method("CID"))
        ms10 = msJ.msfilter(scans, msJ.Level(2));
        @test (ms7.int == ms8.int == ms9.int == ms10.int) == true

    end
end

tests()




