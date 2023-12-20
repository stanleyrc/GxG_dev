

setDTthreads(1)

GENOME="human_g1k_v37.no.extra.chrom.sizes"
Sys.setenv(DEFAULT_GENOME = GENOME)
#file locations for testing interchr
## interchr_dat = "testing_interchr/hic_data_test.rds"
## hic_test_file = "testing_interchr/test.hic"
gmatrix_orig_data = "testing_interchr/gmatrix_before_hic.rds"
hic_test_file = "testing_interchr/test.hic"
mcool_test_file = "testing_interchr/test.mcool"
inter_chr_test_file = "testing_interchr/inter_chr_sim.ecDNA.rds"


test_that("test .hic", {
    suppressWarnings({
        hic.file = system.file('extdata', "test.hic", package = "GxG")

        ## can load in specific chromosomes via character or numeric vector
        ## and the other normalizations
        gm = straw(hic.file, 1:2, res = 5e4)
        gm = straw(hic.file, c(1:3, 'X'), norm = "NONE", res = 1e5)

        ## we can use any GRanges as input to straw and use alternate norms
        ## like KR and VC, though the default is NONE
        ## (which are not provided with the small .hic matrix bundled with
        ## the package but are available to most Juicer outputs
        gm = straw(hic.file, GRanges('1:1-250e6'), norm = 'NONE', res = 5e5)
        expect_true(inherits(gm, "gMatrix"))

        ## check that gTrack can be created
        gt = gm$gtrack(colormap = c('white', 'green', 'blue'), clim = c(0, 50))
        expect_true(inherits(gt, "gTrack"))

        ## test disjoin
        new.ranges = gr.tile(gm$footprint, 2e5)
        gmd = gm$disjoin(new.ranges) ###########this breaks
        expect_true(inherits(gm, "gMatrix"))
        expect_true(length(gmd$gr) >= length(new.ranges))

        ## test whether reading in mcool and .hic matches the original data used to create them
        orig_dat = readRDS(gmatrix_orig_data)
        hic = straw(hic = hic_test_file,gr= c("1","5","15"),res = 1e6)
        cool = cooler(file = mcool_test_file,gr= c("1","5","15"),res = 1e6)
        expect_true(all(c(all(hic$dat==cool$dat),all(hic$dat==orig_dat$dat),all(cool$dat==orig_dat$dat))))
        
        ## test new whether interchr coordinates will be in correct format/order
        dat = readRDS(inter_chr_test_file)
        gr.test = parse.gr("5:1-145138636,15:1-107043718")
        hic = straw(hic = hic_test_file,gr= gr.test,res = 1e6)        
        expect_true(all(hic$dat==dat$dat))
    })
})

