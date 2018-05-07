context("Plotting")

test_that("plot 3 views",{

  data(example_SNP)

  PCs <- cal.pc.linear(simsnp$snp, no.pc = 3)
  res <- plot3views( PCs$PC, sample_labels)

  expect_equal(res, NULL)
  file1 <- file.path(tempdir(),"Rplots.pdf")
  file.remove(file1)
})
