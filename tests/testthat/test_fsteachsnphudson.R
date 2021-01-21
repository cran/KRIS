context("Fst Hudson all SNPs")

test_that("calculate Fst Hudson for all SNPs",{

  data(example_SNP)

  idx1 <- which(sample_labels == 'pop1')
  idx2 <- which(sample_labels == 'pop2')
  fst.pairwise <- fst.each.snp.hudson(simsnp$snp, idx1, idx2)

  expect_equal(fst.pairwise[1], 0.007396462)
  expect_equal(fst.pairwise[2], 0.09957484)
  expect_length(fst.pairwise, 2000)
})
