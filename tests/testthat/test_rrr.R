context("Low rank regression")

test_that("Testing if MOSS gives same results than matrix general LM", {
  set.seed(42)
  X <- matrix(stats::rnorm(1000), 10, 100)
  y <- X[, 1:5] %*% rep(50, 5) + stats::rnorm(10, sd = 0.05)

  b.prod <- crossprod(MASS::ginv(crossprod(X)), crossprod(X, y))

  b.moss <- moss(
    data.blocks = list(y, X), scale.arg = F,
    resp.block = 1,
    method = "lrr",
    K.X = 10,
    K.Y = 1
  )$B

  expect_equal(cor(b.prod, b.moss), matrix(1, 1, 1))
})



test_that("Testing if MOSS gives same results than matrix general LM 
          (multivariate multiple regression)", {
  sim_data <- simulate_data()
  sim_blocks <- sim_data$sim_blocks

  X <- sim_blocks$`Block 1`
  y <- sim_blocks$`Block 3`

  b.prod <- crossprod(MASS::ginv(crossprod(X)), crossprod(X, y))

  b.moss <- moss(
    data.blocks = list(y, X), scale.arg = F, norm.arg = F,
    resp.block = 1,
    method = "lrr",
    K.X = min(dim(X)),
    K.Y = min(dim(X), dim(y))
  )$B
  expect_equal(cor(b.prod[, 1], b.moss[, 1]), 1)
})
