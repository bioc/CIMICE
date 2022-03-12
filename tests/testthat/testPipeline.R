#' Simple test of the CIMICE pipeline
#'
#' Tests the whole pipeline and text (dot) output on a simple test case
#'
#' @return TRUE if the test is successful
#'
#' @examples
#' require(dplyr)
#' require(purrr)
#' test.pipeline()
#'
test.pipeline <- function(){
    require(dplyr)
    require(purrr)
    ans <- to_dot(quick_run(example_dataset()), mode = "geneIDs")
    precom_ans <- "digraph G{\n1[label=\"D\"]\n2[label=\"A\"]\n3[label=\"A, D\"]\n4[label=\"A, C, D\"]\n5[label=\"A, B, D\"]\n6[label=\"Clonal\"]\n6 -> 1[label=\"0.333\"]\n6 -> 2[label=\"0.667\"]\n1 -> 3[label=\"1\"]\n2 -> 3[label=\"1\"]\n3 -> 4[label=\"0.25\"]\n3 -> 5[label=\"0.75\"]\n}"
    test_that("full pipeline", {
        expect_equal(ans, precom_ans)
    })
}

test.pipeline()
