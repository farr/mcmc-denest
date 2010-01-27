#lang scheme

(require schemeunit
         schemeunit/text-ui
         (prefix-in mcmc: "mcmc-test.ss"))

(define tests
  (test-suite
   "all tests"
   mcmc:tests))

(exit (run-tests tests))