#lang scheme

(require schemeunit
         schemeunit/text-ui
         (prefix-in denest: "denest-test.ss")
         (prefix-in mcmc: "mcmc-test.ss"))

(define tests
  (test-suite
   "all tests"
   denest:tests
   mcmc:tests))

(exit (run-tests tests))