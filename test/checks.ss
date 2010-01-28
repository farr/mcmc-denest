#lang scheme

(require schemeunit)

(provide check-close?)

(define-simple-check (check-close? epsabs epsrel x y)
  (let ((dx (abs (- x y))))
    (<= dx (+ epsabs (* epsrel 0.5 (+ (abs x) (abs y)))))))