#lang scheme

(require schemeunit
         schemeunit/text-ui
         scheme/flonum
         scheme/fixnum
         "../denest.ss"
         "../flvector-iter.ss"
         "checks.ss")

(provide tests)

(define ((gaussian mu sigma) x)
  (let ((dx (fl- x mu))
        (sigma2 (fl* sigma sigma)))
    (fl/ (flexp (fl/ (fl* dx dx) (fl* -2.0 sigma2)))
         (fl* (flsqrt (fl* 2.0 pi)) sigma))))

(define (draw-gaussian mu sigma low high)
  (let ((ymax (fl/ 1.0 (fl* (flsqrt (fl* 2.0 pi)) sigma))))
    (let ((x (fl+ low (fl* (fl- high low) (random))))
          (y (fl* ymax (random))))
      (if (fl>= ((gaussian mu sigma) x) y)
          x
          (draw-gaussian mu sigma low high)))))

(define (draw-uniform low high)
  (let* ((n (flvector-length low))
         (result (make-flvector n)))
    (for ((i (in-range n))
          (l (in-flvector low))
          (h (in-flvector high)))
      (flvector-set! result i (fl+ l (fl* (random) (fl- h l)))))
    result))

(define (mean samples)
  (let ((n (length samples)))
    (when (fx= n 0)
      (error 'mean "no samples"))
    (for/fold ((mu (make-flvector (flvector-length (car samples)))))
        ((samp (in-list samples)))
      (for ((i (in-naturals))
            (samp-x (in-flvector samp))
            (mu-x (in-flvector mu)))
        (flvector-set! mu i (fl+ mu-x (/ samp-x n))))
      mu)))

(define (variance samples (mu (mean samples)))
  (let ((n (length samples)))
    (when (fx= n 0)
      (error 'std-dev "no samples"))
    (for/fold ((sigma (make-flvector (flvector-length (car samples)))))
        ((samp (in-list samples)))
      (for ((i (in-naturals))
            (sigma-x (in-flvector sigma))
            (samp-x (in-flvector samp))
            (mu-x (in-flvector mu)))
        (let ((dx (fl- samp-x mu-x)))
          (flvector-set! sigma i (fl+ sigma-x (/ (fl* dx dx) n)))))
      sigma)))

(define tests
  (test-suite
   "denest.ss tests"
   (test-case
    "Re-sampling in two-dimensions."
    (let ((mu0 1.0)
          (mu1 -1.0)
          (sigma0 2.0)
          (sigma1 0.5)
          (N 1000000))
      (let ((samples (for/list ((i (in-range N)))
                       (flvector
                        (draw-gaussian mu0 sigma0
                                       (fl- mu0 (fl* 4.0 sigma0))
                                       (fl+ mu0 (fl* 4.0 sigma0)))
                        (draw-gaussian mu1 sigma1
                                       (fl- mu1 (fl* 4.0 sigma1))
                                       (fl+ mu1 (fl* 4.0 sigma1)))))))
        (let ((tree (points->density-tree samples)))
          (let ((re-samples
                 (for/list ((samp (in-list samples)))
                   (let ((cell (find-cell samp tree)))
                     (draw-uniform (cell-low cell) (cell-high cell))))))
            (let ((mu (mean samples))
                  (re-mu (mean re-samples))
                  (sigma (variance samples))
                  (re-sigma (variance re-samples)))
              (for ((mu (in-flvector mu))
                    (rmu (in-flvector re-mu))
                    (sigma (in-flvector sigma))
                    (rsigma (in-flvector re-sigma)))
                (check-close? 1e-2 1e-2 mu rmu)
                (check-close? 1e-2 1e-2 sigma rsigma))))))))))