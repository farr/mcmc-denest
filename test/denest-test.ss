#lang scheme

(require schemeunit
         schemeunit/text-ui
         scheme/flonum
         scheme/fixnum
         "../denest.ss"
         "../flvector-iter.ss"
         "checks.ss"
         "../mcmc.ss")

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

(define (multi-gaussian-sample mu sigma)
  (let ((result (make-flvector (flvector-length mu))))
    (for ((i (in-naturals))
          (mu (in-flvector mu))
          (sigma (in-flvector sigma)))
      (flvector-set! result i (draw-gaussian mu sigma
                                             (fl- mu (fl* sigma 4.0))
                                             (fl+ mu (fl* sigma 4.0)))))
    result))

(define ((multi-gaussian mu sigma) x)
  (for/fold ((prod 1.0))
      ((mu (in-flvector mu))
       (sigma (in-flvector sigma))
       (x (in-flvector x)))
    (fl* prod ((gaussian mu sigma) x))))

(define (in-bounds? x low high)
  (for/and ((x (in-flvector x))
            (l (in-flvector low))
            (h (in-flvector high)))
    (and (fl<= l x)
         (fl< x h))))

(define tests
  (test-suite
   "denest.ss tests"
   (test-case
    "direct integration of multi-gaussian"
    (let ((mu (flvector 0.0 1.0))
          (sigma (flvector 2.0 3.0)))
      (let ((integral (for*/fold ((sum 0.0))
                          ((x (in-range (fl- (flvector-ref mu 0) (fl* (flvector-ref sigma 0) 3.0))
                                        (fl+ (flvector-ref mu 0) (fl* (flvector-ref sigma 0) 3.0))
                                        (fl/ (flvector-ref sigma 0) 100.0)))
                           (y (in-range (fl- (flvector-ref mu 1) (fl* (flvector-ref sigma 1) 3.0))
                                        (fl+ (flvector-ref mu 1) (fl* (flvector-ref sigma 1) 3.0))
                                        (fl/ (flvector-ref sigma 1) 100.0))))
                        (fl+ sum (fl* (fl/ (flvector-ref sigma 0) 100.0)
                                      (fl* (fl/ (flvector-ref sigma 1) 100.0)
                                           ((multi-gaussian mu sigma) (flvector x y))))))))
        (check-close? 1e-2 1e-2 1.0 integral))))                       
   (test-case
    "Integral of Gaussian MCMC in two dimensions"
    (let ((mu (flvector 0.0 1.0))
          (sigma (flvector 2.0 3.0)))
      (let ((proposal (lambda (x)
                        (let ((x-new (make-flvector 2)))
                          (for ((x (in-flvector x))
                                (i (in-naturals))
                                (sigma (in-flvector sigma)))
                            (flvector-set! x-new i (fl+ (fl- x (fl* 0.5 sigma))
                                                        (fl* sigma (random)))))
                          x-new)))
            (log-jump-probability (lambda args 0.0)) ; Don't care about jumps, since symmetric.
            (log-posterior (lambda (x) (log ((multi-gaussian mu sigma) x))))
            (->coords (lambda (mcmc-state) (car mcmc-state))))
        (let ((mcmc-sample-generate (mcmc-sample proposal log-posterior log-jump-probability mu)))
          (let ((samples (for/list ((i (in-range 100000))) (mcmc-sample-generate))))
            (let ((tree (objects->density-tree ->coords car samples)))
              (let ((integral (fold
                               (lambda (tree sum)
                                 (let ((vol (volume tree))
                                       (obj (object-density-tree-obj tree)))
                                   (fl+ sum (fl* vol (exp (cdr obj))))))
                               0.0
                               tree)))
                (check-close? 1e-1 1e-1 integral 1.0))))))))
   (test-case
    "find-object"
    (let ((low (flvector 0.0 0.0))
          (high (flvector 1.0 1.0)))
      (let ((samples (for/list ((i (in-range 10000))) (draw-uniform low high))))
        (let ((tree (objects->density-tree/bounds values car samples low high)))
          (for ((obj (in-list samples)))
            (let ((obj-tree (find-object values obj tree)))
              (check-true (object-density-tree? obj-tree))
              (check-true (in-bounds? obj
                                      (object-density-tree-low obj-tree)
                                      (object-density-tree-high obj-tree)))))))))))