#lang scheme

(require schemeunit
         schemeunit/text-ui
         scheme/flonum
         "../mcmc.ss")

(define ((log-gaussian mu sigma) x)
  (let ((dx (fl- x mu)))
    (fl- (fl/ (fl* dx dx)
              (fl* -2.0 (fl* sigma sigma)))
         (log sigma))))

(define (random-around x dx)
  (let ((d (fl* 0.5 dx)))
    (fl+ (fl- x d)
         (fl* dx (random)))))

(define (random-around/positive x dx)
  (let ((d (fl* 0.5 dx)))
    (if (fl< x d)
        (fl* dx (random))
        (fl+ (fl- x d)
             (fl* dx (random))))))

(define (biased x dx p-left)
  (let ((left? (< (random) p-left)))
    (if left?
        (fl- x (fl* dx (random)))
        (fl+ x (fl* dx (random))))))

(define-simple-check (check-close? epsabs epsrel x y)
  (let ((dx (abs (- x y))))
    (<= dx (+ epsabs (* epsrel 0.5 (+ (abs x) (abs y)))))))

(define tests
  (test-suite
   "mcmc.ss tests"
   (test-case
    "gaussian posterior"
    (let ((mu-true (random))
          (sigma-true (random))
          (N 1000000)
          (p-left 0.75))
      #;(printf "mu = ~a, sigma = ~a~%" mu-true sigma-true)
      (let ((dx (/ sigma-true 10)))
        (let ((posterior (lambda (x) ((log-gaussian mu-true sigma-true) x)))
              (propose (lambda (x)
                         (biased x dx p-left)))
              (jump-prob (lambda (source target)
                           (if (fl< source target)
                               (log (fl* (fl- 1.0 p-left) (fl/ 1.0 dx)))
                               (log (fl* p-left (fl/ 1.0 dx)))))))
          (let ((samples ((mcmc-sample propose posterior jump-prob mu-true) N)))
            #;(with-output-to-file "mcmc-output.dat"
              (lambda ()
                (for-each
                 (lambda (sample)
                   (printf "~a ~a~%" (car sample) (cdr sample)))
                 samples))
              #:exists 'replace)
            (let* ((sample-x (for/fold ((x-sum 0.0)) ((sample (in-list samples)))
                               (let ((x (car sample)))
                                 (fl+ x-sum (/ x N)))))
                   (sample-sigma (flsqrt 
                                  (for/fold ((s-sum 0.0)) ((sample (in-list samples)))
                                    (let ((x (car sample)))
                                      (let ((dx (fl- x sample-x)))
                                        (fl+ s-sum (/ (fl* dx dx) N))))))))
              (check-close? 1e-2 1e-2 sample-x mu-true)
              (check-close? 1e-2 1e-2 sample-sigma sigma-true)))))))))
            
                          