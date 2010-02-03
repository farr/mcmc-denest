#lang scheme

(require scheme/flonum)

(provide/contract
 (make-mcmc-sampler
  (-> (-> any/c any)
      (-> any/c inexact-real?)
      (-> any/c any/c inexact-real?)
      (-> (cons/c any/c inexact-real?) (cons/c any/c inexact-real?))))
 (in-mcmc
  (-> (-> any/c any)
      (-> any/c inexact-real?)
      (-> any/c any/c inexact-real?)
      any/c
      sequence?)))

(define (make-mcmc-sampler propose log-posterior log-jump-probability)
  (lambda (last)
    (let* ((last-state (car last))
           (last-posterior (cdr last))
           (new-state (propose last-state)))
      (let ((forward (log-jump-probability last-state new-state))
            (backward (log-jump-probability new-state last-state))
            (new-posterior (log-posterior new-state)))
        (let ((log-accept (fl- (fl+ new-posterior backward)
                               (fl+ last-posterior forward))))
          (let ((rand (log (random))))
            (if (fl<= rand log-accept)
                (cons new-state new-posterior)
                last)))))))

(define (in-mcmc propose log-posterior log-jump-probability start)
  (let ((mcmc-sampler (make-mcmc-sampler propose log-posterior log-jump-probability)))
    (make-do-sequence
     (lambda ()
       (values values
               mcmc-sampler
               start
               (lambda args #t)
               (lambda args #t)
               (lambda args #t))))))
  