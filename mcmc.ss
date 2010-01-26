#lang scheme

(require scheme/flonum)

(provide/contract
 (mcmc-sample
  (-> (-> any/c any)
      (-> any/c inexact-real?)
      (-> any/c any/c inexact-real?)
      any/c
      (-> natural-number/c (listof (cons/c any/c inexact-real?))))))

(define ((mcmc-sample propose log-posterior log-jump-probability start) n)
  (reverse
   (for/fold
       ((samples (list (cons start (log-posterior start)))))
       ((i (in-range n)))
     (let ((last (car samples)))
       (let* ((last-state (car last))
              (last-posterior (cdr last))
              (new-state (propose last-state)))
       (let ((forward (log-jump-probability last-state new-state))
             (backward (log-jump-probability new-state last-state))
             (new-posterior (log-posterior new-state)))
         (let ((log-accept (fl- (fl+ new-posterior backward)
                                (fl+ last-posterior forward))))
           (let ((rand (log (random))))
             (if (<= rand log-accept)
                 (cons (cons new-state new-posterior)
                       samples)
                 (cons (cons last-state last-posterior)
                       samples))))))))))