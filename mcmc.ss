#lang scheme

#|  mcmc.ss: Markov-Chain Monte Carlo sampling.
    Copyright (C) 2010 Will M. Farr <wmfarr@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
|#

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
  