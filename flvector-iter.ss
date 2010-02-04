#lang scheme

#|  flvector-iter.ss: Iterators over flvectors.
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

(require scheme/flonum
         scheme/fixnum
         scheme/unsafe/ops)

(provide in-flvector)

(define (*in-flvector v)
  (make-do-sequence
   (lambda ()
     (values
      (lambda (i) (flvector-ref v i))
      add1
      0
      (lambda (i) (fx< i (flvector-length v)))
      (lambda (x) #t)
      (lambda (x i) #t)))))

(define-sequence-syntax in-flvector
  (lambda (stx)
    (syntax-case stx ()
      ((in-flvector body ...)
       (syntax/loc stx (*in-flvector body ...)))))
  (lambda (stx)
    (syntax-case stx ()
      (((id) (in-flvector flvector-expr))
       (syntax/loc stx
         ((id) (:do-in (((the-flvector) flvector-expr))
                       (when (not (flvector? the-flvector))
                         (error 'in-flvector "expected flvector, got: " the-flvector))
                       ((i 0) (N (unsafe-flvector-length the-flvector)))
                       (unsafe-fx< i N)
                       (((id) (unsafe-flvector-ref the-flvector i)))
                       #t
                       #t
                       ((unsafe-fx+ 1 i) N))))))))