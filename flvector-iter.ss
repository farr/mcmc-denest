#lang scheme

(require scheme/flonum
         scheme/fixnum)

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
                       ((i 0) (N (flvector-length the-flvector)))
                       (fx< i N)
                       (((id) (flvector-ref the-flvector i)))
                       #t
                       #t
                       ((add1 i) N))))))))