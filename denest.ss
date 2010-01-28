#lang scheme

(require "flvector-iter.ss"
         scheme/flonum
         scheme/fixnum
         scheme/unsafe/ops)

(define-struct cell
  (low high split-dim left-tree right-tree)
  #:transparent)

(define (density-tree? obj)
  (or (eq? obj #f)
      (cell? obj)))

(provide/contract
 (struct cell ((low flvector?)
               (high flvector?)
               (split-dim natural-number/c)
               (left-tree density-tree?)
               (right-tree density-tree?)))
 (density-tree? (-> any/c boolean?))
 (empty? (-> density-tree? boolean?))
 (points->density-tree (-> (listof flvector?) density-tree?))
 (points->density-tree/bounds (-> (listof flvector?) flvector? flvector? density-tree?))
 (volume (-> density-tree? inexact-real?))
 (find-cell (-> flvector? density-tree? cell?)))

(define (empty? obj) (not obj))

(define (split-bounds low high dim)
  (let ((n (flvector-length low)))
    (let ((new-low (make-flvector n))
          (new-high (make-flvector n)))
    (for ((i (in-range n))
          (l (in-flvector low))
          (h (in-flvector high)))
      (if (fx= i dim)
          (let ((mid (fl* 0.5 (fl+ l h))))
            (unsafe-flvector-set! new-high i mid)
            (unsafe-flvector-set! new-low i mid))
          (begin
            (unsafe-flvector-set! new-low i l)
            (unsafe-flvector-set! new-high i h))))
    (values new-low new-high))))

(define (in-cell? x cell)
  (for/and ((x (in-flvector x))
            (low (in-flvector (cell-low cell)))
            (high (in-flvector (cell-high cell))))
    (and (fl<= low x)
         (fl< x high))))

(define (longest-dimension low high)
  (let-values (((dim extent)
                (for/fold ((dim -1) (extent -inf.0))
                    ((i (in-naturals))
                     (l (in-flvector low))
                     (h (in-flvector high)))
                  (let ((d (unsafe-fl- h l)))
                    (if (unsafe-fl> d extent)
                        (values i d)
                        (values dim extent))))))
    dim))

(define (expand-bounds! low high pt)
  (for ((i (in-naturals))
        (l (in-flvector low))
        (h (in-flvector high))
        (x (in-flvector pt)))
    (when (unsafe-fl< x l)
      (unsafe-flvector-set! low i x))
    (when (unsafe-fl>= x h)
      (unsafe-flvector-set! high i (unsafe-fl+ x (unsafe-fl* 1e-8 (unsafe-fl- x l)))))))

(define (initial-bounds pt)
  (let ((n (flvector-length pt)))
    (let ((low (make-flvector n))
          (high (make-flvector n)))
      (for ((i (in-range n))
            (x (in-flvector pt)))
        (unsafe-flvector-set! low i x)
        (unsafe-flvector-set! high i (unsafe-fl+ x (unsafe-fl* 1e-8 (unsafe-flabs x)))))
      (values low high))))

(define (points->bounds pts)
  (when (null? pts)
    (error 'points->bounds "expected non-null list of points, got " pts))
  (let ((n (flvector-length (car pts))))
    (let-values (((low high) (initial-bounds (car pts))))
      (for/fold ((low low)
                 (high high))
          ((pt (in-list (cdr pts))))
        (expand-bounds! low high pt)
        (values low high)))))

(define (in-bounds? pt low high)
  (for/and ((x (in-flvector pt))
            (l (in-flvector low))
            (h (in-flvector high)))
    (and (unsafe-fl<= l x)
         (unsafe-fl< x h))))

(define (points->density-tree pts)
  (if (null? pts)
      #f
      (let-values (((low high)
                    (points->bounds pts)))
        (points->density-tree/bounds pts low high))))

(define (points->density-tree/bounds pts low high)
  (cond
   ((null? pts) #f)
   ((null? (cdr pts)) (make-cell low high -1 #f #f))
   (else
    (let ((dim (longest-dimension low high)))
      (let-values (((new-low new-high) (split-bounds low high dim)))
        (let-values (((left-pts right-pts) (partition (lambda (pt) (in-bounds? pt low new-high)) pts)))
          (make-cell low high dim
                     (points->density-tree/bounds left-pts low new-high)
                     (points->density-tree/bounds right-pts new-low high))))))))

(define (volume cell)
  (if (empty? cell)
      0.0
      (for/fold ((v 1.0))
          ((l (in-flvector cell-low))
           (h (in-flvector cell-high)))
        (unsafe-fl* v (unsafe-fl- h l)))))

(define (find-cell pt tree)
  (when (empty? tree)
    (error 'find-cell "expected non-empty density-tree, got " tree))
  (let ((left (cell-left-tree tree))
        (right (cell-right-tree tree)))
  (if (and (empty? left)
           (empty? right))
      tree
      (let ((low (cell-low tree))
            (high (cell-high tree))
            (dim (cell-split-dim tree)))
        (let ((mid (unsafe-fl* 0.5 (unsafe-fl+ (unsafe-flvector-ref low dim)
                                               (unsafe-flvector-ref high dim)))))
          (if (unsafe-fl< (unsafe-flvector-ref pt dim) mid)
              (find-cell pt left)
              (find-cell pt right)))))))
                   