#lang scheme

#|  denest.ss: A k-D tree datastructure for partitioning a set of points in R^n.
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

(require "flvector-iter.ss"
         scheme/flonum
         scheme/fixnum
         srfi/67)

(define-struct empty-density-tree
  (low high) #:transparent)

(define-struct object-density-tree
  (low high obj)
  #:transparent)

(define-struct cell-density-tree
  (low high left-tree right-tree)
  #:transparent)

(define (density-tree? obj)
  (or (empty-density-tree? obj)
      (object-density-tree? obj)
      (cell-density-tree? obj)))

(provide/contract
 (struct empty-density-tree
   ((low flvector?)
    (high flvector?)))
 (struct object-density-tree
   ((low flvector?)
    (high flvector?)
    (obj any/c)))
 (struct cell-density-tree
   ((low flvector?)
    (high flvector?)
    (left-tree density-tree?)
    (right-tree density-tree?)))
 (volume (-> density-tree? inexact-real?))
 (objects->density-tree (-> (-> any/c flvector?)
                            (-> (listof any/c) any)
                            (listof any/c)
                            density-tree?))
 (objects->density-tree/bounds (-> (-> any/c flvector?)
                                   (-> (listof any/c) any)
                                   (listof any/c)
                                   flvector?
                                   flvector?
                                   density-tree?))
 (fold (-> (-> object-density-tree? any/c any)
           any/c
           density-tree?
           any))
 (find-object (-> (-> any/c flvector?)
                  any/c
                  density-tree?
                  object-density-tree?)))

(define (bounds-volume low high)
  (for/fold ((prod 1.0))
      ((l (in-flvector low))
       (h (in-flvector high)))
    (fl* prod (fl- h l))))

(define (volume t)
  (cond
   ((empty-density-tree? t) (bounds-volume (empty-density-tree-low t)
                                           (empty-density-tree-high t)))
   ((object-density-tree? t) (bounds-volume (object-density-tree-low t)
                                            (object-density-tree-high t)))
   ((cell-density-tree? t) (bounds-volume (cell-density-tree-low t)
                                          (cell-density-tree-high t)))))

(define (longest-dimension low high)
  (let-values (((i dx)
                (for/fold ((dim -1) (dx-max -inf.0))
                    ((l (in-flvector low))
                     (h (in-flvector high))
                     (i (in-naturals)))
                  (let ((dx (fl- h l)))
                    (if (fl> dx dx-max)
                        (values i dx)
                        (values dim dx-max))))))
    i))

(define (expand-bounds! low high)
  (for ((l (in-flvector low))
        (h (in-flvector high))
        (i (in-naturals)))
    (flvector-set! high i (fl+ h (fl* 1e-8 (fl- h l))))))

(define (in-bounds? x low high)
  (for/and ((x (in-flvector x))
            (l (in-flvector low))
            (h (in-flvector high)))
    (and (fl<= l x)
         (fl< x h))))

(define (in-tree? x t)
  (cond
   ((empty-density-tree? t) (in-bounds? x
                                        (empty-density-tree-low t)
                                        (empty-density-tree-high t)))
   ((object-density-tree? t) (in-bounds? x
                                         (object-density-tree-low t)
                                         (object-density-tree-high t)))
   ((cell-density-tree? t) (in-bounds? x
                                       (cell-density-tree-low t)
                                       (cell-density-tree-high t)))))

(define (objects->bounds ->coords objs)
  (when (null? objs)
    (error 'objects->bounds "no objects"))
  (let* ((coords (->coords (car objs)))
         (ndim (flvector-length coords)))
    (let ((low (make-flvector ndim +inf.0))
          (high (make-flvector ndim -inf.0)))
      (for ((obj (in-list objs)))
        (let ((x (->coords obj)))
          (for ((x (in-flvector x))
                (l (in-flvector low))
                (h (in-flvector high))
                (i (in-naturals)))
            (when (fl< x l)
              (flvector-set! low i x))
            (when (fl> x h)
              (flvector-set! high i x)))))
      (expand-bounds! low high)
      (values low high))))

(define (remove-duplicates/sorted compare list)
  (let loop ((removed '()) (remaining list))
    (cond
     ((null? remaining) (reverse removed))
     ((null? (cdr remaining)) (reverse (cons (car remaining) removed)))
     ((=? compare (car remaining) (cadr remaining))
      (loop removed (cdr remaining)))
     (else
      (loop (cons (car remaining) removed) (cdr remaining))))))

(define (fl-compare x y)
  (cond
   ((fl< x y) -1)
   ((fl> x y) 1)
   (else 0)))

(define (find-split compare objs)
  (let ((sorted-unique-objs (remove-duplicates/sorted compare (sort objs (lambda (o1 o2) (<? compare o1 o2))))))
    (cond
     ((null? sorted-unique-objs) (error 'find-split "no objects"))
     ((null? (cdr sorted-unique-objs)) (error 'find-split "only one object"))
     (else ; At least two unique objects
      (let* ((N/2 (floor (/ (length sorted-unique-objs) 2)))
             (N/2-1 (sub1 N/2)))
        (with-handlers (((lambda (exn) #t)
                         (lambda (exn)
                           (printf "N/2 = ~a~%N/1-1 = ~a~%" N/2 N/2-1)
                           (raise exn))))
          (let-values (((o-left o-right)
                        (let loop ((i 0) (list sorted-unique-objs))
                          (if (fx= i N/2-1)
                              (values (car list) (cadr list))
                              (loop (add1 i) (cdr list))))))
            (values o-left o-right))))))))

(define (sub-bounds ->coords dim o-left o-right low high)
  (let ((N (flvector-length low)))
    (when (not (fx= N (flvector-length high)))
      (error 'sub-bounds "un-equal lengths in initial bounds: ~a vs ~a" N (flvector-length high)))
    (let ((new-low (make-flvector N))
          (new-high (make-flvector N)))
      (for ((i (in-naturals))
            (l (in-flvector low))
            (h (in-flvector high)))
        (flvector-set! new-low i l)
        (flvector-set! new-high i h))
      (let ((mid (fl* 0.5 (fl+ (flvector-ref (->coords o-left) dim)
                               (flvector-ref (->coords o-right) dim)))))
        (flvector-set! new-low dim mid)
        (flvector-set! new-high dim mid)
        (values new-low new-high)))))

(define (flvector=? v1 v2)
  (let ((N (flvector-length v1)))
    (and (fx= N (flvector-length v2))
         (for/and ((x (in-flvector v1))
                   (y (in-flvector v2)))
           (fl= x y)))))

(define (objects->density-tree/bounds ->coords combiner objs low high)
  (cond
   ((null? objs) (make-empty-density-tree low high))
   ((null? (cdr objs)) (make-object-density-tree low high (car objs)))
   ((for/and ((obj (in-list (cdr objs)))) (flvector=? (->coords obj) (->coords (car objs))))
    (make-object-density-tree low high (combiner objs)))
   (else ;; At least two different coordinates
    (let ((dim (call-with-values (lambda () (objects->bounds ->coords objs)) longest-dimension)))
      (define (compare o1 o2)
        (fl-compare (flvector-ref (->coords o1) dim) (flvector-ref (->coords o2) dim)))
      (let-values (((new-low new-high)
                    (call-with-values (lambda () (find-split compare objs))
                      (lambda (o-left o-right)
                        (sub-bounds ->coords dim o-left o-right low high)))))
        (let-values (((left right)
                      (partition (lambda (obj) (in-bounds? (->coords obj) low new-high)) objs)))
          (make-cell-density-tree low high
                                  (objects->density-tree/bounds ->coords combiner left low new-high)
                                  (objects->density-tree/bounds ->coords combiner right new-low high))))))))

(define (objects->density-tree ->coords combiner objs)
  (call-with-values
      (lambda () (objects->bounds ->coords objs))
    (lambda (low high)
      (objects->density-tree/bounds ->coords combiner objs low high))))

(define (fold kons knil tree)
  (cond
   ((empty-density-tree? tree) knil)
   ((object-density-tree? tree) (kons tree knil))
   ((cell-density-tree? tree)
    (let ((result (fold kons knil (cell-density-tree-left-tree tree))))
      (fold kons result (cell-density-tree-right-tree tree))))))

(define (find-object ->coords obj tree)
  (let ((x (->coords obj)))
    (when (not (in-tree? x tree))
      (error 'find-object "object not contained in tree: ~a" obj))
    (let loop ((tree tree))
      (cond
       ((empty-density-tree? tree) (error 'find-object "internal invariant violated: empty tree"))
       ((object-density-tree? tree) tree)
       ((cell-density-tree? tree)
        (let ((left (cell-density-tree-left-tree tree)))
          (if (in-tree? x left)
              (loop left)
              (loop (cell-density-tree-right-tree tree)))))))))