;;;;;; Zhaorong Wang, 2010
;;;;;; In this case, I put the source far away from the grating so
;;;;;; that it becomes plane wave
;;;;;; upon the grating. And the source is also away from the pml so
;;;;;; that it will not disturb the flux collecting.

(define-param sx 5.85)  ; Size of cell in X direction
(define-param sy 50)    ; Size of cell in Y direction
(define-param bottom (* -0.5 sy)); Y coordinate of bottom of the cell
(define-param top (* 0.5 sy))   ; Y coordinate of top of the cell
(define-param left (* -0.5 sx)) ; X coordinate of left of the cell
(define-param right (* 0.5 sx)) ; X coordinate of right of the cell
(define-param pad 8)    ;pad thickness from the top to the source
(define-param lamb 8.0) ; wavelength in vacuum

(set! geometry-lattice (make lattice (size sx sy no-size)))

(define-param h 1.25) ; Height of the grating layer
(define-param w 2.0475) ; width of a grating unit

(define-param no-grating false) ; if true, have no grating structure

(set! geometry
     (if no-grating
         (list )
         (list
          (make block
            (center 0 (+ bottom (* 2 lamb)))
            (size w h infinity)
            (material (make dielectric (epsilon 9)))))))

(set! k-point (vector3 0 0 0))
(set! ensure-periodicity true)

(define-param fcen 0.125) ; pulse center frequency
(define-param df 0.07)    ; pulse width (in frequency)
(set! sources (list
              (make source
                (src (make gaussian-src (frequency fcen) (fwidth df)))
                (component Ez)
                (center 0 (- top pad)) ; pad from the top
)))

(set! pml-layers (list (make pml (thickness 1.0)(direction Y))))
(set-param! resolution 10)

(define-param nfreq 100) ; number of frequencies at which to compute flux
(define trans   ; transmitted flux
     (add-flux fcen df nfreq
                (make flux-region
                  (center 0 (+ bottom (* 1.5 lamb) )) (size w 0))));
0.5 lambda from grating
(define refl    ; reflected flux
     (add-flux fcen df nfreq
                (make flux-region
                  (center 0 (+ bottom (* 2.5 lamb) )) (size w 0))));
0.5 lambda from grating

(if (not no-grating) (load-minus-flux "refl-flux" refl))
(run-sources+
       (stop-when-fields-decayed 50 Ez (vector3 0 (+ bottom (* 1.5 lamb))) 1e-3)
       (at-beginning output-epsilon)
       (to-appended "ez" (at-every 1 output-efield-z)))
(if no-grating (save-flux "refl-flux" refl))

(display-fluxes trans refl)