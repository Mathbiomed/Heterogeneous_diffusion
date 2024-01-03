breed [ cells cell ]
breed [ nuclei nucleus ]

breed [ obstacles obstacle ]
breed [ mRNArepressors mRNArepressor ]
breed [ repressors repressor ]
breed [ phosrepressors phosrepressor ]

breed [ frepressors frepressor ]
breed [ fphosrepressors fphosrepressor ]


globals [

  collision_num
  counting

  ;; Global variables, which are used to set the size of the cell and that of the nucleus
  cellsize
  nuclradi

  ;; Global variables about the number of activators in the local regions for the nucleus (Fig S9A).
  area1
  area2
  anum1
  anum2
  acon

  ;; Global variables, which are used to calcualte the local concentrations of hypophos. PER and hyperphos. PER.


  ;; Global variable, which is used to determine how the cytoplasmic obstacles decrease the speed of the PER protein.
  ;; E.g. slow=0.5 indicates that PER protein speed halved.
  slow
  ;; Global variable, which describe the speed of Per mRNA, that of PER protein and that of obstacle
  speedp
  speedo

  ;; the radius of proteins
;  phosrepressor_radius
  repressor_radius
  ;;
  radius_1st
  radius_2nd
  radius_3rd
  radius_4th
  radius_gap
  ;; output number - the number of proteins in ER region, and output_file_name is for csv file


  protein_num_in_ER
  collision_num_list
  output_file_name
  output_file_name2

  protein-xcor
  protein-ycor

  near_obstacles ;; list of obstacles (agents)
  KK
  maxOverlapObs
  theta_dir
  theta
  ind

;  f_binding
;  f_unbinding

  protein_xcor
  protein_ycor

]

turtles-own [
  stay-time
  is_stuck
  neighbor_er
]

obstacles-own [
  fragment_num
]

extensions [csv matrix]

to load

  clear-all

end

to setup
  set KK 0
  random-seed 31415
  clear-ticks
  clear-turtles
  clear-patches
  clear-drawing
  clear-all-plots
  clear-output

  draw-box

  set counting 0;
  set collision_num 0;
  set cellsize 200
  set nuclradi (cellsize / 6);
;  set phosrepressor_radius 2 / 6 ; small
;  set phosrepressor_radius 0.4
  set repressor_radius 1
  set radius_1st 40
  set radius_2nd 40 + 0.5
  set radius_3rd 40 + 0.5 * 2
  set radius_gap 0.25
;  set radius_4th 40 + 0.5 * 3


  set collision_num 0
  add1 cells 1
  add1 nuclei 1




  ;  add3 mRNArepressors M
  add_obs obstacles 1
  add2 phosrepressors P
;  add2 repressors 0

  set slow 1
;  set speedp 2 / 6 * sqrt(4 / 6 / phosrepressor_radius)
  set speedp 0.016 * sqrt(0.2 / ( phosrepressor_radius / 20))
  set speedo 0

  set protein_num_in_ER []
  set collision_num_list []
;  ask repressors [
;    avoid_obstacles
;  ]



  reset-ticks
end



to go



  ask phosrepressors [
    re_active_transport_binding

  ]

  if (ticks mod 31250 = 0)[
    set output_file_name (word "pic_in_ER_v8_b_" (f_binding * 10) "_u_" (f_unbinding * 100) "_obs_" obstacle_num "_rad_" (100 * ( phosrepressor_radius / 20) ) ticks "_1x_unif.png")
    export-view output_file_name
  ]

  if (ticks mod 50000 = 0)[

    set protein_num_in_ER lput (list count phosrepressors with [(distancexy 0 0) < (nuclradi + 20 + phosrepressor_radius / 2) ]) protein_num_in_ER
    set collision_num_list lput (list collision_num) collision_num_list
    set protein-xcor []
    set protein-ycor []
    get-protein-xcor

;    get-protein-ycor

    set output_file_name (word "xcor_in_ER_v8_b_" (f_binding * 10) "_u_" (f_unbinding * 100) "_obs_" obstacle_num "_rad_" (100 * phosrepressor_radius) ticks "_1x_unif.csv")
    csv:to-file output_file_name protein-xcor
    set output_file_name (word "ycor_in_ER_v8_b_" (f_binding * 10) "_u_" (f_unbinding * 100) "_obs_" obstacle_num "_rad_" (100 * phosrepressor_radius) ticks "_1x_unif.csv")
    csv:to-file output_file_name protein-ycor


  ]


  if (ticks = 5000000)[
    set output_file_name (word "protein_proportion_in_ER_v8_b_" (f_binding * 10) "_u_" (f_unbinding * 100) "_obs_" obstacle_num "_rad_" (100 * phosrepressor_radius) "_1x_unif.csv")
    set output_file_name2 (word "collision_in_ER_v8_b_" (f_binding * 10) "_u_" (f_unbinding * 100) "_obs_" obstacle_num "_rad_" (100 * phosrepressor_radius) "_1x_unif.csv")
    csv:to-file output_file_name protein_num_in_ER
    csv:to-file output_file_name2 collision_num_list

    stop
  ]
  tick
end



to re_active_transport_binding


  let rr -100
  let x1 xcor
  let y1 ycor

;  set f_binding 10
;  set f_unbinding 0.02

  rt random-float 360

  ifelse (is_stuck = 0)[
    fd speedp

    if ((distancexy 0 0) < (nuclradi + 20))[
      if (any? obstacles with [overlap myself > 0])[
        set near_obstacles (obstacles with [distance myself < 0.03])
        if(any? near_obstacles with [overlap myself > 0.01])[
          rt 180
          while [any? near_obstacles with [overlap myself > 0.0]][fd 0.0005]
        ]
         if (random-float 1 > (e ^ (- (f_binding * (1))))) [
        set collision_num (collision_num + 1)
        set is_stuck 1]]
      ]



;        set color [255 54 66]]]
  ][
    if (random-float 1 > (e ^ (- (f_unbinding)))) [
;      set color [0 122 255]
      set is_stuck 0
      fd speedp]

  ]




  let x2 xcor
  let y2 ycor
;  print x2
;  print y2

  let flag1 0
  let flag2 0
  let flag3 0


  let r1 sqrt(x1 * x1 + y1 * y1)
  let r2 sqrt(x2 * x2 + y2 * y2)

  let A 0
  let B 0
  let C 0
  let t -500000


  let x_t 0
  let y_t 0

  ;; Reflecting boundary condition in cell membrane and nuclear membrane
  let xt -100
  let yt -100
  let radiii -100

  if sqrt (x2 * x2 + y2 * y2) - size / 2 < cellsize / 6
    [
      set radiii cellsize / 6 + size / 2
      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 6 + size / 2) ^ 2)

      let temp_discriminant (bv ^ 2 - 4 * av * cv)
      if abs(temp_discriminant) < 1e-7 [set temp_discriminant 0]
      let t1 (- bv + sqrt (temp_discriminant)) / (2 * av)
      let t2 (- bv - sqrt (temp_discriminant)) / (2 * av)



      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
          set xt (t1 * x1 + (1 - t1) * x2)
          set yt (t1 * y1 + (1 - t1) * y2)


      ]
      [
        set xt (t2 * x1 + (1 - t2) * x2)
        set yt (t2 * y1 + (1 - t2) * y2)

      ]

      let xt4 (2 * xt + x2 - 2 * xt / radiii ^ 2 * (x2 * xt + y2 * yt ))
      let yt4 (2 * yt + y2 - 2 * yt / radiii ^ 2 * (x2 * xt + y2 * yt ))

      setxy xt4 yt4


  ]

  if sqrt (x2 * x2 + y2 * y2) + size / 2 > cellsize / 2
    [

    set radiii (cellsize / 2 - size / 2)
    let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
    let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
    let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 2 - size / 2) ^ 2)

    let temp_discriminant (bv ^ 2 - 4 * av * cv)

      let t1 (- bv + sqrt (temp_discriminant)) / (2 * av)
    let t2 (- bv - sqrt (temp_discriminant)) / (2 * av)

    let xt1 (t1 * x1 + (1 - t1) * x2)
    let yt1 (t1 * y1 + (1 - t1) * y2)

    let xt2 (t2 * x1 + (1 - t2) * x2)
    let yt2 (t2 * y1 + (1 - t2) * y2)

    let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
    let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))
      ; 거리로 구분하면 안되고 x1이랑 x2 사이에 있는 점을 골라야 함
      ; 그래서 t2를 쓰는게 항상 맞음;
      ; t1을 고르는 순간 애초에 틀리게 되는거다.

;    [
      set xt (t2 * x1 + (1 - t2) * x2)
      set yt (t2 * y1 + (1 - t2) * y2)
;
;
;    ]

    let xt4 (2 * xt + x2 - 2 * xt / radiii ^ 2 * (x2 * xt + y2 * yt ))
    let yt4 (2 * yt + y2 - 2 * yt / radiii ^ 2 * (x2 * xt + y2 * yt ))

      setxy xt4 yt4

  ]



end

to obstacle_movement

  let x1 xcor
  let y1 ycor

  rt random-float 360
  fd speedo

  let x2 xcor
  let y2 ycor

  if sqrt (x2 * x2 + y2 * y2) - size / 2 < cellsize / 6
    [

      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 6 + size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
          let xt (t1 * x1 + (1 - t1) * x2)
          let yt (t1 * y1 + (1 - t1) * y2)
          setxy xt yt


      ]
      [
        let xt (t2 * x1 + (1 - t2) * x2)
        let yt (t2 * y1 + (1 - t2) * y2)
        setxy xt yt
      ]

  ]

  if sqrt (x2 * x2 + y2 * y2) + size / 2 > cellsize / 2
    [

      let av ((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
      let bv (-2 * (x2 ^ 2 + y2 ^ 2 - x1 * x2 - y1 * y2 ))
      let cv (x2 ^ 2 + y2 ^ 2 - (cellsize / 2 - size / 2) ^ 2)
      let t1 (- bv + sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)
      let t2 (- bv - sqrt (bv ^ 2 - 4 * av * cv)) / (2 * av)

      let xt1 (t1 * x1 + (1 - t1) * x2)
      let yt1 (t1 * y1 + (1 - t1) * y2)

      let xt2 (t2 * x1 + (1 - t2) * x2)
      let yt2 (t2 * y1 + (1 - t2) * y2)

      let dis1 sqrt ((xt1 - x1) * (xt1 - x1) + (yt1 - y1) * (yt1 - y1))
      let dis2 sqrt ((xt2 - x1) * (xt2 - x1) + (yt2 - y1) * (yt2 - y1))

      ifelse dis1 < dis2
        [
          let xt (t1 * x1 + (1 - t1) * x2)
          let yt (t1 * y1 + (1 - t1) * y2)
          setxy xt yt


      ]
      [
        let xt (t2 * x1 + (1 - t2) * x2)
        let yt (t2 * y1 + (1 - t2) * y2)
        setxy xt yt

      ]

  ]

end
to get-protein-xcor

  set protein-xcor []
  set protein-ycor []
  ask phosrepressors [
    set protein-xcor lput (list [xcor] of self) protein-xcor
    set protein-ycor lput (list [ycor] of self) protein-ycor
  ]

end

to-report overlap [ agent ]
  report (size + ( [size] of agent / 20)) / 2 - distance agent
end

to-report fragment_num_report
  report fragment_num
end

to draw-box

  ask patches [set pcolor white]

end

to add1 [kind amount]
  create-turtles amount
  [ set breed kind
    setshape
    setxy 0 0]
end

to add2 [kind amount]
;  set protein_xcor csv:from-file "xcor_in_ER_v8_b_50_u_2_obs_100_rad_40_10x_init.csv"
;  set protein_ycor csv:from-file "ycor_in_ER_v8_b_50_u_2_obs_100_rad_40_10x_init.csv"
;  (foreach range(200) [
;    a -> set ind a
;    create-turtles 1
;    [ set breed kind
;      setshape
;      let thetha random-float 360
;      set is_stuck 0
;      setxy (item 0 (item ind protein_xcor)) (item 0 (item ind protein_ycor))]
;
;  ])


  create-turtles amount
    [ set breed kind
      setshape
;      let radius cellsize - size / 2
      let radius 33 + size / 2 + random-float (cellsize / 2  - 33 - size )
;      let radius cellsize / 6 + size / 2 + random-float (cellsize / 2 - cellsize / 6 - size )
;      let radius 27
      let thetha random-float 360
      set is_stuck 0
      setxy (radius * (cos thetha)) (radius * (sin thetha))]
;  set heading 0]
;  setxy 0 8 ]

end

to add22 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
;      let radius cellsize - size / 2
      let radius 6 + size / 2 + random-float (cellsize / 2  - 6 - size )
;      let radius cellsize / 6 + size / 2 + random-float (cellsize / 2 - cellsize / 6 - size )
;      let radius 27
      let thetha random-float 360
      setxy (radius * (cos thetha)) (radius * (sin thetha))]
;      set heading 0
;setxy 0 8]
end

to add3 [kind amount]
  create-turtles amount
    [ set breed kind
      setshape
      let radius random-float (cellsize / 2) - size / 2
      let thetha random-float 360

      let xvalue (radius * (cos thetha))
      let yvalue (radius * (sin thetha))
      let dist sqrt (xvalue * xvalue + yvalue * yvalue)

      while [dist + size / 2 > cellsize / 6 and dist - size / 2 < cellsize / 6]
      [
        set radius random-float (cellsize / 2) - size / 2
        set thetha random-float 360

        set xvalue (radius * (cos thetha))
        set yvalue (radius * (sin thetha))
        set dist sqrt (xvalue * xvalue + yvalue * yvalue)

      ]
      setxy xvalue yvalue]
end

to add_obs [kind amount]
  if (obstacle_num > 1)[
   set theta 0
  (foreach [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80] [
;      (foreach [ 1 2 3 4 5 6 7 8 9 10 ] [
    a -> set theta (a * 180 / obstacle_num)
;      print(a)
;      print(theta)
    let radius (nuclradi + a * radius_gap)
    create-turtles obstacle_num
    [ set breed kind
      setshape
      let xvalue (radius * (cos theta))
       let yvalue (radius * (sin theta))
      let dist sqrt (xvalue * xvalue + yvalue * yvalue)


      setxy xvalue yvalue



      set theta (theta + (360 / obstacle_num))]
  ])]

end

to add_obs_denser [kind amount]

   set theta 0
;  (foreach [ 1 2 3  6 7 8  11 12 13  16 17 18  21 22 23  26 27 28   31 32 33   36 37 38   41 42 43   46 47 48   51 52 53  56 57 58 ] [
  (foreach [ 1 2 3 4 5 6 7 8 9 10 11 12] [
    a -> set theta 30 * a
    let radius radius_1st
    create-turtles 3
    [ set breed kind
      setshape
      let xvalue (radius * (cos theta))
      let yvalue (radius * (sin theta))
      let dist sqrt (xvalue * xvalue + yvalue * yvalue)


      setxy xvalue yvalue

      set fragment_num 1


      set theta (theta + 10)]
  ])

  set theta 0
;  (foreach [ 1 2 3  6 7 8  11 12 13  16 17 18  21 22 23  26 27 28   31 32 33   36 37 38   41 42 43   46 47 48   51 52 53  56 57 58 ] [
  (foreach [1 2 3 4 5 6 7 8 9 10 ] [
    a -> set theta (36 * a )
    let radius radius_2nd
    create-turtles 3
    [ set breed kind
      setshape
      let xvalue (radius * (cos theta))
      let yvalue (radius * (sin theta))
      let dist sqrt (xvalue * xvalue + yvalue * yvalue)


      setxy xvalue yvalue

      set fragment_num 2


      set theta (theta + 12)]
  ])

   set theta 0
;  (foreach [ 1 2 3  6 7 8  11 12 13  16 17 18  21 22 23  26 27 28   31 32 33   36 37 38   41 42 43   46 47 48   51 52 53  56 57 58 ] [
  (foreach [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15] [
    a -> set theta (24 * a )
    let radius radius_3rd
    create-turtles 1
    [ set breed kind
      setshape
      let xvalue (radius * (cos theta))
      let yvalue (radius * (sin theta))
      let dist sqrt (xvalue * xvalue + yvalue * yvalue)


      setxy xvalue yvalue

      set fragment_num 3


      set theta (theta + 10)]
  ])

  ;  (foreach [ 1 2 3  6 7 8  11 12 13  16 17 18  21 22 23  26 27 28   31 32 33   36 37 38   41 42 43   46 47 48   51 52 53  56 57 58 ] [
  (foreach [ 1 2 3 4 5 6 7 8 9 ] [
    a -> set theta (20 * a )
    let radius radius_4th
    create-turtles 2
    [ set breed kind
      setshape
      let xvalue (radius * (cos theta))
      let yvalue (radius * (sin theta))
      let dist sqrt (xvalue * xvalue + yvalue * yvalue)


      setxy xvalue yvalue

      set fragment_num 4


      set theta (theta + 20)]
  ])
end


to add_obs_maze [kind amount]

   set theta 0
;  (foreach [ 1 2 3  6 7 8  11 12 13  16 17 18  21 22 23  26 27 28   31 32 33   36 37 38   41 42 43   46 47 48   51 52 53  56 57 58 ] [
  (foreach [ 1 3 5 7 9 11  ] [
    a -> set theta 30 * a
    let radius radius_1st
    create-turtles 24
    [ set breed kind
      setshape
      let xvalue (radius * (cos theta))
      let yvalue (radius * (sin theta))
      let dist sqrt (xvalue * xvalue + yvalue * yvalue)


      setxy xvalue yvalue

      set fragment_num floor(a / 5)


      set theta (theta + 1.25)]
  ])

  set theta 0
;  (foreach [ 1 2 3  6 7 8  11 12 13  16 17 18  21 22 23  26 27 28   31 32 33   36 37 38   41 42 43   46 47 48   51 52 53  56 57 58 ] [
  (foreach [ 1  3  5 7  9 11 ] [
    a -> set theta (30 * a + 30)
    let radius radius_2nd
    create-turtles 24
    [ set breed kind
      setshape
      let xvalue (radius * (cos theta))
      let yvalue (radius * (sin theta))
      let dist sqrt (xvalue * xvalue + yvalue * yvalue)


      setxy xvalue yvalue

      set fragment_num floor(a / 5)


      set theta (theta + 1.25)]
  ])

   set theta 0
;  (foreach [ 1 2 3  6 7 8  11 12 13  16 17 18  21 22 23  26 27 28   31 32 33   36 37 38   41 42 43   46 47 48   51 52 53  56 57 58 ] [
  (foreach [ 1   3  5  7  9  11 ] [
    a -> set theta (30 * a )
    let radius radius_3rd
    create-turtles 24
    [ set breed kind
      setshape
      let xvalue (radius * (cos theta))
      let yvalue (radius * (sin theta))
      let dist sqrt (xvalue * xvalue + yvalue * yvalue)


      setxy xvalue yvalue

      set fragment_num floor(a / 5)


      set theta (theta + 1.25)]
  ])
;;;;;;;;;;;;;;;;;;;;

end



;; procedure that assigns a specific shape to a turtle

to setshape

  if breed = cells
  [ set color Grey
    set color lput 20 extract-rgb color
    set shape "circle"
    set size cellsize  ]

  if breed = nuclei
  [set color White
    set shape "circle"

    set size (cellsize / 3)  ]

  if breed = obstacles
  [set color black
;    set color lput 100 extract-rgb color
    set shape "circle"
    set size 0.05]

  if breed = mRNArepressors
  [set color cyan
    set shape "circle"
    set size 0.5]

  if breed = repressors
  [set color [102 145 255] ; blue
    set shape "circle"
    set size 20 * repressor_radius]

  if breed = phosrepressors
  [set color Violet

    set color lput 70 extract-rgb color
;
;    set color [255 54 66] ; red
    set color [102 145 255] ; blue
    set color lput 255 color
    set shape "circle"
    set size 20 * phosrepressor_radius]

  if breed = frepressors
  [set color blue
    set color lput 0 extract-rgb color
    set shape "circle"
    set size 0.5]

  if breed = fphosrepressors
  [set color green
    set color lput 0 extract-rgb color
    set shape "circle"
    set size 0.5]

end
@#$#@#$#@
GRAPHICS-WINDOW
730
94
1392
757
-1
-1
3.224
1
10
1
1
1
0
1
1
1
-101
101
-101
101
1
1
1
ticks
30.0

BUTTON
153
11
259
53
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
270
10
376
52
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
15
176
711
394
Concentration
time
Conc.
0.0
50000.0
0.0
100.0
true
true
"" ""
PENS
"Large" 1.0 0 -8732573 true "" "plot count phosrepressors with [(distancexy 0 0) < (nuclradi + 20 + phosrepressor_radius / 2)]"
"Small" 1.0 0 -817084 true "" "plot count repressors with [(distancexy 0 0) < (5 + repressor_radius / 2)]"

SLIDER
28
101
200
134
P
P
0
500
200.0
5
1
NIL
HORIZONTAL

TEXTBOX
27
79
102
97
Parameter
14
0.0
1

BUTTON
25
11
128
54
load
load
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
227
101
399
134
f_binding
f_binding
0
10
5.0
0.5
1
NIL
HORIZONTAL

SLIDER
411
101
583
134
f_unbinding
f_unbinding
0
0.2
0.02
0.01
1
NIL
HORIZONTAL

SLIDER
227
138
399
171
obstacle_num
obstacle_num
0
1000
1000.0
10
1
NIL
HORIZONTAL

SLIDER
412
139
583
172
phosrepressor_radius
phosrepressor_radius
0
0.5
0.2
0.005
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

This model was extended from the previous mathematical model of the mammalian circadian clock (Kim et al., 2012, MSB). To investigate the effect of the crowdedness of cytoplasmic contents on the circadina clock, this model was developed.

The model simulation suggests that the spatial regulation of PER protein, the repressor of the circadian clock, leads to its sharp switch-like phosphorylation and thus sharp nuclear translocation. This nonlinear nuclear entry plays the critial role in generating the robust circadian rhythms. Please see Beesley et al., for details description of the model.


## HOW TO USE IT

Step 1. Please set the initial number of agents.

Atot: Total number of activators.

M: The initial number of Per mRNA.

HYPOPER: The initial number of hypophosphorylated PER.

HYPERPER: The initial number of hyperphosphorylated PER.

O: Number of cytopalsmic obstacles. Please choose the value of 0 to regulate the crowdedness of the cytoplasmic contents.
Examples
- The model with 150 cytoplasmic obstacles is the normal cell.
- The model with 275 cytoplasmic obstacles is the overcrowded cell.
- The model with more than 300 cytoplasmoc obstacles is the extremely overcrowded cell (i.e. adipocyte). Thus, the model do not simulate the rhythmic PER expression.



Step 2. Please set the parameters, which are described below. The current parameter setting makes the model to simulate the rhythmic PER expression of the normal cell.

pa1: Reaction probability for Per mRNA production for each time step (i.e. tick).
pa2: Reaction probability for PER protein translation for each time step.
pd1: Reaction probability for Per mRNA degradation for each time step.
pd2: Reaction probability for hypophos. PER degradtion for each time step. 
pd3: Reaction probability for hyperphos. PER degradtion for each time step.

Kd: Dissociation constant between hyperphos. PER and activator.
D: Movement step size of PermRNA and PER protein for each time step.
Dobs: Movement step size of cytoplasmic obstacles for each time step.
padvec: Probability that the PER protein is advected to the peri-nucleus by the cytoplasmic flux.
pim: Probability that hyperphos.PER in the cytoplasm is imported to the nucleus for each time step.



Step 3. Load the reation probabilities for hyperphosphorylation and dephosphrylation by pushing "load" buttom.



Step 4. Set the initial condition of the model by pushing "setup" buttom.



Step 5. Run the simulation by pushing "go" buttom. 
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

complex
true
0
Polygon -2674135 true false 76 47 197 150 76 254 257 255 257 47
Polygon -10899396 true false 79 46 198 148 78 254

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

enzyme
true
0
Polygon -2674135 true false 76 47 197 150 76 254 257 255 257 47

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

inhib complex
true
0
Polygon -2674135 true false 76 47 197 150 76 254 257 255 257 47
Polygon -1184463 true false 77 48 198 151 78 253 0 253 0 46

inhibitor
true
0
Polygon -1184463 true false 197 151 60 45 1 45 1 255 60 255

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

substrate
true
5
Polygon -10899396 true true 76 47 197 151 75 256

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.3.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
