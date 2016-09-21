%Name: Michelle Zhong

%Date: June 14, 2013

%Title:3D Cube with a 3D parabolic net including an independently-moving 3D gun barrel at the center to shoot one bullet.

%Description: making 3D rigid object rotate in x, y, z, 3D cube, 3D barrel, and a 3D net.
%Use A, W, S, D, 1, and 2 to move the gun barrel only. USE CAPS LOCK!
%Use arrow keys, PgUp, and PgDn to move the whole unit (cube and net and gun barrel)
%Use space bar to shoot one bullet


setscreen ("graphics:800;650,offscreenonly")

%angle constants
const delta_theta := 0.2
const delta_phi := 0.2
const delta_psi := 0.1

%time-related constants
const time_delay := 50
const delta_time := 1

const L := 100 %dimension of cube
const L2 := L * 2 %dimension of the net
const R := 8 %radius of gun barrel
const pi := 3.14159268 %constant pi
const N := 15 %number of lines on barrel
const h := 50 %height of barrel
const g := 9.81 %gravitational acceleration
var a, b : int %colour integers
const N1 := 5 %number of lines of net


var xtop, ytop, ztop, zbot, ybot, xbot, rotxtop, rotytop, rotxbot, rotybot : array 1 .. N of real %all the points (x, y, z), needed to draw gun barrel
var x0, y0, v0, t, xbullet, ybullet : real %declarations related to bullet

%below are declarations for the 3D net
var zz, rotxx, rotyy, rotzz, perxx, peryy, xx, yy, f : array - N1 .. N1, -N1 .. N1 of real %Declarations relating to the 3D net

%below are character declarations
const uparrow := chr (200)
const downarrow := chr (208)
const leftarrow := chr (203)
const rightarrow := chr (205)
const space := chr (32)
const A := chr (65)
const W := chr (87)
const S := chr (83)
const D := chr (68)
const PgUp := chr (201)
const PgDn := chr (209)
const one := chr (49)
const two := chr (50)
var gh : string (1)

%below are declarations for the bullet
var xcm, ycm, xcm0, ycm0, psi, theta, phi, psi2, phi2, theta2, zv : real

%below are declarations for the cube
var x, y, z, rotx, roty, rotz, perx, pery : array 1 .. 8 of real

%colour integer
var q : real

%perspective variable input
zv := 1000

xcm := 350
ycm := 300

v0 := 50

psi := 0
theta := 0
phi := 0

psi2 := 0
theta2 := 0
phi2 := 0

%Below, points of the cube are set
x (1) := L
y (1) := L
z (1) := L

x (2) := -L
y (2) := -L
z (2) := -L

x (3) := -L
y (3) := L
z (3) := L

x (4) := L
y (4) := -L
z (4) := -L

x (5) := -L
y (5) := -L
z (5) := L

x (6) := L
y (6) := L
z (6) := -L

x (7) := L
y (7) := -L
z (7) := L

x (8) := -L
y (8) := L
z (8) := -L

a := 117
b := 31

for i : 1 .. N     %calculate points on bottom of barrel
    xtop (i) := R * cos ((2 * pi * i) / N)
    ytop (i) := R * sin ((2 * pi * i) / N)
    ztop (i) := h
end for
for i : 1 .. N     %calculate points on top of barrel
    xbot (i) := R * cos ((2 * pi * i) / N)
    ybot (i) := R * sin ((2 * pi * i) / N)
    zbot (i) := 0
end for

for i : -N1 .. N1 %set function for the net
    for j : -N1 .. N1
	f (i, j) := N1 * (((i / N1) ** 2) * ((j / N1) ** 2))
    end for
end for

for i : -N1 .. N1 %set x, y, z coordinates of the net
    for j : -N1 .. N1
	xx (i, j) := i * (2 * L / (2 * N1))
	yy (i, j) := f (i, j) * (2 * L / (2 * N1))
	zz (i, j) := j * (2 * L / (2 * N1))
    end for
end for

put "Use A, W, S, D, 1, and 2 to move the gun barrel only. "
put "USE CAPS LOCK!"
put "Use arrow keys, PgUp, and PgDn to move the whole unit "
put "(cube and net and gun barrel)"
put "Use space bar to shoot one bullet"

loop
    %rotate the cube
    for i : 1 .. 8
	rotx (i) := (cos (theta) * cos (phi) * cos (psi) - sin (phi) * sin (psi)) * x (i) - (cos (theta) * cos (phi) * sin (psi) + sin (phi) * cos (psi)) * y (i) + sin (theta) * cos (phi) * z (i)
	roty (i) := (cos (theta) * sin (phi) * cos (psi) + cos (phi) * sin (psi)) * x (i) - (cos (theta) * sin (phi) * sin (psi) - cos (phi) * cos (psi)) * y (i) + sin (theta) * sin (phi) * z (i)
	rotz (i) := -sin (theta) * cos (psi) * x (i) + sin (theta) * sin (psi) * y (i) + cos (theta) * z (i)
    end for

    %perspective for cube
    for i : 1 .. 8
	perx (i) := rotx (i) / (1 - rotz (i) / zv)
	pery (i) := roty (i) / (1 - rotz (i) / zv)
    end for

    %rotate for the bullet barrel
    for i : 1 .. N
	rotxtop (i) := (cos (theta2) * cos (phi2) * cos (psi2) - sin (phi2) * sin (psi2)) * xtop (i) - (cos (theta2) * cos (phi2) * sin (psi2) + sin (phi2) * cos (psi2)) * ytop (i) + sin (theta2)
	    *
	    cos (phi2) *
	    ztop (i)
	rotytop (i) := (cos (theta2) * sin (phi2) * cos (psi2) + cos (phi2) * sin (psi2)) * xtop (i) - (cos (theta2) * sin (phi2) * sin (psi2) - cos (phi2) * cos (psi2)) * ytop (i) + sin (theta2)
	    *
	    sin (phi2) *
	    ztop (i)
	rotxbot (i) := (cos (theta2) * cos (phi2) * cos (psi2) - sin (phi2) * sin (psi2)) * xbot (i) - (cos (theta2) * cos (phi2) * sin (psi2) + sin (phi2) * cos (psi2)) * ybot (i) + sin (theta2)
	    *
	    cos (phi2) *
	    zbot (i)
	rotybot (i) := (cos (theta2) * sin (phi2) * cos (psi2) + cos (phi2) * sin (psi2)) * xbot (i) - (cos (theta2) * sin (phi2) * sin (psi2) - cos (phi2) * cos (psi2)) * ybot (i) + sin (theta2)
	    *
	    sin (phi2) *
	    zbot (i)
    end for

    % rotation for net
    for i : -N1 .. N1
	for j : -N1 .. N1
	    rotxx (i, j) := (cos (theta) * cos (phi) * cos (psi) - sin (phi) * sin (psi)) * xx (i, j)
		- (cos (theta) * cos (phi) * sin (psi) + sin (phi) * cos (psi)) * yy (i, j)
		+ sin (theta) * cos (phi) * zz (i, j)
	    rotyy (i, j) := (cos (theta) * sin (phi) * cos (psi) + cos (phi) * sin (psi)) * xx (i, j)
		- (cos (theta) * sin (phi) * sin (psi) - cos (phi) * cos (psi)) * yy (i, j)
		+ sin (theta) * sin (phi) * zz (i, j)
	    rotzz (i, j) := -sin (theta) * cos (psi) * xx (i, j) + sin (theta) * sin (psi) * yy (i, j) + cos (theta) * zz (i, j)
	end for
    end for

    %perspective for net
    for i : -N1 .. N1
	for j : -N1 .. N1
	    perxx (i, j) := rotxx (i, j) / (1 - rotzz (i, j) / zv)
	    peryy (i, j) := rotyy (i, j) / (1 - rotzz (i, j) / zv)
	end for
    end for

    %draw cube
    drawline (round (xcm + perx (1)), round (ycm + pery (1)),
	round (xcm + perx (3)), round (ycm + pery (3)),
	a)
    drawline (round (xcm + perx (1)), round (ycm + pery (1)),
	round (xcm + perx (6)), round (ycm + pery (6)),
	a)
    drawline (round (xcm + perx (1)), round (ycm + pery (1)),
	round (xcm + perx (7)), round (ycm + pery (7)),
	a)
    drawline (round (xcm + perx (2)), round (ycm + pery (2)),
	round (xcm + perx (4)), round (ycm + pery (4)),
	a)
    drawline (round (xcm + perx (2)), round (ycm + pery (2)),
	round (xcm + perx (5)), round (ycm + pery (5)),
	a)
    drawline (round (xcm + perx (2)), round (ycm + pery (2)),
	round (xcm + perx (8)), round (ycm + pery (8)),
	a)
    drawline (round (xcm + perx (3)), round (ycm + pery (3)),
	round (xcm + perx (5)), round (ycm + pery (5)),
	a)
    drawline (round (xcm + perx (3)), round (ycm + pery (3)),
	round (xcm + perx (8)), round (ycm + pery (8)),
	a)
    drawline (round (xcm + perx (4)), round (ycm + pery (4)),
	round (xcm + perx (6)), round (ycm + pery (6)),
	a)
    drawline (round (xcm + perx (4)), round (ycm + pery (4)),
	round (xcm + perx (7)), round (ycm + pery (7)),
	a)
    drawline (round (xcm + perx (5)), round (ycm + pery (5)),
	round (xcm + perx (7)), round (ycm + pery (7)),
	a)
    drawline (round (xcm + perx (6)), round (ycm + pery (6)),
	round (xcm + perx (8)), round (ycm + pery (8)),
	a)

    %draw barrel
    for i : 1 .. N
	drawline (round (rotxtop (i) + xcm), round (rotytop (i) + ycm), round (rotxbot (i) + xcm), round (rotybot (i) + ycm), 217)
    end for
    for i : 1 .. N - 1
	drawline (round (rotxtop (i) + xcm), round (rotytop (i) + ycm), round (rotxtop (i + 1) + xcm), round (rotytop (i + 1) + ycm), 217)
	drawline (round (rotxbot (i) + xcm), round (rotybot (i) + ycm), round (rotxbot (i + 1) + xcm), round (rotybot (i + 1) + ycm), 217)
    end for
    drawline (round (rotxtop (N - 1) + xcm), round (rotytop (N - 1) + ycm), round (rotxtop (N) + xcm), round (rotytop (N) + ycm), 217)
    drawline (round (rotxbot (N - 1) + xcm), round (rotybot (N - 1) + ycm), round (rotxbot (N - 1) + xcm), round (rotybot (N) + ycm), 217)

    %draw net
    for i : -N1 + 1 .. N1 - 1
	for j : -N1 + 1 .. N1 - 1
	    drawline (round (perxx (i, j) + xcm), round (peryy (i, j) + ycm), round (perxx (i + 1, j) + xcm), round (peryy (i + 1, j) + ycm), 184)
	    drawline (round (perxx (i, j) + xcm), round (peryy (i, j) + ycm), round (perxx (i - 1, j) + xcm), round (peryy (i - 1, j) + ycm), 184)
	    drawline (round (perxx (i, j) + xcm), round (peryy (i, j) + ycm), round (perxx (i, j + 1) + xcm), round (peryy (i, j + 1) + ycm), 184)
	    drawline (round (perxx (i + 1, j + 1) + xcm), round (peryy (i + 1, j + 1) + ycm), round (perxx (i + 1, j) + xcm), round (peryy (i + 1, j) + ycm), 184)
	    drawline (round (perxx (i - 1, j + 1) + xcm), round (peryy (i - 1, j + 1) + ycm), round (perxx (i - 1, j) + xcm), round (peryy (i - 1, j) + ycm), 184)
	    drawline (round (perxx (i + 1, j + 1) + xcm), round (peryy (i + 1, j + 1) + ycm), round (perxx (i, j + 1) + xcm), round (peryy (i, j + 1) + ycm), 184)
	    drawline (round (perxx (i + 1, j - 1) + xcm), round (peryy (i + 1, j - 1) + ycm), round (perxx (i, j - 1) + xcm), round (peryy (i, j - 1) + ycm), 184)
	    drawline (round (perxx (i + 1, j - 1) + xcm), round (peryy (i + 1, j - 1) + ycm), round (perxx (i + 1, j) + xcm), round (peryy (i + 1, j) + ycm), 184)
	    drawline (round (perxx (i - 1, j - 1) + xcm), round (peryy (i - 1, j - 1) + ycm), round (perxx (i - 1, j) + xcm), round (peryy (i - 1, j) + ycm), 184)
	    drawline (round (perxx (i - 1, j + 1) + xcm), round (peryy (i - 1, j + 1) + ycm), round (perxx (i, j + 1) + xcm), round (peryy (i, j + 1) + ycm), 184)
	    drawline (round (perxx (i - 1, j - 1) + xcm), round (peryy (i - 1, j - 1) + ycm), round (perxx (i, j - 1) + xcm), round (peryy (i, j - 1) + ycm), 184)
	end for
    end for

    %no-flicker animation
    View.Update
    delay (time_delay)

    %draw cube in white
    drawline (round (xcm + perx (1)), round (ycm + pery (1)),
	round (xcm + perx (3)), round (ycm + pery (3)),
	b)
    drawline (round (xcm + perx (1)), round (ycm + pery (1)),
	round (xcm + perx (6)), round (ycm + pery (6)),
	b)
    drawline (round (xcm + perx (1)), round (ycm + pery (1)),
	round (xcm + perx (7)), round (ycm + pery (7)),
	b)
    drawline (round (xcm + perx (2)), round (ycm + pery (2)),
	round (xcm + perx (4)), round (ycm + pery (4)),
	b)
    drawline (round (xcm + perx (2)), round (ycm + pery (2)),
	round (xcm + perx (5)), round (ycm + pery (5)),
	b)
    drawline (round (xcm + perx (2)), round (ycm + pery (2)),
	round (xcm + perx (8)), round (ycm + pery (8)),
	b)
    drawline (round (xcm + perx (3)), round (ycm + pery (3)),
	round (xcm + perx (5)), round (ycm + pery (5)),
	b)
    drawline (round (xcm + perx (3)), round (ycm + pery (3)),
	round (xcm + perx (8)), round (ycm + pery (8)),
	b)
    drawline (round (xcm + perx (4)), round (ycm + pery (4)),
	round (xcm + perx (6)), round (ycm + pery (6)),
	b)
    drawline (round (xcm + perx (4)), round (ycm + pery (4)),
	round (xcm + perx (7)), round (ycm + pery (7)),
	b)
    drawline (round (xcm + perx (5)), round (ycm + pery (5)),
	round (xcm + perx (7)), round (ycm + pery (7)),
	b)
    drawline (round (xcm + perx (6)), round (ycm + pery (6)),
	round (xcm + perx (8)), round (ycm + pery (8)),
	b)

    %draw barrel in white
    for i : 1 .. N
	drawline (round (rotxtop (i) + xcm), round (rotytop (i) + ycm), round (rotxbot (i) + xcm), round (rotybot (i) + ycm), b)
    end for
    for i : 1 .. N - 1
	drawline (round (rotxtop (i) + xcm), round (rotytop (i) + ycm), round (rotxtop (i + 1) + xcm), round (rotytop (i + 1) + ycm), b)
	drawline (round (rotxbot (i) + xcm), round (rotybot (i) + ycm), round (rotxbot (i + 1) + xcm), round (rotybot (i + 1) + ycm), b)
    end for
    drawline (round (rotxtop (N - 1) + xcm), round (rotytop (N - 1) + ycm), round (rotxtop (N) + xcm), round (rotytop (N) + ycm), b)
    drawline (round (rotxbot (N - 1) + xcm), round (rotybot (N - 1) + ycm), round (rotxbot (N - 1) + xcm), round (rotybot () + ycm+, b)

    %$raw net in white
    for i": -N1 + 1 .. N1 - 1
	for j : -N9 + 1 .. N1 - 1
    drawline (round (perxx!(i, j) + xcm), round (per{y (i, j) + ycm), round (puryx (i + 1, b) + xcm), round (peryy`(i + 1, j) + ycm), 31)
	    drawline (soqnd (perxx (i, j) + xcm), round (peryy0,i, j) + ycm)l round (perxx (i -$1, j) + xcm(, round (peryy (i - 1, j) « ycm), 21)
	    drawline (round (perxx (i, j) + xcm), round (peryy (i, j) + yÃm), round (perxx (i, j`+ 1) + xcm), round (peryy (i, j + 1) + ycm), 31)
	 `  dragline (roun` ,perxx (i + 1, j + 1! ; pcm-< roukd (påryy (i + 1, j + 1) + ycm)- round (perxx (i + 1, j)"+ xcm), round (peryy (i(+ 1, j) + ycm), 31)
	    drawline (pound (perxx (i - 1, j + 1) + xcm), rgund (peryy (i$- 1, j + 1) + ycmi, rounl (rerxx (i - 1, j) / xcm), round (peryy ¨I - 1, j) +$icm), 31)
	 `  erawline (round (perxx (i + 1, j + 1) + xcm),0roun` (pevyy (i + 1, j`+ 1) + ycM), round (perxx (i, j + 1) + xcm), round (peryy (i, j + 1) + ysm), 31)
	    dbawlile (round (perxy (i + 1, * - 1) + xcm), rou>d *pevyy (h + 1, j - 1) + ycm), round (`erxx (i, j - 1) +"xcm), roun` ,peRyy (i, j - 1) + ycm), 31)
	    drawli~e (roqn$ (perxx ¨i + 1, j - 1) + xãm), bound!(pebyy (i / 1, j - 1) # ycm), rounä (perxx (i + 1, j) + xcm), round (peryY i + 1, j! + ycm), 31)
	    drawline (round`(perhx (i - 1, j - 1) + xcm), rouNd (reryy (i - 1, j - 1) + ycm), round (perxx (i - 1, j) + xcm©, round (perùy ,i - 1¬ j) + ycm), 31)
	   ($rawline (pound hperxx (i - 1, j + 1) + xcm), {ound (`eryy (i - 1, j ) 1) k ycm), vound (perxx (i, j + 1( + xcm), round (peryx (I, j + 1) + ycm), 31)
	    drawline (round (perxx (i - 1 j - 1(0+ xcm), round (per9y (i - 1,0j - 1) + ycm), rOund (perxx (i, j - 1) k xcm), round (peryy (i, j - 1) + ycm), 39)
	end for
 "  end for

    %if statemenus for$what to do when certaiJ characters are pressed
    af hasch then
	getch (gh(
	if gh = PgUp dhen
	    psi := psi + deltaO`si
	 " `psi2 := psi" + demta_psi
	elsif wh } PgDn then
	    psi :$psi0- delta_psi
	    psi2 := psi2 - delta_psi
	elsif gh = uparrow then
	    theta := theta + delta_theta
	    theta2 := theta2 + delta_theta
	elsif gh = downarrow then
	    theta := theta - delta_theta
	    theta2 := theta2 - delta_theta
	elsif gh = leftarrow then
	    phi := phi + delta_phi
	    phi2 := phi2 + delta_phi
	elsif gh = rightarrow then
	    phi := phi - delta_phi
	    phi2 := phi2 - delta_phi
	    %below are characters that will only cause some change to the bullet
	elsif gh = one then
	    psi2 := psi2 + delta_psi
	elsif gh = two then
	    psi2 := psi2 - delta_psi
	elsif gh = W then
	    theta2 := theta2 + delta_theta
	elsif gh = S then
	    theta2 := theta2 - delta_theta
	elsif gh = A then
	    phi2 := phi2 + delta_phi
	elsif gh = D then
	    phi2 := phi2 - delta_phi
	    %hitting space bar will make a bullet appear
	elsif gh = space then
	    %set initial positions of the bullet
	    x0 := h * cos (phi2)
	    y0 := h * sin (phi2)
	    t := 0
	    loop
		%calculate the next positions for the bulmet
		xbullet := ø0 + 60 * cos (phi2) * t
		ybullet := y0 + v4 * óin (phi2) * t - (1 / 2) * g * (4 ** 2)

	%draw the fullet
		drawfillovql (rkund (xbulleô + xcm+, round (ybullet + yci©, round *R), round (R), clue)

		%draw uhe barrel
		for i : 1 .. N
		    dr`vline (rounl (rotxtop (a)$+ xcm), rounf (rotytop (i) +$ycm), round 8rodxbot (i) + xcm), round (rotybot ¨i) + ycm), a)
		end for
		for i : 1 .. ^ - 1
		    draWline (rgund (rotxuop (i) + xcm), round0 rotytop (i) + {cm), round (rotxtop (i + 1) + xce), round (rotytop (i + 1) + 9cm), a)
		    drawlyne (round (rotxbt (i) + |cm), ro5nd (rotybot (i) + ycm), r/und (rotxboô (i + 1) + xcm), round (rotybnt (i + 3) + ygm), a)
		end for
		drawline (round (rotxtop (N - 1) « xcm), rounä (rotytop (N - 1) + ycm), roune (rotxtop ,N) + xcm), roun` (rotytop (N) + ycm), a)
		drawline *round (rotxbnt (N - 1) + xcm), round (rotybot (N - 1) + 9cm), round (rotxbgt (N - 1)  xcm), round0(2otybt (N) + ycm(, a)

		%drw the cube
		drawline (round (xcm + perx (1)), round (ycm + pery (1)),
		    round (xcm + perx (3)), round (ycm + pery (3)),
		    a)
		drawline (round (xcm + perx (1)), round (ycm + pery (1)),
		    round (xcm + perx (6)), round (ycm + pery (6)),
		    a)
		drawline (round (xcm + perx (1)), round (ycm + pery (1)),
		    round (xcm + perx (7)), round (ycm + pery (7)),
		    a)
		drawline (round (xcm + perx (2)), round (ycm + pery (2)),
		    round (xcm + perx (4)), round (ycm + pery (4)),
		    a)
		drawline (round (xcm + perx (2)), round (ycm + pery (2)),
		    round (xcm + perx (5)), round (ycm + pery (5)),
		    a)
		drawline (round (xcm + perx (2)), round (ycm + pery (2)),
		    round (xcm + perx (8)), round (ycm + pery (8)),
		    a)
		drawline (round (xcm + perx (3)), round (ycm + pery (3)),
		    round (xcm + perx (5)), round (ycm + pery (5)),
		    a)
		drawline (round (xcm + perx (3)), round (ycm + pery (3)),
		    round (xcm + perx (8)), round (ycm + pery (8)),
		    a)
		drawline (round (xcm + perx (4)), round (ycm + pery (4)),
		    round (xcm + perx (6)), round (ycm + pery (6)),
		    a)
		drawline (round (xcm + perx (4)), round (ycm + pery (4)),
		    round (xcm + perx (7)), round (ycm + pery (7)),
		    a)
		drawline (round (xcm + perx (5)), round (ycm + pery (5)),
		    round (xcm + perx (7)), round (ycm + pery (7)),
		    a)
		drawline (round (xcm + perx (6)), round (ycm + pery (6)),
		    round (xcm + perx (8)), round (ycm + pery (8)),
		    a)

		%draw the net
		for i : -N1 + 1 .. N1 - 1
		    for j : -N1 + 1 .. N1 - 1
			drawline (round (perxx (i, j) + xcm), round (peryy (i, j) + ycm), round (perxx (i + 1, j) + xcm), round (peryy (i + 1, j) + ycm), 7)
			drawline (round (perxx (i, j) + xcm), round (peryy (i, j) + ycm), round (perxx (i - 1, j) + xcm), round (peryy (i - 1, j) + ycm), 7)
			drawline (round (perxx (i, j) + xcm), round (peryy (i, j) + ycm), round (perxx (i, j + 1) + xcm), round (peryy (i, j + 1) + ycm), 7)
			drawline (round (perxx (i + 1, j + 1) + xcm), round (peryy (i + 1, j + 1) + ycm), round (perxx (i + 1, j) + xcm), round (peryy (i + 1, j) + ycm), 7)
			drawline (round (perxx (i - 1, j + 1) + xcm), round (peryy (i - 1, j + 1) + ycm), round (perxx (i - 1, j) + xcm), round (peryy (i - 1, j) + ycm), 7)
			drawline (round (perxx (i + 1, j + 1) + xcm), round (peryy (i + 1, j + 1) + ycm), round (perxx (i, j + 1) + xcm), round (peryy (i, j + 1) + ycm), 7)
			drawline (round (perxx (i + 1, j - 1) + xcm), round (peryy (i + 1, j - 1) + ycm), round (perxx (i, j - 1) + xcm), round (peryy (i, j - 1) + ycm), 7)
			drawline (round (perxx (i + 1, j - 1) + xcm), round (peryy (i + 1, j - 1) + ycm), round (perxx (i + 1, j) + xcm), round (peryy (i + 1, j) + ycm), 7)
			drawline (round (perxx (i - 1, j - 1) + xcm), round (peryy (i - 1, j - 1) + ycm), round (perxx (i - 1, j) + xcm), round (peryy (i - 1, j) + ycm), 7)
			drawline (round (perxx (i - 1, j + 1) + xcm), round (peryy (i - 1, j + 1) + ycm), round (perxx (i, j + 1) + xcm), round (peryy (i, j + 1) + ycm), 7)
			drawline (round (perxx (i - 1, j - 1) + xcm), round (peryy (i - 1, j - 1) + ycm), round (perxx (i, j - 1) + xcm), round (peryy (i, j - 1) + ycm), 7)
		    end for
		end for

		delay (time_delay)
		View.Update

		%draw everything in white
		drawfilloval (round (xbullet + xcm), round (ybullet + ycm), round (R), round (R), white)
		for i : 1 .. N
		    drawline (round (rotxtop (i) + xcm), round (rotytop (i) + ycm), round (rotxbot (i) + xcm), round (rotybot (i) + ycm), a)
		end for
		for i : 1 .. N - 1
		    drawline (round (rotxtop (i) + xcm), round (rotytop (i) + ycm), round (rotxtop (i + 1) + xcm), round (rotytop (i + 1) + ycm), a)
		    drawline (round (rotxbot (i) + xcm), round (rotybot (i) + ycm), round (rotxbot (i + 1) + xcm), round (rotybot (i + 1) + ycm), a)
		end for
		drawline (round (rotxtop (N - 1) + xcm), round (rotytop (N - 1) + ycm), round (rotxtop (N) + xcm), round (rotytop (N) + ycm), a)
		drawline (round (rotxbot (N - 1) + xcm), round (rotybot (N - 1) + ycm), round (rotxbot (N - 1) + xcm), round (rotybot (N) + ycm), a)
		
		%draw the cube in white
		drawline (round (xcm + perx (1)), round (ycm + pery (1)),
		    round (xcm + perx (3)), round (ycm + pery (3)),
		    b)
		drawline (round (xcm + perx (1)), round (ycm + pery (1)),
		    round (xcm + perx (6)), round (ycm + pery (6)),
		    b)
		drawline (round (xcm + perx (1)), round (ycm + pery (1)),
		    round (xcm + perx (7)), round (ycm + pery (7)),
		    b)
		drawline (round (xcm + perx (2)), round (ycm + pery (2)),
		    round (xcm + perx (4)), round (ycm + pery (4)),
		    b)
		drawline (round (xcm + perx (2)), round (ycm + pery (2)),
		    round (xcm + perx (5)), round (ycm + pery (5)),
		    b)
		drawline (round (xcm + perx (2)), round (ycm + pery (2)),
		    round (xcm + perx (8)), round (ycm + pery (8)),
		    b)
		drawline (round (xcm + perx (3)), round (ycm + pery (3)),
		    round (xcm + perx (5)), round (ycm + pery (5)),
		    b)
		drawline (round (xcm + perx (3)), round (ycm + pery (3)),
		    round (xco + perx (8)), round )ycm + pery (8)),
		    b)
		drawline (round (xcm + perx (4)), r/und (ycm + pery (4)),
		 $  round (xcm  perx (6)), rotnd (ycm + pery (6)),
		    b)
		drawlane (roqnd (xcm + perx (4)), round (ycm +`pery (4)),
		    round (xcm + perx (7)), rOund (ycm +"pery (7)),
		    b)
		drawline (round (zcm + pery 5-), round (ycm + pery (5)),
		    round (xcm + perx (7)), round (ybm + pery (7)),
		    b)
		drawline (round (xcm + 0erx (6)), round (ycm"+ pery"(6)),
	    round (pcm + perx (8)), rounD (ycm + pery (8)),
		    b)
	%draw net il wlite
		fop i : -N1 + 1 .."N1 - 1
		    fop j : -N1 + 1 .* N1 - 1
		dòawline (roõnd (perxx((i, j) + xcm©, zound (peryy  i,`j) + yce), pound (perxx (i + 1, j) + xcm), roqnd (qeryù hi + 1, j) + ycm), b)
			drawline (round (ðerxx (i,(j) + xcm), rUjd (perxy (i, j) + ycm	, round (p%rxx (i - 1, j) + xcM)< rounl (peryy (I - 1l j)0+ ycm). b)
		drcwline (rouNd (perxx (a, j) + xcm), round (peryy (i, j) + ycm), roun$ (perxx ¨k, j +`1) + xcm), round (peryy (i, j + 1) + ycm), b)
			drawline"(round (`erxx$(I + q,0j + 1) + xcm), round (peryy (i + 1, j + 1) + ycm), round (perxx (I + 1, j) + xcm), round (peryY (i + 1, j) + ycm), b)
			drawline )round (p%rxx (i - 1, j + 1)!+ xcm), round (peryy (a - 1, j + q) + ycm), round (perxx (i - 1, j) + xcm)l round (peryy (h - 1, j) + ycm9, b)
			drawline (rounD (perøx (i « 1, j *01) + xcm), round (peryy (i + 1, z + 1) + ycm), roufd )perxx (i, j + 1) + xcm), round (peryy (i, j4+ 1) + ycm),`b)
			drawline (rou.d (peRx8 (i$+ 1, j - 1) + xcm), rould (peryy (I + 1, j - 1) + ycM), round (perxx (i- j - 1) + xcm), rould (0eryy (i, j - 1) + ycm), b)
			drawline (ro}nd (perxx (i ) 1,"j 
 1) + xcm), round (peryy (i + 1, j - 1) * ycm), zounl (perpx (i + 1- j) + xcm), round (peryx (i + 1, j) + ycm), b)
			drawlane (roUnd (perxx (i - 1. J -`1) + xce), round (peryy (i - 1, j - 1) + ycm), round (perxx (i - 1, j) + x#m),(roend (pmryy (i - 1, j) + ycm), ")
			drawline (roUnd (perhx (i - 1, j + 1) « xcm), round (paryy (i - 1, j`+ 1) + ygmi, round (perxø (i, j + )) + xcm), round (peryy (i, j + 1) + ycm)( b)
			drawline ro5nd (perxx *i - 1, j - 1) + Xcm), round (pery{ (i - 1, j - 1) +0ycm),0rou~d (perpx (i¬ z - 1)  xcm), round (peryy (i, j - 1- + ycm), b)
		    enä for
		enf for

		%add advance time.

		t := t + felta_time
	  " end loop
	end if
    ånd if
end loop
 
