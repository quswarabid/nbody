# nbody
This program calculates and displays the spatial motion of a system of
N point masses bound by gravitational attraction.

## Overview
The program solves the equations of motion with the *leapfrog* method
with the constant step of integration.

To avoid singularities of the gravitational potential in cases of direct
collisions, the force of attraction is artificially set to zero in the
small neighbourhoods of the points. Therefore, the calculated dynamics
resembles that of a system of small spherical shells capable of passing
freely through each other.

## Requirements
The program is written in the standard C99 and should compile on any
POSIX-compliant system. For interactive graphics it requires nothing but
plain Xlib.

## Input
The initial state of the system is read from the standard input. Each
line specifies the position, velocity, and mass of a point:

    x1 y1 z1    v1_x v1_y v1_z    m1
    x2 y2 z2    v2_x v2_y v2_z    m2
    ...

The step of integration is read from the first command line argument.

## Output
The program draws trajectories of the points with colour representing the
current velocity of a point. Bluer colours mean higher speeds.

## Keyboard controls
key | action
----|-------
q   | Quit
d   | Dump state (time, positions, and velocities)
c   | Check conservation laws (energy and momentum)

## Caveats
The step of integration have to be found by trial and error: for
example, by looking at how good the conservation laws are held.

The calculation of the force exerted on a point by the others is the
most time-consuming part of the program. For the system of N points, it
requires O(N*N) operations, and therefore the program is suitable for
simulations of small systems only, about a dozen points or so.

## Customization
One can easily adjust most parameters, such as phase space size,
projection centre and plane, regularization radius, etc. How to do that,
should be evident from the comments in the source.

## Examples
#### 3 bodies, regular motion
A planet is orbiting a star of a binary system.

![d3-stable](https://user-images.githubusercontent.com/29631214/31562119-e9a6d8a4-b062-11e7-849f-773c66905f45.png)

Input:

    0 0 1000     0 10 0     1000000
    0 0 -1000    0 -10 0    1000000
    0 0 1100     0 100 0    1

#### 3 bodies, transient motion
A planet is orbiting a star of a binary system, then the influence of
the second star outweighs, and the planet jumps to an orbit around the
second star. The planet jumps back to the first star at the next
encounter.

![d3-jump](https://user-images.githubusercontent.com/29631214/31562118-e985d23a-b062-11e7-8b1a-d5355ab7e1e3.png)

Input:

    0 0 1000     0 10 0     1000000
    0 0 -1000    0 -10 0    1000000
    0 0 1200     0 60 0     1

#### 3 bodies, irregular motion
Orbital motion is destroyed by a massive body passing nearby.

![d4-destroyed](https://user-images.githubusercontent.com/29631214/31562120-e9c4af6e-b062-11e7-90fa-d83da2bad49b.png)

Input:

    0 800 -500     0 0 0       1000000
    0 1000 -500    0 0 50      100
    0 0 5000       0 0 -100    10000000

#### 9 bodies
No physical meaning, just symmetric motion with sheer beauty.

![d2-symmetry](https://user-images.githubusercontent.com/29631214/31562117-e95994b8-b062-11e7-8cae-6de9ebd4816d.png)

Input:

    0 -500 -500    0 0 0       100000
    0 -500 0       0 0 50      100000
    0 -500 500     0 0 0       100000
    0 0 -500       0 -50 0     100000
    0 0 0          0 0 0       1000000
    0 0 500        0 50 0      100000
    0 500 -500     0 0 0       100000
    0 500 0        0 0 -50     100000
    0 500 500      0 0 0       100000

## License
The program is in the public domain.
