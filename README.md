# nbody
This program calculates and displays the spatial motion of a system of
N point masses bound by gravitational attraction.

## Requirements
The program should compile on any POSIX-compliant system with X11.

## Input
The initial state of the system is read from the `standard input`. Each
line specifies position, velocity, and mass of a point:

    x1 y1 z1    v1_x v1_y v1_z    m1
    x2 y2 z2    v2_x v2_y v2_z    m2
    ...

## Output
The program draws trajectories of the points with color representing the
current velocity of a point. Bluer colors mean higher speeds.

## Keyboard controls
key | action
----|-------
q   | Quit
d   | Dump state (time, positions, and velocities)
c   | Check conservation laws (energy and momentum)

## Caveats
The calculation of the force exerted on a point by the others is the
most time-consuming part of the program. For the system of N points, it
requires O(N*N) operations, and therefore the program is suitable for
simulation of small systems only, about a dozen points or so.

To solve the equations of motion, the program uses the simplest Euler's
integrator. The step of integration is chosen dynamically to achieve
reasonable accuracy, but because of the low order of the method, the
step required often turns out to be very small. This too impedes the
performance of the program.

To avoid singularities of the gravitational potential in cases of direct
collisions, the force of attraction is artificially set to zero in the
small neighbourhoods of the points. Therefore, the calculated dynamics
resembles that of a system of small spherical shells capable of passing
freely through each other.

## Customization
One can easily adjust most parameters, such as phase space size,
projection centre and plane, target accuracy, regularization radius,
etc. How to do that, should be evident from the comments in the source.

## Examples
#### 3 bodies, regular motion
A planet is orbiting a star of a binary system.

![d3-stable](https://user-images.githubusercontent.com/29631214/31294041-1ebbbb6c-aae2-11e7-8f2b-8a18890c2d2d.png)

Input:

    0 0 1000     0 10 0     1000000
    0 0 -1000    0 -10 0    1000000
    0 0 1100     0 100 0    1

#### 3 bodies, transient motion
A planet is orbiting a star of a binary system, then the influence of
the second star outweighs, and the planet jumps to an orbit around the
second star.

![d3-jump](https://user-images.githubusercontent.com/29631214/31294042-1ecf6590-aae2-11e7-9c84-e881c272e036.png)

Input:

    0 0 1000     0 10 0     1000000
    0 0 -1000    0 -10 0    1000000
    0 0 1200     0 60 0     1

#### 3 bodies, irregular motion
Planetary motion is destroyed by a massive body passing nearby.

![d4-destroy](https://user-images.githubusercontent.com/29631214/31294039-1eb6e0ec-aae2-11e7-846b-224dbc699117.png)

Input:

    0 800 -500     0 0 0       1000000
    0 1000 -500    0 0 50      100
    0 0 5000       0 0 -150    10000000

#### 9 bodies
No physical meaning, just symmetric motion with sheer beauty.

![d2-symmetry](https://user-images.githubusercontent.com/29631214/31294043-1ee7f4ac-aae2-11e7-8023-45ab3c446987.png)

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
