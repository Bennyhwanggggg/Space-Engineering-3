# SpaceEngineering
Space Engineering Course Work

Part 1.1 - Satellite Simulation Using TLE Data and Relation Between Orbital Properties and Mission

Introduction

SARAL and O3B FM07 are both satellites that have been launched within the last three years. Both
satellites are used for earth sensing or communication missions but they have very different
inclination which makes them very interesting to research. This section of the report will investigate
how their orbital properties are adjusted to perform their mission.

Satellite in LEO - SARAL

SARAL stands for Satelite with ARgos and ALtika. It is a Low Earth Orbit satellite. The satellite is
operated by Indian Space Research Organisation (ISRO) and CNES. SARAL was launched on February
23rd, 2013 with the mission of studying ocean circulation and sea surface elevation. To perform the
mission, SARAL is equipped with altimeter and data collection system as its payloads. The AL in its
name stands for ALTIKA as AlLTIKA is the altimeter and also the payload of the satellite. The
altimeter operates at Ka band. The altimetric measurements provided by the satellite will have an
impact on the measurement of sea surface height and reduce mapping error. The data collection is
then done by the ARGOS data collection system so the two systems combined resulted in the name
SARAL.

Orbital Parameters Satellite Specifications
Type of Orbit: Sun-Synchronous Launch Mass: 407 kg
Perigee: 791.1km Expected life time: 5 years
Apogee: 791.6 km
Inclination: 98.5 o
Period: 100.5 minutes
Semi Major Axis: 7162 km
Eccentricity: 0.0000401

TLE:
1 39086U 13009A 16062.90226953 +.00000026 +00000-0 +25963-4 0 9997
2 39086 098.5412 251.8101 0000401 050.0426 310.0793 14.32253629157651

Satellite in MEO - O3b FM07
O3b FM07 is part of the O3b satellite constellation. The main purpose of the satellite constellation
is telecommunications. To perform its duties, the satellite is equipped with twelve Ka band antennas
and provides high speed internet and broadcasts to people in remote area, which is actually what
the name of the "Other 3 Billion" satellite network is referencing to.The O3b FM07 was launched on
July 10, 2014 and today the satellite network has more than 12 in the orbit.
Orbital Parameters Satellite Specification
Perigee: 8069.5 km Launch Mass: 650 kg
Apogee: 8076.6 km
Inclination: 0.0359o
Period: 287.9 minutes
Semi Major Axis: 14444 km
Eccentricity: 0.0002445
TLE
1 40081U 14038C 16062.74115278 -.00000022 +00000-0 +00000-0 0 9994
2 40081 000.0359 353.2540 0002445 324.0770 042.6465 05.00115716029950

Methodology

To simulate the orbit of each satellite in 3D over a 24-hour period, the basic Keplerian orbital model
is used. The seven key orbital parameters required to find the location of satellite in any given time
are inclination (i), right ascension of the ascending node (Ω), eccentricity of the orbit (e), argument
of perigee (ω), mean anomaly of the orbit (Mo), time (t) and the semi major axis of the orbit (a).
These parameters can be found by examine the satellite’s TLE.
After these parameters are defined, they are inputted into the Keplerian orbital calculation function.
This function first calculates the mean anomaly at a given time by using
𝑀 = 𝑀𝑜 + 𝑛(𝑡 − 𝑡0)
where 𝑛 = √𝜇/𝑎3 and 𝑡 − 𝑡0 will equal to the time stamp in the Matlab code. The next step is to
find the eccentric anomaly at the given mean anomaly. This can only be done by using iterative
method to solve
𝐸 − 𝑒 𝑠𝑖𝑛(𝐸) − 𝑀 = 0
The iteration is initialized using E0 = M and the next iteration value is found using
𝐸𝑖+1 = 𝐸𝑖 − (𝐸 − 𝑒 𝑠𝑖𝑛(𝐸) − 𝑀 )/(1 − 𝑒 cos(𝐸))
with the eccentric anomaly, E, found, the true anomaly can then be found by rearranging the equation below
tan (𝜃/2) = √((1 + 𝑒)/(1 − 𝑒)) x tan (𝐸/2)
the true anomaly gives the angle the satellite is at in the orbit so now we just need to find the
distance to get the location of the satellite at a certain time in orbit. The distance r is found by using
𝑟 = 𝑝/(1+𝑒 𝑐𝑜𝑠𝜃)
where p is the semilatus rectum and can be found using p = a x (1-e2). With distance and true
anomaly know, we know have the Perifocal frame in polar form and we can convert this to Cartesian
form easily as we have all the parameters required to do so. We are then required to transfer the
coordinates from Perifocal frame to ECI frame which can also be performed by using the transfer
matrix and the result will give us the 3D orbit for satellite. We have also set the time step in the
Matlab code as well so it only plots the orbit over a 24 hour period.
The next step is to plot a ground trace of the orbit. To do this, we simply need to find the satellite's
position in ECEF coordinates. This is simply executed by the eci2ecef function in Matlab which
utilizes the transfer matrix between the two coordinates systems.


Part 1.2 - Orbital Simulation with Perturbation

Introduction

Orbital perturbation can come from many different sources such as Earth oblateness, gravity
harmonics, solar/lunar gravity forces, aerodynamic drag and solar radiation. Perturbations are
external disturbances that affect the ideal orbital dynamics. The magnitude of perturbations of each
different source depends on the altitude of the orbit. Gravitational forces are the most significant
perturbation especially on LEO satellites as they have lower altitude. On the other hand,
perturbations caused by atmospheric drag decreases significantly as altitude increase since there is
less resistance in outer space. Another perturbation is the attraction force of the moon and sun.
These perturbations increase slowly as altitude increase since the higher the altitude, the closer the
satellite will be to them and there will be stronger attraction forces between them. This also stands
true for other planet attraction forces which can cause perturbations as well, but they are less
significant as they are not as close to the satellite compared to the moon or Earth nor do they have
the same attraction force as the sun.

The dominant perturbing source that will be investigated is called J2 zonal perturbations. It is the
most significant perturbation source after gravitational forces especially for LEO satellites. The effect
of J2 perturbation drops with increase in altitude but is still significant for a large range of altitudes
that are close to Earth. J2 perturbation is investigated here and not gravitational forces perturbation
because J2 perturbation is more related to the orbital parameters of the satellite. J2 zonal
perturbation is the result of the effect of Earth oblateness. This is mostly due to the Earth's spin and
results in a gradual change in the ascension of the ascending node and argument of perigee. The
change in ascension of the ascending node occurs because there is more attraction in Earth's
equator since the poles are flatter so the satellite accelerates towards the equator more. As a result,
this extra acceleration also causes the change in argument of perigee since the orbit is no longer a
closed ellipse due to the change in force components.

Methodology

To simulate the modified perturbation orbit model, a starting true anomaly is first found using the
same method as Part 1. With the true anomaly of the starting position found, we can now
compute the equinoctial elements model since we have all the classical orbital parameters required
to calculate the orbit. The conversions equations can be simply applied and are shown in appendix.
This gives us the six parameters required and we need to introduce perturbations into equinoctial
model. To do this we express the J2 perturbation as a perturbing acceleration on the two-body
equation of motion. This acceleration is express as 𝑥̇ and is known as the sate rate equations.
𝑥̇ = 𝑨(𝑥)Δ(𝑥) + 𝒃(𝑥)
and x matrix containing the 6 equinoctial elements that we have found already. The matrix A(x),
Δ(𝑥) and b(x) are expressed in the Matlab code provided in appendix. The Δ(𝑥) matrix is specific for
J2 perturbation and its properties are only affected by the orbital parameters.
The equinoctial model computes the equinoctial elements at each time step and is solved by using
the Runge-Kutta method.

Runge-Kutta method is a type of integration scheme that uses the sate rate equation. To apply it, the
variable k1 is first found by using
𝑘1 = Δ𝑡 𝑓(𝑥𝑘 )
This is equivalent to the integration of the sate rate equation between the time intervals. The
parameters are then updated by adding half k1 and substituted backed into the sate rate equation
to find k2
𝑘2 = Δ𝑡 𝑓 (𝑥𝑘 + 1/2 * 𝑘1)
following the same logic, k3 and k4 are found
𝑘3 = Δ𝑡 𝑓 (𝑥𝑘 + 1/2 * 𝑘2)
𝑘4 = Δ𝑡 𝑓(𝑥𝑘 + 𝑘3)
with these four variables found, we can find the equinoctial elements of the next time step by
𝑥𝑘+1 = 𝑥𝑘 +1/6 * (𝑘1 + 2𝑘2 + 2𝑘3 + 𝑘4)
By repeating this in loops, we can find the equinoctial elements of each time steps. We then convert
them back to ECI position vectors for plotting.

Part 1.3 - Finding the Position of Satellite in Orbit

Introduction

To find a position of a satellite in orbit, the satellite usually sends information about its current
position in orbit to the ground station. This orbital information, such as range, azimuth and elevation
observations, are received by the ground station and used to calculate satellites position. However,
when signals are travelling from satellite to ground station, it has to compete with background noise.
These background noises can come from many different sources. They are either internal or external
source. Internal sources can come from antennas and cannot be entirely eliminated. External
sources include solar radiation, Earth's atmosphere or other natural sources. The amount of noise is
determined by the signal to noise ratio. This ratio is also normally expressed in terms of decibels (dB).
In this section, we will produce a plot of satellite's range, azimuth and elevation throughout an orbit
without perturbation then extract three observations to be the signal data that will be sent to the
ground station. Noise is added onto the signal data before it is sent. In the end, we will obtain
information on the error due to the noise so we can investigate its effect.

Method

First, the ground station was chosen to be Sydney and using orbital parameters from TLE of SARAL
the satellite's ECI position can be found easily using similar method from question one. With
satellite's ECI position found, we now need to find the satellite's relative position to the ground
station if we want to obtain the range, azimuth and elevation. This can be done by using simple
vector subtraction of satellite's ECI position with ground station's ECI position. The relative position
can then be converted from ECI to ECEF then to LG. At last we will obtain the required observation
vector of range, azimuth and elevation by converting from LG Cartesian to LG Polar form.
We are required to identify the period when the satellite can be observed by the ground station in
Earth. This can be estimated by considering that the satellite will be observed when the elevation is
between 0 and 180 degrees. This estimate is not entirely accurate since it considers Earth as flat
which is obviously not the real case.

After we have obtained the observation parameters, White Gaussian noise is added onto those
parameters to indicate that the signal is competing with the noise while it's been transmitted to
ground station. This is done by using the awgn() function in Matlab. (See appendix for more detail)
After it is received, the ground station will use Herrick Gibbs technique to estimate the orbital
parameters. The Herrick Gibbs technique requires the computation of ground station's ECI position
first. This is easily done by converting llhgd to ECEF then ECEF to ECI using the three-time location
where the observation takes place. We then need to convert the data received from polar
coordinates to Cartesian then we can convert it to ECEF then ECI. This gives us the satellite's position
relative to ground station in ECI but we need the satellite's ECI position relative to geodetic center of
Earth so we add the two components together. Therefore, we end up with our r1, r2 and r3 position
vectors.

In Herrick Gibbs technique, we now compute the velocity of the second observation using r1,r2 and
r3. This is given by the formula
𝑣2 = −𝑡32 (1/𝑡21*𝑡31 + 𝜇/12𝑟1^3) 𝒓𝟏 + (𝑡32 − 𝑡21) (1/𝑡21*𝑡32 + 𝜇/12𝑟2^3) 𝒓𝟐 + 𝑡21 (1/𝑡32𝑡31 + 𝜇/12𝑟3^3) 𝒓𝟑

With v2 found, we have all the variables required to find the orbital parameters of tracking and we
simply just need to substitute v2 and r2 in.

To find the error, we subtract the predicted orbital parameters from the original parameters from
TLE. Some errors are measured in terms of percentage and some errors are quantitative. This will be
further explained in discussion. The variation of the error with increase in noise was also tested.
