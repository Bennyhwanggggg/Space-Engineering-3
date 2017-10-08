# SpaceEngineering
Space Engineering Course Work

Part 1
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

Part 2

Introduction

Global Navigation Satellite System (GNSS) utilizes a receiver and a constellation of satellites to help automobiles, aircrafts or just people in general locate themselves on Earth so they can navigate themselves to their desire location. The constellation usually consists of at least 24 space vehicles (SVs) in a semi-synchronous (12 hour) orbits and to locate the position of the receiver. To locate the position of the receiver, which is also the user location, at least 4 pseudorange measurements are required. This means that at least four SVs need to be in the view of the user to provide the pseudorange data. In this question, we have been given the ephemeris data of the space vehicles and the pseudorange data. The ephemeris data provides information on the orbit of the space vehicles, while the pseudorange data provides information about the pseudorange and the time when the pseudorange was taken and which satellite was the pseudorange measurement taken from. First, the orbit of the space vehicles, which are actually satellites, are computed in both ECI and ECEF coordinates using the Keplerian orbital parameters provided ephemeris data and simulated over a 12 hours’ period. The 2D ground track will also be simulated as well. Next, the position of the receiver, which is a UAV in this case, will be estimated using nonlinear least square technique with the pseudorange data provided. Since the true position of the UAV with respect to ground is provided, we will compare our estimated trajectory with the true trajectory by using various plots such as 3D UAV position w.r.t Ground Station, 2D UAV Position w.r.t Ground Station, polar plots of the UAV’s trajectory w.r.t the Ground Station, UAV altitude vs Time w.r.t Ground Station, UAV Clock Bias vs Time w.r.t Ground Station, Number of Visible Satellite vs Time, Range, Azimuth, Elevation vs Time w.r.t Ground Station, various DOP vs Time plots. In part c, we are given a pseudorange file where there is an error in this data and we are required to create an error detection and correction method.

Part 2.1 - Satellite Constellation Simulation

The 3D ECI and ECEF position of each GPS satellite and their 2D ground track are shown in the Appendix Question 1(a) section. Since the Keplerian orbital parameters are provided, the position of the satellite in ECI can be easily calculated by finding the mean anomaly at a time since epoch and use Kepler’s equation to solve for the eccentric anomaly, E, at that time. With the eccentric anomaly known, the true anomaly at that time can be computed so that the distance can be found. The location in orbital frame can then be computed as distance and true anomaly are known. To compute the ECI coordinates and ECEF coordinates, the conversion functions for orbital frame to ECI is used first then the conversion function from ECI to ECEF is used. To generate the 2D ground trace, the ECEF position is converted into longitude, latitude and height using the conversion function previously created. To be able to store multiple satellite positions at multiple time, the Matlab script created uses a 3D matrix where the time is stored in the x axis, satellite number is stored in y and the positions are stores in z. This code allows Matlab to access multiple satellite and their respective positions at each time and be plotted.

Part 2.2 - UAV tracking

the pseudorange data are extracted and stored into their respective variables and the ground station location is stored in longitude, latitude and height coordinates and ECEF coordinates. To convert from ground station longitude, latitude and height coordinates to ECEF, the conversion code created previously for geocentric is used as the true positions of vehicle provided in the true positions file are in LGCV. If ground station coordinate conversion was done using geodetic coordinates, it would cause a translation problem for the estimated trajectory later on. The pseudorange file shows that not all satellite is visible to UAV at all time as the satellite numbers are indicated for each pseudorange data and the time the data is taken is given as well. The first step to estimate the UAV position is to find the position of the satellites that are observing the UAV. This is done using the Matlab code for part A but this time we specify which satellites we are interested in and the time in the orbit. As a result, we can retrieve information on the satellite positions. The Matlab code stores each position with respect to the row number of the pseudorange data. For example, the first row of the pseudorange data shows time at 7.3477373568*10^6 seconds and satellite number is 2, the Matlab code will store the xsv, ysv, zsv values for this row. As a result, the total number of rows the variable xsv, ysv and zsv will have depends on the total number of rows the pseudorange file has. The xsv, ysv and zsv are also in ECEF coordinates and these results are stored in another variable and converted to polar positions w.r.t ground station by first finding the Cartesian coordinates w.r.t ground station for later analysis purposes. The next step is to sort out the pseudorange data so that we can see how many satellites are viewing the UAV at each time. From the pseudorange file, it is obvious that the data with the same time has already been grouped together. Therefore, we only need to detect the index number where the time is not the same as the row before. By detecting those changes, we find that there are 300 data sets. However, for the nonlinear least square method to work, we are required to have four pseudorange as we have four unknowns, which are the x,y,z coordinates and clock bias. As a result, we need each data set to have at least 4 satellite viewing to provide sufficient amount of pseudorange data. We can easily do this by using an if statement in Matlab where if the length of the data set is less than four, that data set would not be valid. The find_t_repeat.m function basically performs the task about but it only extracts the number of satellites of each data set and the time when pseudorange is measured. The reason that this was not done in the pseudorange.m function is because this makes the pseudorange.m function simpler and easier to read. The pseudorange.m file performs the nonlinear least square method to calculate the estimated position of the UAV. In nonlinear least square method, the pseudorange function if first used to compute f(x). The pseudorange function is describe as 
pi = √((Xsv-x)^2 +(Ysv-y)^2 + (Zsv-z)^2) + cb
Where Xsv, Ysv and Zsv are the specific visible satellite position at that time and x,y,z and clock bias (cbu) are the solutions we are trying to solve for so we initialize the starting parameters as [0; 0; 0; 0]. Therefore, to find Δp, the pseudorange for that row is subtracted by the results of the equation shown above. The Jacobian matrix is then found using the relevant formula. Overall, the Matlab code used calculates the Δp value and the H of each row separately and groups them into their respective data set. The Δx value is then found as well using the relevant formula as we already have the required variables to compute it, which are Δp of the data set and the Jacobian Matrix H of the data set. The Δx is then added to the initial iteration value and the whole process is repeated until Δx value is very small (less than 0.0000001 in the Matlab script, as a small Δx means that the result has converged and we will have obtained our UAV estimated position for that time. We then move on to the next data set and repeat the whole process to compute the UAV estimated position for another time. After obtaining the estimated position, they are converted to LGCV so it can be compared with the true positions.

Part 3
Part 3.1 - Rover path find using A* algorithm and Dijkstra algorithm.

Introduction

An attempt was made to plan the path of Mars Rover, Orville, so that it can explore the local terrain of the site it landed on. Motion planning is essential for the Orville rover to operate as it needs to be able to plan a path to move from its initial position to the desired destination. The essential criteria for the optimal path are the time taken and distance travelled and due to the nature of a terrain, some regions are not traversable due to obstacles or would cost a lot of time or energy to travel through due to the environmental conditions. This traversability of a region is often defined by a cost where the higher the cost, the more distance, time or resources would be required to travel in that region. Therefore, the optimal path desired would have the lowest possible cost. In grid based path planning, the terrain is divided into cells and the traversability of each cell is found. Two types of path planning algorithm, the Dijkstra algorithm and A* algorithm, are then implemented and compared. These two algorithm both takes account into the geography of the terrain and uses a similar approach where the path with the minimal cost is found in the end. However, they do have some differences which will be discussed further in a later section of this section of the report.

Methodology

Cell 1: Martian Terrain Map

The first to execute grid based robot motion planning is to give a visualization of the terrain. This image is produced by plotting the data provided in localTerrain.mat
gives an imagery of the terrain point cloud produced by the localTerrainPointCloud.mat file. Point clouds are used to define the surface of the terrain in x,y,z coordinates normally which is more useful compared to just the data provided in the localTerrain.mat as point clouds can be very useful in other calculations later on.

Cell2: Gridded Map

The next step is to create a gridded map that contains cells. Each cell would represent a small area of the terrain in 2D with a top view. The cells were defined to have a size of 20cm by 20cm as this was around an area of a single Orville rover wheel. To generate the cells, the dimension of the map has to be defined first which was a 10m by 10m map. With the map dimension and cell dimension known we create a function called generateGrid. This function will convert coordinate positions into cell positions. Each cell would cover a small area on the map so that coordinate positions that are close together may all belong in the same cell. In order to create the gridded map, the total number of cells in the map is first find by using the map dimension and cell dimension where the center of the cell is used to define the position of the cell. As a result, a 51 by 51 matrix is created. A map of ones which is 51 by 51, as this is the number of cells required, is created. The number one, is the enumeration number that represents the color white. This color enumeration code was already created when the code is provided by extracting RGB values of each corresponding color, which is white, blue, red, black, cyan, orange and yellow from the color map. The starting position and goal position are then found by converting their respective coordinates into cell positions. This operation was done by the pos2cell function and the results were rounded as it does not make sense to have decimal places for cell positions as we only want to know which cell the positions are in
and not where in the cell the position is in. The starting position and goal position are then marked on the map with their corresponding color of blue and red respectively. A key note to take is that the x and y axis are swapped in this case for plotting purposes and the direction of the x axis is reversed as it is increasing when going downwards. Therefore, the x and y points are actually swapped after the cell positions are found. As a result, in the Matlab code, the Y has been defined as the first index of the cell position, where the first index is actually the x coordinate of the cell position. The same theory applied for the second index where X is defined but it is actually the y coordinate of the cell position.

Cell3: Configuration Space

The hazardous areas on the terrain are considered as obstacles and the rovers should avoid it. The coordinates of the obstacles and their size were provided. To generate the obstacles, the rover’s footprint radius has to be considered first and converted into size in terms of cells. Since the cells are always squares, the cell size occupied by the rover’s footprint is found by simply dividing the radius of footprint (0.5m) by one side of the cell dimension (0.2m). The result is 2.5 cells. The next step is to use mesh grid and then the circles formula and then find the indices of the obstacles if it were in a matrix that was the same size of the map. These indices are then imported into our actual map, which is the gridded map created in Cell 2. What the circles formula does is that the xx from mesh grid would be subtracted to the obstacles x coordinate and then squared and added to the y component which is calculated by a similar manner. The result is square rooted to find a radius and if this radius is less than the sum of the rover’s footprint radius and obstacle’s radius, it means that the cell position is considered as hazardous so it would be marked black.

Cell4: The GESTALT Algorithm

The GESTALT Algorithm generates a traversability map that also considers the hazardous area which are not traversable. The algorithm assigns a traversability cost to each cell, which will be between 0 and 255 with the lower score been more traversable. Since the rover is defined by a point and the rover actually has a radius defined by its footprint radius, the rangesearch function is used to find all the point cloud that are covered by the range of the footprint radius. This operation is performed by the getRoverFootPrintPoints function which shows returns the points in the terrainPointCloud.mat file that is covered by the footprint of the rover.
To be able to evaluate the hazard of the terrain, we need to fit a plane onto the set of points covered by the rover’s footprint. This is performed by the fitPlane function which finds the average of the points and the unit normal vector of the plane. With this results, three types of hazard evaluation are performed. The first one evaluates step obstacles which finds the maximum height difference between the points and checks if it above the clearance height divided by 3 and if it is the step hazard for that cell would be computed, otherwise it is zero. The second hazard is from the roughness of the terrain and is computed for every cell. The third hazard is the pitch hazard which indicate how steep the terrain may be by using the fitted plane found previously. The final traversability score is then the maximum of the three hazards mentioned above. Another note is that the cells that were already not traversable, which are the obstacles, or cells that have a higher traversability score than the traversability score will have a traversability of infinite as
they cannot be passed through.

Cell 5: Dijkstra’s Algorithm

To perform path planning, a function that finds the neighboring cells is created first. This function finds all the neighboring cells of the current node and returns them. However, it will not return neighboring cells that were previously already visited (previous nodes), or if the cell is outside the map’s dimension or if the neighboring cell is an obstacle. The last two conditions were applied by simply checking if the indices of the neighboring cells have the corresponding color to obstacles or visited nodes. The first condition is applied by setting the map dimension as a condition when finding neighboring cells. The neighboring cells are all also assigned in index number which is used to reference their position with respect to the visited node.
The Dijkstra’s algorithm finds the shortest path by finding computing the path with the lowest cost after a breadth first search over the traversability map. In the traversability map, the traversability score assigned to each cell is the cost to travel to the center of the cell. The breadth first search is implemented by keep finding the traversable neighbor cells (Not obstacles or cells out of boundary starting with the starting cell) until the goal cell is found. The Dijkstra’s algorithm is similar to Breadth- first algorithm. Each traversable neighboring cells are added to the open list (marked by the cyan color in animations and figure). The cells in open list means that it is under consideration for that the cell may be part of the shortest path to goal cell. However, since there may be many possible paths to travel to each cell, the lowest cost becomes a variable that is compared every time and when another shorter path travelling to a cell in open set is found, that cell is removed from open set and put into closed set (marked by orange color in animations and figure). Cells that are in the closed set would no longer be reviewed again as it is unnecessary since the shortest path to the cell has already been found. In this context, shortest means that the total traversability score is the lowest.
In this algorithm, the starting cell’s cost is set to zero as it is obvious that there is no cost since the rover has not travelled anywhere. The cost to travel to each neighboring cell is computed by finding the sum of the cell’s traversability score divided by two and then multiplied by the proportional distance. For example, the proportional distance of moving to a neighboring cell that is in a diagonal position would be √2 and would simply be 1 if it is not a diagonal movement. As the algorithm computes the cost it takes to travel to the cell from starting position, the total cost would keep growing until the goal position is found and the path with lowest cost would be computed. The parents array stores the parent node of the neighboring cells to help locate new cells and assist in finding the shortest and lowest cost path. The parent array does this by comparing the cost of the neighbor cells with the parent cell and the neighbor cell becomes the parent cell if it has the lower cost. As a result, it would have become the new shorter path which means the costs would be recalculated. [1] The whole process repeats until the open set becomes empty or the goal cell is found. In short, the parent matrix to define the path and in the end the when goal cell is added to open list, the final route will be defined. If the open list becomes empty, there is no possible path for the rover to reach its goal. In Matlab, an if statement that checks whether or not the cost to the goal node is infinite is implemented to check if the goal destination is reached or not. 

Cell 6: A* Algorithm

The A* Algorithm is similar to Dijkstra’s Algorithm where the theory of open set and closed set are still used and still operates in the same way. The main difference between the two
algorithm is the use of two cost function where one is the cost to come cost function (the same cost function implemented in Dijkstra) and the other is the heuristic function. The heuristic function is a function that estimates the cost to go value of each cell and there are many methods that can be used to perform this estimation. The Manhattan and Euclidean distance heuristics Matlab code were provided but I decided to go with the Euclidean method. The Euclidean distance in two dimensional Euclidean space is essentially the same as the distance found by the Pythagorean theorem. With the cost to come values and the cost to go values computed, a matrix, f, is used to represent the sum of the two costs. This sum of the costs becomes the variable that is compared when trying to find the lowest cost path. Just like Dijkstra’s Algorithm, if f is lower, the cell becomes the parent and g and f are both recomputed and the lowest one is placed into the close set. The process repeats until the goal is added into the open set and then the parent matrix is used to find the shortest path. [2] [1] Once again, if the open set is empty at any point in time, it means there is no path. The method used to compute this route and its cost is the same as Dijkstra.

Cell 7: Plotting the Path

This section simply plots the path found on to the terrain map by finding the equivalent cell of the route and convert them into coordinate positions and the z axis value are found by finding the corresponding point clouds
Extension: Geological Survey Expedition
An attempt was made to put the Orville rover into survey expedition by finding sets of efficient and traversable paths between each of the sites of interest. The sites of interest’s coordinates were provided in sitesOfInterest.mat file. The main grid path planner algorithm chosen was the A* algorithm as research showed that it is faster and more efficient compared to Dijkstra. The reason for this would be explained in a later section. A problem that has to be noted here is that the assignment sheet states that there are five sites of interest but the sitesOfInterest.mat file provided actually only have information of four sites of interest. Nevertheless, the Matlab code implemented will work whatever the number of sites of interest is.
The method implemented uses a double for loop where each sites are assigned a number based on their row index and it will take one site of interest as a starting position (a) and another one as a goal position (b). The condition a<b is used here to ensure that no repeated paths are calculated and no starting position and goal positions are the same. The logic behind is this that when a=1, the b loop will go through the rest of the sites by considering the rest as goal positions starting from the second sites of interest. When a = 2, it will then not go through the rest of the calculation if b = 1 or 2 as the site combination and the start position and goal position are the same respectively. Before finding a path, the code implemented checks the location of the sites of interest and if it is located at a hazardous area (obstacles) the rover will not travel there, thus no path to and from that site of interest will be created. The first step taken to find the optimal path sets was to convert the sites’ position from coordinate system to cell coordinates by using the functions previously created to do the same operation for starting position and goal position. After this step, the A* algorithm is implemented by using the same Matlab code from Cell6. However, after the algorithm is finished and the optimal path between the two sites are found, the result is stored in a cell array instead to be used later on. The code will then move on to the next combination of starting position and goal position. The plotting methods is the same as before except the lowest cost route information has to be extracted from the cell array first this time. 

Part 3.2 - Two Burn Orbits Simulation

Introduction

Optimization is a process of finding the optimal parameters input that can maximize a system’s performance by minimizing the cost of operation. In a two-burn orbit transfer operation, optimization techniques are implemented to minimize fuel usage or operation time. Therefore, in the context of two-burn orbit transfer, it would be optimizing parameters such as change in eccentric anomaly at each burn, the required magnitude of change in velocity at each burn, and the elevation and azimuth angles required for each burn. These parameters are optimized in order to help us find the minimal total change in velocity required to complete the two-burn orbit transfer so that fuel usage is minimized. The optimization process uses a nonlinear programming method called Sequential Quadratic Programming, which is essentially an improvement over Quadratic Programming.

Methodology

The relationship between total ΔV minimization and fuel minimization is shown in Tsiolkovosky’s Rocket Equation where ΔV = 𝑣𝑒 𝑙𝑛𝑚𝑜𝑚𝑓. From this equation, it is obvious that the smaller the ΔV is, the less mass difference there is between the initial mass and final mass assuming the exhaust velocity remains constant. Therefore, there is less fuel used when ΔV is smaller since the fuel mass is the only variable changing as rocket mass remains constant.
To minimize total ΔV, the parameters that are related to it are considered. Each ΔV represents the velocity change that will be caused by the burn such that the change in eccentric anomaly where this occurs, the magnitude of change in velocity and its direction will have significant effect on the amount of fuel used. The first step in optimization is to form an objective function and since the aim of this optimization is to minimize total ΔV, the objective function is simply the sum of the velocity changes, F = ΔV1+ ΔV2. Another key function is the constraint function. The constraints are defined by the boundary conditions of the final orbit such that the first constraint is the semi major axis of the target orbit, the second constraint is the velocity required to coast in target orbit, the third ensures circular orbit and the fourth and fifth constraint are used to ensure the desired inclination is achieved.

Constraints: 
1. |𝑟𝑓𝑖𝑛𝑎𝑙|/|𝑟𝑡𝑎𝑟𝑔𝑒𝑡|−1=0
2. |𝑣𝑓𝑖𝑛𝑎𝑙|/|𝑣𝑡𝑎𝑟𝑔𝑒𝑡|−1=0
3. 𝑣𝑓𝑖𝑛𝑎𝑙′𝑟𝑓𝑖𝑛𝑎𝑙/|𝑣𝑓𝑖𝑛𝑎𝑙||𝑟𝑓𝑖𝑛𝑎𝑙|=0
4. 𝑣𝑓𝑖𝑛𝑎𝑙/|𝑣𝑓𝑖𝑛𝑎𝑙|−𝑣𝑡𝑎𝑟𝑔𝑒𝑡/|𝑣𝑡𝑎𝑟𝑔𝑒𝑡|= 0
5. 𝑟𝑓𝑖𝑛𝑎𝑙/|𝑟𝑓𝑖𝑛𝑎𝑙|−𝑟𝑡𝑎𝑟𝑔𝑒𝑡/|𝑟𝑡𝑎𝑟𝑔𝑒𝑡|=0

In the setting up process, users are asked to input the desired optimization technique, such as forward differencing or central differencing and the Hessian update method of either BFGS or second order differencing,
they want to use and the properties of the desired orbit, such as its inclination and right ascension of ascending node. The user can also choose to whether or not to use merit function. However, it is noted that using the merit function is a perquisite for BFGS. After initializing the initial and final orbit parameters, their ECI positions are calculated using Kepler’s equations so that comparisons can be made later on. Furthermore, the final orbit’s semi-major axis and the velocity required are noted using simple orbital mechanics equation. The position and velocity vectors of initial satellite parking orbit are then also found as well to help with initialization of Sequential Quadratic Programming. With these variables found, the user is given an option to choose a very good initial estimate or a very bad one. The effect of a very bad initial estimate will be discussed later. With the parameters initialized, a scaling factor is introduced to assist with numerical stability and plotting purposes.
The initialized parameters are then used in a dynamic model to propagate the satellite trajectory such that the calculation of the objective function, constraint function and the augmented lagrangian occurs here. The user is given a choice to use merit function or not and by choosing to use it, the way λ is calculated would be different as a different formula related to the objective and constraint functions derivative is used. In the dynamic model, the position vector and velocity vector of satellite’s initial coast in the parking orbit before first burn is calculated using universal conic section model with the initial guess of change eccentric anomaly. This guess of change in eccentric anomaly is an approximation of where this burn will occur. With the position and velocity vector before first burn found, the change in velocity required for first burn is determined. Since an initial estimate of change in velocity is made, this velocity is converted to ECI coordinates using 𝜃1and 𝜑1, which are parameters to be optimized, and then added onto the velocity vector determined previously. As a result, the initial velocity vector of the satellite in transfer orbit while the position vector still remains the same assuming the change in velocity is instantaneous. The position vector and velocity vector before transitioning into the target GEO orbit is then found using universal conic section again with the use of the optimization parameter ΔE2. The second change in velocity required is then found using the same method that is used to calculate the first one but with their respective angle parameters. As a result, the velocity vector and position vectors for GEO orbit is obtained. With these position and velocity vectors found, they are put into the constraint function, which is defined by the boundary conditions mentioned previously, objective function, which is described earlier in the report but the scaling factor would also need to be included here for numerical stability, and the augmented lagrangian which is calculated by 𝐿=𝑓−𝜆′∗𝑐+ 𝜌∗𝑐′∗𝑐 where 𝜌 is the penalty parameter that was initialized in the code and f is the objective while c is describes the constraints. With the objective value, augmented lagrangian value and constraints value found, the next step is to calculate their derivatives using either forward differencing or central differencing. The code implemented provides the user to choose either. The major difference between the two scheme is that forward differencing error is linear while central differencing error is quadratic and central differences is also more expensive and sensitive due to the way its calculated. Forward differencing only considers the forward perturbation while central differencing considers perturbation in both directions. All further related discussion made in this report will assume that central differencing scheme was used.
For equality constrained optimization, it is important that the local minimum is calculated and not the local maximum. To perform this check, the second order derivatives should be positive definite and the second derivative is described by the Hessian matrix. After this check is perform, the search direction is determined using the KKT system. The KKT conditions are the necessary condition for a solution in nonlinear programming to be optimal
and utilizes the Lagrange multipliers, the constraint function gradient and objective function gradient. It is a representation of the quadratic objective functions with linear constraint equations. This method of finding update direction is essentially the Newton-method. Another method is steepest descent, which simply uses the opposite direction of the gradient. With search direction found, line search technique is then used to find the size of iterative step size. The step size, α, is initialized as a small value so that f is always decreasing in the search direction. The line search algorithm implemented is a line search algorithm that satisfied the Wolfe conditions. The Wolfe conditions have two components with the first being that there has to be sufficient decrease in objective function and the second being a curvature condition. More on Wolfe condition will be discussed later on. The other backtracking algorithm uses Goldstein condition, which was not implemented.
With the step size found, the optimization parameters are updated and the new Hessian matrix and lambda are found. Two Hessian matrix update method are available for user to choose. The first method is a straight forward method that calculates uses forward differencing to calculate the second derivative which is the definition of the Hessian matrix. The second method is the Quasi-Newton method, which is a low rank approximation method that approximate the Hessian recursively based on old and new gradient and function values. The algorithm implemented for Quasi Newton method is the BFGS method and it preserves symmetry and positive definiteness of Hessian. With the new Hessian matrix calculated, the error between the old optimization parameter and new parameter are found and the whole process is repeated until error is very small indicating that convergence has occurred.
