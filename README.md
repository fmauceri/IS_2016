# IS_2016
###### Grace Handler, Stephen Majercik, Frank Mauceri, Jack Truskowski


![Boids](https://cloud.githubusercontent.com/assets/11000833/12837197/170d39a2-cb8e-11e5-84c1-273e1cb236c9.png)


Instructions for running the simulation can be found in the `ControlPanel...maxpat` file

####Future Work:
- Visually connect boids in same neighborhood (output a matrix from external with each boid and its neighbors) This will make neighborhoods / clusters more visually obvious  
- Implement Quicksort in the external for performance (currently doing a jitter bubble sort in the max patch)  
- Fix visual jump (Improve efficiency in message passing / handling of FlockIDs?)  
- Movable / Multiple Attractors  
- Implement the deletion of boids (age?)  
- Revisit the statistics used to characterize flock behavior  
- Stipulate where boids are born in space when they are added

####Questions:
- Is naive wall deflection okay, or should the velocity change as boids approach walls?
- Make a maximum number of boids? This way we don't have to dynamically make a new array everytime we want to change the number of boids
